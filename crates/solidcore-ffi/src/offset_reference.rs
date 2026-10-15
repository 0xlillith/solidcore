//! Process boundary for the unchanged C++ offset oracle (not kernel code).
//!
//! Every call, including serialization during generation, starts a fresh executable.
//! A timeout or crash is an oracle error, never a fabricated numerical result or an
//! `Unported` result. Only the child created here is killed, and it is always reaped.

use crate::{check, FfiError, SC_ERR_ARG};
use std::io::{Read, Write};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// Upstream did not finish within the per-call oracle deadline.
pub const TIMEOUT: i32 = -1062;
/// The reference child terminated abnormally or its transport failed.
pub const PROCESS_ERROR: i32 = -1063;
const DEADLINE: Duration = Duration::from_secs(5);
const MAX_OUTPUT: u64 = 12 + 8 * 8_388_608;

fn process_error(_: impl std::fmt::Debug) -> FfiError {
    FfiError::Status(PROCESS_ERROR)
}

pub(crate) fn call(op: i32, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    if input.is_empty() || input.len() > 1_048_576 {
        return Err(FfiError::InvalidArgument("offset worker input size".into()));
    }
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&op.to_ne_bytes());
    bytes.extend_from_slice(&(input.len() as u64).to_ne_bytes());
    for value in input {
        bytes.extend_from_slice(&value.to_ne_bytes());
    }
    let mut child = Command::new(env!("SOLIDCORE_OFFSET_WORKER"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(process_error)?;
    let start = Instant::now();
    let result = std::thread::scope(|scope| {
        let stdin = child.stdin.take().ok_or_else(|| process_error("stdin"))?;
        let stdout = child.stdout.take().ok_or_else(|| process_error("stdout"))?;
        let writer = std::thread::Builder::new()
            .name("offset-oracle-input".into())
            .spawn_scoped(scope, move || {
                let mut stdin = stdin;
                stdin.write_all(&bytes)
            })
            .map_err(process_error)?;
        let reader = match std::thread::Builder::new()
            .name("offset-oracle-output".into())
            .spawn_scoped(scope, move || {
                let mut output = Vec::new();
                stdout.take(MAX_OUTPUT + 1).read_to_end(&mut output)?;
                Ok::<_, std::io::Error>(output)
            }) {
            Ok(reader) => reader,
            Err(error) => {
                let _ = child.kill();
                let _ = child.wait();
                return Err(process_error(error));
            }
        };
        let status = loop {
            match child.try_wait() {
                Ok(Some(status)) => break Ok(status),
                Ok(None) if start.elapsed() < DEADLINE => {
                    std::thread::sleep(Duration::from_millis(1));
                }
                Ok(None) => {
                    let _ = child.kill();
                    let _ = child.wait();
                    break Err(FfiError::Status(TIMEOUT));
                }
                Err(error) => {
                    let _ = child.kill();
                    let _ = child.wait();
                    break Err(process_error(error));
                }
            }
        };
        let written = writer.join().map_err(process_error)?;
        let output = reader.join().map_err(process_error)?;
        let status = status?;
        if !status.success() {
            return Err(process_error(status));
        }
        written.map_err(process_error)?;
        output.map_err(process_error)
    });
    // Also covers thread creation/pipe setup errors before the wait loop.
    if result.is_err() {
        let _ = child.kill();
        let _ = child.wait();
    }
    let output = result?;
    decode(&output)
}

fn decode(output: &[u8]) -> Result<Vec<f64>, FfiError> {
    let status = i32::from_ne_bytes(
        output
            .get(..4)
            .ok_or_else(|| process_error("status"))?
            .try_into()
            .map_err(process_error)?,
    );
    if status == SC_ERR_ARG {
        return Err(FfiError::InvalidArgument(
            "offset reference rejected unsafe geometry/state".into(),
        ));
    }
    check(status)?;
    let count = u64::from_ne_bytes(
        output
            .get(4..12)
            .ok_or_else(|| process_error("length"))?
            .try_into()
            .map_err(process_error)?,
    );
    let body = output.get(12..).ok_or_else(|| process_error("body"))?;
    if count > 8_388_608 || count * 8 != body.len() as u64 {
        return Err(process_error("output length"));
    }
    Ok(body
        .as_chunks::<8>()
        .0
        .iter()
        .map(|&bytes| f64::from_ne_bytes(bytes))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_truncated_and_oversized_worker_records() {
        assert!(decode(&[]).is_err());
        let mut header = 0_i32.to_ne_bytes().to_vec();
        header.extend_from_slice(&u64::MAX.to_ne_bytes());
        assert!(decode(&header).is_err());
        header.truncate(4);
        header.extend_from_slice(&1_u64.to_ne_bytes());
        assert!(decode(&header).is_err());
        header.extend_from_slice(&(-0.916_536_277_098_507_5_f64).to_ne_bytes());
        assert_eq!(
            decode(&header).unwrap()[0].to_bits(),
            (-0.916_536_277_098_507_5_f64).to_bits()
        );
    }
}
