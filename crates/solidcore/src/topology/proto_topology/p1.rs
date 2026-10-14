//! Proto topology, part 1: owning-topology lists, assemblies, the copy-B-rep map, the
//! proto-topology containers with their accessors, list management and evaluation,
//! the static helpers of `SmProtoTopology.cpp`, and the API transformations.
//!
//! Upstream: `source/SMLib/inc/SmOwningTopology.h`, `source/SMLib/src/SmOwningTopology.cpp`,
//! `source/SMLib/inc/SmAssembly.h`, `source/SMLib/src/SmAssembly.cpp`,
//! `source/SMLib/inc/SmAssemblyInstance.h`, `source/SMLib/src/SmAssemblyInstance.cpp`,
//! `source/SMLib/inc/SmCopyBrepMap.h`, `source/SMLib/src/SmCopyBrepMap.cpp`,
//! `source/SMLib/inc/SmProtoTopology.h`, `source/SMLib/src/SmProtoTopology.cpp`,
//! `source/SM_API/src/SmApiGeneral.cpp`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`OwningTopology`], [`OwningTopology::new`], [`OwningTopology::destroy`] | `SmOwningTopology`, its constructor and destructor |
//! | [`OwningTopology::size`], [`OwningTopology::list`], [`OwningTopology::set_list`] | `GetSize`, `GetList`, `SetList` |
//! | [`OwningTopology::update_list_size`] | `SmOwningTopology::UpdateListSize` |
//! | [`OwningTopology::insert_after`] | `SmOwningTopology::InsertAfter` |
//! | [`OwningTopology::pre_insert`], [`OwningTopology::post_insert`] | `PreInsert`, `PostInsert` |
//! | [`OwningTopology::remove`] | `SmOwningTopology::Remove` |
//! | [`OwningTopology::is_in_list`] | `SmOwningTopology::IsInList` |
//! | [`OwningTopology::get_all`], [`OwningTopology::get_all_ignoring_status`] | `SmOwningTopology::GetAll` |
//! | [`OwningTopology::reverse_list`] | `SmOwningTopology::ReverseList` |
//! | [`OwningRecords`], [`destroy_topology`] | the `SmOwningTopology*` list owner, `SmTopology::~SmTopology` |
//! | [`owning_topology::TYPE`], [`owning_topology::TYPE_STRING`], [`owning_topology::CLASS_STRING`], [`owning_topology::is_kind_of`] | `SmOwningTopology::GetType`, `GetTypeString`, `GetClassString`, `IsKindOf` |
//! | [`Assembly`], [`Assembly::new`], [`Assembly::destroy`] | `SmAssembly`, its constructor, `~SmAssembly` |
//! | [`Assembly::assembly_instances`] | `SmAssembly::GetAssemblyInstances` |
//! | [`Assembly::add_instance`], [`Assembly::remove_instance`] | `SmAssembly::AddInstance`, `RemoveInstance` |
//! | [`assembly::ASSEMBLY_TYPE`], [`assembly::ASSEMBLY_TYPE_STRING`], [`assembly::ASSEMBLY_CLASS_STRING`], [`assembly::assembly_is_kind_of`] | `SmAssembly::GetType`, `GetTypeString`, `GetClassString`, `IsKindOf` |
//! | [`AssemblyInstance`], [`AssemblyInstance::new`], [`AssemblyInstance::destroy`] | `SmAssemblyInstance`, its constructor, `~SmAssemblyInstance` |
//! | [`AssemblyInstance::parent_assembly`] | `SmAssemblyInstance::GetParentAssembly` |
//! | [`assembly::ASSEMBLY_INSTANCE_TYPE`], [`assembly::ASSEMBLY_INSTANCE_TYPE_STRING`], [`assembly::ASSEMBLY_INSTANCE_CLASS_STRING`], [`assembly::assembly_instance_is_kind_of`] | `SmAssemblyInstance::GetType`, `GetTypeString`, `GetClassString`, `IsKindOf` |
//! | [`CopyBrepMap`], [`CopyBrepMap::new`], [`CopyBrepMap::destroy`] | `SmCopyBrepMap`, its constructor and destructor |
//! | [`CopyBrepMap::set_at`], [`CopyBrepMap::get_at`] | `SmCopyBrepMap::SetAt`, `GetAt` |
//! | [`proto::can_imprint`] | `sm_CanImprint` (static) |
//! | [`proto::topo_dimension`] | `sm_TopoDimension` (static) |
//! | [`proto::find_position_from_pcs`], [`proto::PositionFromPcs`] | `sm_FindPositionFromPCs` (static) |
//! | [`proto::are_same_curve`] | `sm_AreSameCurve` (static) |
//! | [`proto::are_connected`] | the four `sm_AreConnected` overloads (static) |
//! | [`ProtoTopologyManager::which_pv_to_combine`] | `sm_WhichPVToCombine` (static) |
//! | [`proto::is_interior`] | `sm_IsInterior` (static) |
//! | [`proto::are_overlapping`] | `sm_AreOverlapping` (static) |
//! | [`proto::find_ed_on_degen_face`] | `sm_FindEDOnDegenFace` |
//! | [`VertexDefinition`], [`VertexDefinition::new`], [`ProtoTopologyManager::delete_vertex_definition`] | `SmVertexDefinition`, its constructor and destructor |
//! | [`VertexDefinition::position`], [`VertexDefinition::tolerance`], [`VertexDefinition::gap`], [`VertexDefinition::final_vertex`], [`VertexDefinition::proto_vertex`] | `GetPosition`, `GetTolerance`, `GetGap`, `GetFinalVertex`, `GetProtoVertex` |
//! | [`VertexDefinition::brep1_topo`], [`VertexDefinition::brep2_topo`], [`VertexDefinition::brep1_topo_type`], [`VertexDefinition::brep2_topo_type`] | `GetBrep1Topo`, `GetBrep2Topo`, `GetBrep1TopoType`, `GetBrep2TopoType` |
//! | [`VertexDefinition::src_topo`], [`VertexDefinition::src_topo_type`], [`VertexDefinition::brep_t_param`], [`VertexDefinition::brep_uv_param`], [`VertexDefinition::brep_point_classif`] | `GetSrcTopo`, `GetSrcTopo_TYPE`, `GetBrepTParam`, `GetBrepUVParam`, `GetBrepPointClassif` |
//! | [`VertexDefinition::set_position`], [`VertexDefinition::set_proto_vertex`], [`VertexDefinition::set_final_vertex`] | `SetPosition`, `SetProtoVertex`, `SetFinalVertex` |
//! | [`VertexDefinition::set_src_topo`], [`VertexDefinition::set_point_parameters`] | `SmVertexDefinition::SetSrcTopo`, `SetPointParameters` |
//! | [`ProtoTopologyManager::vertex_def_num_proto_edges`] | `SmVertexDefinition::GetNumProtoEdges` |
//! | [`ProtoTopologyManager::resolve_vertex_definition`] | `SmVertexDefinition::Resolve` |
//! | [`EdgeDefinition`], [`EdgeDefinition::new`], [`EdgeDefinition::copy_definition`], [`EdgeDefinition::destroy`] | `SmEdgeDefinition`, its constructors, `~SmEdgeDefinition` |
//! | [`EdgeDefinition::destruct`] | `SmEdgeDefinition::Destruct` |
//! | [`EdgeDefinition::xsect_curve_3d`], [`EdgeDefinition::xsect_interval`], [`EdgeDefinition::tolerance`], [`EdgeDefinition::xsect_uv_curve`] | `GetXSectCurve3d`, `GetXSectInterval`, `GetTolerance`, `GetXSectUVCurve` |
//! | [`EdgeDefinition::src_topo`], [`EdgeDefinition::src_topo_type`], [`EdgeDefinition::proto_edge`], [`EdgeDefinition::final_edge`] | `GetSrcTopo`, `GetSrcTopo_TYPE`, `GetProtoEdge`, `GetFinalEdge` |
//! | [`EdgeDefinition::set_xsect_curve_3d`], [`EdgeDefinition::set_xsect_curve_uv`], [`EdgeDefinition::set_proto_edge`], [`EdgeDefinition::set_src_topo`] | `SetXSectCurve3d`, `SetXSectCurveUV`, `SetProtoEdge`, `SetSrcTopo` |
//! | [`EdgeDefinition::set_xsect_interval_min_max`] | `SmEdgeDefinition::SetXSectInterval(double, double)` |
//! | [`ProtoTopologyManager::edge_def_start_proto_vertex`], [`ProtoTopologyManager::edge_def_end_proto_vertex`] | `SmEdgeDefinition::GetStartProtoVertex`, `GetEndProtoVertex` |
//! | [`ProtoTopologyManager::edge_def_end_position_at`], [`EdgeDefinition::end_position`] | the two `SmEdgeDefinition::GetEndPosition` |
//! | [`EdgeDefinition::evaluate`], [`EdgeDefinition::evaluate_mid`] | `SmEdgeDefinition::Evaluate`, `EvaluateMid` |
//! | [`EdgeDefinition::reverse_orientation`], [`EdgeDefinition::reparameterize`] | `SmEdgeDefinition::ReverseOrientation`, `Reparameterize` |
//! | [`ProtoTopologyManager::is_same_edge_definition`] | `SmEdgeDefinition::IsSameEdgeDefinition` |
//! | [`EdgeDefinition::is_interior_point`] | `SmEdgeDefinition::IsInteriorPoint` |
//! | [`ProtoVertex`], [`ProtoTopologyManager::new_proto_vertex`], [`ProtoTopologyManager::new_proto_vertex_with_pcs`], [`ProtoTopologyManager::delete_proto_vertex_object`] | `SmProtoVertex`, two of its constructors, `~SmProtoVertex` |
//! | [`ProtoVertex::position`], [`ProtoVertex::tolerance`], [`ProtoVertex::num_vertex_definitions`], [`ProtoVertex::num_proto_edges`], [`ProtoVertex::vertex_definitions`], [`ProtoVertex::proto_edges`] | `GetPosition`, `GetTolerance`, `GetNumVertexDefinitions`, `GetNumProtoEdges`, `GetVertexDefinitions`, `GetProtoEdges` |
//! | [`ProtoVertex::clear_vertex_definitions`], [`ProtoVertex::add_proto_edge`] | `ClearVertexDefinitions`, `AddProtoEdge` |
//! | [`ProtoVertex::remove_vertex_definition`], [`ProtoVertex::remove_proto_edge`] | `SmProtoVertex::RemoveVertexDefinition`, `RemoveProtoEdge` |
//! | [`ProtoTopologyManager::add_vertex_definition`] | `SmProtoVertex::AddVertexDefinition(SmVertexDefinition*)` |
//! | [`ProtoTopologyManager::find_matching_vertex_definition`] | `SmProtoVertex::FindMatchingVertexDefinition` |
//! | [`ProtoTopologyManager::closest_vertex_definition`] | `SmProtoVertex::GetClosestVertexDefinition` |
//! | [`ProtoTopologyManager::proto_vertex_final_vertex`] | `SmProtoVertex::GetFinalVertex` |
//! | [`ProtoTopologyManager::combine_proto_vertex`] | `SmProtoVertex::CombineProtoVertex` |
//! | [`ProtoEdge`], [`ProtoTopologyManager::new_proto_edge`], [`ProtoTopologyManager::copy_proto_edge`], [`ProtoTopologyManager::delete_proto_edge_object`] | `SmProtoEdge`, its constructors, `~SmProtoEdge` |
//! | [`ProtoEdge::xsect_interval`], [`ProtoEdge::tolerance`], [`ProtoEdge::start_proto_vertex`], [`ProtoEdge::end_proto_vertex`], [`ProtoEdge::overlapping_pes`] | `GetXSectInterval`, `GetTolerance`, `GetStartProtoVertex`, `GetEndProtoVertex`, both `GetOverlappingPEs` |
//! | [`ProtoEdge::num_edge_definitions`], [`ProtoEdge::edge_definitions`] | `GetNumEdgeDefinitions`, `GetEdgeDefinitions` |
//! | [`ProtoEdge::set_start_proto_vertex`], [`ProtoEdge::set_end_proto_vertex`], [`ProtoEdge::add_overlapping_pe`], [`ProtoEdge::append_overlapping_pes`], [`ProtoEdge::clear_edge_definitions`] | `SetStartProtoVertex`, `SetEndProtoVertex`, `AddOverlappingPE`, `AppendOverlappingPEs`, `ClearEdgeDefinitions` |
//! | [`ProtoEdge::replace_proto_vertex`], [`ProtoEdge::remove_edge_definition`], [`ProtoEdge::remove_proto_vertex`] | `SmProtoEdge::ReplaceProtoVertex`, `RemoveEdgeDefinition`, `RemoveProtoVertex` |
//! | [`ProtoTopologyManager::overlapping_pe_group`] | `SmProtoEdge::GetOverlappingPEGroup` |
//! | [`ProtoTopologyManager::proto_edge_end_position_at`], [`ProtoTopologyManager::proto_edge_end_position`] | the two `SmProtoEdge::GetEndPosition` |
//! | [`ProtoTopologyManager::relate_coincident_pf`] | `SmProtoEdge::RelateCoincidentPF` |
//! | [`ProtoTopologyManager::proto_edge_evaluate`], [`ProtoTopologyManager::proto_edge_evaluate_mid`] | `SmProtoEdge::Evaluate`, `EvaluateMid` |
//! | [`ProtoTopologyManager::is_same_proto_edge`] | `SmProtoEdge::IsSameProtoEdge` |
//! | [`ProtoFace`], [`ProtoTopologyManager::new_proto_face`] | `SmProtoFace`, its constructor and destructor |
//! | [`ProtoFace::src_topo`], [`ProtoFace::coinc_proto_edges`], [`ProtoFace::set_src_topo`] | `SmProtoFace::GetSrcTopo`, both `GetCoincProtoEdges`, `SetSrcTopo` |
//! | [`ProtoTopologyManager`], [`ProtoTopologyManager::new`], [`ProtoTopologyManager::destroy`], [`ProtoTopologyManager::reset`] | `SmProtoTopologyManager`, its constructor, destructor, `Reset` |
//! | [`ProtoTopologyManager::which_brep`], [`ProtoTopologyManager::tolerance`], [`ProtoTopologyManager::proto_vertices`], [`ProtoTopologyManager::proto_edges`], [`ProtoTopologyManager::proto_faces`] | `GetWhichBrep`, `GetTolerance`, `GetProtoVertices`, `GetProtoEdges`, `GetProtoFaces` |
//! | [`ProtoTopologyManager::add_proto_vertex`], [`ProtoTopologyManager::add_proto_edge`], [`ProtoTopologyManager::add_proto_face`] | `AddProtoVertex`, `AddProtoEdge`, `AddProtoFace` |
//! | [`ProtoTopologyManager::find_proto_vertex`], [`ProtoTopologyManager::find_edge_definition_from_edge`] | `FindProtoVertex`, `FindEdgeDefinitionFromEdge` |
//! | [`ProtoTopologyManager::all_vertex_definitions`], [`ProtoTopologyManager::all_edge_definitions`] | `GetAllVertexDefinitions`, `GetAllEdgeDefinitions` |
//! | [`ProtoTopologyManager::remove_vertex_definition`], [`ProtoTopologyManager::remove_proto_vertex`], [`ProtoTopologyManager::remove_proto_edge`] | `SmProtoTopologyManager::RemoveVertexDefinition`, `RemoveProtoVertex`, `RemoveProtoEdge` |
//! | [`ProtoTopologyManager::check_overlapping_proto_edges`] | `SmProtoTopologyManager::CheckOverlappingProtoEdges` |
//! | [`ProtoTopologyManager::create_proto_face`] | `SmProtoTopologyManager::CreateProtoFace` |
//! | [`ProtoTopologyManager::can_combine_edge_definitions`] | `SmProtoTopologyManager::CanCombineEdgeDefinitions` |
//! | [`ProtoTopologyManager::remove_coincident_no_edge_proto_vertex`] | `SmProtoTopologyManager::RemoveCoincidentNoEdgeProtoVertex` |
//! | [`proto::MANAGER_TYPE`], [`proto::MANAGER_TYPE_STRING`], [`proto::MANAGER_CLASS_STRING`], [`proto::manager_is_kind_of`] | `SmProtoTopologyManager::GetType`, `GetTypeString`, `GetClassString`, `IsKindOf` |
//! | [`api_transform::transform_independent_object`], [`api_transform::TransformObject`] | `sm_TransformIndependentObject` (anonymous namespace) |
//! | [`api_transform::api_transform`], [`api_transform::api_scale`], [`api_transform::api_scale_by_pt`], [`api_transform::api_translate`], [`api_transform::api_rotate`] | `SmApiTransform`, `SmApiScale`, `SmApiScaleByPt`, `SmApiTranslate`, `SmApiRotate` |
//!
//! Deviations:
//! - **Lists.** An owner's `m_pList` / `m_lListSize` is an [`OwningTopology`]; the
//!   members' links stay in their [`Topology`](crate::topology::topology::Topology)
//!   records, edited through a [`TopologyRecords`](crate::topology::topology::TopologyRecords)
//!   store as the friend class edits them. Operations that dereference a null link or
//!   walk a list that never returns to its head (undefined or non-terminating
//!   upstream) return [`SolidError::InvalidInput`](crate::SolidError) instead; the
//!   walks keep a visited set for that. `GetAll`'s partial array on its `SM_ERR` path
//!   is [`OwningTopology::get_all_ignoring_status`]. `FixBackPointers`, `UpdateList`,
//!   `ConcatenateList`, `TypicalListOwnerUse`, `AssertValid` and `Dump` are not part of
//!   this unit (unused in the call graph).
//! - **Destructors.** `~SmTopology`'s removal from the owner's list is
//!   [`destroy_topology`] (the owner's list data comes from an [`OwningRecords`] store).
//!   The assembly destructors and constructors take the store and the object's key.
//!   Names are plain strings (the `SM_TBLOCK_SIZE` buffer is not modelled).
//! - **Proto objects** live in the [`ProtoTopologyManager`]'s store, addressed by
//!   typed ids; upstream allocates them with `new (pOwner)`. Members that reach the
//!   manager through a back pointer are manager methods taking the object's id. A
//!   deleted object's slot is empty; dereferencing a deleted or null object (undefined
//!   upstream) is [`SolidError::InvalidInput`](crate::SolidError), except in the
//!   destructors, where such references (a use after free upstream) are skipped.
//! - **Handles.** Topology objects, curves, B-reps and the topology intersector are
//!   caller-chosen handles `O` ([`ClassObject`](crate::intersect::point_class::ClassObject)),
//!   as in the point classification; the curves behind them come from a
//!   [`ClassGeometry`](crate::intersect::point_class::ClassGeometry), or a
//!   [`ProtoGeometry`] for the curve edits. Destructors that delete curves return the
//!   handles deleted. `SM_CAST_PTR(SmTopology, ..)` of a point object (undefined
//!   upstream) treats it as non-topology.
//! - **Constructor tolerance.** `SmProtoTopologyManager` only reads the topology
//!   intersector's `GetThisApproxTol3d`, passed as a value; the override tolerance is
//!   overwritten by it, as upstream.
//! - **Uninitialised members** of objects upstream allocates with `smos_Calloc` (the
//!   final edges of an edge definition) start as `None`.
//! - **Out-params** become return values: `IsInList`'s index, `FindProtoVertex`'s
//!   distance, `IsSameEdgeDefinition`'s identity flag, `IsInteriorPoint`'s parameter
//!   (meaningful when the result is true), `CanCombineEdgeDefinitions`' choice, the
//!   positions and tangents of the evaluations (`None` for what was not asked),
//!   `CheckOverlappingProtoEdges`' proto edges (upstream appends them).
//! - **Ignored statuses** of evaluations leave the point `SM_UNDEF` (upstream's
//!   `SetUninitialized`), and are otherwise dropped, except [`SolidError::Unported`](crate::SolidError).
//! - `sm_FindPositionFromPCs` does not compute `SmTol::GetXSectTol3d` of the two
//!   objects: it only feeds a diagnostic message.
//! - `CanCombineEdgeDefinitions` keeps upstream's comparison of the second curve over
//!   the first definition's interval. Its degenerate-face branch is unreachable (two
//!   different edges are never connected) but is ported.
//! - **Callees in units not ported yet** return [`SolidError::Unported`](crate::SolidError):
//!   the B-rep topology walks of `sm_AreConnected` (`topology.vertex`,
//!   `topology.edge`), `SmVertex::GetPoint`, `SmEdge::GetBrep`, `SmBrep::GetTolerance`
//!   (`topology.brep_p1`) in the manager's fallback tolerance,
//!   `SmProtoVertex::ContainsPoint` (`query.face_p2`) in `FindProtoVertex` and
//!   `sm_WhichPVToCombine`, `SmBrep::Transform` (`topology.brep_p5`) and
//!   `SmPolyBrep::Transform` (`query.poly`).
//! - `Dump`, `Draw` and the `SM_DEBUG_CODE` blocks are not ported.
//! - **Recursion.** None.

pub mod api_transform;
pub mod assembly;
pub mod copy_brep_map;
pub mod owning_topology;
pub mod proto;

pub use assembly::{Assembly, AssemblyInstance};
pub use copy_brep_map::CopyBrepMap;
pub use owning_topology::{destroy_topology, OwningRecords, OwningTopology};
pub use proto::{
    EdgeDefId, EdgeDefinition, ProtoEdge, ProtoEdgeId, ProtoFace, ProtoFaceId, ProtoGeometry,
    ProtoTopologyManager, ProtoVertex, ProtoVertexId, VertexDefId, VertexDefinition,
};
