bin_PROGRAMS = ryuto

ryuto_SOURCES=Logger/logger.cc Controller/bam_controller.cc Controller/base_controller.cc Reader/Chromosome/interval.cc Reader/Chromosome/read_collection.cc Reader/Chromosome/chromosome.cc Reader/Chromosome/exon.cc Reader/Chromosome/connection_iterator.cc Reader/Chromosome/read.cc Reader/Chromosome/raw_series_counts.cc Reader/Chromosome/raw_atom.cc Reader/Chromosome/connected.cc Reader/gffReader.cc Reader/reader_base.cc Reader/bam_guided_reader.cc Reader/bam_reader.cc Graph/flow_manager/mincost_flow_square_hcost.cc Graph/flow_manager/path_finder/path_finder.cc Graph/flow_manager/path_finder/single_path_heuristic_evidence_length_penalty.cc Graph/flow_manager/path_finder/single_path_heuristic.cc Graph/flow_manager/mincost_flow_base.cc Graph/flow_manager/base_manager_buildgraph.cc Graph/flow_manager/base_manager_flow_handler.cc Graph/flow_manager/mincost_flow_square.cc Graph/output/transcript.cc Graph/output/alternative_transcript_collection.cc Graph/output/transcript_unsecurity.cc Graph/output/unsecurity_id.cc Graph/overlap_graph/range_helper.cc Graph/overlap_graph/contained_node.cc Graph/overlap_graph/overlap_node.cc Graph/overlap_graph/arc_range.cc Graph/flow_graph/edge_length.cc Graph/flow_graph/build/region.cc Graph/flow_graph/build/count_raw_node.cc Graph/flow_graph/build/count_raw_edge.cc Graph/flow_graph/arc_bridge.cc Graph/flow_graph/arc_back_bridge.cc Graph/flow_graph/coverage/capacity_mean.cc Graph/flow_graph/coverage/flow_series.cc Graph/flow_graph/arc_identifier.cc Graph/flow_graph/path_evidence.cc Graph/flow_graph/exon_edge.cc Graph/flow_graph/path.cc Graph/pre_graph/exon_meta.cc Graph/pre_graph/exon_group.cc Graph/pre_graph/pre_graph.cc Graph/pre_graph/exon_group_count.cc Options/options.cc main.cc 

# samtools includes
AM_CPPFLAGS = $(HTS_CPPFLAGS)
AM_LDFLAGS = $(HTS_LDFLAGS)
ryuto_LDADD = $(HTS_LIBS)

# boost includes
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_LDFLAGS += $(BOOST_LDFLAGS)
ryuto_LDADD += $(BOOST_PROGRAM_OPTIONS_LIB)

# lemon includes
AM_CPPFLAGS += $(LEMON_CPPFLAGS)
AM_LDFLAGS += $(LEMON_LDFLAGS)
ryuto_LDADD += $(LEMON_LIBS)

#OpenMP
AM_CPPFLAGS += $(OPENMP_CXXFLAGS)
AM_LDFLAGS += $(OPENMP_CXXFLAGS)

#CLP
AM_CPPFLAGS += $(CLP_CPPFLAGS)
AM_LDFLAGS += $(CLP_LDFLAGS)
ryuto_LDADD += $(CLP_LIBS)

#Debug Flag
#CPPFLAGS = -DALLOW_DEBUG
