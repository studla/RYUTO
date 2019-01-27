/* 
 * File:   path_finder.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on June 8, 2016, 2:22 PM
 */

#ifndef PATH_FINDER_H
#define	PATH_FINDER_H

#include <lemon/list_graph.h>
#include <lemon/dfs.h>
#include <vector>
#include "../../output/alternative_transcript_collection.h"
#include "../../flow_graph/exon_edge.h"
#include "../../flow_graph/edge_types.h"
#include "../../flow_graph/edge_length.h"
#include "../../flow_graph/coverage/capacity_mean.h"
#include "../../flow_graph/arc_identifier.h"
#include "../../flow_graph/arc_bridge.h"
#include "../../flow_graph/arc_back_bridge.h"
#include "../../../Datatype_Templates/topsort.h"
#include "../../pre_graph/exon_group.h"
#include "../../output/unsecurity_id.h"

using namespace lemon;

class path_finder {
public:
    path_finder(ListDigraph& wc,
        ListDigraph::Node& s,
        ListDigraph::Node& t,
        ListDigraph::ArcMap<flow_series>& fc,   
        ListDigraph::ArcMap<arc_identifier> &ai,
        ListDigraph::NodeMap<unsigned int>& cni,
        ListDigraph::ArcMap<arc_bridge>& kp,
        ListDigraph::ArcMap<arc_back_bridge>& kbp,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id> &unsecurityId,
        std::set<int>& input_ids,
        std::map<int, alternative_transcript_collection>& transcripts,    
        unsigned int size);
    virtual ~path_finder();
    
    virtual void extract_transcripts( alternative_transcript_collection &results, int guiding) = 0; 
    
    void extract_guided_transcripts(alternative_transcript_collection &results, graph_list<exon_group *> guided, int guiding_id);
    void extract_guided_transcripts_linear_order(alternative_transcript_collection& results, graph_list<exon_group*> guided, int guiding_id);
    
    static path_finder *create_path_finder(ListDigraph& wc,
        ListDigraph::Node& s,
        ListDigraph::Node& t,
        ListDigraph::ArcMap<flow_series>& fc,   
        ListDigraph::ArcMap<arc_identifier> &ai,
        ListDigraph::NodeMap<unsigned int>& cni,
        ListDigraph::ArcMap<arc_bridge>& kp,
        ListDigraph::ArcMap<arc_back_bridge>& kbp,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id> &unsecurityId,
        std::set<int>& input_ids,
        std::map<int, alternative_transcript_collection>& transcripts,
        unsigned int size, exon_meta* meta);
     
protected:
    
    void add_path_to_collection(std::deque<ListDigraph::Arc> &stack, alternative_transcript_collection &cc,
        ListDigraph &wc, ListDigraph::ArcMap<arc_identifier> &ai, ListDigraph::ArcMap<flow_series> &fs,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,  std::set<int> &input_ids);
    flow_series remove_full_path(std::deque<ListDigraph::Arc> &stack, 
        gmap<int, capacity_type> &capacities, int guiding, ListDigraph &wc,
        ListDigraph::ArcMap<arc_identifier> &ai,ListDigraph::ArcMap<flow_series> &fs,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id> &unsecurityId, std::set<int> &input_ids, alternative_transcript_collection &trans);
    
    void reverse_arc_cross_map(ListDigraph& g, ListDigraph& wc, ListDigraph::ArcMap<ListDigraph::Arc> & arc_ref, ListDigraph::ArcMap<ListDigraph::Arc> & arc_ref_rev);
    
    void update_known(ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph &g , ListDigraph &wc, ListDigraph::ArcMap<ListDigraph::Arc> &ref);
    void update_known_back(ListDigraph::ArcMap<arc_back_bridge>& know_back_paths, ListDigraph &g, ListDigraph& wc, ListDigraph::ArcMap<ListDigraph::Arc>& ref);
    
    bool recursive_arc_backtracing(exon_edge& goal, unsigned int next_start, unsigned int goal_end_index, ListDigraph::Node n, ListDigraph::Arc la, std::deque<ListDigraph::Arc> &path, bool exon) ;
    
    ListDigraph& wc;
    ListDigraph::Node& s;
    ListDigraph::Node& t;
    ListDigraph::ArcMap<flow_series>& fc;
    ListDigraph::ArcMap<arc_identifier> &ai; 
    ListDigraph::NodeMap<unsigned int>& cni;
    ListDigraph::ArcMap<arc_bridge>& kp;
    ListDigraph::ArcMap<arc_back_bridge>& kbp;
    ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > >& unsecurityArc;
    ListDigraph::NodeMap< unsecurity_id >& unsecurityId;
    
    std::set<int>& input_ids;
    std::map<int, alternative_transcript_collection>& transcripts;
    
    unsigned int size;
};

#endif	/* PATH_FINDER_H */

