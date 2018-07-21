/* 
 * File:   exhaustive_enum.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on June 8, 2016, 3:13 PM
 */

#ifndef EXHAUSTIVE_ENUM_H
#define	EXHAUSTIVE_ENUM_H

#include "path_finder.h"

class exhaustive_enum : public path_finder{
public:
    exhaustive_enum(ListDigraph& wc,
        ListDigraph::Node& s,
        ListDigraph::Node& t,
        ListDigraph::ArcMap<capacity_type>& cfc,
        ListDigraph::ArcMap<capacity_mean> &mc,     
        ListDigraph::ArcMap<exon_edge>& ces,
        ListDigraph::ArcMap<edge_types::edge_type>& cet,
        ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int>& cni,
        ListDigraph::ArcMap<arc_bridge>& kp,
        ListDigraph::ArcMap<arc_back_bridge>& kbp,
        ListDigraph::ArcMap<unsigned int>& cycle_id_in,
        ListDigraph::ArcMap<unsigned int>& cycle_id_out,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
        unsigned int size);
    
    virtual ~exhaustive_enum();
    
    virtual void extract_transcripts( alternative_transcript_collection &results);
    
private:
    
    void full_enum_recursion(ListDigraph &wc, ListDigraph::Node &s, ListDigraph::Node &t,
        ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int> &cni,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        ListDigraph::ArcMap<unsigned int>& cycle_id_in, ListDigraph::ArcMap<unsigned int>& cycle_id_out,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
        alternative_transcript_collection &collection,
        alternative_transcript_collection& best, unsigned int* min_count);

    bool advance_stack(std::deque<ListDigraph::OutArcIt> &stack);
    
    void advance_to_end(std::deque<ListDigraph::OutArcIt> &stack, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::Node &s, ListDigraph::Node &t);
    
};

#endif	/* EXHAUSTIVE_ENUM_H */

