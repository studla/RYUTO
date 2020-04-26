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
    
    virtual ~exhaustive_enum();
    
    virtual void extract_transcripts( alternative_transcript_collection &results, int guiding);

private:
    
    void full_enum_recursion(ListDigraph& wc,
        ListDigraph::Node& s,
        ListDigraph::Node& t,
        ListDigraph::ArcMap<flow_series>& fc,   
        ListDigraph::ArcMap<arc_identifier> &ai,
        ListDigraph::NodeMap<unsigned int>& cni,
        ListDigraph::ArcMap<arc_bridge>& kp,
        ListDigraph::ArcMap<arc_back_bridge>& kbp,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id> &unsecurityId,
        
        alternative_transcript_collection &collection,
        alternative_transcript_collection& best, unsigned int* min_count,
        
        int guiding);

    bool advance_stack(std::deque<ListDigraph::OutArcIt> &stack);
    
    void advance_to_end(std::deque<ListDigraph::OutArcIt> &stack, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::Node &s, ListDigraph::Node &t);
    
};

#endif	/* EXHAUSTIVE_ENUM_H */

