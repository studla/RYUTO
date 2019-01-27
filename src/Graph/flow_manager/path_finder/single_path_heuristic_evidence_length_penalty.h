/* 
 * File:   single_path_heuristic.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on June 8, 2016, 3:03 PM
 */

#ifndef SINGLE_PATH_HEURISTIC_EVIDENCE_LENGTH_PENALTY_H
#define	SINGLE_PATH_HEURISTIC_EVIDENCE_LENGTH_PENALTY_H

#include "path_finder.h"

class single_path_heuristic_evidence_length_penalty : public path_finder {
public:
    single_path_heuristic_evidence_length_penalty(ListDigraph& wc,
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

    virtual ~single_path_heuristic_evidence_length_penalty();
    
    virtual void extract_transcripts( alternative_transcript_collection &results, int guiding);
    
private:

    void find_min_max(std::deque<ListDigraph::Arc> &path,  int guiding);
};

#endif	/* SINGLE_PATH_HEURISTIC_EVIDENCE_PENALTY_H */

