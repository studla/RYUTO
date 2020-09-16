/* 
 * File:   exhaustive_enum.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on June 8, 2016, 3:13 PM
 */

#include "exhaustive_enum.h"
#include <list>
#include <limits>

#include <lemon/lgf_writer.h>

exhaustive_enum::exhaustive_enum(ListDigraph& wc,
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
        unsigned int size)
: path_finder(wc, s, t, fc, ai, cni,kp, kbp, unsecurityArc, unsecurityId, input_ids, transcripts, size) {
}


exhaustive_enum::~exhaustive_enum() {
}

void exhaustive_enum::extract_transcripts(alternative_transcript_collection& results, int guiding) {

    ListDigraph::ArcIt has_arc(wc);
    if ( has_arc == INVALID ) {
        // no more arcs left in the graph, we are at the end!  
        return;
    }
    
    alternative_transcript_collection first;
    unsigned int min_count = std::numeric_limits<unsigned int>::max();
    
    bool found = true;
    while (found) {
        found = false;
        for (ListDigraph::OutArcIt a(wc, s); a!=INVALID; ++a) {
            if (wc.target(a) == t) {
               std::deque<ListDigraph::Arc> path;
               path.push_back(a);

               add_path_to_collection(path, first, wc, ai, fc, kp, unsecurityArc, unsecurityId, input_ids);

               gmap<int, capacity_type> capacities;
               for (gmap<int, transcript::series_struct>::iterator iss = first.transcripts.back()->series.begin(); iss != first.transcripts.back()->series.end(); ++iss) { // contains all inputs!
                   capacities[iss->first] = iss->second.flow;
               }

               // now we need to remove this path from the copy
               flow_series flows = remove_full_path(path, capacities, guiding, wc, ai, fc, kp, kbp, unsecurityArc, unsecurityId, input_ids, transcripts[guiding]);

               for (std::set<int>::iterator iii =  input_ids.begin(); iii != input_ids.end(); ++iii) {    
                   first.transcripts.back()->series[*iii].mean = flows.get_mean(*iii).mean;
                   first.transcripts.back()->series[*iii].effective_length = flows.get_mean(*iii).weight;
                   first.transcripts.back()->series[*iii].score = flows.get_mean(*iii).compute_score();
               }
               found = true;
            }
        }
    }
    
    full_enum_recursion(wc, s, t, fc, ai, cni, kp, kbp, unsecurityArc, unsecurityId, first, results, &min_count, guiding); 
}


void exhaustive_enum::full_enum_recursion(ListDigraph& wc,
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
        
        int guiding) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Enter recursion " + std::to_string(collection.size()) + " " + std::to_string(*min_count) + " " + std::to_string(guiding) +  "\n");
    #endif

//           digraphWriter(wc, std::cout)
//            .arcMap("ai", ai)
//            .arcMap("flow", fc)
//            .node("source", s)
//            .node("drain", t)
//            .run(); 
    
    if (collection.size() >= *min_count) { // branch and bound condition
        return;
    }
   
    ListDigraph::ArcIt has_arc(wc);
    if ( has_arc == INVALID ) {
        // no more acrs left in the graph, we are at the end!
        if (*min_count > collection.size()) {
            *min_count = collection.size();
            best = collection; // this only works thanks to deep copy
        }
        
        return;
    }
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Advance to end.\n");
    #endif

    std::deque<ListDigraph::OutArcIt> stack; // we keep a stack for path enumeration
    advance_to_end(stack, wc, kp, s, t); // we get the first arc
    
    bool more_path = true;
    while (more_path) {
             
        // we can iterator from here, first we copy the copy...
        ListDigraph wcc; // working copy for modification
        ListDigraph::ArcMap<flow_series> fcc(wcc);
        ListDigraph::ArcMap<arc_identifier> cai(wcc);
        ListDigraph::NodeMap<unsigned int> ccni(wcc);
        ListDigraph::ArcMap<arc_bridge> ckp(wcc);
        ListDigraph::ArcMap<arc_back_bridge> ckbp(wcc);
       
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > cua(wcc);
        ListDigraph::NodeMap< unsecurity_id > cui(wcc);
        
        ListDigraph::Node ct, cs;
        
        ListDigraph::ArcMap<ListDigraph::Arc> arc_ref(wcc);
        ListDigraph::ArcMap<ListDigraph::Arc> arc_ref_rev(wc);

        DigraphCopy<ListDigraph, ListDigraph> copy(wc,wcc);
        copy.arcMap(fc, fcc).arcMap(ai, cai)
            .arcMap(unsecurityArc,cua).nodeMap(unsecurityId, cui)
            .arcMap(kp, ckp).arcMap(kbp, ckbp)
            .nodeMap(cni, ccni).node(t,ct).node(s,cs)
            .arcCrossRef(arc_ref);
       
        std::deque<ListDigraph::Arc> path;
        path.resize(stack.size());
        
        int i = 0;
        for ( std::deque<ListDigraph::OutArcIt>::iterator a = stack.begin(); a!= stack.end(); ++a, ++i) {
            copy.arc(*a, path[i]);
        }
        copy.run();

        
        reverse_arc_cross_map(wc, wcc, arc_ref, arc_ref_rev);
        // unfortunately we now also need to update all known and back_known, since IDs are not consistent
        update_known(ckp, wc, wcc, arc_ref_rev);
        update_known_back(ckbp, wc, wcc, arc_ref_rev);
        
        // we need to add new path to a collection copy
        alternative_transcript_collection cc = collection;
        
//        logger::Instance()->debug("Copy.\n");
//             digraphWriter(wcc, std::cout)
//            .arcMap("ai", cai)
//            .arcMap("flow", fcc)
//            .node("source", cs)
//            .node("drain", ct)
//            .run();
        
        add_path_to_collection(path, cc, wcc, cai, fcc, ckp, cua, cui, input_ids);       
        gmap<int, capacity_type> capacities;
        for (gmap<int, transcript::series_struct>::iterator iss = cc.transcripts.back()->series.begin(); iss != cc.transcripts.back()->series.end(); ++iss) { // contains all inputs!
            capacities[iss->first] = iss->second.flow;
        }
        
        // now we need to remove this path from the copy
        flow_series flows = remove_full_path(path, capacities, guiding, wcc, cai, fcc, ckp, ckbp, cua, cui, input_ids, cc);

        for (std::set<int>::iterator iii =  input_ids.begin(); iii != input_ids.end(); ++iii) {    
            cc.transcripts.back()->series[*iii].mean = flows.get_mean(*iii).mean;
            cc.transcripts.back()->series[*iii].score = flows.get_mean(*iii).compute_score();
        }
        
        // start the recursion
        full_enum_recursion(wcc, cs, ct, fcc, cai, cni, ckp, ckbp, cua, cui, cc, best, min_count, guiding);
        
        if (advance_stack(stack)) {
            advance_to_end(stack, wc, kp, s, t);
        } else {
            more_path = false;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("exit\n");
        #endif
        
    }

}


void exhaustive_enum::advance_to_end(std::deque<ListDigraph::OutArcIt> &stack,
        ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::Node &s, ListDigraph::Node &t) {
    
    if (!stack.empty() && wc.target(stack.back()) == t) { // we hit drain, we are done here
        return;
    }
    
    ListDigraph::Node next_node;
    if (stack.empty()) {
        next_node = s;
    } else {
        next_node = wc.target(stack.back());
    }
    
    for (ListDigraph::OutArcIt a(wc, next_node); a!=INVALID; ++a) {
        
        if ( stack.empty() || !know_paths[stack.back()].is_evidenced_path() || know_paths[stack.back()].has_path_evidence(wc.id(a)) )  {
            stack.push_back(a);
            advance_to_end(stack, wc, know_paths, s, t);
            return;
        }
    }
    
    // this should never happen for correct flow graphs by decomposition lemma
}

bool exhaustive_enum::advance_stack(std::deque<ListDigraph::OutArcIt> &stack) {
    
    ++stack.back(); // we advance the current iterator
    if (stack.back() == INVALID) {
        // this means we track back one step
        stack.pop_back();
        if (!stack.empty()) {
            
            return advance_stack(stack);
        } else {
            return false;
        }
    }
    return true;
}
