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
        unsigned int size) 
: path_finder(wc, s, t, cfc,mc, ces, cet, cel, cni,kp, kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size) {
}


exhaustive_enum::~exhaustive_enum() {
}

void exhaustive_enum::extract_transcripts(alternative_transcript_collection& results) {

    ListDigraph::ArcIt has_arc(wc);
    if ( has_arc == INVALID ) {
        // no more arcs left in the graph, we are at the end!  
        return;
    }
    
    alternative_transcript_collection first;
    unsigned int min_count = std::numeric_limits<unsigned int>::max();
    
    full_enum_recursion(wc, s, t, cfc, ces, cet, cel, cni, kp, kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, first, results, &min_count); 
}


void exhaustive_enum::full_enum_recursion(ListDigraph &wc, ListDigraph::Node &s, ListDigraph::Node &t,
        ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int> &cni,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        ListDigraph::ArcMap<unsigned int>& cycle_id_in, ListDigraph::ArcMap<unsigned int>& cycle_id_out,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
        alternative_transcript_collection &collection,
        alternative_transcript_collection& best, unsigned int* min_count) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Enter recursion " + std::to_string(collection.size()) + " " + std::to_string(*min_count) + "\n");
    #endif

//           digraphWriter(wc, std::cout)
//            .arcMap("edge_specifier", ces)
//            .arcMap("edge_type", cet)
//            .arcMap("flow", fc)
//            .arcMap("length", cel)
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
    advance_to_end(stack, wc, know_paths, s, t); // we get the first arc
    
    bool more_path = true;
    while (more_path) {
             
        // we can iterator from here, first we copy the copy...
        ListDigraph wcc; // working copy for modification
        ListDigraph::ArcMap<capacity_type> fcc(wcc);
        ListDigraph::ArcMap<exon_edge> cces(wcc); 
        ListDigraph::ArcMap<edge_types::edge_type> ccet(wcc);
        ListDigraph::ArcMap<edge_length> ccel(wcc);
        ListDigraph::NodeMap<unsigned int> ccni(wcc);
        ListDigraph::ArcMap<arc_bridge> ckp(wcc);
        ListDigraph::ArcMap<arc_back_bridge> ckbp(wcc);
        ListDigraph::ArcMap<unsigned int> ccidin(wcc);
        ListDigraph::ArcMap<unsigned int> ccidout(wcc);
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > cua(wcc);
        ListDigraph::NodeMap< unsecurity_id > cui(wcc);
        
        ListDigraph::Node ct, cs;
        
        ListDigraph::ArcMap<ListDigraph::Arc> arc_ref(wcc);
        ListDigraph::ArcMap<ListDigraph::Arc> arc_ref_rev(wc);

        DigraphCopy<ListDigraph, ListDigraph> copy(wc,wcc);
        copy.arcMap(fc, fcc).arcMap(cel, ccel)
            .arcMap(ces, cces).arcMap(cet, ccet)
            .arcMap(cycle_id_in,ccidin).arcMap(cycle_id_out,ccidout)
            .arcMap(unsecurityArc,cua).nodeMap(unsecurityId, cui)
            .arcMap(know_paths, ckp).arcMap(know_back_paths, ckbp)
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
      //MC_FIX  capacity_type cap = add_path_to_collection(path, cc, wcc, fcc, cces, ccet, ckp, ccidin, ccidout, cua, cui);
               
        // now we need to remove this path from the copy
     //MC_FIX   remove_full_path(path, cap, wcc, fcc, cces, ccet, ccel, ckp, ckbp, ccidin, ccidout, cua, cui); // problem copying in KNOWN path and backpaths

        // start the recursion
        full_enum_recursion(wcc, cs, ct, fcc, cces, ccet, ccel, ccni, ckp, ckbp, ccidin, ccidout, cua, cui, cc, best, min_count);
        
        if (advance_stack(stack)) {
            advance_to_end(stack, wc, know_paths, s, t);
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
