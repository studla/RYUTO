/* 
 * File:   path_finder.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on June 8, 2016, 2:22 PM
 */

#include "path_finder.h"
#include "single_path_heuristic.h"
#include "exhaustive_enum.h"
#include "multi_path_heuristic.h"
#include "single_path_heuristic_evidence_penalty.h"
#include "single_path_heuristic_evidence_length_penalty.h"
#include "single_path_heuristic_evidence_offset.h"
#include "single_path_heuristic_kmax_length.h"
#include "single_path_heuristic_longest_max.h"

#include "../base_manager.h"
#include "single_path_heuristic_no_neighbour.h"
#include <unordered_map>

path_finder::path_finder(ListDigraph& wc,
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
        ListDigraph::ArcMap< lazy<std::set< transcript_unsecurity > > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
        unsigned int size) 
 : wc(wc), s(s), t(t), cfc(cfc), mc(mc), ces(ces), cet(cet), cel(cel), cni(cni), kp(kp), kbp(kbp), cycle_id_in(cycle_id_in), cycle_id_out(cycle_id_out), unsecurityArc(unsecurityArc), unsecurityId(unsecurityId), size(size)
{

}

path_finder::~path_finder() {
}

path_finder* path_finder::create_path_finder(ListDigraph& wc, ListDigraph::Node& s, ListDigraph::Node& t, ListDigraph::ArcMap<capacity_type>& cfc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge>& ces, ListDigraph::ArcMap<edge_types::edge_type>& cet, ListDigraph::ArcMap<edge_length>& cel, ListDigraph::NodeMap<unsigned int>& cni, ListDigraph::ArcMap<arc_bridge>& kp, ListDigraph::ArcMap<arc_back_bridge>& kbp, ListDigraph::ArcMap<unsigned int>& cycle_id_in, ListDigraph::ArcMap<unsigned int>& cycle_id_out, ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc, ListDigraph::NodeMap< unsecurity_id> &unsecurityId, unsigned int size, exon_meta* meta) {
    
//    unsigned int complexity = 1;
//    for (ListDigraph::NodeIt n(wc); n != INVALID; ++n) {
//        if (n == t) {
//            continue;
//        }
//                
//        unsigned int count = 0;
//        for (ListDigraph::OutArcIt o(wc, n); o != INVALID; ++o) {
//            ++count;
//        }
//        if (count != 0) complexity = complexity * count;
//        if (complexity > 11) {
//            break;
//        }
//    }
//    
//    if (complexity < 11) {
//        logger::Instance()->info("Exhaustive " + std::to_string(complexity)+ "\n"); 
//        return new exhaustive_enum(wc,s,t,cfc,ces,cet,cel,cni,kp,kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size);
//    }
    
    return new single_path_heuristic_evidence_length_penalty(wc,s,t,cfc, mc,ces,cet,cel,cni,kp,kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size);

    
    
    // for now a simple return of one version is good enough
   // return new single_path_heuristic(wc,s,t,cfc, mc,ces,cet,cel,cni,kp,kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size);
    // return new single_path_heuristic_evidence_penalty(wc,s,t,cfc,ces,cet,cel,cni,kp,kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size);
    // return new single_path_heuristic_evidence_offset(wc,s,t,cfc,ces,cet,cel,cni,kp,kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size);
   //   return new exhaustive_enum(wc,s,t,cfc,ces,cet,cel,cni,kp,kbp,cycle_id_in, cycle_id_out, size);
    // return new multi_path_heuristic(wc,s,t,cfc,ces,cet,cel,cni,kp,kbp,cycle_id_in, cycle_id_out, size);  
     
    //return new single_path_heuristic_kmax_length(wc,s,t,cfc,ces,cet,cel,cni,kp,kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size);
    //return new single_path_heuristic_no_neighbour(wc,s,t,cfc,ces,cet,cel,cni,kp,kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size, meta);
    //return new single_path_heuristic_longest_max(wc,s,t,cfc,ces,cet,cel,cni,kp,kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size);
}



capacity_type path_finder::add_path_to_collection(std::deque<ListDigraph::Arc> &stack, alternative_transcript_collection &cc,
        ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<exon_edge> &ces, ListDigraph::ArcMap<edge_types::edge_type> &cet,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<unsigned int>& cycle_id_in,
        ListDigraph::ArcMap<unsigned int>& cycle_id_out, ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId) {
     
    // we collect every info over the path
    cc.transcripts.push_back(lazy<transcript>());
    
    std::deque<ListDigraph::Arc>::iterator a = stack.begin();
    // there is always at least one element in the stack
    cc.transcripts.back()->found_edge = ces[*a];
    cc.transcripts.back()->flow = fc[*a];  
    cc.transcripts.back()->cycle_id_in = cycle_id_in[*a];
 
    std::deque<ListDigraph::Arc>::iterator first = a;
    ++a;
    
    if (cet[*first] == edge_types::HELPER && a != stack.end()) {
        cc.transcripts.back()->found_edge = ces[*a];
        
        if (cet[*a] == edge_types::HELPER && ces[*a].node_index == -1) {
            cc.transcripts.back()->found_edge = ces[*first];
        }
        
    } else {
         for (std::set<transcript_unsecurity>::iterator unsecs = unsecurityArc[*first]->begin(); unsecs != unsecurityArc[*first]->end(); ++unsecs) {
            cc.transcripts.back()->unsecurity_id.push_back(*unsecs);
        }
    }
    
    for(; a!= stack.end(); ++a) {
        
        std::deque<ListDigraph::Arc>::iterator prev = a;
        --prev;
        
        
        transcript_unsecurity::evidence ev;
        if (unsecurityId[wc.source(*a)].resolvable) {
            if (know_paths[*prev].has_path_evidence(wc.id(*a))) {
                ev = transcript_unsecurity::EVIDENCED;
            } else {
                ev = transcript_unsecurity::UNEVIDENCED;
            }
        } else {
            ev = transcript_unsecurity::GUESSED;
        }
        
        cc.transcripts.back()->unsecurity_id.push_back( transcript_unsecurity(unsecurityId[wc.source(*a)].id, ev));
        
        for (std::set<transcript_unsecurity>::iterator unsecs = unsecurityArc[*a]->begin(); unsecs != unsecurityArc[*a]->end(); ++unsecs) {
            cc.transcripts.back()->unsecurity_id.push_back(*unsecs);
        }
        
        if (cet[*a] == edge_types::EXON) {
             cc.transcripts.back()->found_edge.id |= ces[*a].id;
        }
        
        if (fc[*a] < cc.transcripts.back()->flow )  {
            cc.transcripts.back()->flow = fc[*a];
        }
    }
    --a;
    cc.transcripts.back()->cycle_id_out = cycle_id_out[*a];

    // this is a case if we have single exon flow extraction
    if (cc.transcripts.back()->found_edge.id.empty()) {
        // we have to just create the foundedge to the id
       int node_index = cc.transcripts.back()->found_edge.node_index;
       
       #ifdef ALLOW_DEBUG
       logger::Instance()->debug("node_index " + std::to_string(node_index)+"\n"); 
       #endif
       cc.transcripts.back()->found_edge = exon_edge(size);
       cc.transcripts.back()->found_edge.set(node_index, true); 
        
    }
    
    return cc.transcripts.back()->flow;
}


capacity_mean path_finder::remove_full_path(std::deque<ListDigraph::Arc> &stack, capacity_type cap,
        ListDigraph &wc,
        ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        ListDigraph::ArcMap<unsigned int>& cycle_id_in, ListDigraph::ArcMap<unsigned int>& cycle_id_out, 
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id> &unsecurityId) {
    
    std::deque<ListDigraph::Arc>::iterator a = stack.begin(); 
    ListDigraph::Arc left = *a;
    ++a;
    for (;a != stack.end(); ++a) {
        
        // logging stuff to remove
        
//        for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
//            
//           logger::Instance()->debug("Arc " + std::to_string(wc.id(a))); 
//           for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
//               logger::Instance()->debug(" f " + std::to_string(ab->first)); 
//           }
//           for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
//               logger::Instance()->debug(" b " + std::to_string(*ab)); 
//           }
//           logger::Instance()->debug("\n" );  
//        }
        
        // end logging stuff
        
        left = base_manager::extract_path(left, *a, cap,  know_paths, know_back_paths, wc, fc, mc, ces, cet , cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
        
//        logger::Instance()->debug("Removed " + std::to_string( wc.id(*a) ) +" " + std::to_string(wc.id(left))+ "\n");
//        digraphWriter(wc, std::cout)
//            .arcMap("edge_specifier", ces)
//            .arcMap("edge_type", cet)
//            .arcMap("flow", fc)
//            .run();  
        
    }
      
    capacity_mean res = mc[left];
    wc.erase(left); // erased source target arc of path we extracted stepwise
      
    return res;
}


void path_finder::update_known(ListDigraph::ArcMap<arc_bridge>& know_paths, ListDigraph &g, ListDigraph& wc, ListDigraph::ArcMap<ListDigraph::Arc>& ref) {
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        
        lazy<path_evidence_map<int, rcount> > new_bridge;
        
        for (arc_bridge::iterator it = know_paths[a].begin(); it != know_paths[a].end(); ++it) {
//             logger::Instance()->debug("KF Trans " + std::to_string(wc.id(a))  + " ; " + std::to_string(it->first) + " to " + std::to_string(wc.id( ref[ g.arcFromId(it->first) ] )) + "\n");
            new_bridge->insert(std::make_pair( wc.id( ref[ g.arcFromId(it->first) ] ) , it->second)); 
        }
        
        know_paths[a].bridges = new_bridge;
    }
}


void path_finder::update_known_back(ListDigraph::ArcMap<arc_back_bridge>& know_back_paths, ListDigraph &g, ListDigraph& wc, ListDigraph::ArcMap<ListDigraph::Arc>& ref) {
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        
        lazy<path_evidence_set<int>  > new_bridge;
        
        for (arc_back_bridge::iterator it = know_back_paths[a].begin(); it != know_back_paths[a].end(); ++it) {
//            logger::Instance()->debug("KB Trans "  + std::to_string(wc.id(a))  + " ; " + std::to_string(*it) + " to " + std::to_string(wc.id( ref[ g.arcFromId(*it) ])) + "\n");
            new_bridge->insert(wc.id( ref[ g.arcFromId(*it) ] ));
        }
        
        know_back_paths[a].bridges = new_bridge;
    }
}

void path_finder::reverse_arc_cross_map(ListDigraph& g, ListDigraph& wc, ListDigraph::ArcMap<ListDigraph::Arc> & arc_ref, ListDigraph::ArcMap<ListDigraph::Arc> & arc_ref_rev) {
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        ListDigraph::Arc tmp = arc_ref[a];
        arc_ref_rev[tmp] = a;
        
    }
    
}


void path_finder::extract_guided_transcripts(alternative_transcript_collection& results, graph_list<exon_group*> guided) {
 
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("-------------- Establish Guided Transcript order.\n");
    #endif

    std::deque<std::deque<ListDigraph::Arc> > matched_guided;
    graph_list<exon_group *> matched_exons;
    
    for (graph_list<exon_group *>::iterator eg_it = guided.begin(); eg_it != guided.end(); ++eg_it) { 
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Use Guide " + (*eg_it)->bin_mask.to_string() +".\n");
        #endif
        
        exon_edge * transc = &(*eg_it)->bin_mask;
        
        // we search if there is a fully evidenced path around this
        
        matched_guided.push_back(std::deque<ListDigraph::Arc>());
        
        unsigned int start_index = transc->id.find_first();
        boost::dynamic_bitset<>::size_type index = start_index;
        unsigned int end_index = start_index;
        while(index != boost::dynamic_bitset<>::npos) {
            end_index = index;
            index = transc->id.find_next(index);
        } 
        
        ListDigraph::Arc null_arc = ListDigraph::ArcIt(INVALID);
        if (recursive_arc_backtracing(*transc, start_index, end_index, s, null_arc, matched_guided.back(), false) ) {
            matched_exons.push_back(*eg_it);
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Found \n");
            #endif
        } else {
            matched_guided.pop_back();
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("GUIDE ERROR 1\n");
            #endif
        }
    }
    
    std::unordered_map<int,int> counter;
    // we try to find a guide that has a unique arc
    for (std::deque<std::deque<ListDigraph::Arc> >::iterator fi = matched_guided.begin(); fi != matched_guided.end(); ++fi) {
        for (std::deque<ListDigraph::Arc>::iterator inner_it = fi->begin(); inner_it != fi->end(); ++inner_it) {
            ++counter[wc.id(*inner_it)];
        }
    }
    
    graph_list<exon_group *> ordered_guides;
    
    while (!matched_guided.empty()) {

        bool found = false;
        
        graph_list<exon_group *>::iterator me_it = matched_exons.begin();
        std::deque<std::deque<ListDigraph::Arc> >::iterator fi = matched_guided.begin();
        for (; fi != matched_guided.end() && !found; fi++, ++me_it) {
            for (std::deque<ListDigraph::Arc>::iterator inner_it = fi->begin(); inner_it != fi->end() && !found; ++inner_it) {
                if ( counter[wc.id(*inner_it)] == 1 ) {

                    for (std::deque<ListDigraph::Arc>::iterator cup_it = fi->begin(); cup_it != fi->end(); ++cup_it) {
                        --counter[wc.id(*cup_it)];
                    }
                    
                    ordered_guides.push_back(*me_it);
                    
                    matched_exons.erase(me_it);
                    matched_guided.erase(fi);
                    
                    found = true; // this breaks out the fors
                }
            }
        }
        
        if (!found) { // this means ALL transcripts
                // TODO: MAJOR, in this cause we have a full system of possible solutions that cannot be simply opened
                // can we just take a non-max junction? that could risk double transcripts...
                // for now we just choose a random one
                
            for (std::deque<ListDigraph::Arc>::iterator cup_it = matched_guided.front().begin(); cup_it != matched_guided.front().end(); ++cup_it) {
                --counter[wc.id(*cup_it)];
            }
            ordered_guides.push_back(matched_exons.front());
            
            matched_guided.pop_front();
            matched_exons.pop_front();
        } 
    }
   
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("-------------- Done.\n");
    #endif

    extract_guided_transcripts_linear_order(results, ordered_guides);
}


// OLD LINEAR VERSION, better if we find a corrected order among them first
void path_finder::extract_guided_transcripts_linear_order(alternative_transcript_collection& results, graph_list<exon_group*> guided) { 
 
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("-------------- Extract Guided Transcripts.\n");
    #endif

    int offset = results.transcripts.size();
    for (graph_list<exon_group *>::iterator eg_it = guided.begin(); eg_it != guided.end(); ++eg_it) { 
        
         exon_edge * transc = &(*eg_it)->bin_mask;
        
        // we search if there is a fully evidenced path around this
        std::deque<ListDigraph::Arc> path;
        unsigned int start_index = transc->id.find_first();
        boost::dynamic_bitset<>::size_type index = start_index;
        unsigned int end_index = start_index;
        while(index != boost::dynamic_bitset<>::npos) {
            end_index = index;
            index = transc->id.find_next(index);
        } 
        
        ListDigraph::Arc null_arc = ListDigraph::ArcIt(INVALID);
        if ( recursive_arc_backtracing(*transc, start_index, end_index, s, null_arc, path, false) ) {
           add_path_to_collection(path, results, wc, cfc, ces, cet, kp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
           results.transcripts.back()->guided = true;
           results.transcripts.back()->guide_reference = (*eg_it)->reference_name;    
        }
    }
    
    alternative_transcript_collection ac;
    int i = 0;
    for (graph_list<exon_group *>::iterator eg_it = guided.begin(); eg_it != guided.end(); ++eg_it, ++i) { 
          
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Use Guide " + (*eg_it)->bin_mask.to_string() +".\n");
        #endif
        
        exon_edge * transc = &(*eg_it)->bin_mask;
        
        // we search if there is a fully evidenced path around this
        std::deque<ListDigraph::Arc> path;
        unsigned int start_index = transc->id.find_first();
        boost::dynamic_bitset<>::size_type index = start_index;
        unsigned int end_index = start_index;
        while(index != boost::dynamic_bitset<>::npos) {
            end_index = index;
            index = transc->id.find_next(index);
        } 
        
        ListDigraph::Arc null_arc = ListDigraph::ArcIt(INVALID);
        if ( recursive_arc_backtracing(*transc, start_index, end_index, s, null_arc, path, false) ) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Path found " + std::to_string( path.size() )+ " ");
            for (std::deque<ListDigraph::Arc>::iterator it = path.begin(); it!= path.end(); ++it) {
                  logger::Instance()->debug(std::to_string( wc.id(*it) )+ ", " + ces[*it].to_string() + "\n");
            }
            logger::Instance()->debug("\n");
            #endif
            
            capacity_type cap = add_path_to_collection(path, ac, wc, cfc, ces, cet, kp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
            // now we need to remove this path from the copy
            capacity_mean mean = remove_full_path(path, cap, wc, cfc, mc,  ces, cet, cel, kp, kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId); 

            results.transcripts[offset+i]->flow = cap;
            results.transcripts[offset+i]->mean = mean.mean;
            results.transcripts[offset+i]->score = mean.compute_score();
        } else {

            results.transcripts[offset+i]->flow = 1;
            results.transcripts[offset+i]->mean = 1;
            results.transcripts[offset+i]->score = 1;
          //  logger::Instance()->error("GUIDE ERROR 2\n");
        }
    }
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("-------------- Done.\n");
    #endif
}

bool path_finder::recursive_arc_backtracing(exon_edge& goal, unsigned int next_start, unsigned int goal_end_index, ListDigraph::Node n, ListDigraph::Arc la, std::deque<ListDigraph::Arc> &path, bool exon) {

    if (n == t) {
        // we did it, perhaps
        return exon && next_start == goal_end_index;
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Recursion " + std::to_string(next_start) + " n " + std::to_string(wc.id(n))+ " la " + std::to_string(wc.id(la)) +".\n");
     #endif
    
    // we check all outward arcs to find the correct path through the maze
    // we are linear in # edges, but in reality much better
    for (ListDigraph::OutArcIt a(wc, n) ; a!=INVALID; ++a) {
        
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Look a " + std::to_string(wc.id(a))+".\n");
        #endif
        
        if (cet[a] == edge_types::BACKLINK) { // ignoring circular references, otherwise we have a BAD time
            continue;
        }
        
        bool correct_edge = false;
        if (cet[a] != edge_types::EXON) { // all edges except for circular and labeled exon edges
            // just jump add
            bool correct = exon;
            if (ces[a].node_index == next_start && ces[a].node_index == goal_end_index) {
                correct = true;
            }
            
            correct_edge = recursive_arc_backtracing(goal, next_start, goal_end_index, wc.target(a), a, path, correct);
        } else {
            
            unsigned int start_index = ces[a].id.find_first();
            
            // this is an EXON here
            if (start_index != next_start) {
                continue;
            }
            
            // for guides path evidences do not matter!
            
            unsigned int end_index; // this is slow but meh, doesn't happen too often
            boost::dynamic_bitset<>::size_type index = start_index;
            while(index != boost::dynamic_bitset<>::npos) {
                end_index = index;
                index = ces[a].id.find_next(index);
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test Cont " + ces[a].to_string() + " vs " + goal.to_string() + " ; " + std::to_string(start_index) + "-"+  std::to_string(end_index) +".\n");
            #endif
            // exon edge, we only follow this if it is compliant to the goal
            if ( ces[a].is_contained_in(goal, start_index, end_index) ) {
                
                // it is contained, so we take the route
                correct_edge = recursive_arc_backtracing(goal, end_index, goal_end_index, wc.target(a), a, path, true);
            }
        }
        if (correct_edge) { // we found the correct edge and propagate the solution down
            
            path.push_front(a);
            return true;
        }
    }
    
    return false;
}
