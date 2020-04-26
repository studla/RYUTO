/* 
 * File:   base_manager_flow_handler.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on 28.04.2016, 4:08 PM
 */

#include "base_manager.h"
#include <tuple>
#include <limits>
#include <unordered_set>
#include <deque>
#include "../../Options/options.h"
#include "../../Datatype_Templates/sets.h"
#include "path_finder/path_finder.h"
#include <unordered_map>
#include <lemon/lp.h>

#include <coin/ClpSimplex.hpp>
#include <coin/CoinHelperFunctions.hpp>
#include <coin/CoinBuild.hpp>

//#include <iostream>
//#include <fstream>

//void base_manager::extract_transcripts_from_flow(std::ostream &gs) {

void base_manager::extract_transcripts_from_flow() {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract Transcripts. " + std::to_string(input_count) + " \n");
    #endif
        
    if (meta->size == 1) { // TODO: REAL handling
        return;
    }
   
    #ifdef ALLOW_DEBUG
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif

//    std::ofstream pre_denoise;
//    pre_denoise.open ("pre_denoise.graphs", std::fstream::app);
//    digraphWriter(g, pre_denoise)
//            .arcMap("Identifier", ai)
//            .arcMap("Flow/Capacity", fs)
//            .node("Source", s)
//            .node("Drain", t)
//            .run();
//    pre_denoise << "--" << std::endl;
//    pre_denoise.close();

    // downsize neighbouring exon junctions, as they are overrepresented often
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::EXON 
                && ai[a].edge_specifier.left_consecutive && ai[a].edge_specifier.right_consecutive 
                && ai[a].edge_specifier.id.count()) {
            
            // downsize this one!
            ListDigraph::InArcIt left_node(g, g.source(a));
            ListDigraph::OutArcIt right_node(g, g.target(a));
            
            for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin(); ssi != fs[a].series.end(); ++ssi) {
                
                int id = ssi->first;
                
                float this_ave = ssi->second.region_average;
                float left_ave = fs[left_node].series[id].region_average;
                float right_ave = fs[right_node].series[id].region_average;

                float sum_left = 0;
                for (ListDigraph::OutArcIt io(g, g.target(left_node)); io != INVALID; ++io) {
                    sum_left += fs[io].series[id].region_average;
                }

                float sum_right = 0;
                for (ListDigraph::InArcIt io(g, g.source(right_node)); io != INVALID; ++io) {
                    sum_right += fs[io].series[id].region_average;
                }

                if (sum_left == 0 || sum_right == 0) {
                    continue;
                }
              
                float min_ave = std::min(left_ave*this_ave/sum_left, right_ave*this_ave/sum_right);

                float cap = std::min(min_ave, (float) ssi->second.capacity);
            
                ssi->second.capacity = cap;
                ssi->second.mean = capacity_mean(cap, ssi->second.length); // reset mean also
            }
        }
    }  
      
    #ifdef ALLOW_DEBUG
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif
    
    for (ListDigraph::OutArcIt o(g,s); o != INVALID; ++o) {
        
        ListDigraph::OutArcIt pe(g, g.target(o));
        
        ListDigraph::InArcIt tmp_in(g, g.source(pe));
        if (tmp_in == o ) { ++tmp_in;}
        if (tmp_in != INVALID ) { 
            continue;
        }
        
        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[pe].series.begin(); ssi != fs[pe].series.end(); ++ssi) { 

            if (ssi->second.length == 0) {
                continue; // skip fully uninitialized
            }

            float cap = ssi->second.average_to_first_zero_from_right;
            cap = std::max(cap, 1.0f);
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Reset Source " + std::to_string(g.id(pe)) + " at "+ std::to_string(ssi->first) + " " + std::to_string(cap)+ "-"+ std::to_string(ssi->second.region_average) + ".\n");
            #endif
            
            ssi->second.capacity = cap;
            ssi->second.mean = capacity_mean(cap, ssi->second.length); // reset mean also
        }
           
    }
    
    for (ListDigraph::InArcIt i(g,t); i != INVALID; ++i) {
        
        ListDigraph::InArcIt pe(g, g.source(i));
        
        ListDigraph::OutArcIt tmp_out(g, g.target(pe));
        if (tmp_out == i ) { ++tmp_out;}
        if (tmp_out != INVALID ) { 
            continue;
        }
        
        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[pe].series.begin(); ssi != fs[pe].series.end(); ++ssi) { 
            if (ssi->second.length == 0) {
                continue; // skip fully uninitialized
            }

            float cap = ssi->second.average_to_first_zero_from_left;
            cap = std::max(cap, 1.0f);
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Reset Drain " + std::to_string(g.id(pe)) + " at "+ std::to_string(ssi->first) + " " + std::to_string(cap)+ "-"+ std::to_string(ssi->second.region_average) + ".\n");
            #endif
            
            ssi->second.capacity = cap;
            ssi->second.mean = capacity_mean(cap, ssi->second.length); // reset mean also
        }
           
    }
    
    // now combine means of all edges  
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        fs[a].combined_capacity = 0;
        gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin();
        if (ssi == fs[a].series.end()) {
            fs[a].combined_mean = capacity_mean();
            continue;
        }
        
//        fs[a].combined_mean = capacity_mean();
//        for(; ssi != fs[a].series.end(); ++ssi) {
//            fs[a].combined_capacity += ssi->second.capacity;
//            if (ssi->second.mean.mean > fs[a].combined_mean.mean) {
//                fs[a].combined_mean = ssi->second.mean;
//            }
//        }       
        
        fs[a].combined_mean = ssi->second.mean;
        fs[a].combined_capacity += ssi->second.capacity;
        unsigned int index = 2;
        ++ssi;
        for(; ssi != fs[a].series.end(); ++ssi) {
            fs[a].combined_capacity += ssi->second.capacity;
            fs[a].combined_mean.sum_up(ssi->second.mean, index);
            ++index;
        }       
        
        fs[a].combined_mean.scores.clear();
        fs[a].combined_mean.scores.push_back(fs[a].combined_mean.mean);
    }
    
    #ifdef ALLOW_DEBUG
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif
         
   ListDigraph::ArcMap<bool> guided_saves(g);
   std::deque<std::deque<ListDigraph::Arc> > guided_paths;
   find_guide_sources_and_drains(guided_saves, guided_paths);
    
   // Guide Stats
   unsigned int total_edges = 0;
   unsigned int guided_edges = 0;
   for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
       total_edges++;
       if(guided_saves[a]) {
            guided_edges++;
       }
   }


   // this is a set of exons directly neighbouring other exons
   // we only mark those which are either filled introns or in between sources or drains
   ListDigraph::ArcMap<bool> consecutive_s_t(g); // inits as false
   ListDigraph::ArcMap<bool> marked_source(g);
   ListDigraph::ArcMap<bool> marked_drain(g);
   find_s_t_exons(guided_saves, marked_source, marked_drain, consecutive_s_t);
   
    #ifdef ALLOW_DEBUG
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::NODE && !marked_source[a] && !marked_drain[a]) {
           for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin(); ssi != fs[a].series.end(); ++ssi) { 
               ssi->second.mean.scores.clear();
           }
           fs[a].combined_mean.scores.clear();
        }
    }
   
   prune_sources_and_drains_adjacent_starts(guided_saves); // directly adjacent/starts at exactly inner exons 

   #ifdef ALLOW_DEBUG
   print_graph_debug_copy(std::cout, g, fs, ai, s, t);
   #endif   
   
   if (!options::Instance()->is_keep_introns()) {
        filter_introns(guided_saves, marked_source, marked_drain);
        filter_broken_introns(guided_saves, marked_source, marked_drain);
   }

   #ifdef ALLOW_DEBUG
   print_graph_debug_copy(std::cout, g, fs, ai, s, t);
   #endif
   
    float low_mark = options::Instance()->get_low_edge_mark();
    remove_too_short_source_drain(guided_saves, marked_source, marked_drain, consecutive_s_t); 
    while (true) {
                
        remove_dead_ends();         
        ListDigraph::ArcIt a(g);
        if(a == INVALID) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Emergency Out.\n");
            #endif
            return;
        }
        
        #ifdef ALLOW_DEBUG
        print_graph_debug_copy(std::cout, g, fs, ai, s, t);
        #endif
        
        if (erase_low_boundary_st(guided_saves, marked_source, marked_drain, consecutive_s_t)) { //1
            continue;
        }
        
        if (erase_low_deviation_st(guided_saves, marked_source, marked_drain, consecutive_s_t)) { //2
            continue;
        }
        
//        if (remove_small_low_st(guided_saves, marked_source, marked_drain, consecutive_s_t)) { //3
//            continue;
//        }
        
        if (prune_small_junctions(guided_saves, marked_source, marked_drain, consecutive_s_t)) { //4
            continue;
        }

        if (threshold_filter(guided_saves, marked_source, marked_drain, consecutive_s_t, low_mark)) { //5 
            continue;
        }
        
        break;
    }
   
    remove_dead_ends(); 
    ListDigraph::ArcIt a(g);
    if(a == INVALID) {
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Locus destroyed.\n");
        #endif
        return;
    }

    
    for (std::set<int>::iterator idi = input_ids.begin(); idi != input_ids.end(); ++idi){
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Denoise/Flow " + std::to_string(*idi) + ".\n"); 
        print_graph_debug_copy(std::cout, g, fs, ai, s, t);
        #endif

        for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {  // we get rid of 0s in multi sets

            if(fs[a].series[*idi].capacity == 0 && ai[a].edge_type != edge_types::HELPER) {
                //logger::Instance()->debug("Correct " + std::to_string(g.id(a)) + " " + std::to_string(*idi) + ".\n");

                fs[a].series[*idi].capacity = 1;
                fs[a].series[*idi].mean = capacity_mean(1, std::max(fs[a].series[*idi].length, 1ul));
            }
        }
        
        if (options::Instance()->is_graph_denoising()) {

            if (guided_paths.empty()) {
                denoise_graph(guided_saves, marked_source, marked_drain, *idi);
            } else {
                bool double_denoise = (guided_edges *100 / total_edges) < 90;

                if (double_denoise) { 
                    denoise_graph(guided_saves, marked_source, marked_drain, *idi);
                }
                denoise_graph_guided(guided_paths, *idi, double_denoise);
            }
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Before Flow.\n");
        print_graph_debug_copy(std::cout, g, fs, ai, s, t);
        #endif
        
        finalize_flow_graph(*idi);
        compute_flow(*idi);
    }
    
    // join individual flow for joint work
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        fs[a].combined_flow = 0;
        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin(); ssi != fs[a].series.end(); ++ssi) {
            fs[a].combined_flow += ssi->second.flow;
            ssi->second.capacity = 0;
        }
        fs[a].combined_capacity = 0;        
    }
    
//    std::ofstream post_denoise;
//    post_denoise.open ("post_denoise.graphs", std::fstream::app);
//    digraphWriter(g, post_denoise)
//            .arcMap("Identifier", ai)
//            .arcMap("Flow/Capacity", fs)
//            .node("Source", s)
//            .node("Drain", t)
//            .run();
//    post_denoise << "--" << std::endl;
//    post_denoise.close();
//    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("After Flow.\n");
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif
    
    // from here on out, we treat the combined data as a unique input from here on out, if we have more than one id!
    if (input_ids.size() > 1) input_ids.insert(-1); 
    
    // in some cases we might have lost guides though due to poor performance of flow!
    for(std::deque<std::deque<ListDigraph::Arc> >::iterator it = guided_paths.begin(); it != guided_paths.end(); ++it) {
        for(std::deque<ListDigraph::Arc>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
            for (std::set<int>::iterator idi = input_ids.begin(); idi != input_ids.end(); ++idi){
                ++fs[*it2].get_flow(*idi);
            }
        }
    }
    
    // remove all possible 0 flow edges by vote!
    for (ListDigraph::ArcIt a(g); a != INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        
        bool all_zero = true;
        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[arc].series.begin(); ssi != fs[arc].series.end(); ++ssi) {
            if (ssi->second.flow != 0) { // online erase edges
                all_zero = false;
                break;
            } 
        }
        if (!all_zero) {
            continue;
        }
    
        ListDigraph::Node target = g.target(arc);
        ListDigraph::Node source = g.source(arc);
        g.erase(arc);
        if (target != t) {
            ListDigraph::OutArcIt o(g,target);
            ListDigraph::InArcIt i(g,target);
            if (i == INVALID && o == INVALID) {
                g.erase(target);
            }
        }
        if (source != s) {
            ListDigraph::OutArcIt o(g,source);
            ListDigraph::InArcIt i(g,source);
            if (i == INVALID && o == INVALID) {
                g.erase(source);
            }
        }   
    }
    
    // we need degree counts from here on:
    InDegMap<ListDigraph> in_deg(g);
    OutDegMap<ListDigraph> out_deg(g);
    
    ListDigraph::ArcMap<arc_bridge> know_paths(g);
    ListDigraph::ArcMap<arc_back_bridge> know_back_paths(g);

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Cleaned.\n");
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif
    
    // we remove path node that are one to one, one to many, or many to one
    // this keeps flow rules intact
    while ( contract_composites(g, fs, ai, in_deg, out_deg) ); // of this runs too often we could make a top. sorting

//    std::ofstream contracted;
//    contracted.open ("contracted.graphs", std::fstream::app);
//    digraphWriter(g, contracted)
//            .arcMap("Identifier", ai)
//            .arcMap("Flow/Capacity", fs)
//            .node("Source", s)
//            .node("Drain", t)
//            .run();
//    contracted << "--" << std::endl;
//    contracted.close();
//    return;
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Contracting done.\n");
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif
    
    struct raw_unsec {
        raw_unsec() {}
        raw_unsec(rpos r, rpos l, unsigned int id) : left(l), right(r), id(id) { 
        }
        
        rpos left;
        rpos right;
        unsigned int id;
    };
    
    ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > unsecurityArc(g);
    ListDigraph::NodeMap< unsecurity_id > unsecurityId(g);
    gmap<unsigned int, raw_unsec> unsecIdToNodeIndex;
    // we keep this for the transcripts
    
    unsigned int nid = 1;
    for (ListDigraph::NodeIt n(g); n != INVALID;++n) {
        unsecurityId[n] = unsecurity_id(nid, true);
        unsecIdToNodeIndex.insert(std::make_pair(nid, raw_unsec(meta->exons[node_index[n]].left, meta->exons[node_index[n]].right, node_index[n])));
        ++nid;
    }
    
    
     #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract all IDs in order.\n");
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif
    for (std::set<int>::iterator idi = input_ids.begin(); idi != input_ids.end(); ++idi) {
        
        if (input_ids.size() > 1 && !options::Instance()->is_compute_all_singles() && *idi != -1) {
            continue;
        }
        
        int guiding = *idi;
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Guide ID " + std::to_string(guiding) + ".\n");
        #endif
        
        ListDigraph gc;
        ListDigraph::ArcMap<arc_identifier> aic(gc);
        ListDigraph::ArcMap<flow_series> fsc(gc);
        ListDigraph::ArcMap<arc_bridge> kpc(gc);
        ListDigraph::ArcMap<arc_back_bridge> kbpc(gc);
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > uac(gc);
        ListDigraph::NodeMap< unsecurity_id > uic(gc);
        ListDigraph::NodeMap<unsigned int> nic(gc);
        ListDigraph::Node sc, tc;

        // we need a working copy for each one
        DigraphCopy<ListDigraph, ListDigraph> copy(g, gc);
        copy.arcMap(ai, aic)
            .arcMap(fs, fsc)
            .arcMap(know_paths, kpc)
            .arcMap(know_back_paths, kbpc)
            .arcMap(unsecurityArc, uac)
            .nodeMap(node_index, nic)
            .nodeMap(unsecurityId, uic)
            .node(s,sc)
            .node(t,tc);
        copy.run();
//        // self made copy for better comparison to old tool!
//        ListDigraph::NodeMap<ListDigraph::Node> c_ref(g);
//        for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
//            ListDigraph::Node nn = gc.addNode();
//            c_ref[n] = nn;
//            nic[nn] = node_index[n];
//            uic[nn] = unsecurityId[n];
//            
//             #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Copy N " + std::to_string(g.id(n)) + " to " + std::to_string(g.id(nn)) + ".\n");
//            #endif
//        }
//        std::deque<std::pair<ListDigraph::Arc, capacity_type> > edges_cs;
//        for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//            edges_cs.push_back(std::make_pair(a, fs[a].get_flow(guiding)));
//        }
//        std::sort(edges_cs.begin(), edges_cs.end(), []( std::pair<ListDigraph::Arc, capacity_type> & a, std::pair<ListDigraph::Arc, capacity_type> & b) -> bool {return a.second > b.second;});
//        
//        for (std::deque<std::pair<ListDigraph::Arc, capacity_type> >::iterator ecsi = edges_cs.begin(); ecsi != edges_cs.end(); ++ecsi) {
//            ListDigraph::Arc a = ecsi->first;
//            ListDigraph::Arc na = gc.addArc(c_ref[g.source(a)], c_ref[g.target(a)]);
//            aic[na] = ai[a];
//            fsc[na] = fs[a];
//            kpc[na] = know_paths[a];
//            kbpc[na] = know_back_paths[a];
//            uac[na] = unsecurityArc[a];
//            
//             #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Copy A " + std::to_string(g.id(a)) + " to " + std::to_string(g.id(na)) + ".\n");
//            #endif
//        }
//        sc = c_ref[s];
//        tc = c_ref[t];
//        
        // remove possible zero edges for this individual
        for (ListDigraph::ArcIt a(gc); a != INVALID; ) {
            ListDigraph::Arc arc(a);
            ++a;

            if (fsc[arc].get_flow(guiding) > 0) {
                continue;
            }

            ListDigraph::Node target = gc.target(arc);
            ListDigraph::Node source = gc.source(arc);
            erase_arc(gc, arc, fsc, input_ids, transcripts[guiding]);
            if (target != t) {
                ListDigraph::OutArcIt o(gc,target);
                ListDigraph::InArcIt i(gc,target);
                if (i == INVALID && o == INVALID) {
                    gc.erase(target);
                }
            }
            if (source != s) {
                ListDigraph::OutArcIt o(gc,source);
                ListDigraph::InArcIt i(gc,source);
                if (i == INVALID && o == INVALID) {
                    gc.erase(source);
                }
            }   
        }
        
        ListDigraph::ArcIt a(gc);
        if(a == INVALID) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Locus destroyed.\n");
            #endif
            return;
        }
        
        InDegMap<ListDigraph> in_deg_c(gc);
        OutDegMap<ListDigraph> out_deg_c(gc);
        while ( contract_composites(gc, fsc, aic, in_deg_c, out_deg_c) ); // of this runs too often we could make a top. sorting
        
       #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Guide Before.\n");
        print_graph_debug_copy(std::cout, gc, fsc, aic, sc, tc);
        #endif
        
        path_finder* pf = path_finder::create_path_finder(gc, sc, tc, fsc, aic, nic, kpc, kbpc, uac, uic, input_ids, transcripts, raw->size, meta);
        pf->extract_guided_transcripts(transcripts[guiding], raw->guide_transcripts, guiding);    // test and extract evidenced transcripts 
        delete pf; 

        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Guide after.\n");
        print_graph_debug_copy(std::cout, gc, fsc, aic, sc, tc);
        #endif

        ListDigraph::ArcMap<bool> barred(gc);
        while(simplify_ambiguous_ILP(gc, sc, tc, fsc, aic, nic, kpc, kbpc, uac, uic, guiding, barred, input_ids, transcripts[guiding]));

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Simplify done.\n");
        print_graph_debug_copy(std::cout, gc, fsc, aic, sc, tc);
        #endif

        // we have a simplified graph that we can now extract transcripts from
        // at this point we get the maximum amount of transcripts possible out of the existing DAG
        pf = path_finder::create_path_finder(gc, sc, tc, fsc, aic, nic, kpc, kbpc, uac, uic, input_ids, transcripts, raw->size, meta);
        pf->extract_transcripts(transcripts[guiding], guiding); // normal path extraction
        delete pf; 

        // set exon positions for unsecurities
        for (graph_list<lazy<transcript> >::iterator t = transcripts[guiding].transcripts.begin(); t != transcripts[guiding].transcripts.end(); ++t) {

            if (!(*t).has_ref()) {
                continue;
            }

            for ( std::deque<transcript_unsecurity>::iterator u = (*t)->unsecurity_id.begin(); u!= (*t)->unsecurity_id.end(); ++u) {
                u->position = unsecIdToNodeIndex[u->id].id;
                u->left = unsecIdToNodeIndex[u->id].left;
                u->right = unsecIdToNodeIndex[u->id].right;
            }
        }
    }
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Pathfinder done " + std::to_string(transcripts.size())  + ".\n");
    print_graph_debug_copy(std::cout, g, fs, ai, s, t);
    #endif
    
}

void base_manager::find_guide_sources_and_drains(ListDigraph::ArcMap<bool> &guided_saves, std::deque<std::deque<ListDigraph::Arc> > &paths) {
    
        logger::Instance()->debug("Save Guide Arcset " +std::to_string(raw->guide_transcripts.size()) + " .\n");

    for (graph_list<exon_group *>::iterator eg_it = raw->guide_transcripts.begin(); eg_it != raw->guide_transcripts.end(); ++eg_it) {
        
        logger::Instance()->debug("Attempt save Guide Arcset " + (*eg_it)->bin_mask.to_string() + ".\n");

        exon_edge * transc = &(*eg_it)->bin_mask;
        
        std::deque<ListDigraph::Arc> path;
        unsigned int start_index = transc->id.find_first();
        boost::dynamic_bitset<>::size_type index = start_index;
        unsigned int end_index = start_index;
        while(index != boost::dynamic_bitset<>::npos) {
            end_index = index;
            index = transc->id.find_next(index);
        } 
        ListDigraph::Arc null_arc = ListDigraph::ArcIt(INVALID);

        if (recursive_arc_backtracing(*transc, start_index, end_index, s, null_arc, path, false) ) {
            
            paths.push_back(std::deque<ListDigraph::Arc>());
            for (std::deque<ListDigraph::Arc>::iterator pi = path.begin(); pi!= path.end(); ++pi) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Save Arc: " + std::to_string(g.id(*pi)) + "\n");
                #endif
                
                guided_saves[*pi] = true;
                paths.back().push_back(*pi);
            }
        }
    }
 
}


void base_manager::prune_sources_and_drains_adjacent_starts(ListDigraph::ArcMap<bool> &guided_saves) {
    
    for (ListDigraph::OutArcIt o(g,s); o != INVALID; ) {
        
        ListDigraph::Arc arc(o);
        ++o;
        
        if (guided_saves[arc]) { // no need to test, keep it!
            continue;
        }
        
        bool source_without_own_arc = false;
        ListDigraph::InArcIt ni(g, g.target(arc));
        for (;ni != INVALID; ++ni) {
            if (ni != arc) { 
                source_without_own_arc = true;
            }
        }
        
        if ( source_without_own_arc) {
            //erase
            #ifdef ALLOW_DEBUG
	    logger::Instance()->debug("Erase Source " + std::to_string(g.id(arc)) + ".\n");
            #endif

            g.erase(arc);
        }
    }

    for (ListDigraph::InArcIt i(g,t); i != INVALID; ) {
        
        ListDigraph::Arc arc(i);
        ++i;
        
        if (guided_saves[arc]) { // no need to test, keep it!
            continue;
        }
        
        bool drain_without_own_arc = false;
	ListDigraph::OutArcIt no(g, g.source(arc)); // first node this belongs to!
        for (;no != INVALID; ++no) {
            if (no != arc) { 
                drain_without_own_arc = true;
            }
        }

        if ( drain_without_own_arc) {
            //erase
            #ifdef ALLOW_DEBUG
	    logger::Instance()->debug("Erase Drain " + std::to_string(g.id(arc)) + ".\n");
            #endif
            
            g.erase(arc);
                       
        }
    }
}


void base_manager::find_s_t_exons(ListDigraph::ArcMap<bool> &guided_saves,
        ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Mark Sources and Drains \n");
    #endif
    
    // we mark all sources and drain connections first
    std::set<int> sources_to_test; 
    for (ListDigraph::OutArcIt o(g,s); o != INVALID; ++o) {
    
        ListDigraph::Arc follow(o);
        marked_source[follow] = true;
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Mark S " + std::to_string(g.id(follow)) +"\n");
        #endif
        
        while (true) {
            ListDigraph::InArcIt tmp_in(g, g.target(follow));
            if (tmp_in == follow ) { ++tmp_in;}
            ListDigraph::OutArcIt tmp_out(g, g.target(follow));
            if (tmp_out == INVALID ) { // next has out 
                break;
            }
            
            ListDigraph::Arc next(tmp_out);
            ++tmp_out;
            
            if (tmp_out != INVALID || tmp_in != INVALID) { // as long as we have only one connection, but no in edge!, we follow!
                sources_to_test.insert(g.id(next)); // has to be a node!
                break;
            }
            
            follow = next;
            marked_source[follow] = true;
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Mark S " + std::to_string(g.id(follow)) +"\n");
            #endif
        }
    }
    std::set<int> drains_to_test; 
    for (ListDigraph::InArcIt i(g,t); i != INVALID; ++i) {
    
        ListDigraph::Arc follow(i);
        marked_drain[follow] = true;
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Mark T " + std::to_string(g.id(follow)) +"\n");
        #endif
        
        while (true) {
            ListDigraph::InArcIt tmp_in(g, g.source(follow));
            if (tmp_in == INVALID ) { // next has arc 
                break;
            }
            ListDigraph::OutArcIt tmp_out(g, g.source(follow));
            if (tmp_out == follow ) { ++tmp_out;}
            
            ListDigraph::Arc next(tmp_in);
            ++tmp_in;
            
            if (tmp_out != INVALID || tmp_in != INVALID) { // as long as we have only one connection, but no in edge!, we follow!
                drains_to_test.insert(g.id(next)); // has to be a node!
                break;
            }
            
            follow = next;
            marked_drain[follow] = true;
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Mark T " + std::to_string(g.id(follow)) +"\n");
            #endif
        }
    }
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Find Consecutive \n");
    #endif
    // neighbouring sources
    
    std::vector<std::deque<ListDigraph::Arc> > same_node(meta->size);
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::NODE) { //with opened middle node
            int ni = ai[a].edge_specifier.node_index;
            same_node[ni].push_back(a); 
        }
    }
    
    
    for (std::set<int>::iterator src = sources_to_test.begin(); src != sources_to_test.end(); ++src ) {
        
        ListDigraph::Arc node2 = g.arcFromId(*src);
        
        if (ai[node2].edge_type != edge_types::NODE) { // cannot happen if we end
            continue;
        }
        
        // do we have a non-source adjourning node?
        bool non_source_adjourning = false;
        for (std::deque<ListDigraph::Arc>::iterator nr = same_node[ai[node2].edge_specifier.node_index].begin(); nr != same_node[ai[node2].edge_specifier.node_index].end(); ++nr) {
            for (ListDigraph::InArcIt i(g, g.source(*nr) ); i != INVALID; ++i) {
                if (!marked_source[i]) {
                    non_source_adjourning = true;
                    break;
                }
            }
        }
        
        if (!non_source_adjourning) {
            continue;
        }
        
        // this means we could have a problematic inner intron end source!
        for (std::deque<ListDigraph::Arc>::iterator nr = same_node[ai[node2].edge_specifier.node_index].begin(); nr != same_node[ai[node2].edge_specifier.node_index].end(); ++nr) {
            for (ListDigraph::InArcIt i(g, g.source(*nr) ); i != INVALID; ++i) {
                if (marked_source[i] && ai[i].edge_type == edge_types::EXON) {

                    ListDigraph::Arc edge1(i);
                    ListDigraph::InArcIt node1(g, g.source(edge1));

                    if (guided_saves[edge1]) {
                        //don't block guides!
                        continue;
                    }

    //                if (regions[node1].get_average() < 3) {
    //                    // we erase all nodes below this threshold, no matter of overlap or not
    //                    #ifdef ALLOW_DEBUG
    //                    logger::Instance()->debug("Erase S by Coverage" + std::to_string(g.id(edge1)) + "\n");
    //                    #endif
    //                    g.erase(node1);
    //                    g.erase(edge1);
    //                }

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Source Test Complex " + std::to_string(g.id(edge1)) + " " + std::to_string(g.id(node2)) + "\n");
                    #endif
    //                rcount node2_count = regions[node2].get_average();

                    if (ai[edge1].edge_specifier.right_consecutive) {
                        consecutive_s_t[edge1] = true;
                        consecutive_s_t[node1] = true;
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Block Push S " + std::to_string(g.id(edge1)) + "\n");
                        #endif
                    }
                }
            } 
        }
    }
    // neighbouring drains
    for (std::set<int>::iterator drn = drains_to_test.begin(); drn != drains_to_test.end(); ++drn ) {
        
        ListDigraph::Arc node2 = g.arcFromId(*drn);
        
        if (ai[node2].edge_type != edge_types::NODE) { // cannot happen if we end
            continue;
        }
        
        // do we have a non-source adjourning node?
        bool non_drain_adjourning = false;
        for (std::deque<ListDigraph::Arc>::iterator nr = same_node[ai[node2].edge_specifier.node_index].begin(); nr != same_node[ai[node2].edge_specifier.node_index].end(); ++nr) {
            for (ListDigraph::OutArcIt o(g, g.target(*nr) ); o != INVALID; ++o) {
                if (!marked_drain[o]) {
                    non_drain_adjourning = true;
                    break;
                }
            }
        }
        
        if (!non_drain_adjourning) {
            continue;
        }
        
        // this means we could have a problematic inner intron end source!
        for (std::deque<ListDigraph::Arc>::iterator nr = same_node[ai[node2].edge_specifier.node_index].begin(); nr != same_node[ai[node2].edge_specifier.node_index].end(); ++nr) {
            for (ListDigraph::OutArcIt o(g, g.target(*nr) ); o != INVALID; ++o) {
                if (marked_drain[o] && ai[o].edge_type == edge_types::EXON) {

                    ListDigraph::Arc edge1(o);
                    ListDigraph::OutArcIt node1(g, g.target(edge1));

                     if (guided_saves[edge1]) {
                        //don't block guides!
                        continue;
                    }


    //                if (regions[node1].get_average() < 3) {
    //                    // we erase all nodes below this threshold, no matter of overlap or not
    //                    #ifdef ALLOW_DEBUG
    //                    logger::Instance()->debug("Erase T by Coverage" + std::to_string(g.id(edge1)) + "\n");
    //                    #endif
    //                    g.erase(node1);
    //                    g.erase(edge1);
    //                }

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Drain Test Complex " + std::to_string(g.id(edge1)) + " " + std::to_string(g.id(node2)) + "\n");
                    #endif
    //                rcount node2_count = regions[node2].get_average();

                    if (ai[edge1].edge_specifier.left_consecutive) {               
                        consecutive_s_t[edge1] =true;
                        consecutive_s_t[node1] =true;
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Block Push T " + std::to_string(g.id(edge1)) + "\n");
                        #endif
                    }
                }
            } 
        }
    }
}


void base_manager::filter_introns(ListDigraph::ArcMap<bool> &guided_saves,
        ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain) {
    
    
    std::vector<std::deque<ListDigraph::Node> > same_node_r(meta->size);
    std::vector<std::deque<ListDigraph::Node> > same_node_l(meta->size);
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::NODE) { //with opened middle node
            
            int ni = ai[a].edge_specifier.node_index;
            same_node_r[ni].push_back(g.target(a));
            same_node_l[ni].push_back(g.source(a));
            
        }
    }
    
    // test for in between neighbours 
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        if (guided_saves[a]) {
            continue;
        }
        
        if (ai[a].edge_type == edge_types::NODE) { //with opened middle node

            ListDigraph::InArcIt ni(g, g.source(a));
            ListDigraph::OutArcIt no(g, g.target(a));

            ListDigraph::Arc left_edge(ni);
            ListDigraph::Arc rigth_edge(no);

            ++ni;
            ++no;
            
            if (ni == INVALID && no == INVALID && ai[left_edge].edge_type == edge_types::EXON && ai[rigth_edge].edge_type == edge_types::EXON) {
                // first condition: only one connecting edge in each direction, otherwise this can't be an intron

                // each edge can have but one node
                ListDigraph::InArcIt left_node(g, g.source(left_edge));
                ListDigraph::OutArcIt right_node(g, g.target(rigth_edge));

                // check if this is really neighbouring and if bot exons are actually connected!
                if (meta->exons[ai[left_node].edge_specifier.node_index].right + 1 == meta->exons[ai[a].edge_specifier.node_index].left 
                        && meta->exons[ai[a].edge_specifier.node_index].right + 1 == meta->exons[ai[right_node].edge_specifier.node_index].left) {
                    
                    // potential intron here!
                    // check all inputs

                    bool breaked = false;
                    // look for direct connection, a single edge has to exist or intron would be further split, which is legally impossible
                    for (std::deque<ListDigraph::Node>::iterator nr = same_node_r[ai[left_node].edge_specifier.node_index].begin(); nr != same_node_r[ai[left_node].edge_specifier.node_index].end() && !breaked; ++nr) {
                        ListDigraph::OutArcIt search(g, *nr); // search is legally an edge and con only have one node output!
                        for (;search != INVALID && !breaked; ++search) {
                            
                            bool hit_target = false;
                            for (std::deque<ListDigraph::Node>::iterator nl = same_node_l[ai[right_node].edge_specifier.node_index].begin(); nl != same_node_l[ai[right_node].edge_specifier.node_index].end(); ++nl) {
                                
                                #ifdef ALLOW_DEBUG
                                logger::Instance()->debug("A " + std::to_string(g.id(search)) + " " + std::to_string(g.id(*nl)) + " " + std::to_string(g.id(*nr)) + "\n");
                                #endif  
                                
                                if (g.target(search) == *nl) {
                                    hit_target = true;
                                    break;
                                }
                            }
                            
                            if (hit_target) {

                                #ifdef ALLOW_DEBUG
                                logger::Instance()->debug("Test " + std::to_string(g.id(a)) + "\n");
                                #endif  
                                
                                unsigned int vote_count = input_count - fs[a].series.size();
                                for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin(); ssi != fs[a].series.end(); ++ssi) {
                                    if (fs[search].series[ssi->first].capacity != 0 && ( (ssi->second.capacity < 45 && ssi->second.region_average *100 / fs[search].series[ssi->first].capacity < options::Instance()->get_intron_retention_threshold()) || (ssi->second.capacity <= 2) ) )
                                    {
                                        ++vote_count;
                                    }
                                }
                               
                                if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_low)) {
                                    #ifdef ALLOW_DEBUG
                                    logger::Instance()->debug("Block Intron Edge " + std::to_string(g.id(a)) + " " + std::to_string(g.id(left_edge)) + " " + std::to_string(g.id(rigth_edge)) + "\n");
                                    #endif                                

                                    g.erase(a);
                                    g.erase(left_edge);
                                    g.erase(rigth_edge);
                                    breaked = true;
                                }
                            }
                        } 
                    }
                }
            }
        } else if (ai[a].edge_type == edge_types::EXON) {
            if (ai[a].edge_specifier.id.count() == 3) { // only then can middle node be the intron!
                
                unsigned int left = ai[a].edge_specifier.id.find_first();
                unsigned int middle = ai[a].edge_specifier.id.find_next(left);
                unsigned int right = ai[a].edge_specifier.id.find_next(middle);
                               
                if (meta->exons[left].right + 1 == meta->exons[middle].left 
                        && meta->exons[middle].right + 1 == meta->exons[right].left) {
                    
                    
                    ListDigraph::OutArcIt search(g, g.source(a));
                    for (;search != INVALID; ++search) {
                         
                        if (g.target(search) == g.target(a) && search != a) {
                    
                            ListDigraph::InArcIt left_node(g, g.source(a));
                            ListDigraph::OutArcIt right_node(g, g.target(a));

                            unsigned int vote_count = input_count - fs[a].series.size();
                            for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin(); ssi != fs[a].series.end(); ++ssi) {
                                if (fs[search].series[ssi->first].capacity != 0 && ( (ssi->second.capacity < 45 && ssi->second.region_average *100 / fs[search].series[ssi->first].capacity < options::Instance()->get_intron_retention_threshold()) || (ssi->second.capacity <= 2) ) )
                                {
                                    ++vote_count;
                                }
                            }
                            
                            if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_low)) {
                                #ifdef ALLOW_DEBUG
                                logger::Instance()->debug("Block Intron Edge " + std::to_string(g.id(a)) + "\n");
                                #endif                           

                                g.erase(a);

                            }
                        }
                    }
                }
            }
        }
    }
}

void base_manager::filter_broken_introns(ListDigraph::ArcMap<bool> &guided_saves,
        ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain) {
    
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Filter Broken Intron\n");
    #endif
    
    std::vector<std::deque<ListDigraph::Node> > same_node_r(meta->size);
    std::vector<std::deque<ListDigraph::Node> > same_node_l(meta->size);
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::NODE) { //with opened middle node
            
            int ni = ai[a].edge_specifier.node_index;
            same_node_r[ni].push_back(g.target(a));
            same_node_l[ni].push_back(g.source(a));
            
        }
    }    
        
    ListDigraph::ArcMap<bool> kill_status(g);
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::EXON && !marked_drain[a] && !marked_source[a]) {
            
            ListDigraph::InArcIt left_node(g, g.source(a));
            ListDigraph::OutArcIt right_node(g, g.target(a));
            
            
            if (left_node != INVALID && right_node != INVALID && ai[left_node].edge_type == edge_types::NODE && ai[right_node].edge_type == edge_types::NODE ) {
//                && edge_specifier[right_node].node_index == edge_specifier[left_node].node_index + 3 ) {
            //        && (!marked_drain[left_node] && !marked_drain[right_node]) && (!marked_source[left_node] && !marked_source[right_node]) {
            
                bool has_source = false;
                ListDigraph::Arc source, source_base;

                bool has_drain = false;
                ListDigraph::Arc drain, drain_base;

                for (std::deque<ListDigraph::Node>::iterator nl = same_node_l[ai[right_node].edge_specifier.node_index].begin(); nl != same_node_l[ai[right_node].edge_specifier.node_index].end(); ++nl) {
                    for (ListDigraph::InArcIt oi(g, *nl); oi != INVALID ;++oi) {
                        ListDigraph::InArcIt lr(g, g.source(oi));

                        if (marked_source[oi] && lr != INVALID && meta->exons[ai[lr].edge_specifier.node_index].right + 1 == meta->exons[ai[right_node].edge_specifier.node_index].left) {
                            has_source = true;
                            source = oi;
                            source_base = lr;
                        }
                    }
                }
                for (std::deque<ListDigraph::Node>::iterator nr = same_node_r[ai[left_node].edge_specifier.node_index].begin(); nr != same_node_r[ai[left_node].edge_specifier.node_index].end(); ++nr) {
                    for (ListDigraph::OutArcIt oi(g, *nr); oi != INVALID ;++oi) {
                        ListDigraph::OutArcIt rl(g, g.target(oi));

                        if (marked_drain[oi] && rl != INVALID && meta->exons[ai[left_node].edge_specifier.node_index].right + 1 == meta->exons[ai[rl].edge_specifier.node_index].left) {
                            has_drain = true;
                            drain = oi;
                            drain_base = rl;
                        }
                    }
                }
                
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Broken Intron Test " + std::to_string(g.id(a)) +" s " + std::to_string(g.id(source)) +" d "+ std::to_string(g.id(drain)) + "\n");
                #endif
                unsigned int vote_left_count = 0;
                unsigned int vote_right_count = 0;
                for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin(); ssi != fs[a].series.end(); ++ssi) {
                    
                    int id = ssi->first;
                    if (has_drain && has_source && fs[source].series[id].capacity < ssi->second.capacity && fs[drain].series[id].capacity < ssi->second.capacity) {
                    // likely we can kill both

                        rcount source_cap = fs[source].series[id].capacity;
                        rcount drain_cap = fs[drain].series[id].capacity;

                        rcount source_cap_base = fs[source_base].series[id].mean.mean;
                        rcount drain_cap_base = fs[drain_base].series[id].mean.mean;
                        rcount min = std::min(source_cap_base, drain_cap_base);
                        rcount max = std::max(source_cap_base, drain_cap_base);

                        if (source_cap_base > 30 && drain_cap_base > 30) {
                            continue;
                        }

                        if ( (source_cap_base == drain_cap_base || (max - min)*100/(max) < 30 || (max < 15 && min > 5)) && ssi->second.capacity !=0 && source_cap_base *100 / ssi->second.capacity < options::Instance()->get_broken_intron_retention_threshold() 
                                && drain_cap_base*100 / ssi->second.capacity < options::Instance()->get_broken_intron_retention_threshold()) {

                            ++vote_right_count;
                            ++vote_left_count;

                        } else if (source_cap_base > drain_cap_base && drain_cap_base*100 / ssi->second.capacity < options::Instance()->get_broken_intron_retention_threshold()) {
                            ++vote_right_count;
                        } else if (source_cap_base < drain_cap_base && source_cap_base *100 / ssi->second.capacity < options::Instance()->get_broken_intron_retention_threshold()) {
                            ++vote_left_count;
                        }
                    }
                }
                
                if (options::Instance()->vote(vote_left_count, input_count, options::delete_on::group_vote_low)) {
                    kill_status[source] = true;
                }
                if (options::Instance()->vote(vote_right_count, input_count, options::delete_on::group_vote_low)) {
                    kill_status[drain] = true;
                }
            }
        }
    }
    
   for (ListDigraph::ArcIt a(g); a != INVALID; ) {
        if (kill_status[a] == true && guided_saves[a] == false) {
            ListDigraph::Arc c(a);
            ++a;
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Erase Intron " + std::to_string(g.id(c)) +"\n");
            #endif
            g.erase(c);
        } else {
            ++a;
        }
    }
}

void base_manager::remove_too_short_source_drain(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Pre-Filter\n");
    #endif
    
    for (ListDigraph::OutArcIt o(g,s); o != INVALID; ) {
        
        ListDigraph::Arc arc(o);
        ++o;
        
        // first node this belongs to!
        ListDigraph::OutArcIt no(g, g.target(arc)); // first node
        
        if (no == INVALID || guided_saves[no] || guided_saves[arc] || !consecutive_s_t[no] || !marked_source[no]) {
            continue;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Test " + std::to_string(g.id(no)) + "\n");
        #endif
        
        unsigned int vote_count = input_count - fs[no].series.size();
        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[no].series.begin(); ssi != fs[no].series.end(); ++ssi) {
            if (ssi->second.length_to_first_zero_right <= 20 || ssi->second.region_average < 2) {
                ++vote_count;
            }    
        }   
        if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_high)) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Erase Short Source " + std::to_string(g.id(arc)) + ".\n");
            #endif
            g.erase(arc);
            g.erase(no);
        }
    }
    
    for (ListDigraph::InArcIt i(g,t); i != INVALID; ) {
        
        ListDigraph::Arc arc(i);
        ++i;
        
        ListDigraph::InArcIt ni(g, g.source(arc));
        
        if (ni == INVALID || guided_saves[ni] || guided_saves[arc] || !consecutive_s_t[ni] || !marked_drain[ni]) { // no need to test, keep it!
            continue;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Test " + std::to_string(g.id(ni)) + "\n");
        #endif
        
        unsigned int vote_count = input_count - fs[ni].series.size();;
        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[ni].series.begin(); ssi != fs[ni].series.end(); ++ssi) {
            if (ssi->second.length_to_first_zero_left <= 20 || ssi->second.region_average < 2) {
                ++vote_count;
            }    
        }   
        if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_high)) {
            #ifdef ALLOW_DEBUG
	    logger::Instance()->debug("Erase Short Drain " + std::to_string(g.id(arc)) + ".\n");
            #endif
            g.erase(arc);
            g.erase(ni);  
        }
    }
    
}
  

bool base_manager::remove_small_low_st(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t) {

        #ifdef ALLOW_DEBUG
    logger::Instance()->debug("small low st\n");
    #endif
    
    bool res = false;
    for (ListDigraph::ArcIt a(g); a != INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        
        if (ai[arc].edge_type == edge_types::NODE) {

            bool left_neighbour = false;
            bool left_source = false;
            bool left_other = false;
            for(ListDigraph::InArcIt ei(g, g.source(arc)); ei != INVALID; ++ei) {
               if (ai[ei].edge_specifier.right_consecutive) {
                   left_neighbour = true;
               } else if (ai[ei].edge_type == edge_types::HELPER) {
                   left_source = true;
               }
            }
           
            bool right_neighbour = false;
            bool right_drain = false;
            for(ListDigraph::OutArcIt eo(g, g.target(arc)); eo != INVALID; ++eo) {
                if (ai[eo].edge_specifier.left_consecutive ) {
                    right_neighbour = true;
                } else if (ai[eo].edge_type == edge_types::HELPER) {
                   right_drain = true;
                }
            }
           
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test Arc " +  std::to_string(g.id(arc))  + " " + std::to_string(left_source) + " " + std::to_string(right_drain) + "\n");
            #endif

            unsigned int vote_count = input_count - fs[arc].series.size();
            for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[arc].series.begin(); ssi != fs[arc].series.end(); ++ssi) {
                            
                rpos length = ssi->second.length;
                float average = ssi->second.region_average;

                if ( (left_source || right_drain) && (!left_neighbour || !right_neighbour) && (length <= 15 || average <= 1.5) && !guided_saves[arc]) {
                     ++vote_count;
                }
            }
            
            if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_high)) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Erase Empty " + std::to_string(g.id(arc)) + ".\n");
                #endif
                
                g.erase(arc);
                ListDigraph::InArcIt ei(g, g.source(arc));
                if (ei != INVALID) g.erase(ei);
                ListDigraph::OutArcIt eo(g, g.target(arc));
                if (eo != INVALID) g.erase(eo);
                
                res = true;
            } 
        }
    }
    return res;
}
    

bool base_manager::erase_low_deviation_st(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t) {

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("low dev st\n");
    #endif

    std::vector<std::deque<ListDigraph::Node> > same_node_r(meta->size);
    std::vector<std::deque<ListDigraph::Node> > same_node_l(meta->size);
    std::vector<std::deque<ListDigraph::Arc> > same_node(meta->size);
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::NODE) { //with opened middle node
            
            int ni = ai[a].edge_specifier.node_index;
            same_node_r[ni].push_back(g.target(a));
            same_node_l[ni].push_back(g.source(a));
            same_node[ni].push_back(a);
            
        }
    }
    
    bool res = false;
    for (ListDigraph::ArcIt a(g); a != INVALID;) {
        
        ListDigraph::Arc arc(a);
        ++a;
        
        if (guided_saves[arc]) {
            continue;
        }
        
        if (ai[arc].edge_type == edge_types::NODE ) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Check Dev"  + std::to_string(g.id(arc)) + "\n");
            #endif
            
            
            ListDigraph::InArcIt arc_left(g, g.source(arc));
            ListDigraph::OutArcIt arc_right(g, g.target(arc));
            
            if (arc_left == INVALID || arc_right == INVALID) {
                continue;
            }
            ListDigraph::InArcIt cl(arc_left);
            ++cl;
            ListDigraph::OutArcIt cr(arc_right);
            ++cr;
            
            if (cl != INVALID || cr != INVALID) {
                continue;
            }
            
            if (ai[arc_left].edge_type != edge_types::HELPER && ai[arc_right].edge_type != edge_types::HELPER) {
                continue;
            }
            
            if (ai[arc_left].edge_type == edge_types::HELPER) {
                // we test if this is unique!
                unsigned int counter = 0;
                for (ListDigraph::InArcIt oi(g, g.target(arc_right)); oi != INVALID ;++oi) {
                    ++counter;
                }
                if (counter == 1) {
                    continue;
                }
            }
            
            if (ai[arc_right].edge_type == edge_types::HELPER) {
                // we test if this is unique!
                unsigned int counter = 0;
                for (ListDigraph::OutArcIt oi(g, g.source(arc_left)); oi != INVALID ;++oi) {
                    ++counter;
                }
                if (counter == 1) {
                    continue;
                }
            }
               
            unsigned int vote_count = input_count - fs[arc].series.size();
            for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[arc].series.begin(); ssi != fs[arc].series.end(); ++ssi) {
                if (ssi->second.deviation <= 0.01) {
                    ++vote_count;
                }
            }

            if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_high)) {
                res = true;

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Erase by Dev " + std::to_string(g.id(arc)) + ".\n");
                #endif
                
                g.erase(arc);
                g.erase(arc_left);
                g.erase(arc_right);
            }
        }
    }
    return res;
}
    

bool base_manager::erase_low_boundary_st(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t) {

    std::vector<std::deque<ListDigraph::Node> > same_node_r(meta->size);
    std::vector<std::deque<ListDigraph::Node> > same_node_l(meta->size);
    std::vector<std::deque<ListDigraph::Arc> > same_node(meta->size);
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::NODE) { //with opened middle node
            
            int ni = ai[a].edge_specifier.node_index;
            same_node_r[ni].push_back(g.target(a));
            same_node_l[ni].push_back(g.source(a));
            same_node[ni].push_back(a);
            
        }
    }
        
    bool res = false;
    for (ListDigraph::ArcIt a(g); a != INVALID;) {
        
        ListDigraph::Arc arc(a);
        ++a;
        
        if (guided_saves[arc]) {
            continue;
        }
        
        if (ai[arc].edge_type == edge_types::EXON && g.source(arc) != s && g.target(arc) != t) {
            
            // we go over every node an check for possible removal of arcs

            ListDigraph::InArcIt node_left(g, g.source(arc));
            ListDigraph::OutArcIt node_right(g, g.target(arc));
            
            if (node_left == INVALID || node_right == INVALID) {
                    continue;
            }
            
            unsigned int out_degree = 0;
            bool drain_exists = false;
            for (std::deque<ListDigraph::Node>::iterator nr = same_node_r[ai[node_left].edge_specifier.node_index].begin(); nr != same_node_r[ai[node_left].edge_specifier.node_index].end(); ++nr) {
                for (ListDigraph::OutArcIt oi(g, *nr); oi != INVALID ;++oi) {
                    if (ai[oi].edge_type == edge_types::HELPER) {
                        drain_exists = true;
                        continue;
                    }
                    ++out_degree;
                }
            }
            
            unsigned int in_degree = 0;
            bool source_exists = false;
            for (std::deque<ListDigraph::Node>::iterator nl = same_node_l[ai[node_right].edge_specifier.node_index].begin(); nl != same_node_l[ai[node_right].edge_specifier.node_index].end(); ++nl) {
                for (ListDigraph::InArcIt oi(g, *nl); oi != INVALID ;++oi) {
                    if (ai[oi].edge_type == edge_types::HELPER) {
                        source_exists = true;
                        continue;
                    }
                    ++in_degree;
                }
            }


            if (meta->exons[ai[node_left].edge_specifier.node_index].right + 1 == meta->exons[ai[node_right].edge_specifier.node_index].left) {
                continue; // disregard neighbouring here!
            }
            
            unsigned int vote_count = 0;
            for (std::set<int>::iterator idi = input_ids.begin(); idi != input_ids.end(); ++idi) {
                
                gfmap<int, flow_series::series_struct>::iterator ssi = fs[arc].series.find(*idi);
                
                int id = *idi;

                rcount left = 0;
                for (std::deque<ListDigraph::Arc>::iterator ni = same_node[ai[node_left].edge_specifier.node_index].begin(); ni != same_node[ai[node_left].edge_specifier.node_index].end(); ++ni) {
                    left += fs[*ni].series[id].region_average;
                }

                rcount right = 0;
                for (std::deque<ListDigraph::Arc>::iterator ni = same_node[ai[node_right].edge_specifier.node_index].begin(); ni != same_node[ai[node_right].edge_specifier.node_index].end(); ++ni) {
                    right += fs[*ni].series[id].region_average;
                }

                if ( ssi != fs[arc].series.end() ) {
                    
                    rcount we = ssi->second.capacity;
                    if ((out_degree == 1 && left>= 10.0 * we * we + 10.0) || (in_degree == 1 && right>= 10.0 * we * we + 10.0)) {
                        ++vote_count;
                    }
                    
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Test N Scallop " + std::to_string(g.id(arc)) + " : " + std::to_string(we) + " " + std::to_string(left) + " " + std::to_string(right) + " : " + std::to_string(meta->exons[ai[node_left].edge_specifier.node_index].right) + " " + std::to_string(meta->exons[ai[node_right].edge_specifier.node_index].left)  + ".\n");
                    #endif
                    
                } else if ((out_degree == 1 && left>= 10.0) || (in_degree == 1 && right>= 10.0)) {

                    ++vote_count;
                    
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Test N Scallop " + std::to_string(g.id(arc)) + " empty vote.\n");
                    #endif
                }
                
            }
            if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_high)) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Erase by boundary " + std::to_string(g.id(arc)) + ".\n");
                #endif

                res = true;

                g.erase(arc);
                if (!drain_exists && out_degree == 1) {
                    ListDigraph::Arc narc = g.addArc(g.target(node_left), t);
                    ai[narc].edge_type = edge_types::HELPER;
                    initialize_source_drain_arc(narc);
                }
                if (!source_exists && in_degree == 1)  {
                    ListDigraph::Arc narc = g.addArc(s, g.source(node_right));
                    ai[narc].edge_type = edge_types::HELPER;
                    initialize_source_drain_arc(narc);
                }
            }
        } 
    } 
    return res;
}


bool base_manager::prune_small_junctions(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Small Junctions\n");
    #endif
    
    std::vector<std::deque<ListDigraph::Arc> > same_node(meta->size);
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::NODE) { //with opened middle node
            int ni = ai[a].edge_specifier.node_index;
            same_node[ni].push_back(a);  
        }
    }
    
    bool res = false;
    for (ListDigraph::ArcIt a(g); a != INVALID;) {
        
        ListDigraph::Arc arc(a);
        ++a;
        
        if (ai[arc].edge_type == edge_types::NODE) {
            
            std::map<int, unsigned int> vote_counts; // if current has no coverage, no evidence for removal of other edges
            for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[arc].series.begin(); ssi != fs[arc].series.end(); ++ssi) {
                int id = ssi->first;
                            
                capacity_type current = ssi->second.region_average;
                capacity_type max_in = 0;

                for (std::deque<ListDigraph::Arc>::iterator ni = same_node[ai[arc].edge_specifier.node_index].begin(); ni != same_node[ai[arc].edge_specifier.node_index].end(); ++ni) {
                    for (ListDigraph::InArcIt oi(g, g.source(*ni)); oi != INVALID ;++oi) {

                        ListDigraph::InArcIt node_left(g, g.source(oi));
                        if (node_left == INVALID) continue;

                        if (!ai[oi].edge_specifier.right_consecutive) continue;

                        capacity_type nla = fs[node_left].series[id].region_average;
                        capacity_type wl = std::max(nla , (capacity_type) fs[oi].series[id].capacity);

                        if (wl > max_in) {
                            max_in = wl;
                        }
                    }
                }
                for (std::deque<ListDigraph::Arc>::iterator ni = same_node[ai[arc].edge_specifier.node_index].begin(); ni != same_node[ai[arc].edge_specifier.node_index].end(); ++ni) {
                    for (ListDigraph::InArcIt oi(g, g.source(*ni)); oi != INVALID ;) {

                        ListDigraph::Arc oik(oi);
                        ++oi;

                        if (ai[oik].edge_specifier.right_consecutive) continue;

                        capacity_type w = fs[oik].series[id].capacity;

                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Test scallop L " + std::to_string(g.id(oik)) + " to " + std::to_string(g.id(arc)) + " : " + std::to_string(w) + " " + std::to_string(max_in) + " " + std::to_string(current) + ".\n");
                        #endif

                         if ( max_in >= 2.0 * w * w + 18.0 && current >= 2.0 * w * w + 18.0 && !guided_saves[oik]) {
                            ++vote_counts[g.id(oik)];
                        }

                    }
                }

                capacity_type max_out = 0;
                for (std::deque<ListDigraph::Arc>::iterator ni = same_node[ai[arc].edge_specifier.node_index].begin(); ni != same_node[ai[arc].edge_specifier.node_index].end(); ++ni) {
                    for (ListDigraph::OutArcIt oi(g, g.target(*ni)); oi != INVALID ;++oi) {

                        ListDigraph::OutArcIt node_right(g, g.target(oi));
                        if (node_right == INVALID) continue;
                        if (!ai[oi].edge_specifier.left_consecutive) continue;

                        capacity_type nra = fs[node_right].series[id].region_average;
                        capacity_type wr = std::max(nra , (capacity_type) fs[oi].series[id].capacity);
                        if (wr > max_out) {
                            max_out = wr;
                        }
                    }
                }

                for (std::deque<ListDigraph::Arc>::iterator ni = same_node[ai[arc].edge_specifier.node_index].begin(); ni != same_node[ai[arc].edge_specifier.node_index].end(); ++ni) {
                    for (ListDigraph::OutArcIt oi(g, g.target(*ni)); oi != INVALID ;) {

                        ListDigraph::Arc oik(oi);
                        ++oi;

                        if (ai[oik].edge_specifier.left_consecutive) continue;

                        capacity_type w = fs[oik].series[id].capacity;

                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Test scallop R " + std::to_string(g.id(oik)) + " to " + std::to_string(g.id(arc)) + " : " + std::to_string(w) + " " + std::to_string(max_out) + " " + std::to_string(current) + ".\n");
                        #endif

                        if ( max_out >= 2.0 * w * w + 18.0 && current >= 2.0 * w * w + 18.0 && !guided_saves[oik]) {
                            ++vote_counts[g.id(oik)];
                        }
                    }
                }
            }
            
            for(std::map<int, unsigned int>::iterator vi = vote_counts.begin(); vi != vote_counts.end(); ++vi) {
                if (options::Instance()->vote(vi->second, input_count, options::delete_on::group_vote_high)) {
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Kill scallop " + std::to_string(vi->first) + ".\n");
                    #endif
                    g.erase(g.arcFromId(vi->first));

                    res = true;
                }
            }
        } 
    }
    return res;
}
    
bool base_manager::threshold_filter(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t, float low_mark) {
    
    #ifdef ALLOW_DEBUG
	logger::Instance()->debug("Threshold " + std::to_string(low_mark) + ".\n");
    #endif
    
    std::vector<std::deque<ListDigraph::Arc> > same_node(meta->size);

    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::NODE) { //with opened middle node
            int ni = ai[a].edge_specifier.node_index;
            same_node[ni].push_back(a);  
        }
    }
    
    ListDigraph::ArcMap<bool> filtered_sd(g);
    ListDigraph::ArcMap<bool> kill_status(g); // potentially killed edges!
    for (ListDigraph::OutArcIt o(g,s); o != INVALID; ) {
        
        ListDigraph::Arc arc(o);
        ++o;
        
        // first node this belongs to!
        ListDigraph::OutArcIt no(g, g.target(arc));
        if (no == INVALID) continue;
        ListDigraph::OutArcIt no2(g, g.target(no));
        
        if (no2 == INVALID || guided_saves[no] || !marked_source[no]) {
            continue;
        }
        
        #ifdef ALLOW_DEBUG
	logger::Instance()->debug("Test low Source " + std::to_string(g.id(no)) + ".\n");
        #endif
        
        unsigned int vote_count = input_count - fs[no].series.size();
        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[no].series.begin(); ssi != fs[no].series.end(); ++ssi) {
            int id = ssi->first;  
            if (ssi->second.average_to_first_zero_from_right <= low_mark || ssi->second.length <= 20 ) {
                //erase
                ++vote_count;
            }
        }
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Vote Count S " + std::to_string(g.id(arc)) + " : " + std::to_string(vote_count) + "\n");
        #endif
        if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_low)) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Mark Short Source " + std::to_string(g.id(arc)) + ".\n");
            #endif
            kill_status[no] = true;
            filtered_sd[no] = true;
            for (ListDigraph::OutArcIt no3(g, g.target(no)); no3 != INVALID; ++no3) {
                kill_status[no3] = true;
                filtered_sd[no3] = true;
                #ifdef ALLOW_DEBUG
	        logger::Instance()->debug("Mark " + std::to_string(g.id(no3)) + ".\n");
                #endif
            } 
        }
    }
    for (ListDigraph::InArcIt i(g,t); i != INVALID; ) {
        
        ListDigraph::Arc arc(i);
        ++i;
        
        ListDigraph::InArcIt ni(g, g.source(arc));
        if (ni == INVALID) continue; 
        ListDigraph::InArcIt ni2(g, g.source(ni));
        
        if (ni2 == INVALID || guided_saves[ni] || !marked_drain[ni]) { // no need to test, keep it!
            continue;
        }
        
        #ifdef ALLOW_DEBUG
	logger::Instance()->debug("Test low Drain " + std::to_string(g.id(ni)) + ".\n");
        #endif
        unsigned int vote_count = input_count - fs[ni].series.size();
        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[ni].series.begin(); ssi != fs[ni].series.end(); ++ssi) {
            int id = ssi->first;  
            if (ssi->second.average_to_first_zero_from_left <= low_mark  || ssi->second.length <= 20 ) {
                //erase
                ++vote_count;
            }
        }
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Vote Count D " + std::to_string(g.id(arc)) + " : " + std::to_string(vote_count) + "\n");
        #endif
        if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_low)) {
            #ifdef ALLOW_DEBUG
	    logger::Instance()->debug("Mark Short Drain " + std::to_string(g.id(arc)) + ".\n");
            #endif
            kill_status[ni] = true;
            filtered_sd[ni] = true;
            for (ListDigraph::InArcIt ni3(g, g.source(ni)); ni3 != INVALID; ++ni3) {
                kill_status[ni3] = true;
                filtered_sd[ni3] = true;
                #ifdef ALLOW_DEBUG
	        logger::Instance()->debug("Mark " + std::to_string(g.id(ni3)) + ".\n");
                #endif
            }  
        }
    }    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (ai[a].edge_type == edge_types::EXON) {
            unsigned int vote_count = input_count - fs[a].series.size();
            for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin(); ssi != fs[a].series.end(); ++ssi) {
                if (ssi->second.capacity <= low_mark) {
                    ++vote_count;
                }
            }
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Vote Count E " + std::to_string(g.id(a)) + " : " + std::to_string(vote_count) + "\n");
            #endif
            if (options::Instance()->vote(vote_count, input_count, options::delete_on::group_vote_high)) {
                kill_status[a] = true;
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Mark Exon by threshold " + std::to_string(g.id(a)) + ".\n");
                 #endif
            }
        }
    }


    bool node_leftover = false;
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        if (ai[a].edge_type == edge_types::NODE) { 
      
            bool has_in_arc = false;
            bool has_in_helper = false;
            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                if (ai[i].edge_type == edge_types::HELPER) {
                    has_in_helper = true;
                } else if (kill_status[i] != true) {
                    has_in_arc = true;
                }
            }

            bool has_out_arc = false;
            bool has_out_helper = false;
            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                if (ai[o].edge_type == edge_types::HELPER) {
                    has_out_helper = true;
                } else if(kill_status[o] != true) {
                    has_out_arc = true;
                }
            }

            bool has_in = has_in_arc || has_in_helper;
            bool has_out = has_out_arc || has_out_helper;

            if (!has_in_arc && !has_out_arc) {
                if (has_in_helper || has_out_helper ) {
                    // keep it!
                } else {
                    kill_status[a] = true;
                     #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Mark Node by threshold " + std::to_string(g.id(a)) + ".\n");
                    #endif
                }
            } else {
                node_leftover = true;
            }
        }
    }
    if (!node_leftover) {
        std::map<int, int> max_vote;        
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
            if (ai[a].edge_type == edge_types::NODE) {
                for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[a].series.begin(); ssi != fs[a].series.end(); ++ssi) {
                    int id = ssi->first;
                    if (max_vote.find(id) == max_vote.end() || ssi->second.capacity > fs[g.arcFromId(max_vote[id])].series[id].capacity) {
                        max_vote[id] = g.id(a);
                    }
                }
            }
        }
        if ( !max_vote.empty() ) {
            kill_status[g.arcFromId(most_common_mapped(max_vote.begin(), max_vote.end()))] = false;
        }
    }

    // rescue mission!    
    while (true) {
        bool changed = false;
        ListDigraph::ArcMap<bool> next_rescue(g); // mark all rescues here!
                
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
            
            if (kill_status[a] || ai[a].edge_type == edge_types::HELPER) {
                continue;
            }
            
            bool has_in_arc = false;
            bool has_in_helper = false;
            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                if (ai[i].edge_type == edge_types::HELPER) {
                    has_in_helper = true;
                } else if (kill_status[i] != true) {
                    has_in_arc = true;
                }
            }

            bool has_out_arc = false;
            bool has_out_helper = false;
            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                if (ai[o].edge_type == edge_types::HELPER) {
                    has_out_helper = true;
                } else if(kill_status[o] != true) {
                    has_out_arc = true;
                }
            }
            bool has_in = has_in_arc || has_in_helper;
            bool has_out = has_out_arc || has_out_helper; 
            
            bool has_non_sd_filtered_in = false;
            std::map<int, float> max_count_in;
            bool has_non_sd_filtered_out = false;
            std::map<int, float> max_count_out;
            if (ai[a].edge_type == edge_types::NODE) {
                for (std::deque<ListDigraph::Arc>::iterator ni = same_node[ai[a].edge_specifier.node_index].begin(); ni != same_node[ai[a].edge_specifier.node_index].end(); ++ni) {
                    for (ListDigraph::InArcIt i(g, g.source(*ni)); i!=INVALID; ++i) {
                        if (!filtered_sd[i]) {
                            has_non_sd_filtered_in = true;
                        }
                        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[i].series.begin(); ssi != fs[i].series.end(); ++ssi) {
                            int id = ssi->first;
                            float cap = ssi->second.mean.mean;
                            if (max_count_in[id] < cap) {
                                max_count_in[id] = cap;
                            }
                        }
                    }
                    for (ListDigraph::OutArcIt o(g, g.target(*ni)); o!=INVALID; ++o) {
                        if (!filtered_sd[o]) {
                            has_non_sd_filtered_out = true;
                        }
                        for(gfmap<int, flow_series::series_struct>::iterator ssi = fs[o].series.begin(); ssi != fs[o].series.end(); ++ssi) {
                            int id = ssi->first;
                            float cap = ssi->second.mean.mean;
                            if (max_count_out[id] < cap) {
                                max_count_out[id] = cap;
                            }
                        }
                    }
                }
            }
            
            if (!has_in) {
                std::deque<int> vote;
                std::deque<int> secondary_vote;
                for (std::set<int>::iterator idi = input_ids.begin(); idi != input_ids.end(); ++idi) {
                    bool has_max = false;
                    float max_count = 0;                
                    std::deque<ListDigraph::Arc> max_arc_all;
                    for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                        float cap = fs[i].series[*idi].mean.mean;

                        if (cap > max_count || !has_max) {
                            max_arc_all.clear();
                            max_arc_all.push_back(i);
                            max_count = cap;
                            has_max = true;
                        } else if (cap == max_count) {
                            max_arc_all.push_back(i);
                        }  
                    }
                    for (std::deque<ListDigraph::Arc>::iterator it = max_arc_all.begin(); it != max_arc_all.end(); ++it) {
                        if (!filtered_sd[*it] || !has_non_sd_filtered_in || max_count_in[*idi] == max_count) { // arc is not low drain or source or only filtered ones exist
                            if (max_count > 4) {
                                vote.push_back(g.id(*it));
                            }
                            secondary_vote.push_back(g.id(*it));
                        }
                    }
                }
                if (vote.empty()) { // so all are below 4
                    vote = secondary_vote;
                } 
                if (!vote.empty()) {
                    std::set<int> vr;
                    most_common_multi(vote.begin(), vote.end(), vr);
                    for (std::set<int>::iterator ivr = vr.begin(); ivr != vr.end(); ++ivr) {
                        next_rescue[g.arcFromId(*ivr)] = true;
                        changed = true;       
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Rescue In " + std::to_string(g.id(a)) + " " + std::to_string(*ivr) +".\n"); 
                        #endif
                    }
                }
            }
            

            if (!has_out) {
                std::deque<int> vote;
                std::deque<int> secondary_vote;
                for (std::set<int>::iterator idi = input_ids.begin(); idi != input_ids.end(); ++idi) {                
                    bool has_max = false;
                    float max_count = 0;                
                    std::deque<ListDigraph::Arc> max_arc_all;

                    for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {

                        float cap = fs[o].series[*idi].mean.mean;

                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Cap " + std::to_string(g.id(o)) + " " + std::to_string(cap) +".\n");
                        #endif

                        if (cap > max_count || !has_max) {
                            max_arc_all.clear();
                            max_arc_all.push_back(o);
                            max_count = cap;
                            has_max = true;
                        } else if (cap == max_count) {
                            max_arc_all.push_back(o);
                        }  
                    }
                    for (std::deque<ListDigraph::Arc>::iterator it = max_arc_all.begin(); it != max_arc_all.end(); ++it) {
                         if (!filtered_sd[*it] || !has_non_sd_filtered_out || max_count_out[*idi] == max_count) { // arc is not low drain or source or only filtered ones exist
                             if (max_count > 4) {
                                vote.push_back(g.id(*it));
                            }
                            secondary_vote.push_back(g.id(*it));
                         }
                    }
                }
                if (vote.empty()) { // so all are below 4 
                    vote = secondary_vote;
                }
                if (!vote.empty()) {
                    std::set<int> vr;
                    most_common_multi(vote.begin(), vote.end(), vr);
                    for (std::set<int>::iterator ivr = vr.begin(); ivr != vr.end(); ++ivr) {
                        next_rescue[g.arcFromId(*ivr)] = true;
                        changed = true;       
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Rescue Out " + std::to_string(g.id(a)) + " " + std::to_string(*ivr) +".\n"); 
                        #endif
                    }
                }
            }
        }
        if (!changed) {
            break;
        }
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
            kill_status[a] = (kill_status[a] && !next_rescue[a]);
        }
    }            
         
    bool changed = false;
    for (ListDigraph::ArcIt a(g); a != INVALID; ) {
        if (kill_status[a] == true && guided_saves[a] == false) {
            ListDigraph::Arc c(a);
            ++a;
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Erase " + std::to_string(g.id(c)) +"\n");
            #endif
            g.erase(c);
            changed = true;
        } else {
            ++a;
        }
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Initial Erased Done \n");
    #endif
    return changed;
}

void base_manager::remove_dead_ends() {
    
   ListDigraph::ArcIt a(g);
   if(a == INVALID) {
       return;
   }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Remove Dead End \n");
    #endif
    
    // create a topological ordering
    std::deque<ListDigraph::Node> top_order;
    pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > visitor(top_order);
    DfsVisit<ListDigraph, pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > > dfs(g, visitor);
    dfs.init();
    dfs.addSource(s);
    dfs.start();

    
    for (std::deque<ListDigraph::Node>::reverse_iterator n = top_order.rbegin(); n!= top_order.rend(); ) {
        ListDigraph::OutArcIt o(g, *n);
        if (*n != t && o == INVALID) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Erase at node " + std::to_string(g.id(*n)) + " \n");
            #endif
            
            for (ListDigraph::InArcIt a(g, *n);a!=INVALID;) {
                ListDigraph::Arc sa(a);
                ++a;
                g.erase(sa);
            }
            ListDigraph::Node sn(*n);
            ++n;
            g.erase(sn);
            continue;
        }
        ++n;
    }
    
    top_order.clear();
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
        if (!dfs.reached(n)) {
            dfs.addSource(n);
            dfs.start();
        }
    }

    for (std::deque<ListDigraph::Node>::iterator n = top_order.begin(); n!= top_order.end(); ) {
        ListDigraph::InArcIt i(g, *n);
        if (*n != s && i == INVALID) {
            
            #ifdef ALLOW_DEBUG
             logger::Instance()->debug("Erase at node " + std::to_string(g.id(*n)) + " \n");
            #endif
            
            for (ListDigraph::OutArcIt a(g, *n);a!=INVALID; ) {
                ListDigraph::Arc sa(a);
                ++a;
                g.erase(sa);
            }
            ListDigraph::Node sn(*n);
            ++n;
            g.erase(sn);
            continue;
        }
        ++n;
    }
    
    // fully erased
}


void base_manager::denoise_graph(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, int id) {
    
    // first we get a topological sorting of all nodes
    std::deque<ListDigraph::Node> top_order;
    pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > visitor(top_order);
    DfsVisit<ListDigraph, pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > > dfs(g, visitor);
    dfs.init();
    dfs.addSource(s);
    dfs.start();

    
//    ListDigraph::ArcMap<rcount> average_push_fwd(g); // init as 0, we set a value
//    ListDigraph::ArcMap<rcount> average_push_bwd(g);
    
    ListDigraph::ArcMap<rcount> average_push(g); // init as 0, we set a value

    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) { // init as the average!
        
        average_push[a]= fs[a].series[id].mean.mean;
         
        if (marked_source[a]) {
            average_push[a] = rec_score_backward(a, id);
        }
        if (marked_drain[a]) {
            average_push[a] = rec_score_forward(a, id);
        }
    }
   
    // forward potential!
    ListDigraph::ArcMap<rcount> cp_fwd(g); // init as 0
    push_potential_forward(cp_fwd, top_order, average_push, id);

    // backward potential!
    ListDigraph::ArcMap<rcount> cp_bw(g); // init as 0
    push_potential_backward(cp_bw, top_order, average_push, id);

    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("add " + std::to_string(g.id(a)) + " " + std::to_string(cp_fwd[a]) + " " + std::to_string(cp_bw[a]) + "\n");
        #endif
        
        capacity_type max = std::max(cp_fwd[a], cp_bw[a]);
        if (ai[a].edge_type != edge_types::HELPER) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Correct Capacity " + std::to_string(g.id(a)) + " from " + std::to_string(fs[a].series[id].capacity) + " to " + std::to_string(fs[a].series[id].capacity+max) + ".\n");
            #endif
            fs[a].series[id].capacity += max;
        }
    }
}

float base_manager::rec_score_forward(ListDigraph::Arc a, int id) {
    
    unsigned int count = 0;
    for (ListDigraph::InArcIt i(g, g.source(a)); i != INVALID ;++i) {
         ++count;
    }
    if (count > 1) {
        return 0;
    }
    
    float score_left = fs[a].series[id].mean.compute_score();
    
    float score_right = 0;
    for (ListDigraph::OutArcIt i(g, g.target(a)); i != INVALID ;++i) {
           score_right+=rec_score_forward(i, id);
    }
    if (score_right == 0) {
        return score_left;
    }
    
    return sqrt(score_left * score_right);
}


float base_manager::rec_score_backward(ListDigraph::Arc a, int id) {
    
    unsigned int count = 0;
    for (ListDigraph::OutArcIt i(g, g.target(a)); i != INVALID ;++i) {
         ++count;
    }
    if (count > 1) {
        return 0;
    }
    
    float score_right = fs[a].series[id].mean.compute_score();
    
    float score_left = 0;
    for (ListDigraph::InArcIt i(g, g.source(a)); i != INVALID ;++i) {
           score_left+=rec_score_backward(i, id);
    }
    if (score_left == 0) {
        return score_right;
    }
    
    return sqrt(score_left * score_right);
}

void base_manager::push_potential_forward(ListDigraph::ArcMap<rcount> &cp_fwd, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, int id) {
        
     ListDigraph::ArcMap<rcount> cs(g);  

     for (std::deque<ListDigraph::Node>::iterator n = top_order.begin(); n!= top_order.end(); ++n) {
        
//          logger::Instance()->debug("Push_Node " + std::to_string(g.id(*n)) + "\n");
                  
        // first test if we can propagate
        ListDigraph::InArcIt innt(g, *n);
        bool node_left= (innt != INVALID && ai[innt].edge_type == edge_types::NODE);
        if (node_left) {
            
            // this means we can only have one valid left, to of course multiple outs edges
            ListDigraph::Arc node(innt);
            
            rcount in_push = 0;
            for (ListDigraph::InArcIt inpi(g, g.source(node)); inpi != INVALID ;++inpi) {
                in_push += cp_fwd[inpi];
            }
            
            rcount left = 0;
            for (ListDigraph::InArcIt inpi(g, g.source(node)); inpi != INVALID ;++inpi) {
                left += cs[inpi];
            }
                        
            // compute potentials
            rcount max = fs[node].series[id].region_max;
            rcount right = fs[node].series[id].right;
            left = std::min(fs[node].series[id].left, left);
            
            rcount left_potential = 0;
            if (max > left){
                left_potential = max - left;
            }
            
            rcount right_potential = 0;
            if (max > right){
                right_potential = max - right;
            }
            
//            logger::Instance()->debug("Pots " + std::to_string(g.id(*n)) + ": " + std::to_string(left_potential)  + " " + std::to_string(right_potential) + " PUSH " + std::to_string(in_push) + "\n");
            
            // find increases
            if (in_push > left_potential) {
                in_push -= left_potential;
            } else {
                in_push = 0;
            }
            // we set the push
            cp_fwd[node] = in_push; // the increase value
//            logger::Instance()->debug("SetNode " + std::to_string(g.id(node)) + " " + std::to_string(cp_fwd[node]) + "\n");
            
            right_potential += in_push; // we transport excess potential to the right
            
            if ( left < 15 && right < 15) {
                right_potential = in_push; // negate this
            }
            
            // we push the potential to the next adjourning edges
            rcount total_out = 0;
            bool helper = false;
            for (ListDigraph::OutArcIt oi(g, *n); oi != INVALID ;++oi) {
                total_out += average_push[oi];
                if (ai[oi].edge_type == edge_types::HELPER) {
                    helper = true;
                }
            }
            
            for (ListDigraph::OutArcIt oi(g, *n); oi != INVALID ;++oi) {
                if (total_out != 0) {
                    cp_fwd[oi] = right_potential * average_push[oi]/float(total_out); // upscaling according to own size
                    cs[oi] = right * average_push[oi]/float(total_out);
                } else {
                    cp_fwd[oi] = right_potential; // upscaling according to own size
                    cs[oi] = right;
                }

//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Set " + std::to_string(g.id(oi)) + " " + std::to_string(cp_fwd[oi]) + "\n");
//                #endif
            }
        }
    }
}


void base_manager::push_potential_backward(ListDigraph::ArcMap<rcount> &cp_bw, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, int id) {
   
     ListDigraph::ArcMap<rcount> cs(g);  

     for (std::deque<ListDigraph::Node>::reverse_iterator n = top_order.rbegin(); n!= top_order.rend(); ++n) {
         // first test if we can propagate
        ListDigraph::OutArcIt onnt(g, *n);
        bool node_right= (onnt != INVALID && ai[onnt].edge_type == edge_types::NODE);
        if (node_right) {
            
//            logger::Instance()->debug("Push_Node " + std::to_string(g.id(*n)) + "\n");
            
            // this means we can only have one valid left, to of course multiple outs edges
            ListDigraph::Arc node(onnt);
            
            rcount in_push = 0;
            for (ListDigraph::OutArcIt inpi(g, g.target(node)); inpi != INVALID ;++inpi) {
                in_push += cp_bw[inpi];
            }
            
            rcount right = 0;
            for (ListDigraph::OutArcIt inpi(g, g.target(node)); inpi != INVALID ;++inpi) {
                right += cs[inpi];
            }
            
            // compute potentials
            rcount max = fs[node].series[id].region_max;
            right = std::min(fs[node].series[id].right, right);
            rcount left = fs[node].series[id].left;
            
//            logger::Instance()->debug("BASE " + std::to_string(max) + " " + std::to_string(left) + " " +  std::to_string(right) + "  Inpush " + std::to_string(in_push) + "\n");
            
            rcount left_potential = 0;
            if (max > left){
                left_potential = max - left;
            }
            
            rcount right_potential = 0;
            if (max > right){
                right_potential = max - right;
            }
            
//           logger::Instance()->debug("Pots " + std::to_string(g.id(*n)) + ": " + std::to_string(left_potential)  + " " + std::to_string(right_potential) + " PUSH " + std::to_string(in_push) + "\n");
            
            // find increases
            if (in_push > right_potential) {
                in_push -= right_potential;
            } else {
                in_push = 0;
            }
            // we set the push
            cp_bw[node] = in_push; // the increase value
//            logger::Instance()->debug("SetNode " + std::to_string(g.id(node)) + " " + std::to_string(cp_bw[node]) + "\n");
            
            left_potential += in_push; // we transport excess potential to the right
            
            if ( left < 15 && right < 15) {
                left_potential = in_push; // negate this
            }
            
            // we push the potential to the next adjourning edges
            rcount total_out = 0;
            bool helper = false;
            for (ListDigraph::InArcIt oi(g, *n); oi != INVALID ;++oi) {
                total_out += average_push[oi];
                if (ai[oi].edge_type == edge_types::HELPER) {
                    helper = true;
                }
            }
            
            for (ListDigraph::InArcIt oi(g, *n); oi != INVALID ;++oi) {
                if (total_out != 0) {
                    cp_bw[oi] = left_potential * average_push[oi]/float(total_out); // upscaling according to own size
                    cs[oi] = left * average_push[oi]/float(total_out);
                } else {
                    cp_bw[oi] = left_potential; // upscaling according to own size
                    cs[oi] = left;                        
                }
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Set " + std::to_string(g.id(oi)) + " " + std::to_string(cp_bw[oi]) + "\n");
//                #endif
            }
        }
    }
}

void base_manager::denoise_graph_guided(std::deque<std::deque<ListDigraph::Arc> > &guides, int id, bool double_denoise) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Graph Denoise Guided \n");
    #endif
    
    ListDigraph::NodeMap<float> ratios(g);
    ListDigraph::ArcMap<rcount> guided_capacity(g);
    ListDigraph::ArcMap<float> guided_percentage(g);
    
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
        
        bool in_helper = false;
        bool out_helper = false;
		
        rcount in_push = 0;
        for (ListDigraph::InArcIt inpi(g, n); inpi != INVALID ;++inpi) {
            if (ai[inpi].edge_type == edge_types::HELPER) {
                in_helper = true;
            }
            in_push += fs[inpi].series[id].capacity;
        }
        
        rcount out_push = 0;
        for (ListDigraph::OutArcIt onpi(g, n); onpi != INVALID ;++onpi) {
            if (ai[onpi].edge_type == edge_types::HELPER) {
                out_helper = true;
            }
            out_push += fs[onpi].series[id].capacity;
        }

    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Graph Factor " + std::to_string(g.id(n)) + " : in " + std::to_string(in_push) + " out " + std::to_string(out_push)  + "\n");
    #endif

        if (out_push == 0 || in_push == 0) {
			ratios[n] = 1;
		} else if (in_helper || out_helper) {
            // we use the exon borders to estimate how much is actually lost by the helpers
            if (!double_denoise) {
			    if (out_helper) {
                    ListDigraph::InArcIt exon_arc(g, n); // there can only be one here as we are not compacted
                    out_push = fs[exon_arc].series[id].right;
                    if (out_push == 0) {
                         out_push = 1;
                    }
                //    logger::Instance()->debug("Correct Helper Out " + std::to_string(out_push)  + "\n");

			    } else {  // cannot be both
                    ListDigraph::OutArcIt exon_arc(g, n); // there can only be one here as we are not compacted
                    in_push = fs[exon_arc].series[id].left;
                    if (in_push == 0) {
                         in_push = 1;
                    }
                //    logger::Instance()->debug("Correct Helper In " + std::to_string(in_push)  + "\n");

			    }
            }
            ratios[n] = out_push/float(in_push);
        } else {
            ratios[n] = out_push/float(in_push);
        }
    }
    
    for (std::deque<std::deque<ListDigraph::Arc> >::iterator g_it = guides.begin(); g_it !=  guides.end(); ++g_it) {
        // for each guide
        
        float way_factor = 1;
        rcount guide_cap = 0;
        
        // loop first time
        std::deque<ListDigraph::Arc>::iterator a_it = g_it->begin();
        ++a_it;
        
        for (; a_it != g_it->end(); ++a_it) {
            
            if (ai[*a_it].edge_type == edge_types::HELPER) {
                continue;
            }
            
            float fac = ratios[g.source(*a_it)];
            way_factor = way_factor * fac;
            rcount local_limit = std::ceil(fs[*a_it].series[id].capacity * 1/way_factor);
            
            logger::Instance()->debug("Cap " + std::to_string(g.id(*a_it)) + " " + std::to_string(fs[*a_it].series[id].capacity) + "\n");
            logger::Instance()->debug("Overshoot Test " + std::to_string(g.id(*a_it)) + " " + std::to_string(local_limit) + " " + std::to_string(fac) + " " + std::to_string(way_factor) + "\n");
            if (local_limit < guide_cap || guide_cap == 0) {
                guide_cap = local_limit ;
            }
        }
        
        logger::Instance()->debug("Guide Cap " + std::to_string(guide_cap) + "\n");
        
        // set guide_increase in second and third run!
        rcount mean = 0;
        rcount guide_current = guide_cap;
        a_it = g_it->begin();
        mean += guide_current;
        //rcount max = guide_current;
        
        if (ai[*a_it].edge_type == edge_types::HELPER) {
            guided_percentage[*a_it] = 1;
        } else {
            if (fs[*a_it].series[id].capacity > 0) {
                guided_percentage[*a_it] += guide_current/float(fs[*a_it].series[id].capacity); 
            } else {
                guided_percentage[*a_it] = 1;
            }
        }
        
        ++a_it;
        for (; a_it != g_it->end(); ++a_it) {
            
            float fac = ratios[g.source(*a_it)];
            guide_current = std::ceil(guide_current * fac);
            mean += guide_current;
            //if (guide_current > max) {
            //    max = guide_current;
            //}
            
            logger::Instance()->debug("Current " + std::to_string(g.id(*a_it)) + " " + std::to_string(guide_current) + " " + std::to_string(fac) + "\n");
            
            if (ai[*a_it].edge_type == edge_types::HELPER) {
                guided_percentage[*a_it] = 1;
            } else {
                if (fs[*a_it].series[id].capacity > 0) {
                    guided_percentage[*a_it] += guide_current/float(fs[*a_it].series[id].capacity); 
                } else {
                    guided_percentage[*a_it] = 1;
                }
            }
        }
        mean = mean / (g_it->size() - 1);
        
        //logger::Instance()->debug("Guide Mean " + std::to_string(mean) +  "\n");
        
        for (std::deque<ListDigraph::Arc>::iterator f_it = g_it->begin(); f_it != g_it->end(); ++f_it) {
            guided_capacity[*f_it] += mean;
        }
    }
    
    float truth_limit = 0.8; //1 - options::Instance()->get_guide_trust_level();
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Set " + std::to_string(g.id(a)) + " " + std::to_string(guided_percentage[a]) + " " + std::to_string(guided_capacity[a]) + "\n");
        #endif
        
        if (guided_percentage[a] >= truth_limit) {
            fs[a].series[id].capacity = guided_capacity[a];
        } else if (guided_percentage[a] > 0) {
            fs[a].series[id].capacity = guided_capacity[a]  + fs[a].series[id].capacity * ( 1 - guided_percentage[a] );
        }
    } 
}


bool base_manager::recursive_arc_backtracing(exon_edge& goal, unsigned int next_start, unsigned int goal_end_index, ListDigraph::Node n, ListDigraph::Arc la, std::deque<ListDigraph::Arc> &path, bool exon) {

    if (n == t) {
        // we did it, perhaps
        #ifdef ALLOW_DEBUG
         logger::Instance()->debug("end\n");
        #endif
        return exon && next_start == goal_end_index;
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Recursion " + std::to_string(next_start) + " n " + std::to_string(g.id(n))+ " la " + std::to_string(g.id(la)) +".\n");
     #endif
    
    // we check all outward arcs to find the correct path through the maze
    // we are linear in # edges, but in reality much better
    for (ListDigraph::OutArcIt a(g, n) ; a!=INVALID; ++a) {
        
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Look a " + std::to_string(g.id(a))+".\n");
        #endif
        
        if (ai[a].edge_type == edge_types::BACKLINK) { // ignoring circular references, otherwise we have a BAD time
            continue;
        }
        
        bool correct_edge = false;
        if (ai[a].edge_type != edge_types::EXON) { // all edges except for circular and labeled exon edges
            // just jump add
            
            bool correct = exon;
            if (ai[a].edge_specifier.node_index == next_start && ai[a].edge_specifier.node_index == goal_end_index) {
                correct = true;
            }
            
            correct_edge = recursive_arc_backtracing(goal, next_start, goal_end_index, g.target(a), a, path, correct);
        } else {
            
            unsigned int start_index = ai[a].edge_specifier.id.find_first();
            
            // this is an EXON here
            if (start_index != next_start) {
                continue;
            }
            
            // for guides path evidences do not matter!
            
            unsigned int end_index; // this is slow but meh, doesn't happen too often
            boost::dynamic_bitset<>::size_type index = start_index;
            while(index != boost::dynamic_bitset<>::npos) {
                end_index = index;
                index = ai[a].edge_specifier.id.find_next(index);
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test Cont " + ai[a].edge_specifier.to_string() + " vs " + goal.to_string() + " ; " + std::to_string(start_index) + "-"+  std::to_string(end_index) +".\n");
            #endif
            // exon edge, we only follow this if it is compliant to the goal
            if ( ai[a].edge_specifier.is_contained_in(goal, start_index, end_index) ) {
                
                // it is contained, so we take the route
                correct_edge = recursive_arc_backtracing(goal, end_index, goal_end_index, g.target(a), a, path, true);
            }
        }
        if (correct_edge) { // we found the correct edge and propagate the solution down
            
            path.push_front(a);
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Correct Out\n");
            #endif

            return true;
        }
    }
    
    #ifdef ALLOW_DEBUG
     logger::Instance()->debug("False Out\n");
    #endif
    return false;
}



bool base_manager::contract_composites(ListDigraph& gc,
        ListDigraph::ArcMap<flow_series>& fsc,   
        ListDigraph::ArcMap<arc_identifier> &aic,
        InDegMap<ListDigraph> &in_deg,
        OutDegMap<ListDigraph> &out_deg) {
    
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("contract_composites.\n");
    #endif
    
    bool changed = false;
    
     // loop over all node that can be a potential start
    for (ListDigraph::NodeIt n(gc); n != INVALID; ++n) {
        // we can erase while looping, but be careful
        
        ListDigraph::Node node(n);
                
        if ( out_deg[node] != 1 || in_deg[node] != 1) { // not an in between node from in perspective
            // test backwards arcs
            
            for (ListDigraph::OutArcIt a(gc, node); a!=INVALID; ) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Iterator " + std::to_string(g.id(a)) + "\n");
                #endif
    
                ListDigraph::Node next = gc.target(a);

                if ( (in_deg[next] == 1 || out_deg[next] == 1) && out_deg[next] != 0) { // not an in between node from out perspective

                    // this is the perfect seed to get a path

                    changed = true;

                    // we do the merge
                    // this means we build a composite edge from all encountered

                    flow_series path_fs = fsc[a]; 
                    
                    exon_edge path_ident = aic[a].edge_specifier;
                    edge_length path_length = aic[a].edge_lengths;
                    unsigned int cidin = aic[a].cycle_id_in ;
                    
                    if (aic[a].edge_type == edge_types::NODE) {
                        path_length.first_exon = path_length.middle;
                        path_length.middle = 0;
                    }
                    
                    // if in-degree is > 1 we can not yet delete nodes and arcs, only update them
                    bool del = in_deg[next] <= 1;
                    
                    ListDigraph::Arc arc_tmp(a);
                    ++a;
                    gc.erase(arc_tmp);
                                       
                    // follow contraction
                    follow_contraction(gc, fsc, aic, in_deg, out_deg, path_fs, path_ident, path_length, cidin, node, next, del);
                            
                } else {
                    ++a;
                }
            }
            
        }
    }
    
    return changed;
}


void base_manager::follow_contraction(ListDigraph& gc,
        ListDigraph::ArcMap<flow_series>& fsc,   
        ListDigraph::ArcMap<arc_identifier> &aic,
        InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg,
        flow_series path_fs, exon_edge path_ident, edge_length path_length,
        unsigned int cidin, ListDigraph::Node first, ListDigraph::Node current, bool del) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Follow Contraction " + path_ident.to_string() + " " + std::to_string(del) + "; "+  std::to_string(g.id(first)) + " " + std::to_string(g.id(current)) +  "\n");
    #endif
    
    for (ListDigraph::OutArcIt a(gc, current); a!=INVALID; ) {
        
        ListDigraph::Node next = gc.target(a);
    
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("iterator " + aic[a].edge_specifier.to_string() + " " +  std::to_string(g.id(a)) +  " ; "+ std::to_string(g.id(next)) + "\n");
            for (std::set<int>::iterator idi = input_ids.begin(); idi != input_ids.end(); ++idi) {
                logger::Instance()->debug("means at id " + std::to_string(*idi) + ": " + std::to_string(path_fs.get_mean(*idi).mean) + " " +  std::to_string(path_fs.get_mean(*idi).weight) +  " ; "+ std::to_string(fsc[a].get_flow(*idi)) + " " + std::to_string(path_fs.get_flow(*idi))+ "\n");
            }
        #endif
        
        flow_series next_path_fs = path_fs;
        for (std::set<int>::iterator idi = input_ids.begin(); idi != input_ids.end(); ++idi) {
            
            capacity_type path_flow = next_path_fs.get_flow(*idi);
            capacity_type flow = fsc[a].get_flow(*idi);
            
            logger::Instance()->debug("FLOW " + std::to_string(*idi) +  ": "+ std::to_string(path_flow) +  " "+ std::to_string(flow) + "\n");
            
            if (flow >= path_flow) { // this means we first enter with a second incoming edge !
                capacity_mean new_mean = fsc[a].get_mean(*idi);
                new_mean.reduce(path_flow/(float) flow);
                
                logger::Instance()->debug("NM " + std::to_string(new_mean.mean) +  " "+ std::to_string(new_mean.weight) + "\n");
                
                next_path_fs.get_mean(*idi).update(new_mean); 

            } else { //one to many
                next_path_fs.get_mean(*idi).reduce(flow/(float) path_flow);            
                next_path_fs.get_mean(*idi).update(fsc[a].get_mean(*idi)); 
            }

            if (flow < path_flow) {
                next_path_fs.get_flow(*idi) = flow;
            }
        }
        
        exon_edge up_path_ident = path_ident;
        if (!aic[a].edge_specifier.id.empty()) {
            if (path_ident.id.empty()) {
                up_path_ident = aic[a].edge_specifier;
            } else {
                up_path_ident.join_edge(aic[a].edge_specifier);
            }
        } else {
            if (path_ident.id.empty() && path_ident.node_index == -1) {
                up_path_ident = aic[a].edge_specifier;
            }
        }   
        
        edge_length up_path_length;
        if (aic[a].edge_type == edge_types::EXON) {
            up_path_length = aic[a].edge_lengths;
            up_path_length.middle += path_length.middle + path_length.last_exon;
            up_path_length.first_exon = path_length.first_exon;
        } else {
            up_path_length = path_length;
            if (up_path_length.first_exon == 0) {
                up_path_length.first_exon = aic[a].edge_lengths.middle;
            }
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Next degrees " + std::to_string(in_deg[next]) + " " + std::to_string(out_deg[next]) +  "\n");
        #endif
        
        if ( del && ( in_deg[next] == 1 || out_deg[next] == 1) && out_deg[next] != 0) {
    
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Next 1 \n");
            #endif
            
            // can continue to merge this one :)
            // we can even already remove this arc
         
            bool up_del = in_deg[next] == 1;
            
            ListDigraph::Arc arc_tmp(a);
            ++a;
            gc.erase(arc_tmp);
            if (out_deg[current] == 0 ) {
                gc.erase(current);
            }
            
            follow_contraction(gc, fsc, aic, in_deg, out_deg, next_path_fs, up_path_ident, up_path_length, cidin, first, next, up_del);
            
        } else if ( !del &&  out_deg[next] == 1 ) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Next 2 \n");
            #endif
            
            // can continue to merge this one :)
            // however, another contraction blocks this one
            
            ++a;
            follow_contraction(gc, fsc, aic, in_deg, out_deg, next_path_fs, up_path_ident, up_path_length, cidin, first, next, del);
            
        } else {
 
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Next 3 Arc: " + up_path_ident.to_string() + " " + std::to_string(gc.id(first)) + " " + std::to_string(gc.id(next)) + "\n");
            #endif
            
            // I won’t give up, no I won’t give in
            // Till I reach the end
            // And then I’ll start again
            
            // we have found the end, create new edges if we deleted to up here and have fun :)
            
            ListDigraph::Arc arc_tmp(a);
                ++a;
            
            ListDigraph::Arc narc = gc.addArc(first, next);
            fsc[narc] = next_path_fs;
            aic[narc].edge_specifier = up_path_ident;
            if (up_path_ident.id.empty()) {
                aic[narc].edge_type = edge_types::HELPER; // in rare cases helper are only contracted into the node 
            } else {
                aic[narc].edge_type = edge_types::EXON;  
            }
            aic[narc].edge_lengths = up_path_length;
            aic[narc].cycle_id_in = cidin;
            aic[narc].cycle_id_out = aic[arc_tmp].cycle_id_out;
            
            if (del) {    
                 gc.erase(arc_tmp);
                 
                 if (a==INVALID) {
                     gc.erase(current);
                 }
            }
            
        }    
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Exit \n");
    digraphWriter(gc, std::cout)
        .arcMap("ai", aic)
        .arcMap("flow", fsc)
        .run();
    #endif
    
}

bool base_manager::simplify_ambiguous_ILP(ListDigraph& gc,
        ListDigraph::Node& sc,
        ListDigraph::Node& tc,
        ListDigraph::ArcMap<flow_series>& fsc,   
        ListDigraph::ArcMap<arc_identifier> &aic,
        ListDigraph::NodeMap<unsigned int>& nic,
        ListDigraph::ArcMap<arc_bridge>& kpc,
        ListDigraph::ArcMap<arc_back_bridge>& kbpc,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &uac,
        ListDigraph::NodeMap< unsecurity_id> &uic,
        int guiding_id,
        ListDigraph::ArcMap<bool> &barred, std::set<int>& ip_ids, ATC &trans) {
        
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Simplify ambiguous ILP \n");
    #endif
     
    for (ListDigraph::ArcIt a(gc); a != INVALID;++a) {
        for (std::set<int>::iterator iii = ip_ids.begin(); iii != ip_ids.end(); ++iii) {
            fsc[a].get_mean(*iii).reduce_score();
        }
    }
    
    ListDigraph::ArcMap<arc_bridge> temp_know_paths(gc);
    ListDigraph::ArcMap<arc_back_bridge> temp_know_back_paths(gc);
    ListDigraph::ArcMap<arc_bridge> temp_know_paths_add_max(gc);
    ListDigraph::ArcMap<arc_back_bridge> temp_know_back_paths_add_max(gc);
    ListDigraph::NodeMap<lazy<classify_component> > classify(gc);
    ListDigraph::NodeMap<std::unordered_set<int> > single_to_single(gc);
    ListDigraph::NodeMap<std::unordered_set<int> > block_delete_high_res(gc);
    ListDigraph::NodeMap<std::unordered_set<int> > block_delete(gc);
    
    ListDigraph::NodeMap<std::map<int, evidence_group> > left_groups_ids(gc);
    ListDigraph::NodeMap<std::map<int, evidence_group> > right_groups_ids(gc);
    
    // ######### Classify Each Node ######### //
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Classify ##################### \n");
    #endif
  
    for (ListDigraph::NodeIt n(gc); n != INVALID;++n) {
        // we can erase while looping, but be careful

        ListDigraph::Node node(n);       
        
        // out_deg and in_deg unfunction because lemon are fucking stupip fucks
        
        unsigned int pre_size = 0;
        unsigned int post_size = 0;
        
        for (ListDigraph::InArcIt a(gc, node); a!=INVALID; ++a) {
                if (barred[a]) { // do NOT include barred edges
                    continue;
                }
                ++pre_size;
            }
            
            for (ListDigraph::OutArcIt a(gc, node); a!=INVALID; ++a) {
                if (barred[a]) { // do NOT include barred edges
                    continue;
                }
                ++post_size;
            }
        
        if ( pre_size > 1 && post_size > 1) {
            // this is an unresolvable node
                        
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Classify Node " + std::to_string(gc.id(n)) + "\n");
            #endif
                        
            std::deque<path> pre_path;
            std::deque<path> post_path;
            
            bool helper = false;
            
            // we init pre and post with directly adjacent edges only
            for (ListDigraph::InArcIt a(gc, node); a!=INVALID; ++a) {
                if (aic[a].edge_type == edge_types::HELPER) {
                    // helper can't possibly be disproven/proven
                    helper = true;
                    continue;
                } else if (barred[a] || aic[a].edge_type == edge_types::RESOLVE_HELPER) { // do NOT include barred edges unless for final resolutions!
                    continue;
                }
                pre_path.push_back( path(&aic[a].edge_specifier, nic[gc.source(a)], a, gc.source(a)) );
            }
            
            for (ListDigraph::OutArcIt a(gc, node); a!=INVALID; ++a) {
                if (aic[a].edge_type == edge_types::HELPER) {
                    // helper can't possibly be disproven/proven
                    helper = true;
                    continue;
                } else if (barred[a] || aic[a].edge_type == edge_types::RESOLVE_HELPER) { // do NOT include barred edges unless for final resolutions!
                    continue;
                }
                post_path.push_back( path(&aic[a].edge_specifier, nic[gc.target(a)], a, gc.target(a)) );
            }
            
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("l-Path: ");
            for ( std::deque<path>::iterator it = pre_path.begin(); it!= pre_path.end(); ++it) {
                logger::Instance()->debug("(" + it->identifier.to_string()+")" );
            }
            logger::Instance()->debug("\n");
            logger::Instance()->debug("r-Path: ");
            for ( std::deque<path>::iterator it = post_path.begin(); it!= post_path.end(); ++it) {
                logger::Instance()->debug("(" + it->identifier.to_string()+")" );
            }
            logger::Instance()->debug("\n");
            #endif
            
            // we need to get the id of the exon we bridge here
            ListDigraph::Arc ei;
            ListDigraph::OutArcIt a(gc, node);
            while (a != INVALID && (aic[a].edge_type == edge_types::HELPER || aic[a].edge_type == edge_types::RESOLVE_HELPER)) {
                ++a;
            }
            ei = a;
            if (a==INVALID) {
                ListDigraph::InArcIt a(gc, node);
                while (aic[a].edge_type == edge_types::HELPER || aic[a].edge_type == edge_types::RESOLVE_HELPER) {
                    ++a;
                }
                ei = a;
            }
            
            unsigned int exon_index = aic[ei].edge_specifier.id.find_first();
                 
            unsigned int pre_size_paths = pre_path.size();
            unsigned int post_size_paths = post_path.size();
            
            // we are only interested in the direct arcs
            // all paths going further cannot be resolved, because more path through these nodes might exist
            std::set<int> pre_path_match;
            std::set<int> post_path_match;
            
            // we do the full UNBLOCKED count into this map!
            std::map<std::set<int> , std::map<std::set<int>,  rcount > > hit_reference;
            
            // +++++++++++++++++++++  Overarching  +++++++++++++++++++++
            // we start with fragments overarching this node before dealing with paired data
            exon_edge last;
            for (graph_list<std::pair<exon_group *, gmap<int, rcount> > >::iterator si = raw->single_for_exon[exon_index].begin(); si != raw->single_for_exon[exon_index].end() ;++si) {
                
                exon_edge left, right;
                si->first->bin_mask.left_split(exon_index, left);
                si->first->bin_mask.right_split(exon_index, right);
                
                if (guiding_id != -1 && si->second.find(guiding_id) == si->second.end()) {
                    continue;
                }
                
                rcount count = 0;
                if (guiding_id == -1) {
                    for(gmap<int, rcount>::iterator ci = si->second.begin(); ci != si->second.end(); ++ci) {
                        count += ci->second;
                    }
                } else {
                    count = si->second[guiding_id];
                }
                if (count == 0) {
                    continue;
                }
                count = count * 2;
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Test Overarch:  " + left.to_string() + " " + right.to_string() + " " + std::to_string(count) + "\n");
                #endif
                
                unsigned int leftmost = si->first->range_start;
                unsigned int rightmost = si->first->range_end;
                
                if (last != left) {
                    // finalize last set :)
                    
                    pre_path_match.clear();

                    // go over all left paths first
                    for ( unsigned i = 0; i < pre_path.size(); ++i) {

                        path *pp = &pre_path[i];
                        
                        // can this match? path long enough?
                        if (pp->border_index > leftmost) {
                            // extendpath
                            extend_path_left(pre_path, pp, leftmost, gc, aic, nic, kpc);
                        }
                        
                        if (left.is_contained_in(pp->identifier, leftmost, exon_index)) {
                            pre_path_match.insert( gc.id(pp->starting_arc) );
//                            #ifdef ALLOW_DEBUG
//                            logger::Instance()->debug("Match L: " + std::to_string(gc.id(pp->starting_arc)) + "\n");
//                            #endif  
                        }
                    }
                    last = left; 
                }
                    
                post_path_match.clear();
                   
                // then we go over right paths
                for ( rpos i = 0; i < post_path.size(); ++i ) {
                    
                    path *pp = &post_path[i];
                    
                    // can this match? path long enough?
                    if (pp->border_index < rightmost) {
                        // extendpath
                        extend_path_right(post_path, pp, rightmost, gc, aic, nic, kpc);
                    }
                    
                    if (right.is_contained_in(pp->identifier, exon_index, rightmost)) {
                        post_path_match.insert(gc.id(pp->starting_arc));
//                        #ifdef ALLOW_DEBUG
//                        logger::Instance()->debug("Match R: " + std::to_string(gc.id(pp->starting_arc)) + "\n");
//                        #endif 
                    }
                }
                
                if (!pre_path_match.empty() && !post_path_match.empty()) { 
                    hit_reference[pre_path_match][post_path_match] += count; 
                }
            }
            
            pre_path_match.clear();
            post_path_match.clear();
            
            // +++++++++++++++++++++  Treat Pairs as Overarching when possible  +++++++++++++++++++++
            // we start with fragments overarching this node before dealing with paired data
            last = exon_edge();
            std::unordered_set<paired_exon_group *> opairs;
            for (graph_list<paired_exon_group *>::iterator pi = raw->pairs_for_exon[exon_index].begin(); pi != raw->pairs_for_exon[exon_index].end() ; ++pi) {
                

                if ( (*pi)->left_read->range_end < exon_index && (*pi)->right_read->range_start > exon_index) { 
                    continue;
                } 
                
                exon_edge join = (*pi)->left_read->bin_mask;
                join.join_edge((*pi)->right_read->bin_mask);
                
                exon_edge left, right;
                join.left_split(exon_index, left);
                join.right_split(exon_index, right);

                if (guiding_id != -1 && (*pi)->count.find(guiding_id) == (*pi)->count.end()) {
                    continue;
                }
                
                rcount count = 0;
                if (guiding_id == -1) {
                    for(gmap<int, rcount>::iterator ci = (*pi)->count.begin(); ci != (*pi)->count.end(); ++ci) {
                        count += ci->second;
                    }
                } else {
                    count = (*pi)->count[guiding_id];
                }
                if (count == 0) {
                    continue;
                }
    
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Test OverarchB: " + left.to_string() + " " + right.to_string() + " " + std::to_string(count) + "\n");
                logger::Instance()->debug("       BasedOn: " + (*pi)->left_read->bin_mask.to_string() + " " + (*pi)->right_read->bin_mask.to_string() + "\n");
                #endif
                
                unsigned int leftmost = (*pi)->left_read->range_start;
                unsigned int rightmost = (*pi)->right_read->range_end;
                
                if (last != left) {
                    // finalize last set :)
                    
                    pre_path_match.clear();

                    // go over all left paths first
                    for ( unsigned i = 0; i < pre_path.size(); ++i) {

                        path *pp = &pre_path[i];
                        
                        // can this match? path long enough?
                        if (pp->border_index > leftmost) {
                            // extendpath
                            extend_path_left(pre_path, pp, leftmost, gc, aic, nic, kpc);
                        }
                        
                        if (left.is_contained_in(pp->identifier, leftmost, exon_index)) {
                            pre_path_match.insert( gc.id(pp->starting_arc) );
//                            #ifdef ALLOW_DEBUG
//                            logger::Instance()->debug("Match L: " + std::to_string(gc.id(pp->starting_arc)) + "\n");
//                            #endif  
                        }
                    }
                    last = left; 
                }
                    
                post_path_match.clear();
                   
                // then we go over right paths
                for ( rpos i = 0; i < post_path.size(); ++i ) {
                    
                    path *pp = &post_path[i];
                    
                    // can this match? path long enough?
                    if (pp->border_index < rightmost) {
                        // extendpath
                        extend_path_right(post_path, pp, rightmost, gc, aic, nic, kpc);
                    }
                    
                    if (right.is_contained_in(pp->identifier, exon_index, rightmost)) {
                        post_path_match.insert(gc.id(pp->starting_arc));
//                        #ifdef ALLOW_DEBUG
//                        logger::Instance()->debug("Match R: " + std::to_string(gc.id(pp->starting_arc)) + "\n");
//                        #endif 
                    }
                }
                
                if (!pre_path_match.empty() && !post_path_match.empty()) { 
                    hit_reference[pre_path_match][post_path_match] += count; 
                    opairs.insert(*pi);
                }
            }
            
            pre_path_match.clear();
            post_path_match.clear();
            
            // +++++++++++++++++++++  Paired  +++++++++++++++++++++
            // in the second step we apply paired end data
            exon_group* lastp = NULL;
            for (graph_list<paired_exon_group *>::iterator pi = raw->pairs_for_exon[exon_index].begin(); pi != raw->pairs_for_exon[exon_index].end() ; ++pi) {
                
                if (opairs.find(*pi) != opairs.end()) {
                    continue;
                }
                
                exon_edge left_a, left_b, right_a, right_b;
                rpos la1, la2, lb1, lb2, ra1, ra2, rb1, rb2;
                
                if ( (*pi)->left_read->range_end > exon_index ) { 
                    (*pi)->left_read->bin_mask.left_split(exon_index, left_a);
                    (*pi)->left_read->bin_mask.right_split(exon_index, left_b);
                    
                    la1 = (*pi)->left_read->range_start;
                    la2 = exon_index;
                    lb1 = exon_index;
                    lb2 = (*pi)->left_read->range_end;
                } else {
                    left_a = (*pi)->left_read->bin_mask;
                    la1 = (*pi)->left_read->range_start;
                    la2 = (*pi)->left_read->range_end;
                    lb1 = 0;
                    lb2 = 0;
                }
                
                if ( (*pi)->right_read->range_start < exon_index ) { 
                    
                    (*pi)->right_read->bin_mask.left_split(exon_index, right_a);
                    (*pi)->right_read->bin_mask.right_split(exon_index, right_b);

                    ra1 = (*pi)->right_read->range_start;
                    ra2 = exon_index;
                    rb1 = exon_index;
                    rb2 = (*pi)->right_read->range_end;
                } else {
                    right_b = (*pi)->right_read->bin_mask;
                    ra1 = 0;
                    ra2 = 0;
                    rb1 = (*pi)->right_read->range_start;
                    rb2 = (*pi)->right_read->range_end;
                }
                
                unsigned int leftmost = (*pi)->left_read->range_start;
                unsigned int rightmost = (*pi)->right_read->range_end;
                                
                if (guiding_id != -1 && (*pi)->count.find(guiding_id) == (*pi)->count.end()) {
                    continue;
                }
                
                rcount count = 0;
                if (guiding_id == -1) {
                    for(gmap<int, rcount>::iterator ci = (*pi)->count.begin(); ci != (*pi)->count.end(); ++ci) {
                        count += ci->second;
                    }
                } else {
                    count = (*pi)->count[guiding_id];
                }
                if (count == 0) {
                    continue;
                }
    
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Test Paired:    " + (*pi)->left_read->bin_mask.to_string() + " " + (*pi)->right_read->bin_mask.to_string() + " " + std::to_string(count) + "\n");
                #endif
                
                if (lastp != (*pi)->left_read ) {
                    // finalize last set :)
                    
                    pre_path_match.clear();

                    // go over all left paths first
                    for ( unsigned i = 0; i < pre_path.size(); ++i) {

                        path *pp = &pre_path[i];
                        
                        // can this match? path long enough?
                        if (pp->border_index > leftmost) {
                            // extendpath
                            extend_path_left(pre_path, pp, leftmost, gc, aic, nic, kpc);
                        }
                        
                        if ( left_a.is_contained_in(pp->identifier, la1, la2 ) && (ra1 == 0 || right_a.is_contained_in(pp->identifier, ra1, ra2)) ) {
                            pre_path_match.insert( gc.id(pp->starting_arc) );
//                            #ifdef ALLOW_DEBUG
//                            logger::Instance()->debug("Match L: " + std::to_string(gc.id(pp->starting_arc)) + "\n");
//                            #endif 
                        }
                    }
                    lastp = (*pi)->left_read;
                }
                
                post_path_match.clear();
                                
                // then we go over right paths
                for ( unsigned j = 0; j < post_path.size(); ++j ) {
                    
                    path *pp = &post_path[j];
                    
                    // can this match? path long enough?
                    if (pp->border_index < rightmost) {
                        // extendpath
                        extend_path_right(post_path, pp, rightmost, gc, aic, nic, kpc);
                    }
                    
                    if ( right_b.is_contained_in(pp->identifier, rb1, rb2) && (lb1 == 0 || left_b.is_contained_in(pp->identifier, lb1, lb2 )) ) {
                        post_path_match.insert(gc.id(pp->starting_arc));
//                        #ifdef ALLOW_DEBUG
//                        logger::Instance()->debug("Match R " + std::to_string(gc.id(pp->starting_arc)) + "\n");
//                        #endif 
                    }
                }
                
                bool left_unique = pre_path_match.size() == 1;
                bool right_unique = post_path_match.size() == 1;
                if (!left_unique && !pre_path_match.empty()) {
                    std::set<int>::iterator first = pre_path_match.begin();
                    std::set<int>::iterator pi = first;
                    ++pi;
                    bool breaked = false;
                    for (; pi !=  pre_path_match.end(); ++pi) {
                        rpos pos = std::max(nic[gc.source(gc.arcFromId(*first))], nic[gc.source(gc.arcFromId(*pi))]);
                        if (!aic[gc.arcFromId(*first)].edge_specifier.is_contained_in(aic[gc.arcFromId(*pi)].edge_specifier, std::max(la1, pos), exon_index)) {
                            breaked = true;
                            break;
                        }
                    }
                    if (!breaked) {
                        left_unique = true;
                    }
                }
                if (!right_unique && !post_path_match.empty()) {
                    std::set<int>::iterator first = post_path_match.begin();
                    std::set<int>::iterator pi = first;
                    ++pi;
                    bool breaked = false;
                    for (; pi !=  post_path_match.end(); ++pi) {
                        rpos pos = std::max(nic[gc.target(gc.arcFromId(*first))], nic[gc.target(gc.arcFromId(*pi))]);
                        if (!aic[gc.arcFromId(*first)].edge_specifier.is_contained_in(aic[gc.arcFromId(*pi)].edge_specifier, std::min(pos, rb2), exon_index)) {
                            breaked = true;
                            break;
                        }
                    }
                    if (!breaked) {
                        right_unique = true;
                    }
                }
                
                if (left_unique && right_unique) { 
                    hit_reference[pre_path_match][post_path_match] += count;
                }
            }

            pre_path_match.clear();
            post_path_match.clear();
            
            // +++++++++++++++++++++  Get Simplified Block Info  +++++++++++++++++++++
            
            // we check for delete-groups as well!
            std::set<std::set<int> > left_raw_groups;
            std::set<std::set<int> > right_raw_groups;
            
            // valid!
            std::map<std::set<int> , std::map<std::set<int>,  bool > > valid_hits;
            
            // transfer reference to sorted list
            typedef std::deque<std::pair< std::set<int>, rcount> > ISRT;
            typedef std::deque<std::pair<std::set<int>, ISRT > > SRT;
            
            unsigned int max_mult = 0;
            
            SRT s_reference;
            for (std::map<std::set<int> , std::map<std::set<int>,  rcount > >::iterator i = hit_reference.begin(); i != hit_reference.end(); ++i) {
                s_reference.push_back(std::make_pair(i->first, std::deque<std::pair< std::set<int>, rcount> >()));
                
                if (i->first.size() > max_mult) {
                    max_mult = i->first.size();
                }
                
                for (std::map<std::set<int>,  rcount >::iterator i2 = i->second.begin(); i2 != i->second.end(); ++i2) {
                    
                    if (i2->first.size() > max_mult) {
                        max_mult = i2->first.size();
                    }
                    
                    s_reference.back().second.push_back(std::make_pair(i2->first, i2->second));
                    
                }    
                std::sort(s_reference.back().second.begin(), s_reference.back().second.end(), []( std::pair< std::set<int>, rcount > & a, std::pair< std::set<int>, rcount > & b) -> bool {return a.first.size() < b.first.size();});
            }
            std::sort(s_reference.begin(), s_reference.end(), []( std::pair<std::set<int>,  std::deque<std::pair< std::set<int>, rcount> > > & a,std::pair<std::set<int>,  std::deque<std::pair< std::set<int>, rcount> > > & b) -> bool {return a.first.size() < b.first.size();});
            
            // we go over the evidences!
            path_evidence_map<int, path_evidence > evidences;
            bool has_evidence = false;     
            
            rcount threshold = 10;
//            std::map<int, std::pair<unsigned int, rcount> > individual_count_left;
//            std::map<int, std::pair<unsigned int, rcount> > individual_count_right;
                
//            for (unsigned int max_mult_counter = 1; max_mult_counter <= max_mult; ++max_mult_counter) {
//                
//                std::map<int, std::pair<unsigned int, rcount> > step_count_left;
//                std::map<int, std::pair<unsigned int, rcount> > step_count_right;
                
            for (SRT::iterator i = s_reference.begin(); i != s_reference.end(); ++i) { 

                if (i->first.size() == pre_size_paths) {
                    continue;
                }

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("-----\n");
                for (std::set<int>::iterator li = i->first.begin(); li != i->first.end(); ++li) {
                    logger::Instance()->debug("GL:  " + std::to_string(*li) + "\n");
                }
                #endif

                path_evidence_map<int, rcount> all_found;
                path_evidence_set<int> pre_left, pre_right;
                rcount total_count = 0;
                for (ISRT::iterator i2 = i->second.begin(); i2 != i->second.end(); ++i2) {

//                        if (max_mult_counter != i->first.size() && max_mult_counter != i2->first.size()) {
//                            continue;
//                        }

                    if (i2->first.size() == post_size_paths) {
                        continue;
                    }

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("> " + std::to_string(i2->second) + "\n");
                    for (std::set<int>::iterator it = i2->first.begin(); it!=i2->first.end(); ++it) {
                        logger::Instance()->debug("GR:  " + std::to_string(*it)  + "\n");
                    }
                    #endif

                    if (i->first.size() == 1 && i2->first.size() == 1) {

                        single_to_single[n].insert(*i->first.begin());
                        single_to_single[n].insert(*i2->first.begin());

                        block_delete_high_res[n].insert(*i->first.begin());
                        block_delete_high_res[n].insert(*i2->first.begin());
                    }

                    
//                    bool valid = true;
//                    for (std::set<int>::iterator ri = i2->first.begin(); ri != i2->first.end(); ++ri) {
//                        if (all_found.find(*ri) != all_found.end()) {
//                            valid = false;
//                            break;
//                        }
//                    }
//                    if (!valid) {
//                        #ifdef ALLOW_DEBUG
//                        logger::Instance()->debug("invalid 1\n");
//                        #endif
//                        continue;
//                    }
//                    
//                    for (std::set<int>::iterator li = i->first.begin(); li != i->first.end(); ++li) {
//                        for (std::set<int>::iterator ri = i2->first.begin(); ri != i2->first.end(); ++ri) {
//                            if (evidences[*li].is_evidence(*ri)) {
//                                valid = false;
//                                pre_left.insert(*li);
//                                pre_right.insert(*ri);
//                            }
//                        }
//                    }
                    
                    rcount sub_count = 0;
                    for (SRT::iterator i_p = s_reference.begin(); i_p != i+1; ++i_p) { // we thest if we have a smaller version of this

                        bool left_fit = true;
                        for (std::set<int>::iterator li = i_p->first.begin(); li != i_p->first.end(); ++li) {  // can ALL be found in the current left?
                            if (i->first.find(*li) == i->first.end()) {
                                left_fit = false;
                                break;
                            }
                        }
                        if (!left_fit) {
                            continue;
                        }

                        for (ISRT::iterator i2_p = i_p->second.begin(); i2_p != i_p->second.end(); ++i2_p) {

                            if (i_p == i && i2_p >= i2) {
                                continue;
                            }

                            bool right_fit = true;
                            for (std::set<int>::iterator ri = i2_p->first.begin(); ri != i2_p->first.end(); ++ri) {

                                if (i2->first.find(*ri) == i2->first.end()) {
                                    right_fit = false;
                                    break;
                                }
                            }
                            if (!right_fit) {
                                continue;
                            }
                                sub_count += i2_p->second;
                        }
                    }

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Sub Count " + std::to_string(sub_count) + "\n");
                    #endif
                    if (sub_count > 4) {
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("invalid\n");
                        #endif
                        continue;
                    }

                    for (std::set<int>::iterator it = i2->first.begin(); it!=i2->first.end(); ++it) {
                        all_found[*it] += i2->second; 
                        total_count += i2->second;
                    }

                    right_raw_groups.insert(i2->first);   
                    valid_hits[i->first][i2->first] = true;
                }

                if (total_count == 0) {
                        continue;
                }

                left_raw_groups.insert(i->first);

                if (total_count > threshold) { // we have true evidence!
                    has_evidence = true;

                    // collect block-list
                    path_evidence_set<int> block_set;
                    for (ListDigraph::OutArcIt a(gc, n); a!=INVALID; ++a) {
                        if (aic[a].edge_type == edge_types::HELPER) {
                            continue;  // helper can't possibly be disproven/proven
                        } 
                        int id = gc.id(a);
                        path_evidence_map<int, rcount>::iterator fi = all_found.find(id);
                        if ( fi == all_found.end()  && pre_right.find(id) == pre_right.end()) {
                            block_set.insert(id);
                        }
                    }                                           

                    for (std::set<int>::iterator li = i->first.begin(); li != i->first.end(); ++li) {

                        for (path_evidence_map<int, rcount>::iterator ri = all_found.begin(); ri != all_found.end(); ++ri) {
                            if (!evidences[*li].is_blocked(ri->first)) evidences[*li].add_evidence(ri->first, ri->second);
                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Add:  " + std::to_string(*li) + " " + std::to_string(ri->first) + "\n");
                            #endif
                        }
                        for (path_evidence_set<int>::iterator ri = block_set.begin(); ri != block_set.end(); ++ri) {
                             if (!evidences[*li].is_evidence(*ri)) evidences[*li].add_blocked(*ri);
                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Block:  " + std::to_string(*li) + " " + std::to_string(*ri) + "\n");
                            #endif
                        }
                    }
                } else { 
                    for (std::set<int>::iterator li = i->first.begin(); li != i->first.end(); ++li) {
                        for (path_evidence_map<int, rcount>::iterator ri = all_found.begin(); ri != all_found.end(); ++ri) {

                            evidences[*li].add_evidence(ri->first, ri->second);

//                                step_count_left[*li].first += 1;
//                                step_count_left[*li].second += ri->second;
//                                step_count_right[ri->first].first += 1;
//                                step_count_right[ri->first].second += ri->second;

                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Add Weak:  " + std::to_string(*li) + " " + std::to_string(ri->first) + "\n");
                            #endif
                        }
                    }   
                }
            }
//            }
            
            // +++++++++++++++++++++  Add in Info  +++++++++++++++++++++
            
            // now resolve what can be resolved
            // we use know_paths and know_back_paths
            
            // we iterate over all possible paths and now mark the found on direct arcs
 
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Add In evidences\n");
            #endif
//
//            path_evidence_map<int, rcount> count_left;
//            path_evidence_map<int, rcount> count_right;
//            for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
//                int pp = gc.id(a);
//                for (path_evidence_map<int, rcount>::iterator ei = evidences[pp].begin(); ei!= evidences[pp].end(); ++ei) {
//                    if (!evidences[pp].is_blocked(ei->first)) {
//                        count_left[pp] += ei->second;
//                        count_right[ei->first] += ei->second;
//                    }
//                }
//            }
//            
            
            for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
                int pp = gc.id(a);
                for (path_evidence_map<int, rcount>::iterator ei = evidences[pp].begin(); ei!= evidences[pp].end(); ++ei) {
                    if (!evidences[pp].is_blocked(ei->first) ) {//&& (ei->second * 100 / count_left[pp] >= 5 || ei->second * 100 / count_right[ei->first] >= 5) ) {
                        temp_know_paths[gc.arcFromId(pp)][ei->first] = ei->second;
                        temp_know_back_paths[gc.arcFromId(ei->first)].bridges->insert(pp);
                        
                        block_delete[n].insert(pp);
                        block_delete[n].insert(ei->first);
                    }
                }
            }
                 
            uic[n].resolvable = has_evidence;
            
            // we compute the groups first 
            std::deque<std::set<int> > left_groups;
            std::deque<std::set<int> > right_groups;
            compute_edge_groups(left_raw_groups, left_groups);
            compute_edge_groups(right_raw_groups, right_groups);
            
            for (std::deque<std::set<int> >::iterator lgi = left_groups.begin(); lgi != left_groups.end(); ++lgi) {
                std::set<int>::iterator s = lgi->begin();
                std::set<int>::iterator s2 = s;
                int max = *s;
                float max_score = fsc[gc.arcFromId(*s)].get_mean(guiding_id).compute_score();
                ++s;
                for(; s != lgi->end(); ++s) {
                    float score = fsc[gc.arcFromId(*s)].get_mean(guiding_id).compute_score();
                    if (score > max_score) {
                        max_score = score;
                        max = *s;
                    }
                }
//                for(; s2 != lgi->end(); ++s2) {
//                    float score = means[gc.arcFromId(*s2)].compute_score();
//                    if ( max_score * 0.5 <= score && block_delete[n].find(*s2) != block_delete[n].end()) block_delete_high_res[n].insert(*s2);
//                }
                if (block_delete[n].find(max) != block_delete[n].end()) block_delete_high_res[n].insert(max);
            }
            
            for (std::deque<std::set<int> >::iterator lgi = right_groups.begin(); lgi != right_groups.end(); ++lgi) {
                std::set<int>::iterator s = lgi->begin();
                std::set<int>::iterator s2 = s;
                int max = *s;
                float max_score = fsc[gc.arcFromId(*s)].get_mean(guiding_id).compute_score();
                ++s;
                for(; s != lgi->end(); ++s) {
                    float score = fsc[gc.arcFromId(*s)].get_mean(guiding_id).compute_score();
                    if (score > max_score) {
                        max_score = score;
                        max = *s;
                    }
                }
//                for(; s2 != lgi->end(); ++s2) {
//                    float score = means[gc.arcFromId(*s2)].compute_score();
//                    if ( max_score * 0.5 <= score && block_delete[n].find(*s2) != block_delete[n].end()) block_delete_high_res[n].insert(*s2);
//                }
                if (block_delete[n].find(max) != block_delete[n].end()) block_delete_high_res[n].insert(max);
            }
            
            int id_iterator = 0;
            for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
                left_groups_ids[n][gc.id(a)] = evidence_group(id_iterator, false);
                ++id_iterator;
            }
            for (ListDigraph::OutArcIt a(gc, n); a!=INVALID; ++a) {
                right_groups_ids[n][gc.id(a)] = evidence_group(id_iterator, false);
                ++id_iterator;
            }
            
            // we will not resolve those in big groups
            if (left_groups.size() == 1) {
                left_groups.clear();
            }
            for (std::deque<std::set<int> >::iterator lgi = left_groups.begin(); lgi != left_groups.end(); ) {
                if (lgi->size() == 1) {
                    lgi = left_groups.erase(lgi);
                } else {
                    for(std::set<int>::iterator s = lgi->begin(); s != lgi->end(); ++s) {
                        left_groups_ids[n][*s] = evidence_group(id_iterator, true);
                        #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Group L Mark"+ std::to_string(id_iterator) +"\n");
                        #endif  
//                        for(std::set<int>::iterator s2 = lgi->begin(); s2 != lgi->end(); ++s2) {
//                            for (arc_bridge::iterator ab = temp_know_paths[gc.arcFromId(*s)].begin(); ab != temp_know_paths[gc.arcFromId(*s)].end(); ++ab) {
//                                 temp_know_back_paths[gc.arcFromId(ab->first)].bridges->insert(*s2);
//                                 temp_know_paths[gc.arcFromId(*s2)][ab->first] += 1;
//                            }
//                        }
                        //if (block_delete[n].find(*s) != block_delete[n].end()) block_delete_high_res[n].insert(*s);
                    }
                    ++id_iterator;
                    ++lgi;
                }
            }
            if (right_groups.size() == 1) {
                right_groups.clear();
            }
            for (std::deque<std::set<int> >::iterator rgi = right_groups.begin(); rgi != right_groups.end(); ) {
                if (rgi->size() == 1) {
                    rgi = right_groups.erase(rgi);
                } else {
                    for(std::set<int>::iterator s = rgi->begin(); s != rgi->end(); ++s) {
                        right_groups_ids[n][*s] = evidence_group(id_iterator, true);
                        #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Group R Mark"+ std::to_string(id_iterator) +"\n");
                        #endif  
//                        for(std::set<int>::iterator s2 = rgi->begin(); s2 != rgi->end(); ++s2) {
//                            for (arc_back_bridge::iterator ab = temp_know_back_paths[gc.arcFromId(*s)].begin(); ab != temp_know_back_paths[gc.arcFromId(*s)].end(); ++ab) {
//                                 temp_know_back_paths[gc.arcFromId(*s2)].bridges->insert(*ab);
//                                 temp_know_paths[gc.arcFromId(*ab)][*s2] += 1;
//                            }
//                        }
                        //if (block_delete[n].find(*s) != block_delete[n].end()) block_delete_high_res[n].insert(*s);
                    }
                    ++id_iterator;
                    ++rgi;
                }
            }
            
            // identify the components
            for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
                
                if(barred[a]) {
                    continue;
                }
                
                bool matched = false;
                for (std::deque<component>::iterator ci = classify[n]->components.begin(); ci !=  classify[n]->components.end(); ++ci) {
                    if ( ci->nodes.find(gc.id(a)) != ci->nodes.end()) {
                        matched = true;
                        break;
                    }
                }
                
                if (!matched) {
                    classify[n]->components.push_back(component());
                    find_component_left( classify[n]->components.back().nodes, gc.id(a), gc, temp_know_paths, temp_know_back_paths);
                }
            }
            // now set up the right numbers!
            for (std::deque<component >::iterator ci = classify[n]->components.begin(); ci !=  classify[n]->components.end(); ) {
                        
                if (ci->nodes.size() == 1) {
//                     #ifdef ALLOW_DEBUG
//                    logger::Instance()->debug("CI Single"+ std::to_string(gc.id(n)) + " as " + std::to_string(*ci->nodes.begin())+"\n");
//                    #endif
                    ci = classify[n]->components.erase(ci);
                    continue;
                }
                
                unsigned int in = 0;
                unsigned int out = 0;
                unsigned int edges = 0;
                for ( std::set<int>::iterator it = ci->nodes.begin(); it != ci->nodes.end(); ++it) {                   
                    ListDigraph::Arc a = gc.arcFromId(*it);
                    if (gc.target(a) == node ) {
                        ++in;
                        edges += temp_know_paths[a].size();
                    } else {
                        ++out;
                        edges += temp_know_back_paths[a].size();
                    }
                }

                ci->in_nodes = in;
                ci->out_nodes = out;
                ci->edges = edges;
                
                classify[n]->assigned_in += in;
                classify[n]->assigned_out += out;
                classify[n]->assigned_edges += edges;
//                
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("CI "+ std::to_string(gc.id(n)) + " as " + std::to_string(in)+"-"+std::to_string(out) + " " + std::to_string(edges) +"\n");
//                #endif
                
                if (in * out == edges / 2 && in > 1 && out > 1) {
                    
                    // do remove?
                    bool kill_this = true;
                    for ( std::set<int>::iterator it = ci->nodes.begin(); it != ci->nodes.end(); ++it) { 
                        if (single_to_single[n].find(*it) != single_to_single[n].end()) {
                            kill_this = false;
                            break;
                        }
                    }
                    
                    if (kill_this) {
                    // remove from knows!
                        for ( std::set<int>::iterator it = ci->nodes.begin(); it != ci->nodes.end(); ++it) {                 
                            ListDigraph::Arc a = gc.arcFromId(*it);
                            if (gc.target(a) == node ) {
                               temp_know_paths[a].clear();
                            } else {
                               temp_know_back_paths[a].clear();
                            }
                            block_delete_high_res[n].erase(*it);
                            block_delete[n].erase(*it);
                        }     
                        // rewind
                        classify[n]->assigned_in -= in;
                        classify[n]->assigned_out -= out;
                        classify[n]->assigned_edges -= edges;
                        
                        ci->in_nodes = 0;
                        ci->out_nodes = 0;
                        ci->edges = 0;
                        
                        ci = classify[n]->components.erase(ci);
                        continue;
                        
                    } else {    
                        // we are fully connected!
                        classify[n]->fully_connected_in += in;
                        classify[n]->fully_connected_out += out;
                        classify[n]->fully_connected_edges += edges;
                    }
                }
                
                // save up the fragments we used to come to this decision
                for (std::map<std::set<int> , std::map<std::set<int>,  rcount > >::iterator i1 = hit_reference.begin(); i1 != hit_reference.end(); ++i1) {
                    if (ci->nodes.find(*i1->first.begin()) == ci->nodes.end()) {
                        continue;
                    }
                    for (std::map<std::set<int>,  rcount >::iterator i2 = i1->second.begin(); i2 != i1->second.end(); ++i2) {
                        
                        if (!valid_hits[i1->first][i2->first]) {
                            continue;
                        }
                        
                        ci->fragments.push_back(resolve_count_set());
                        ci->fragments.back().left = i1->first;
                        ci->fragments.back().right = i2->first;
                        ci->fragments.back().count = i2->second;
                    }
                }
                
                ++ci;
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("CLASSIFY SET On "+ std::to_string(gc.id(n)) + " with " + std::to_string(pre_size)  + "-" + std::to_string(classify[n]->assigned_in) + "  " + std::to_string(post_size)+ "-" + std::to_string(classify[n]->assigned_out)+"\n");
            #endif
            
            classify[n]->unassigned_in = pre_size - classify[n]->assigned_in;
            classify[n]->unassigned_out = post_size - classify[n]->assigned_out;   
         
            if (classify[n]->components.size() == 1) {
                
                std::deque<component >::iterator ci = classify[n]->components.begin();
                if (classify[n]->unassigned_in == 0) {
                    // join out to max in!
                    capacity_type max = 0;
                    ListDigraph::Arc marc;
                    for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
                        if (ci->nodes.find(gc.id(a)) != ci->nodes.end()) {
                            if (fsc[a].get_flow(guiding_id) > max) {
                                max = fsc[a].get_flow(guiding_id);
                                marc = a;
                            }
                        }
                    }
                    for (ListDigraph::OutArcIt a(gc, n); a!=INVALID; ++a) {
                        if (ci->nodes.find(gc.id(a)) == ci->nodes.end() && !barred[a]) {

                            temp_know_paths_add_max[marc][gc.id(a)] += 1;
                            temp_know_back_paths_add_max[a].bridges->insert(gc.id(marc));
                            ci->fragments.push_back(resolve_count_set());
                            ci->fragments.back().left.insert(gc.id(marc));
                            ci->fragments.back().right.insert(gc.id(a));
                            ci->fragments.back().count = 1;
                            ci->nodes.insert(gc.id(a));
                        }
                    }
                } else if (classify[n]->unassigned_out == 0) {
                     // join in to max out!
                    capacity_type max = 0;
                    ListDigraph::Arc marc;
                    for (ListDigraph::OutArcIt a(gc, n); a!=INVALID; ++a) {
                        if (ci->nodes.find(gc.id(a)) != ci->nodes.end()) {
                            if (fsc[a].get_flow(guiding_id) > max) {
                                max = fsc[a].get_flow(guiding_id);
                                marc = a;
                            }
                        }
                    }
                    for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
                        if (ci->nodes.find(gc.id(a)) == ci->nodes.end() && !barred[a]) {

                            temp_know_paths_add_max[a][gc.id(marc)] += 1;
                            temp_know_back_paths_add_max[marc].bridges->insert(gc.id(a));
                            ci->fragments.push_back(resolve_count_set());
                            ci->fragments.back().left.insert(gc.id(a));
                            ci->fragments.back().right.insert(gc.id(marc));
                            ci->fragments.back().count = 1;
                            ci->nodes.insert(gc.id(a));
                        }
                    }
                }
            }

            
            unsigned int edges_in_groups = 0;
            unsigned int nodes_offset = 0;
            std::set<int> lgs;
            std::set<int> rgs;
            std::set<std::pair<int, int> > combinations;
            for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
                ListDigraph::Arc l = a;
                for (arc_bridge::iterator bi = kpc[l].begin(); bi != kpc[l].end(); ++bi) {
                    ListDigraph::Arc r = gc.arcFromId(bi->first);

                    // get the original group combos!
                    int i1 = 0;
                    for ( ;i1 < left_groups.size(); ++i1) {
                        if(left_groups[i1].find(gc.id(l)) != left_groups[i1].end()) {
                            lgs.insert(i1);
                            break;
                        }
                    }
                    int i2 = 0;
                    for ( ;i2 < right_groups.size(); ++i2) {
                        if(right_groups[i2].find(gc.id(r)) != right_groups[i2].end()) {
                            rgs.insert(i2);
                            break;
                        }
                    }
                    
                    if (i1 != left_groups.size() || i2 != right_groups.size()) {    
                        ++edges_in_groups;
                    }
                    
                    if (i1 != left_groups.size() && i2 != right_groups.size()) {
                        combinations.insert(std::make_pair(i1, i2));
                    } else if (i1 == left_groups.size() && i2 != right_groups.size()) {
                        combinations.insert(std::make_pair(i1 + gc.id(l), i2));
                    } else if (i1 != left_groups.size() && i2 == right_groups.size()) {
                        combinations.insert(std::make_pair(i1, i2 + gc.id(r)));
                    }
                    
                }
            }
            for ( std::set<int>::iterator ig = lgs.begin(); ig != lgs.end(); ++ig) {
                nodes_offset += left_groups[*ig].size() - 1;
            }
            for ( std::set<int>::iterator ig = rgs.begin(); ig != rgs.end(); ++ig) {
                nodes_offset += right_groups[*ig].size() - 1;
            }
            
            classify[n]->edge_group_offset = edges_in_groups - combinations.size();
            classify[n]->node_group_offset = nodes_offset;
            
//                         save the maximal member of each group!
//            for(std::set<std::set<int> >::iterator group = left_raw_groups.begin(); group != left_raw_groups.end(); ++group) {
//                float max = 0;
//                ListDigraph::Arc arc; 
//                
//                for (std::set<int>::iterator gi = group->begin(); gi != group->end(); ++gi) {    
//                    ListDigraph::Arc a = gc.arcFromId(*gi);
//                    
//                    float score = means[a].compute_score();
//                    if (score > max) {
//                        max = score;
//                        arc = a;
//                    }
//                }
//                
//                if (block_delete[n].find(gc.id(arc)) != block_delete[n].end()) block_delete_high_res[n].insert(gc.id(arc));
//            }
//                    
//            for(std::set<std::set<int> >::iterator group = right_raw_groups.begin(); group != right_raw_groups.end(); ++group) {
//                float max = 0;
//                ListDigraph::Arc arc; 
//                
//                for (std::set<int>::iterator gi = group->begin(); gi != group->end(); ++gi) {    
//                    ListDigraph::Arc a = gc.arcFromId(*gi);
//                    
//                    float score = means[a].compute_score();
//                    if (score > max) {
//                        max = score;
//                        arc = a;
//                    }
//                }
//                
//                if (block_delete[n].find(gc.id(arc)) != block_delete[n].end()) block_delete_high_res[n].insert(gc.id(arc));
//            }
            
        }
    }
    
    
    // ######### Test Actions! ######### //
    #ifdef ALLOW_DEBUG
    if (options::Instance()->is_debug()) {
    logger::Instance()->debug("Actions ##################### \n");
    digraphWriter(gc, std::cout)
                .arcMap("ai", aic)
                .arcMap("fs", fsc)
                .arcMap("barred", barred)
                .node("source", sc)
                .node("drain", tc)
                .run();  

//    logger::Instance()->debug("Has Block Mark \n");
//    for (ListDigraph::ArcIt a(gc); a != INVALID; ++a) {
//        logger::Instance()->debug("Arc " + std::to_string(gc.id(a)) + "\n"); 
//        bool barred = false;
//        for(std::set<transcript_unsecurity>::iterator u_it = uac[a]->begin(); u_it != uac[a]->end() ; ++u_it) {
//
//            logger::Instance()->debug(std::to_string(u_it->id) + "-" + std::to_string(u_it->evidenced) + ", " );
//        }
//        logger::Instance()->debug("\n"); 
//    }
//    logger::Instance()->debug("Nodes \n");
//    for (ListDigraph::NodeIt n(gc); n != INVALID; ++n) {
//        logger::Instance()->debug("Node " + std::to_string(gc.id(n)) + " : " + std::to_string(uic[n].id) + "\n"); 
//    }

    logger::Instance()->debug("Classify done \n");
    for (ListDigraph::ArcIt a(gc); a != INVALID; ++a) {
       logger::Instance()->debug("Arc " + std::to_string(gc.id(a))); 
       for (arc_bridge::iterator ab = temp_know_paths[a].begin(); ab != temp_know_paths[a].end(); ++ab) {
           logger::Instance()->debug(" f " + std::to_string(ab->first)); 
       }
       for (arc_bridge::iterator ab = temp_know_paths_add_max[a].begin(); ab != temp_know_paths_add_max[a].end(); ++ab) {
           logger::Instance()->debug(" fm " + std::to_string(ab->first)); 
       }
       for (arc_back_bridge::iterator ab = temp_know_back_paths[a].begin(); ab != temp_know_back_paths[a].end(); ++ab) {
           logger::Instance()->debug(" b " + std::to_string(*ab)); 
       }
       for (arc_back_bridge::iterator ab = temp_know_back_paths_add_max[a].begin(); ab != temp_know_back_paths_add_max[a].end(); ++ab) {
           logger::Instance()->debug(" bm " + std::to_string(*ab)); 
       }
       logger::Instance()->debug("\n" );  
    }
    
//    logger::Instance()->debug("-------\n" ); 
//    for (ListDigraph::ArcIt a(gc); a != INVALID; ++a) {
//       logger::Instance()->debug("Arc " + std::to_string(gc.id(a))); 
//       for (arc_bridge::iterator ab = kpc[a].begin(); ab != kpc[a].end(); ++ab) {
//           logger::Instance()->debug(" cf " + std::to_string(ab->first)); 
//       }
//       for (arc_back_bridge::iterator ab = kbpc[a].begin(); ab != kbpc[a].end(); ++ab) {
//           logger::Instance()->debug(" cb " + std::to_string(*ab)); 
//       }
//       logger::Instance()->debug("\n" );  
//    }
    
    for (ListDigraph::NodeIt n(gc); n != INVALID;++n) {
        logger::Instance()->debug("For n " + std::to_string(gc.id(n)) +": \n");
        logger::Instance()->debug("Left Groups\n");
        for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
            logger::Instance()->debug("Arc " + std::to_string(gc.id(a)) + " G " + std::to_string(left_groups_ids[n][gc.id(a)].id)  + "\n"); 
        }
        logger::Instance()->debug("Right Groups\n");
        for (ListDigraph::OutArcIt a(gc, n); a!=INVALID; ++a) {
            logger::Instance()->debug("Arc " + std::to_string(gc.id(a)) + " G " + std::to_string(right_groups_ids[n][gc.id(a)].id)  + "\n"); 
        }
    }
    }
    #endif

    
//    std::unordered_set<int> delete_block;
//    for (ListDigraph::ArcIt a(g); a != INVALID;++a) {  
//        
//        ListDigraph::Node ns = gc.source(a);
//        ListDigraph::Node nt = gc.target(a);
//        
//        bool lblock = block_delete_high_res[ns].find(gc.id(a)) !=  block_delete_high_res[ns].end();
//        bool rblock = block_delete_high_res[nt].find(gc.id(a)) !=  block_delete_high_res[nt].end();
//       
//        if ( rblock || lblock ) {
//            delete_block.insert(gc.id(a));
//        }
//    }
    
    
    std::deque<ListDigraph::Node> top_order;

    pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > visitor(top_order);
    DfsVisit<ListDigraph, pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > > dfs(gc, visitor);
    dfs.init();
    dfs.addSource(sc);
    dfs.start();
    
    std::deque<std::pair<ListDigraph::Node, float> > ratio_order;
    for (ListDigraph::NodeIt n(gc); n != INVALID;++n) {  
        float error = 0;
        capacity_type f_count = 0;
        for ( std::deque<component>::iterator ci = classify[n]->components.begin(); ci != classify[n]->components.end(); ++ci) {
            compute_ratio(n, temp_know_paths, temp_know_back_paths, ci->nodes, gc, fsc, guiding_id, error, f_count);
            
                for (std::set<int>::iterator cii = ci->nodes.begin(); cii != ci->nodes.end(); ++cii) {
                    logger::Instance()->debug(std::to_string(*cii) + ", ");
                }            
        }
        ratio_order.push_back(std::make_pair(n, error / (float) f_count));
    }
    std::sort(ratio_order.begin(), ratio_order.end(), [](const std::pair<ListDigraph::Node, float> &a1, const std::pair<ListDigraph::Node, float> &a2){return a1.second < a2.second;});
    
    
    const float low_barr_ratio = 0.40;
        
    // we have a perfect match
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
        
        if ( classify[n]->unassigned_in + classify[n]->assigned_in < 2 || classify[n]->unassigned_out + classify[n]->assigned_out < 2 ) {
            continue;
        }
        
        unsigned int max = classify[n]->assigned_in;
        if (classify[n]->assigned_out > max) {
            max = classify[n]->assigned_out;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Testing "+ std::to_string(gc.id(n)) + " un " + std::to_string(classify[n]->unassigned_in)+"-"+std::to_string(classify[n]->unassigned_out) + " as " + std::to_string(classify[n]->assigned_in)+"-"+std::to_string(classify[n]->assigned_out) + " " + std::to_string(classify[n]->assigned_edges) +"\n");
        #endif
        
        if (classify[n]->unassigned_in == 0 && classify[n]->unassigned_out == 0 && classify[n]->assigned_edges / 2 == max && classify[n]->components.size() == 1) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option A "+ std::to_string(gc.id(n)) +"\n");
            #endif
            //this is perfectly resolved! use ALL
            insert_local_known(n, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
            add_local_known(n, gc, temp_know_paths_add_max, temp_know_back_paths_add_max, kpc, kbpc);

            capacity_type max = 0;
            for ( std::deque<component>::iterator ci = classify[n]->components.begin(); ci != classify[n]->components.end(); ++ci) {
                capacity_type count = unravel_evidences_ILP(n, guiding_id, left_groups_ids[n], right_groups_ids[n], barred, ci->fragments, ci->nodes, gc, aic, fsc, nic, kpc, kbpc, uac, uic, meta->size, input_ids, trans);
                if (count > max) {
                    max = count;
                }
            }
            clean_ILP_leftovers(n, guiding_id, barred, max * low_barr_ratio, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);
            return true;
        }  
    }

    // we have a near perfect match
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
    
         if ( classify[n]->unassigned_in + classify[n]->assigned_in < 2 || classify[n]->unassigned_out + classify[n]->assigned_out < 2 ) {
            continue;
        }
        
        unsigned int edge_count = classify[n]->assigned_edges / 2 - classify[n]->edge_group_offset;
        unsigned int nodes = classify[n]->assigned_in + classify[n]->assigned_out - classify[n]->node_group_offset;
        
        if (edge_count + 2 <= 1 + nodes
                && classify[n]->unassigned_in == 0 && classify[n]->unassigned_out == 0 && classify[n]->components.size() <= 2) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option B "+ std::to_string(gc.id(n)) +"\n");
            #endif
            // one component and easy degree, also use ILP!
            insert_local_known(n, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
            add_local_known(n, gc, temp_know_paths_add_max, temp_know_back_paths_add_max, kpc, kbpc);
            
            capacity_type max = 0;
            for ( std::deque<component>::iterator ci = classify[n]->components.begin(); ci != classify[n]->components.end(); ++ci) {
                capacity_type count = unravel_evidences_ILP(n, guiding_id, left_groups_ids[n], right_groups_ids[n], barred, ci->fragments, ci->nodes, gc, aic, fsc, nic, kpc, kbpc, uac, uic, meta->size, input_ids, trans);
                if (count > max) {
                    max = count;
                }
            }
            clean_ILP_leftovers(n, guiding_id, barred, max * low_barr_ratio, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);
            return true;
        }
        
    }
      
    { // we set these off for local variables, ugly but mehhhhh
        bool barr_edge = false;
        ListDigraph::Node node;
        ListDigraph::Arc barc;
        float dr = 0;
        for (std::deque<ListDigraph::Node>::iterator tn = top_order.begin(); tn != top_order.end(); ++tn) {
            ListDigraph::Node n(*tn);
            
            
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("CHECK BARR NODE"+ std::to_string(gc.id(n)) +"\n");
                #endif
            
            if ( classify[n]->unassigned_in + classify[n]->assigned_in < 2 || classify[n]->unassigned_out + classify[n]->assigned_out < 2 ) {
                continue;
            }
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("IN "+ std::to_string(gc.id(n)) +"\n");
            #endif
                
            ListDigraph::Arc nbarc;
            float ndr = 0;
            if (bar_smallest_edge(n, gc, fsc, guiding_id, barred, 0.33, nbarc, ndr, block_delete_high_res[n])) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Report "+ std::to_string(gc.id(nbarc)) + " " + std::to_string(ndr) +"\n");
                #endif
                if ( (ndr < dr || !barr_edge) ) {
                    dr = ndr;
                    barc = nbarc;
                    barr_edge = true;
                    node = n;
                }
            }
        }
        if (barr_edge) {
                ListDigraph::Node node2;
                if (gc.target(barc) == node) {
                    node2 = gc.source(barc);
                } else {
                    node2 = gc.target(barc);
                }
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Option C.1 "+ std::to_string(gc.id(node)) + " Arc " + std::to_string(gc.id(barc)) +"\n");
                #endif
                barred[barc] = true;
                uac[barc]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
                insert_local_known(node, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
                clean_barred_leftovers(node, barred, guiding_id, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);
                clean_barred_leftovers(node2, barred, guiding_id, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);

                return true;
        }
    }
             
    { // we set these off for local variables, ugly but mehhhhh
        bool barr_edge = false;
        ListDigraph::Node node;
        ListDigraph::Arc barc;
        float dr = 0;
        for (std::deque<ListDigraph::Node>::iterator tn = top_order.begin(); tn != top_order.end(); ++tn) {
            ListDigraph::Node n(*tn);
            
            if ( classify[n]->unassigned_in + classify[n]->assigned_in < 2 || classify[n]->unassigned_out + classify[n]->assigned_out < 2 ) {
                continue;
            }
            
            ListDigraph::Arc nbarc;
            float ndr = 0;
            if (bar_negligible_edge(n, gc, fsc, guiding_id, barred, 0.05, nbarc, ndr, block_delete_high_res[n])) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Report "+ std::to_string(gc.id(nbarc)) + " " + std::to_string(ndr) +"\n");
                #endif
                if ( (ndr < dr || !barr_edge) ) {
                    dr = ndr;
                    barc = nbarc;
                    barr_edge = true;
                    node = n;
                }
            }
        }
        if (barr_edge) {
                ListDigraph::Node node2;
                if (gc.target(barc) == node) {
                    node2 = gc.source(barc);
                } else {
                    node2 = gc.target(barc);
                }
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Option N.1 "+ std::to_string(gc.id(node)) + " Arc " + std::to_string(gc.id(barc)) +"\n");
                #endif
                barred[barc] = true;
                uac[barc]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
                insert_local_known(node, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
                clean_barred_leftovers(node, barred, guiding_id, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);
                clean_barred_leftovers(node2, barred, guiding_id, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);

                return true;
        }
    }
    
    
        // include near perfect multi component matches
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
        
        if ( classify[n]->unassigned_in + classify[n]->assigned_in < 2 || classify[n]->unassigned_out + classify[n]->assigned_out < 2 ) {
            continue;
        }
        
        unsigned int components = classify[n]->unassigned_in + classify[n]->unassigned_out + classify[n]->components.size();
        
        unsigned int edge_count = classify[n]->assigned_edges / 2 - classify[n]->edge_group_offset;
        unsigned int nodes = classify[n]->assigned_in + classify[n]->assigned_out - classify[n]->node_group_offset;

        if (edge_count + 2 * components <= 3 + classify[n]->unassigned_in + classify[n]->unassigned_out + nodes
            && components > 1 && classify[n]->assigned_edges > 0) {
            // multi component and easy degree, also use ILP!
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option D "+ std::to_string(gc.id(n)) +"\n");
            #endif
            insert_local_known(n, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
            add_local_known(n, gc, temp_know_paths_add_max, temp_know_back_paths_add_max, kpc, kbpc);
            
            capacity_type max = 0;
            for ( std::deque<component>::iterator ci = classify[n]->components.begin(); ci != classify[n]->components.end(); ++ci) {             
                capacity_type count = unravel_evidences_ILP(n, guiding_id, left_groups_ids[n], right_groups_ids[n], barred, ci->fragments, ci->nodes, gc, aic, fsc, nic, kpc, kbpc, uac, uic, meta->size, input_ids, trans);
                if (count > max) {
                    max = count;
                }
            }
            clean_ILP_leftovers(n, guiding_id, barred, max * low_barr_ratio, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);
            return true;
        }
    }

        // clean up bad multi unique
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
        
        if ( classify[n]->unassigned_in + classify[n]->assigned_in < 2 || classify[n]->unassigned_out + classify[n]->assigned_out < 2 ) {
            continue;
        }
        
        if (classify[n]->assigned_in > classify[n]->fully_connected_in && classify[n]->assigned_out > classify[n]->fully_connected_out && classify[n]->assigned_edges > classify[n]->fully_connected_edges ) {
            
            bool has_loose_ends = false;
            
            std::set<int> in_groups, out_groups;
            std::map<int, std::set<int> > rev_in_groups;
            std::map<int, std::set<int> > rev_out_groups;
            for (ListDigraph::InArcIt a(gc, n); a!=INVALID; ++a) {
                if (aic[a].edge_type == edge_types::HELPER) continue;
                int id = gc.id(a);
                in_groups.insert(left_groups_ids[n][id].id);
                rev_in_groups[left_groups_ids[n][id].id].insert(id);
            }
            for (ListDigraph::OutArcIt a(gc, n); a!=INVALID; ++a) {
                if (aic[a].edge_type == edge_types::HELPER) continue;
                int id = gc.id(a);
                out_groups.insert(right_groups_ids[n][id].id);
                rev_out_groups[right_groups_ids[n][id].id].insert(id);
            }
            for (std::set<int>::iterator si = in_groups.begin(); si != in_groups.end() && !has_loose_ends; ++si) {
                std::set<int> kg;
                for (std::set<int>::iterator fi = rev_in_groups[*si].begin(); fi != rev_in_groups[*si].end(); ++fi) {
                    for (arc_bridge::iterator bi = kpc[gc.arcFromId(*fi)].begin(); bi != kpc[gc.arcFromId(*fi)].end(); ++bi) {
                        kg.insert(right_groups_ids[n][bi->first].id);
                    }
                }
                if (kg.size() == 1) {
                    has_loose_ends = true;
                }
            }
            for (std::set<int>::iterator si = out_groups.begin(); si != out_groups.end() && !has_loose_ends; ++si) {
                std::set<int> kg;
                for (std::set<int>::iterator fi = rev_out_groups[*si].begin(); fi != rev_out_groups[*si].end(); ++fi) {
                    for (arc_back_bridge::iterator bi = kbpc[gc.arcFromId(*fi)].begin(); bi != kbpc[gc.arcFromId(*fi)].end(); ++bi) {
                        kg.insert(left_groups_ids[n][*bi].id);
                    }
                }
                if (kg.size() == 1) {
                    has_loose_ends = true;
                }
            }
            
            if (has_loose_ends) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Option F "+ std::to_string(gc.id(n)) +"\n");
                #endif
                insert_local_known(n, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
                
                unravel_evidences_groups(n, guiding_id, left_groups_ids[n], right_groups_ids[n], barred, gc, aic, fsc, nic, kpc, kbpc, uac, uic, meta->size, ip_ids, trans);

//                unravel_single(n, kpc, kbpc, gc, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, uac, uic);
                clean_barred_leftovers(n, barred, guiding_id, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);

                return true;
            }
        }
    }
        
    // we have a match although bad
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
        
         if ( classify[n]->unassigned_in + classify[n]->assigned_in < 2 || classify[n]->unassigned_out + classify[n]->assigned_out < 2 ) {
            continue;
        }
        
         #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option G TEST "+  std::to_string(gc.id(n)) + ": "+ std::to_string(classify[n]->unassigned_in) + " " + std::to_string(classify[n]->unassigned_out) + " " + std::to_string(classify[n]->components.size())  +"\n");
            #endif
         
        if (classify[n]->unassigned_in + classify[n]->unassigned_out + classify[n]->components.size() == 1) { // one component, although bad edgecount
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option G "+ std::to_string(gc.id(n)) +"\n");
            #endif
            // one component and easy degree, also use ILP!
            insert_local_known(n, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
            add_local_known(n, gc, temp_know_paths_add_max, temp_know_back_paths_add_max, kpc, kbpc);
            
            capacity_type max = 0;
            for ( std::deque<component>::iterator ci = classify[n]->components.begin(); ci != classify[n]->components.end(); ++ci) {
                capacity_type count = unravel_evidences_ILP(n, guiding_id, left_groups_ids[n], right_groups_ids[n], barred, ci->fragments, ci->nodes, gc, aic, fsc, nic, kpc, kbpc, uac, uic, meta->size, input_ids, trans);
                if (count > max) {
                    max = count;
                }
            }
            clean_ILP_leftovers(n, guiding_id, barred, max * low_barr_ratio, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);
            return true;
        }
    }
        
    { // we set these off for local variables, ugly but mehhhhh
        bool barr_edge = false;
        ListDigraph::Node node;
        ListDigraph::Arc barc;
        float dr = 0;
        for (std::deque<ListDigraph::Node>::iterator tn = top_order.begin(); tn != top_order.end(); ++tn) {
            ListDigraph::Node n(*tn);
        
            if ( classify[n]->unassigned_in + classify[n]->assigned_in < 2 || classify[n]->unassigned_out + classify[n]->assigned_out < 2 ) {
                continue;
            }
            
            ListDigraph::Arc nbarc;
            float ndr = 0;
            if (bar_negligible_edge(n, gc, fsc, guiding_id, barred, 0.75, nbarc, ndr, block_delete_high_res[n])) {
                if ((ndr < dr || !barr_edge)) { //&& !temp_know_paths[nbarc].is_evidenced_path() && !temp_know_back_paths[nbarc].is_evidenced_path() ) {
                    dr = ndr;
                    barc = nbarc;
                    barr_edge = true;
                    node = n;
                }
            }
        }
        if (barr_edge) {
            
                ListDigraph::Node node2;
                if (gc.target(barc) == node) {
                    node2 = gc.source(barc);
                } else {
                    node2 = gc.target(barc);
                }
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Option H "+ std::to_string(gc.id(node)) + " Arc " + std::to_string(gc.id(barc)) +"\n");
                #endif
                barred[barc] = true;
                uac[barc]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
                insert_local_known(node, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
                clean_barred_leftovers(node, barred, guiding_id, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);
                clean_barred_leftovers(node2, barred, guiding_id, gc, aic, fsc, kpc, kbpc, uac, uic, ip_ids, trans);
                return true;
        }
    } 
    
    for (ListDigraph::NodeIt n(gc); n != INVALID;++n) {
         insert_local_known(n, gc, temp_know_paths, temp_know_back_paths, kpc, kbpc);
    }

    return false;
}


void base_manager::compute_edge_groups(std::set<std::set<int> > &rawhits, std::deque<std::set<int> > &groups) {
    

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Compute Edgegroup\n");
    #endif

    std::deque<std::set<int> > hits;
    for(std::set<std::set<int> >::iterator ri = rawhits.begin(); ri != rawhits.end(); ++ri) {
        hits.push_back(*ri);
    }

    std::sort(hits.begin(), hits.end(), []( std::set<int> & a, std::set<int> & b) -> bool {return a.size() < b.size();});

    for (std::deque<std::set<int> >::iterator hi = hits.begin(); hi != hits.end(); ) {

        bool global_contained = false;
        for (std::deque<std::set<int> >::iterator hi2 = hi+1 ; hi2 != hits.end(); ++hi2) {
            bool contained = true;
            for (std::set<int>::iterator si = hi->begin(); si != hi->end(); ++si) {
                
                if ( hi2->find(*si) == hi2->end() ) {
                    contained = false;
                    break;
                }
            }
            if (contained) {
                global_contained = true;
                break;
            }
        }
        
        if (global_contained) {
            hi = hits.erase(hi); 
        } else {
            ++hi;
        }
    }
    
    unsigned int size = hits.size();
    
    for (unsigned int i = 0; i< size; ++i) {
        
        std::set<int> *it = &hits[i];
        
//        #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("-------- " + std::to_string(it->size()) + "\n");
//            for (std::set<int>::iterator g = it->begin(); g != it->end(); ++g) {  
//                logger::Instance()->debug(std::to_string(*g) + ", ");
//            }   
//            logger::Instance()->debug("\n");
//        #endif
        
        bool match = false;
        for(std::deque<std::set<int> >::iterator lgi = groups.begin(); lgi != groups.end();++lgi) {

            std::set<int> intersection;
            std::set_intersection(lgi->begin(), lgi->end(), it->begin(), it->end(), std::inserter(intersection,intersection.begin()));

//            #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("CEG " + std::to_string(match) + "\n");
//                for (std::set<int>::iterator g = lgi->begin(); g!= lgi->end(); ++g) {
//                    logger::Instance()->debug(std::to_string(*g) + ", ");
//                }   
//                logger::Instance()->debug("\nI");
//                for (std::set<int>::iterator g = intersection.begin(); g!= intersection.end(); ++g) {
//                    logger::Instance()->debug(std::to_string(*g) + ", ");
//                }   
//                logger::Instance()->debug("\n");
//            #endif
                
            if (intersection.empty()) {
//                logger::Instance()->debug("Empty\n");
                continue;
            } 

            // so there is an intersection
            if (intersection.size() == it->size() && intersection.size() == lgi->size()) {
                // they are identical!
                match = true;
                continue;
            }

            std::set<int> A;
            std::set_difference(lgi->begin(), lgi->end(), intersection.begin(), intersection.end(), std::inserter(A,A.begin()));

            std::set<int> B;
            std::set_difference(it->begin(), it->end(), intersection.begin(), intersection.end(), std::inserter(B,B.begin()));

            groups.erase(lgi);
            if (!A.empty()) {
                hits.push_back(A);
                ++size;
            }
            if (!B.empty()) {
                hits.push_back(B);
                ++size;
            }
            hits.push_back(intersection);
            ++size;

            match = true;
            break;
        }
        
        if (!match) {
//            logger::Instance()->debug("ADD\n");
            groups.push_back(*it);
        }
   
    }
}

bool base_manager::bar_negligible_edge(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<flow_series> &fs, int guiding_id, ListDigraph::ArcMap<bool> &barred, float ratio, ListDigraph::Arc &barc, float &value, std::unordered_set<int> &block_delete) {
        
    float max_left = 0;
    float sum_left = 0;
    float min_left = -1;
    ListDigraph::Arc larc;
    
    float max_right = 0;
    float sum_right = 0;
    float min_right = -1;
    ListDigraph::Arc rarc;
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (barred[a] ) {
            continue;
        }
        
        float score = fs[a].get_mean(guiding_id).compute_score();
        
        if (score > max_left) {
            max_left = score;
        }
        if ((score < min_left || min_left < 0) && block_delete.find(wc.id(a)) == block_delete.end()) {
            min_left = score;
            larc = a;
        }
        sum_left += score;
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (barred[a]) {
            continue;
        }

        float score = fs[a].get_mean(guiding_id).compute_score();
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("RA " + std::to_string(wc.id(a)) + " " + std::to_string(score) + " " + std::to_string(block_delete.find(wc.id(a)) == block_delete.end()) + "\n");
        #endif
        
        if (score > max_right) {
            max_right = score;
        }
        if ((score < min_right || min_right < 0) && block_delete.find(wc.id(a)) == block_delete.end()) {
            min_right = score;
            rarc = a;
        }
        sum_right += score;
    }
    
    bool left_bar = max_left != 0 && min_left > 0 && (float) min_left <= ratio * (float) max_left;
    bool right_bar = max_right != 0 && min_right > 0 && (float) min_right <= ratio * (float) max_right;  
    
    float lr = 1;
    if (sum_left != 0 && min_left > 0) lr = min_left / sum_left;
    float rr = 1;
    if (sum_right != 0 && min_right > 0) rr = min_right / sum_right;
    
//    bool left_bar = lr <= ratio;
//    bool right_bar = rr <= ratio; 
            
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Barr Arcs " + std::to_string(wc.id(node)) + " " + std::to_string(wc.id(larc)) + " " + std::to_string(wc.id(rarc)) + "\n");
    logger::Instance()->debug("Left Barr " + std::to_string(left_bar) + " " + std::to_string(min_left) + " " + std::to_string(lr) + "\n");
    logger::Instance()->debug("Right Barr " + std::to_string(right_bar) + " " + std::to_string(min_right) + " " + std::to_string(rr) + "\n");
    #endif
    
    if (left_bar && right_bar) {
        if (lr < rr) {
            barc = larc;
            value = lr;
        } else {
            barc = rarc;
            value = rr;
        }
        return true;
    }
    if (left_bar) {
        barc = larc;
        value = lr;
        return true;
    }
    if (right_bar) {
        barc = rarc;
        value = rr;
        return true;
    }
        
    return false;
}


bool base_manager::bar_smallest_edge(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<flow_series> &fs, int guiding_id, ListDigraph::ArcMap<bool> &barred, float ratio, ListDigraph::Arc &barc, float &value, std::unordered_set<int> &block_delete) {
        
    float max_left = 0;
    float sum_left = 0;
    float min_left = -1;
    ListDigraph::Arc larc;
    
    float max_right = 0;
    float sum_right = 0;
    float min_right = -1;
    ListDigraph::Arc rarc;
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (barred[a] ) {
            continue;
        }
        
        float score = fs[a].get_mean(guiding_id).compute_score();
        
        if (score > max_left) {
            max_left = score;
        }
        if ((score < min_left || min_left < 0)) {
            min_left = score;
            larc = a;
        }
        sum_left += score;
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (barred[a]) {
            continue;
        }

        float score = fs[a].get_mean(guiding_id).compute_score();
        
        if (score > max_right) {
            max_right = score;
        }
        if ((score < min_right || min_right < 0)) {
            min_right = score;
            rarc = a;
        }
        sum_right += score;
    }
    
    bool left_bar = max_left != 0 && min_left > 0 && (float) min_left <= ratio * (float) sum_left;
    bool right_bar = max_right != 0 && min_right > 0 && (float) min_right <= ratio * (float) sum_right;  
    
    float lr = 1;
    if (sum_left != 0 && min_left > 0) lr = min_left / max_left;
    float rr = 1;
    if (sum_right != 0 && min_right > 0) rr = min_right / max_right;
    
//    bool left_bar = lr <= ratio;
//    bool right_bar = rr <= ratio; 
            
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Barr Arcs " + std::to_string(wc.id(node)) + " " + std::to_string(wc.id(larc)) + " " + std::to_string(wc.id(rarc)) + "\n");
    logger::Instance()->debug("Left Barr " + std::to_string(left_bar) + " " + std::to_string(min_left) + " " + std::to_string(lr) + "\n");
    logger::Instance()->debug("Right Barr " + std::to_string(right_bar) + " " + std::to_string(min_right) + " " + std::to_string(rr) + "\n");
    #endif
    
    if (left_bar && right_bar) {
        if (lr < rr) {
            if (block_delete.find(wc.id(larc)) != block_delete.end()) {
                return false;
            }
            barc = larc;
            value = lr;
        } else {
            if (block_delete.find(wc.id(rarc)) != block_delete.end()) {
                return false;
            }
            barc = rarc;
            value = rr;
        }
        return true;
    }
    if (left_bar) {
        if (block_delete.find(wc.id(larc)) != block_delete.end()) {
            return false;
        }
        barc = larc;
        value = lr;
        return true;
    }
    if (right_bar) {
        if (block_delete.find(wc.id(rarc)) != block_delete.end()) {
            return false;
        }
        barc = rarc;
        value = rr;
        return true;
    }
        
    return false;
}

bool base_manager::bar_negligible_edge_vote(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<flow_series> &fs, std::set<int>& ip_ids, int guiding_id, ListDigraph::ArcMap<bool> &barred, float ratio, ListDigraph::Arc &barc, float &value, std::unordered_set<int> &block_delete) {
            
    std::map<int, int> vote_min_left;
    std::map<int, int> vote_min_right;
    for (std::set<int>::iterator iii = ip_ids.begin(); iii != ip_ids.end(); ++iii) {
        int id = *iii;
        
        float min_left = -1;
        float min_right = -1;
        ListDigraph::Arc larc;
        ListDigraph::Arc rarc;
        
        for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {

            if (barred[a] ) {
                continue;
            }

            float score = fs[a].get_mean(guiding_id).compute_score();
            if ((score < min_left || min_left < 0) && block_delete.find(wc.id(a)) == block_delete.end()) {
                min_left = score;
                larc = a;
            }
        }
        for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {

            if (barred[a]) {
                continue;
            }

            float score = fs[a].get_mean(guiding_id).compute_score();
            if ((score < min_right || min_right < 0) && block_delete.find(wc.id(a)) == block_delete.end()) {
                min_right = score;
                rarc = a;
            }
        }
        vote_min_left[id] = wc.id(larc);
        vote_min_right[id] = wc.id(rarc);
    }
    
    float max_left = 0;
    float sum_left = 0;
    float min_left = -1;
    ListDigraph::Arc larc = wc.arcFromId(most_common_mapped(vote_min_left.begin(), vote_min_left.end()));
    
    float max_right = 0;
    float sum_right = 0;
    float min_right = -1;
    ListDigraph::Arc rarc = wc.arcFromId(most_common_mapped(vote_min_right.begin(), vote_min_right.end()));
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (barred[a] ) {
            continue;
        }
        
        float score = fs[a].get_mean(guiding_id).compute_score();
        
        if (score > max_left) {
            max_left = score;
        }
        if (a == larc) {
            min_left = score;
        }
        sum_left += score;
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (barred[a]) {
            continue;
        }

        float score = fs[a].get_mean(guiding_id).compute_score();
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("RA " + std::to_string(wc.id(a)) + " " + std::to_string(score) + " " + std::to_string(block_delete.find(wc.id(a)) == block_delete.end()) + "\n");
        #endif
        
        if (score > max_right) {
            max_right = score;
        }
        if (a == rarc) {
            min_right = score;
        }
        sum_right += score;
    }
    
    bool left_bar = max_left != 0 && min_left > 0 && (float) min_left <= ratio * (float) max_left;
    bool right_bar = max_right != 0 && min_right > 0 && (float) min_right <= ratio * (float) max_right;  
    
    float lr = 1;
    if (sum_left != 0 && min_left > 0) lr = min_left / sum_left;
    float rr = 1;
    if (sum_right != 0 && min_right > 0) rr = min_right / sum_right;
    
//    bool left_bar = lr <= ratio;
//    bool right_bar = rr <= ratio; 
            
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Barr Arcs " + std::to_string(wc.id(node)) + " " + std::to_string(wc.id(larc)) + " " + std::to_string(wc.id(rarc)) + "\n");
    logger::Instance()->debug("Left Barr " + std::to_string(left_bar) + " " + std::to_string(min_left) + " " + std::to_string(lr) + "\n");
    logger::Instance()->debug("Right Barr " + std::to_string(right_bar) + " " + std::to_string(min_right) + " " + std::to_string(rr) + "\n");
    #endif
    
    if (left_bar && right_bar) {
        if (lr < rr) {
            barc = larc;
            value = lr;
        } else {
            barc = rarc;
            value = rr;
        }
        return true;
    }
    if (left_bar) {
        barc = larc;
        value = lr;
        return true;
    }
    if (right_bar) {
        barc = rarc;
        value = rr;
        return true;
    }
        
    return false;
}


bool base_manager::bar_smallest_edge_vote(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<flow_series> &fs, std::set<int>& ip_ids, int guiding_id, ListDigraph::ArcMap<bool> &barred, float ratio, ListDigraph::Arc &barc, float &value, std::unordered_set<int> &block_delete) {
       
    std::map<int, int> vote_min_left;
    std::map<int, int> vote_min_right;
    for (std::set<int>::iterator iii = ip_ids.begin(); iii != ip_ids.end(); ++iii) {
        int id = *iii;
        
        float min_left = -1;
        float min_right = -1;
        ListDigraph::Arc larc;
        ListDigraph::Arc rarc;
        
        for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {

            if (barred[a] ) {
                continue;
            }

            float score = fs[a].get_mean(guiding_id).compute_score();
            if ((score < min_left || min_left < 0)) {
                min_left = score;
                larc = a;
            }
        }
        for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {

            if (barred[a]) {
                continue;
            }

            float score = fs[a].get_mean(guiding_id).compute_score();
            if ((score < min_right || min_right < 0)) {
                min_right = score;
                rarc = a;
            }
        }
        vote_min_left[id] = wc.id(larc);
        vote_min_right[id] = wc.id(rarc);
    }
    
    float max_left = 0;
    float sum_left = 0;
    float min_left = -1;
    ListDigraph::Arc larc = wc.arcFromId(most_common_mapped(vote_min_left.begin(), vote_min_left.end()));;
    
    float max_right = 0;
    float sum_right = 0;
    float min_right = -1;
    ListDigraph::Arc rarc = wc.arcFromId(most_common_mapped(vote_min_right.begin(), vote_min_right.end()));;
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (barred[a] ) {
            continue;
        }
        
        float score = fs[a].get_mean(guiding_id).compute_score();
        
        if (score > max_left) {
            max_left = score;
        }
        if (a == larc) {
            min_left = score;
        }
        sum_left += score;
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (barred[a]) {
            continue;
        }

        float score = fs[a].get_mean(guiding_id).compute_score();
        
        if (score > max_right) {
            max_right = score;
        }
        if (a == rarc) {
            min_right = score;
        }
        sum_right += score;
    }
    
    bool left_bar = max_left != 0 && min_left > 0 && (float) min_left <= ratio * (float) sum_left;
    bool right_bar = max_right != 0 && min_right > 0 && (float) min_right <= ratio * (float) sum_right;  
    
    float lr = 1;
    if (sum_left != 0 && min_left > 0) lr = min_left / max_left;
    float rr = 1;
    if (sum_right != 0 && min_right > 0) rr = min_right / max_right;
    
//    bool left_bar = lr <= ratio;
//    bool right_bar = rr <= ratio; 
            
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Barr Arcs " + std::to_string(wc.id(node)) + " " + std::to_string(wc.id(larc)) + " " + std::to_string(wc.id(rarc)) + "\n");
    logger::Instance()->debug("Left Barr " + std::to_string(left_bar) + " " + std::to_string(min_left) + " " + std::to_string(lr) + "\n");
    logger::Instance()->debug("Right Barr " + std::to_string(right_bar) + " " + std::to_string(min_right) + " " + std::to_string(rr) + "\n");
    #endif
    
    if (left_bar && right_bar) {
        if (lr < rr) {
            if (block_delete.find(wc.id(larc)) != block_delete.end()) {
                return false;
            }
            barc = larc;
            value = lr;
        } else {
            if (block_delete.find(wc.id(rarc)) != block_delete.end()) {
                return false;
            }
            barc = rarc;
            value = rr;
        }
        return true;
    }
    if (left_bar) {
        if (block_delete.find(wc.id(larc)) != block_delete.end()) {
            return false;
        }
        barc = larc;
        value = lr;
        return true;
    }
    if (right_bar) {
        if (block_delete.find(wc.id(rarc)) != block_delete.end()) {
            return false;
        }
        barc = rarc;
        value = rr;
        return true;
    }
        
    return false;
}

 void base_manager::find_component_left(std::set<int> &group, int ia, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths) {
     
    group.insert(ia);
    ListDigraph::Arc a = wc.arcFromId(ia);
        
    for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
        int na = ab->first;
        if (group.find(na) == group.end()) {
            find_component_right(group, na, wc, know_paths, know_back_paths);
        }
    } 
 }
 
 void base_manager::find_component_right(std::set<int> &group, int ia, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths) {
    group.insert(ia);
    ListDigraph::Arc a = wc.arcFromId(ia);
     
    for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
        int na = *ab;
        if (group.find(na) == group.end()) {
            find_component_left(group, na, wc, know_paths, know_back_paths);
        }
    }
 }
    
ListDigraph::Arc base_manager::extract_path(ListDigraph::Arc left_arc, ListDigraph::Arc right_arc,
            gmap<int, capacity_type> &capacities, int guiding, 
            ListDigraph &wc,
            ListDigraph::ArcMap<arc_identifier> &ai, ListDigraph::ArcMap<flow_series> &fs,
             ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, std::set<int>& ip_ids, ATC &trans) {
    
    // the idea is to extract the flow until it is a full separate egde that can be removed without any further problems
    // this is one call on two edges
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract Path. " + std::to_string(wc.id(left_arc)) + " " + std::to_string(wc.id(right_arc))  + " \n");
    #endif

    // new we create a new edge out of the two (this is the extraction part, in the end we have a full s-t arc)
    ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));
    ai[new_arc].edge_type = edge_types::EXON;
        
    for (gmap<int, capacity_type>::iterator ci = capacities.begin(); ci != capacities.end(); ++ci) {
        int id = ci->first;
        capacity_type capacity = ci->second;
        
        capacity_type flow_left = fs[left_arc].get_flow(id);
        capacity_type flow_right = fs[right_arc].get_flow(id);
        
        capacity_mean c1 = fs[left_arc].get_mean(id);
        capacity_mean c2 = fs[right_arc].get_mean(id);

        if (flow_left != 0) {
            c1.reduce(capacity/ float(flow_left));
            fs[left_arc].get_mean(id).reduce((flow_left - capacity)/ float(flow_left));
        } else {
            c1.reduce(0);
        }
        if (flow_right != 0) {
            c2.reduce(capacity/ float(flow_right));
            fs[right_arc].get_mean(id).reduce((flow_right - capacity)/ float(flow_right));
        } else {
            c2.reduce(0);
        }
        
        c1.update(c2);

        // we know we can extract the capacity freely
        fs[left_arc].get_flow(id) -= capacity;
        fs[right_arc].get_flow(id) -= capacity;
        
        fs[new_arc].get_flow(id) = capacity;
        fs[new_arc].get_mean(id) = c1;
    
    }
    
    // we now need to transfer the know paths of the right and left side
    know_paths[new_arc].bridges.ref() = know_paths[right_arc].bridges.ref(); // deep copy 
   // know_back_paths[new_arc] = know_back_paths[left_arc]; // deep copy, if called right this is always empty!
    
    if (fs[left_arc].get_flow(guiding) != 0 && fs[right_arc].get_flow(guiding) != 0) {
        // we are done already, because nothing changed in the evidence of the edges!
        return new_arc;
    }
    
    // if we are here, one of the edges has been removed
    // first we remove all invalidated
    
    if (fs[left_arc].get_flow(guiding) == 0) {
        for (arc_bridge::iterator bi = know_paths[left_arc].begin(); bi != know_paths[left_arc].end(); ++bi) {
            ListDigraph::Arc a = wc.arcFromId(bi->first);
            know_back_paths[a].remove_id(wc.id(left_arc));
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("REMOVE BACK " + std::to_string(wc.id(a)) + " " + std::to_string(wc.id(left_arc))  + " \n");
            #endif
        }
    } 
    if (fs[right_arc].get_flow(guiding) == 0) {
        
        // we update for all right neighbouring, this is the trick why this extraction can work
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].replace_evidence(wc.id(right_arc), wc.id(new_arc));
        }
        
        know_paths[left_arc].remove_id(wc.id(right_arc));
        for (arc_back_bridge::iterator bi = know_back_paths[right_arc].begin(); bi != know_back_paths[right_arc].end(); ++bi) {
            ListDigraph::Arc a = wc.arcFromId(*bi);
            know_paths[a].remove_id(wc.id(right_arc));
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("REMOVE FORW " + std::to_string(wc.id(a)) + " " + std::to_string(wc.id(right_arc))  + " \n");
            #endif
           
        }
    } else {
        
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].add_evidence_if(wc.id(right_arc), wc.id(new_arc));
        }
    }
            
    // only now we can test for extraction
   
    unravel_evidences(wc.target(left_arc), guiding, wc, ai, fs, know_paths, know_back_paths, unsecurityArc, unsecurityId, ip_ids, trans);

    return new_arc;
}


void base_manager::insert_local_known(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &temp_know_paths, ListDigraph::ArcMap<arc_back_bridge> &temp_know_back_paths,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths) {
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        know_paths[a] = temp_know_paths[a]; // we explicitly overwrite!
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        know_back_paths[a] = temp_know_back_paths[a]; // we explicitly overwrite!
    }
    
}

void base_manager::add_local_known(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<arc_bridge> &temp_know_paths, ListDigraph::ArcMap<arc_back_bridge> &temp_know_back_paths,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths) {
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        for (arc_bridge::iterator bi = temp_know_paths[a].begin(); bi != temp_know_paths[a].end(); ++bi) {
            know_paths[a][bi->first] += bi->second;
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        for (arc_back_bridge::iterator bi = temp_know_back_paths[a].begin(); bi != temp_know_back_paths[a].end(); ++bi) {
            know_back_paths[a].bridges->insert(*bi);
        }
    }
}

void base_manager::compute_ratio(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            std::set<int> &component,
            ListDigraph &wc,
            ListDigraph::ArcMap<flow_series>& fsc, int guiding_id, float &error, capacity_type &flow) {
    
    std::deque<ListDigraph::Arc> inArc, outArc;

    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        if (component.find(wc.id(a)) != component.end()) {
            flow += fsc[a].get_flow(guiding_id);
            inArc.push_back(a);
        }
    }

    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        if (component.find(wc.id(a)) != component.end()) {
            flow += fsc[a].get_flow(guiding_id);
            outArc.push_back(a);
        }
    }
    
    Lp lp;
    std::map<std::pair< int, int > , Lp::Col> bi_edges;
    std::deque<Lp::Col> left_z;
    std::deque<Lp::Col> right_z;
    Lp::Expr optimize;
    
    for (std::deque<ListDigraph::Arc>::iterator a_it = inArc.begin(); a_it != inArc.end(); ++a_it) {
        ListDigraph::Arc a = *a_it;
        
        for (arc_bridge::iterator bi = know_paths[a].begin(); bi != know_paths[a].end(); ++bi) {
            ListDigraph::Arc right = wc.arcFromId(bi->first);
            bi_edges.insert(std::make_pair(std::make_pair(wc.id(a), wc.id(right)), lp.addCol()));
            
            lp.colLowerBound(bi_edges[std::make_pair(wc.id(a), wc.id(right))], 0);
            lp.colUpperBound(bi_edges[std::make_pair(wc.id(a), wc.id(right))], Lp::INF);
            
        }
    }

    for (std::deque<ListDigraph::Arc>::iterator a_it = inArc.begin(); a_it != inArc.end(); ++a_it) {
        ListDigraph::Arc left = *a_it;
        
        Lp::Expr abs1, abs2;
        abs1 += fsc[left].get_flow(guiding_id);
        for (arc_bridge::iterator bi = know_paths[left].begin(); bi != know_paths[left].end(); ++bi) {
            ListDigraph::Arc right = wc.arcFromId(bi->first);
            
            abs1 -= bi_edges[std::make_pair(wc.id(left), wc.id(right))];
            abs2 += bi_edges[std::make_pair(wc.id(left), wc.id(right))]; 
        }
        
        left_z.push_back(lp.addCol());
        lp.addRow(left_z.back() >= abs1);
        lp.addRow( abs2 <= fsc[left].get_flow(guiding_id)); 
        optimize += left_z.back();// * fc[left];
    }

    for (std::deque<ListDigraph::Arc>::iterator a_it = outArc.begin(); a_it != outArc.end(); ++a_it) {
        ListDigraph::Arc right = *a_it;
        
        Lp::Expr abs1, abs2;
        abs1 += fsc[right].get_flow(guiding_id);
        for (arc_back_bridge::iterator bi = know_back_paths[right].begin(); bi != know_back_paths[right].end(); ++bi) {
            ListDigraph::Arc left = wc.arcFromId(*bi);
            
            abs1 -= bi_edges[std::make_pair(wc.id(left), wc.id(right))];
            abs2 += bi_edges[std::make_pair(wc.id(left), wc.id(right))]; 
        }
        
        right_z.push_back(lp.addCol());
        lp.addRow(right_z.back() >= abs1);
        lp.addRow(abs2 <= fsc[right].get_flow(guiding_id)); 
        optimize += right_z.back();// * fc[right];
    }
        
    // Specify the objective function
    lp.min();
    lp.obj(optimize);
    
    // Solve the problem using the underlying LP solver
    lp.solve();
    // Print the results
    if (lp.primalType() != Lp::OPTIMAL) {
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("LP Error A " + std::to_string(wc.id(node)) + "\n");
        #endif
        error += flow;
        return;
    }
    
    for (std::deque<ListDigraph::Arc>::iterator a_it = inArc.begin(); a_it != inArc.end(); ++a_it) {
        ListDigraph::Arc left = *a_it;

        int xee_sum = 0; 
        for (arc_bridge::iterator bi = know_paths[left].begin(); bi != know_paths[left].end(); ++bi) {
            ListDigraph::Arc right = wc.arcFromId(bi->first);
            
            xee_sum += std::round(lp.primal(bi_edges[std::make_pair(wc.id(left), wc.id(right))]));
          
        }
        error += std::sqrt(fsc[left].get_flow(guiding_id) - xee_sum);
    }
    
    for (std::deque<ListDigraph::Arc>::iterator a_it = outArc.begin(); a_it != outArc.end(); ++a_it) {
        ListDigraph::Arc right = *a_it;
        
        int xee_sum = 0; 
        for (arc_back_bridge::iterator bi = know_back_paths[right].begin(); bi != know_back_paths[right].end(); ++bi) {
            ListDigraph::Arc left = wc.arcFromId(*bi);
            
            xee_sum += std::round(lp.primal(bi_edges[std::make_pair(wc.id(left), wc.id(right))]));      
        }
        error += std::sqrt(fsc[right].get_flow(guiding_id) - xee_sum);
    }  
     
}


void base_manager::compute_flow_two_step(ListDigraph::Node node, int guiding_id,
            std::map<std::pair< int, int >, capacity_type > & group_connections,
            std::map<int, evidence_group > & left_groups, std::map<int, evidence_group > & right_groups,
            std::map<int, std::set<int> > & rev_in_groups, std::map<int, std::set<int> > & rev_out_groups,
            std::deque< resolve_count_set > & hit_counter_all,  
            ListDigraph &wc,
            ListDigraph::ArcMap<flow_series>& fsc,  
            ListDigraph::ArcMap<arc_identifier> &aic,
            ListDigraph::NodeMap<unsigned int>& nic,
            ListDigraph::ArcMap<arc_bridge>& kpc,
            ListDigraph::ArcMap<arc_back_bridge>& kbpc) {
      
    
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Flow 2 Step " + std::to_string(guiding_id) + "\n");
    #endif
    
    capacity_type flow_reference_in = 0;
    capacity_type flow_reference_out = 0;
    capacity_type flow_evidence = 0;
    
    std::set<int> in_groups, out_groups;
    std::map<int, std::set<int> > fwd_in_groups, fwd_out_groups;
    
    for (std::map<std::pair< int, int >, capacity_type >::iterator ggi = group_connections.begin(); ggi != group_connections.end(); ++ggi) {
        int lg = ggi->first.first;
        int rg = ggi->first.second; 
        
        in_groups.insert(lg);
        fwd_in_groups[lg].insert(rg);
        out_groups.insert(rg);
        fwd_out_groups[rg].insert(lg);
    }
    
    for(std::set<int>::iterator gi = in_groups.begin(); gi != in_groups.end(); ++gi) {
        for (std::set<int>::iterator ri = rev_in_groups[*gi].begin(); ri != rev_in_groups[*gi].end(); ++ri) {
            flow_reference_in += fsc[wc.arcFromId(*ri)].get_flow(guiding_id);
        }
    }
    
    for(std::set<int>::iterator gi = out_groups.begin(); gi != out_groups.end(); ++gi) {
        for (std::set<int>::iterator ri = rev_out_groups[*gi].begin(); ri != rev_out_groups[*gi].end(); ++ri) {
            flow_reference_out += fsc[wc.arcFromId(*ri)].get_flow(guiding_id);
        }
    }
    
    capacity_type flow_reference = std::min(flow_reference_in, flow_reference_out);

     #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Flow Reference " + std::to_string(flow_reference) + "\n");
    #endif

    for (std::deque< resolve_count_set >::iterator ha = hit_counter_all.begin(); ha != hit_counter_all.end(); ++ha) {
        flow_evidence += ha->count;
    }
    
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Flow Evidence " + std::to_string(flow_evidence) + "\n");
    #endif
        
        // now vamp up the evidences to match!
    for (std::deque< resolve_count_set >::iterator ha = hit_counter_all.begin(); ha != hit_counter_all.end(); ++ha) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Previous " + std::to_string(ha->count) + "\n");
        logger::Instance()->debug("Left: ");
        for (std::set<int>::iterator it = ha->left.begin(); it != ha->left.end(); ++it) {
            logger::Instance()->debug(std::to_string(*it) + ", ");
        }
        logger::Instance()->debug("\n");
        logger::Instance()->debug("Right: ");
        for (std::set<int>::iterator it = ha->right.begin(); it != ha->right.end(); ++it) {
            logger::Instance()->debug(std::to_string(*it) + ", ");
        }
        logger::Instance()->debug("\n");
        #endif

        ha->corrected_count = ha->count * flow_reference/ (float) flow_evidence;
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Correct to " + std::to_string(ha->corrected_count) + "\n");
        #endif
    }
    
    Lp lp;
    std::map<std::pair< int, int > , Lp::Col> bi_edges;
    std::deque<Lp::Col> left_z;
    std::deque<Lp::Col> right_z;
    Lp::Expr optimize;
    
    for (std::map<std::pair< int, int >, capacity_type >::iterator ggi = group_connections.begin(); ggi != group_connections.end(); ++ggi) {
        int lg = ggi->first.first;
        int rg = ggi->first.second; 

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Bi " + std::to_string(lg) + " " + std::to_string(rg) + "\n");
        #endif

        bi_edges.insert(std::make_pair(std::make_pair(lg, rg), lp.addCol()));

        lp.colLowerBound(bi_edges[std::make_pair(lg, rg)], 0);
        lp.colUpperBound(bi_edges[std::make_pair(lg, rg)], Lp::INF);

    }

    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        
        capacity_type group_cap = 0;
        std::set<int> know_groups;        
        for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
            ListDigraph::Arc left = wc.arcFromId(*ri);
            group_cap += fsc[left].get_flow(guiding_id);
        }
        
        Lp::Expr sum, sub;
        sub += group_cap;
        for ( std::set<int>::iterator fig = fwd_in_groups[lg].begin(); fig != fwd_in_groups[lg].end(); ++fig) {
            int rg = *fig;
            sum += bi_edges[std::make_pair(lg, rg)]; 
            sub -= bi_edges[std::make_pair(lg, rg)];
        }
        
        left_z.push_back(lp.addCol());
        lp.addRow(left_z.back() >= sub);
        lp.addRow( sum <= group_cap);  
        optimize += left_z.back();
    }
    
    for (std::set<int>::iterator rg_it = out_groups.begin(); rg_it != out_groups.end(); ++rg_it) {
        int rg = *rg_it;
        
        capacity_type group_cap = 0;
        std::set<int> know_groups;
        for (std::set<int>::iterator ri = rev_out_groups[rg].begin(); ri != rev_out_groups[rg].end(); ++ri) {
            ListDigraph::Arc right = wc.arcFromId(*ri);
            group_cap += fsc[right].get_flow(guiding_id);
        }
        
        Lp::Expr sum, sub;
        sub += group_cap;
        for ( std::set<int>::iterator fig = fwd_out_groups[rg].begin(); fig != fwd_out_groups[rg].end(); ++fig) {
            int lg = *fig;
            sum += bi_edges[std::make_pair(lg, rg)]; 
            sub -= bi_edges[std::make_pair(lg, rg)];
        }
        
        left_z.push_back(lp.addCol());
        lp.addRow(left_z.back() >= sub);
        lp.addRow( sum <= group_cap);  
        optimize += left_z.back();
    }
    
    Lp::Expr optimize2;
    std::deque<Lp::Col> opts;
    std::deque<Lp::Col> constraints;

    std::map<std::pair< int, int > , Lp::Expr> constraint_map_plus;
    std::map<std::pair< int, int > , Lp::Expr> constraint_map_minus;
    for (std::deque< resolve_count_set >::iterator ha = hit_counter_all.begin(); ha != hit_counter_all.end(); ++ha) {
        
        std::set<int> lg_set, rg_set;
        for (std::set<int>::iterator li = ha->left.begin(); li != ha->left.end(); ++li) {
            lg_set.insert(left_groups[*li].id);
        }
        for (std::set<int>::iterator ri = ha->right.begin(); ri != ha->right.end(); ++ri) {
            rg_set.insert(right_groups[*ri].id);
        }        
        
//        logger::Instance()->debug("AllCombo START\n");
        Lp::Expr all_combos = 0;
        bool existing_edge = false;
        for (std::set<int>::iterator li = lg_set.begin(); li != lg_set.end(); ++li) {
            for (std::set<int>::iterator ri = rg_set.begin(); ri != rg_set.end(); ++ri) {
                
                if (bi_edges.find(std::make_pair(*li,*ri)) == bi_edges.end()) { // we only add in ACTUAL valid edges!
                    continue;
                }
                
                existing_edge = true;
                
                constraints.push_back(lp.addCol());
                lp.colLowerBound(constraints.back(), 0);
                lp.colUpperBound(constraints.back(), ha->corrected_count);
                all_combos += constraints.back();
                
                if (constraint_map_plus.find(std::make_pair(*li,*ri)) ==  constraint_map_plus.end()) {
                    constraint_map_plus[std::make_pair(*li,*ri)] = 0;
                }
                if (constraint_map_minus.find(std::make_pair(*li,*ri)) ==  constraint_map_minus.end()) {
                    constraint_map_minus[std::make_pair(*li,*ri)] = bi_edges[std::make_pair(*li,*ri)];
                }
                    
                constraint_map_plus[std::make_pair(*li,*ri)] += constraints.back();
                constraint_map_minus[std::make_pair(*li,*ri)] -= constraints.back(); 

//                logger::Instance()->debug("C " + std::to_string(*li) + " " + std::to_string(*ri) + "\n");
            }
        }
        
//        logger::Instance()->debug("AllCombo " + std::to_string(ha->count) + "\n");
        if (existing_edge) {
            lp.addRow(all_combos == ha->corrected_count);
        }
    }
    
    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;

       for (std::set<int>::iterator fig = fwd_in_groups[lg].begin(); fig != fwd_in_groups[lg].end(); ++fig) {
            int rg = *fig;
            
            opts.push_back(lp.addCol());
            lp.addRow(opts.back() >= constraint_map_plus[std::make_pair(lg, rg)] - bi_edges[std::make_pair(lg, rg)]);
            lp.addRow(opts.back() >= constraint_map_minus[std::make_pair(lg, rg)]); 
            optimize2 += opts.back(); 
        }
    }
    
    // compute second solution
    lp.min();
    lp.obj(optimize2);
    
    // Solve the problem using the underlying LP solver
    lp.solve();
    // Print the results
    if (lp.primalType() != Lp::OPTIMAL) {
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("LP Error A " + std::to_string(wc.id(node)) + "\n");
        #endif
        return;
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Costs " + std::to_string(lp.primal()) + " as " + std::to_string(lp.primal(optimize2)) + "\n");
    for (std::deque<Lp::Col>::iterator it = constraints.begin(); it != constraints.end(); ++it) {
         logger::Instance()->debug("Con " + std::to_string(lp.primal(*it)) + "\n");
    }
    std::deque<Lp::Col>::iterator opi = opts.begin();
    for (std::map<std::pair< int, int > , Lp::Col>::iterator it = bi_edges.begin(); it != bi_edges.end(); ++it, ++opi) {
         logger::Instance()->debug("XX1 " + std::to_string(it->first.first) + "-" + std::to_string(it->first.second) + " " + std::to_string(lp.primal(it->second)) + " W " + std::to_string(lp.primal(constraint_map_plus[std::make_pair(it->first.first, it->first.second)]))  + "\n");
    }
    logger::Instance()->debug("LP A Done\n");
    #endif

    // after optimizing constraints we optimize flow in a second step
    
    
    Lp lp2;
    std::map<std::pair< int, int > , Lp::Col> bi_edges2;
    std::deque<Lp::Col> left_z2;
    std::deque<Lp::Col> right_z2;
    Lp::Expr optimize3;
    
    for (std::map<std::pair< int, int >, capacity_type >::iterator ggi = group_connections.begin(); ggi != group_connections.end(); ++ggi) {
        int lg = ggi->first.first;
        int rg = ggi->first.second; 
            
        bi_edges2.insert(std::make_pair(std::make_pair(lg, rg), lp2.addCol()));

        lp2.colLowerBound(bi_edges2[std::make_pair(lg, rg)], lp.primal(bi_edges[std::make_pair(lg, rg)]));
        lp2.colUpperBound(bi_edges2[std::make_pair(lg, rg)], Lp::INF);

    }

    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        
        capacity_type group_cap = 0;
        for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
            ListDigraph::Arc left = wc.arcFromId(*ri);
            group_cap += fsc[left].get_flow(guiding_id);
        }
        
        Lp::Expr abs1, abs2;
        abs1 += group_cap;
        for ( std::set<int>::iterator fig = fwd_in_groups[lg].begin(); fig != fwd_in_groups[lg].end(); ++fig) {
            int rg = *fig;
            
            abs1 -= bi_edges2[std::make_pair(lg, rg)];
            abs2 += bi_edges2[std::make_pair(lg, rg)]; 
        }
        
        left_z2.push_back(lp2.addCol());
        lp2.addRow(left_z2.back() >= abs1);
        lp2.addRow( abs2 <= group_cap); 
        optimize3 += left_z2.back();
    }
    
    for (std::set<int>::iterator rg_it = out_groups.begin(); rg_it != out_groups.end(); ++rg_it) {
        int rg = *rg_it;
        
        capacity_type group_cap = 0;
        for (std::set<int>::iterator ri = rev_out_groups[rg].begin(); ri != rev_out_groups[rg].end(); ++ri) {
            ListDigraph::Arc right = wc.arcFromId(*ri);
            group_cap += fsc[right].get_flow(guiding_id);
        }
        
        Lp::Expr abs1, abs2;
        abs1 += group_cap;
        for ( std::set<int>::iterator fig = fwd_out_groups[rg].begin(); fig != fwd_out_groups[rg].end(); ++fig) {
            int lg = *fig;
            
            abs1 -= bi_edges2[std::make_pair(lg, rg)];
            abs2 += bi_edges2[std::make_pair(lg, rg)]; 
        }
        
        left_z2.push_back(lp2.addCol());
        lp2.addRow(left_z2.back() >= abs1);
        lp2.addRow( abs2 <= group_cap); 
        optimize3 += left_z2.back();
    }
    
        
    lp2.min();
    lp2.obj(optimize3);
    
    // Solve the problem using the underlying LP solver
    lp2.solve();
    // Print the results
    if (lp2.primalType() != Lp::OPTIMAL) {
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("LP Error B " + std::to_string(wc.id(node)) + "\n");
        #endif
        return;
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Costs " + std::to_string(lp2.primal()) + " as " + std::to_string(lp2.primal(optimize)) + "\n");
    for (std::map<std::pair< int, int > , Lp::Col>::iterator it = bi_edges2.begin(); it != bi_edges2.end(); ++it) {
         logger::Instance()->debug("XX2 " + std::to_string(it->first.first) + "-" + std::to_string(it->first.second) + " " + std::to_string(lp2.primal(it->second)) + "\n");
    }
    logger::Instance()->debug("LP B Done\n");
    #endif

    for (std::map<std::pair< int, int >, capacity_type >::iterator ggi = group_connections.begin(); ggi != group_connections.end(); ++ggi) {
        int lg = ggi->first.first;
        int rg = ggi->first.second; 
        
        ggi->second = lp2.primal(bi_edges2[std::make_pair(lg, rg)]) + 0.001;
    }
    
}



capacity_type base_manager::unravel_evidences_ILP(ListDigraph::Node node, int guiding_id,
            std::map<int, evidence_group> & left_groups, std::map<int, evidence_group> & right_groups, ListDigraph::ArcMap<bool> &barred,
            std::deque< resolve_count_set > &hit_counter_all,
            std::set<int> &component,
            ListDigraph &wc,
            ListDigraph::ArcMap<arc_identifier> &aic,
            ListDigraph::ArcMap<flow_series>& fsc,  
            ListDigraph::NodeMap<unsigned int>& nic,
            ListDigraph::ArcMap<arc_bridge>& kpc,
            ListDigraph::ArcMap<arc_back_bridge>& kbpc,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, unsigned int size,
            std::set<int> &input_ids, ATC &transcripts) {

    #ifdef ALLOW_DEBUG
      logger::Instance()->debug("Unravel Evidences ILP " + std::to_string(guiding_id) + "\n");
      logger::Instance()->debug("Component ");
      for (std::set<int>::iterator ci = component.begin(); ci != component.end(); ++ci) {
          logger::Instance()->debug(std::to_string(*ci) + ", ");
      }
      logger::Instance()->debug("\n");
    #endif
    
    using namespace lemon;

    std::deque<ListDigraph::Arc> inArc;
    
    std::set<int> in_groups, out_groups;
    std::map<int, std::set<int> > rev_in_groups;
    std::map<int, std::set<int> > rev_out_groups;
    std::map<int, bool > rev_in_groups_ev;
    std::map<int, bool > rev_out_groups_ev;

    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        int id = wc.id(a);
        if (component.find(id) != component.end()) {
            inArc.push_back(a);
            in_groups.insert(left_groups[id].id);
            rev_in_groups[left_groups[id].id].insert(id);
            rev_in_groups_ev[left_groups[id].id] = left_groups[id].block;
        }
    }

    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        int id = wc.id(a);
        if (component.find(id) != component.end()) { 
            out_groups.insert(right_groups[id].id);
            rev_out_groups[right_groups[id].id].insert(id);
            rev_out_groups_ev[right_groups[id].id] = right_groups[id].block;
        }
    }
    
    typedef std::map<std::pair< int, int >, capacity_type > GCM;
    typedef std::map<int, GCM> IGCM;

    IGCM group_connections; 
    for (std::deque<ListDigraph::Arc>::iterator a_it = inArc.begin(); a_it != inArc.end(); ++a_it) {
        ListDigraph::Arc left = *a_it;
        
        for (arc_bridge::iterator bi = kpc[left].begin(); bi != kpc[left].end(); ++bi) {
            ListDigraph::Arc right = wc.arcFromId(bi->first);
            
            int lg = left_groups[wc.id(left)].id;
            int rg = right_groups[wc.id(right)].id;
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Bi " + std::to_string(wc.id(left)) + " " + std::to_string(wc.id(right)) + " : " + std::to_string(lg) + " " + std::to_string(rg) + "\n");
             #endif

            if ( group_connections[guiding_id].find(std::make_pair(lg, rg)) != group_connections[guiding_id].end() ) {
                continue; // only set up this variable once
            }
            
            group_connections[guiding_id].insert(std::make_pair(std::make_pair(lg, rg), 0));
        }
    }

   
    /// -------------
    
    compute_flow_two_step(node, guiding_id,
            group_connections[guiding_id],
            left_groups, right_groups,
            rev_in_groups, rev_out_groups,
            hit_counter_all,  
            wc,
            fsc,  
            aic,
            nic,
            kpc,
            kbpc); 
    
     /// -------------
    
    
    GCM group_connection_others;
    std::deque< resolve_count_set > hit_counter_others;
    for(typename GCM::iterator gci = group_connections[guiding_id].begin();  gci != group_connections[guiding_id].end(); ++gci) {
        if (gci->second > 0) {
            group_connection_others.insert(std::make_pair(gci->first, 0));
            hit_counter_others.push_back(resolve_count_set());
            hit_counter_others.back().count = gci->second;
            hit_counter_others.back().left.insert(*rev_in_groups[gci->first.first].begin());
            hit_counter_others.back().right.insert(*rev_out_groups[gci->first.second].begin());
        }
    }
    for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
        if (*iii != guiding_id) {
            group_connections[*iii] = group_connection_others;
            compute_flow_two_step(node, *iii,
                group_connections[*iii],
                left_groups, right_groups,
                rev_in_groups, rev_out_groups,
                hit_counter_others,  
                wc, fsc, aic, nic, kpc, kbpc);
        }
    }    

    // we extract exactly with them!
    
    // move around blocked nodes
    std::map<int, ListDigraph::Node > new_group_nodes_left;
    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        if (rev_in_groups_ev[lg]) {
            new_group_nodes_left[lg] = wc.addNode();
            nic[new_group_nodes_left[lg]] = nic[node];
            unsecurityId[new_group_nodes_left[lg]] = unsecurityId[node];
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Migrate Group " + std::to_string(lg) + " to new Group " + std::to_string(wc.id(new_group_nodes_left[lg])) + "\n");
            #endif

            //move over the edges
            for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
                wc.changeTarget(wc.arcFromId(*ri), new_group_nodes_left[lg]);
            }
        }
    }
    std::map<int, ListDigraph::Node > new_group_nodes_right;
    for (std::set<int>::iterator rg_it = out_groups.begin(); rg_it != out_groups.end(); ++rg_it) {
        int rg = *rg_it;
        if (rev_out_groups_ev[rg]) {
            new_group_nodes_right[rg] = wc.addNode();
            nic[new_group_nodes_right[rg]] = nic[node];
            unsecurityId[new_group_nodes_right[rg]] = unsecurityId[node];
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Migrate Group " + std::to_string(rg) + " to new Group " + std::to_string(wc.id(new_group_nodes_right[rg])) + "\n");
            #endif
            
            //move over the edges
            for (std::set<int>::iterator ri = rev_out_groups[rg].begin(); ri != rev_out_groups[rg].end(); ++ri) {
                wc.changeSource(wc.arcFromId(*ri), new_group_nodes_right[rg]);
            }
        }
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract Groups, guiding " + std::to_string(guiding_id) + "\n");
    #endif

    capacity_type max_cap = 0;
    for(typename GCM::iterator gci = group_connections[guiding_id].begin();  gci != group_connections[guiding_id].end(); ++gci) {
        int lg = gci->first.first;
        int rg = gci->first.second;
           
        capacity_type cap = gci->second;

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Group: " + std::to_string(lg) + "-" + std::to_string(rg) + ":"+ std::to_string(cap) + "\n");
        #endif
        
        if (cap == 0) {
            continue;
        }

        if (cap > max_cap) {
            max_cap = cap;
        }

         #ifdef ALLOW_DEBUG
         logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "\n");
         #endif

        // which groups do we resolve?
        if (rev_in_groups_ev[lg] && rev_out_groups_ev[rg]) {

            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Double Block \n");
            #endif

            ListDigraph::Arc new_arc = wc.addArc(new_group_nodes_left[lg], new_group_nodes_right[rg]);
            aic[new_arc].edge_type = edge_types::RESOLVE_HELPER;
            aic[new_arc].edge_specifier = exon_edge(size);

            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;

                fsc[new_arc].get_flow(id) = group_connections[id][gci->first];

                float mean_left = 0;
                float score_left = 0;
                capacity_type weight_left = 0;
                for (ListDigraph::InArcIt a(wc, new_group_nodes_left[lg]); a!=INVALID; ++a) {
                    if (!barred[a]) {
                        mean_left += fsc[a].get_flow(id) * fsc[a].get_mean(id).mean;
                        score_left += fsc[a].get_flow(id) * fsc[a].get_mean(id).compute_score();
                        weight_left += fsc[a].get_flow(id);
                    }
                }
                mean_left = mean_left / (float) weight_left;
                score_left = score_left / (float) weight_left;

                float mean_right = 0;
                float score_right = 0;
                capacity_type weight_right = 0;
                for (ListDigraph::OutArcIt a(wc, new_group_nodes_right[rg]); a!=INVALID; ++a) {
                    if (!barred[a]) {
                        mean_right += fsc[a].get_flow(id) * fsc[a].get_mean(id).mean;
                        score_right += fsc[a].get_flow(id) * fsc[a].get_mean(id).compute_score();
                        weight_right += fsc[a].get_flow(id);
                    }
                }
                mean_right = mean_right / (float) weight_right;
                score_right = score_right / (float) weight_right;

                float mean;
                float score;
                capacity_type weight;
                if (score_left > score_right) {
                    mean = mean_left;
                    score = score_left;
                    weight = weight_left;
                } else {
                    mean = mean_right;
                    score = score_right;
                    weight = weight_right;
                }

                float ratio = (float) group_connections[id][gci->first] / (float) weight;
                fsc[new_arc].get_mean(id).hidden_score = score * ratio;
                fsc[new_arc].get_mean(id).weight = 0;
                fsc[new_arc].get_mean(id).assign_mean(mean * ratio);
            }

        } else if (rev_in_groups_ev[lg] || rev_out_groups_ev[rg]) {

            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Single Block \n");
            #endif

            ListDigraph::Arc base_arc, new_arc;

            if (rev_in_groups_ev[lg]) {
                base_arc = wc.arcFromId(*rev_out_groups[rg].begin());
                new_arc = wc.addArc(new_group_nodes_left[lg], wc.target(base_arc));
            } else { // rev_out_groups_ev[rg].block
                base_arc = wc.arcFromId(*rev_in_groups[lg].begin());
                new_arc = wc.addArc(wc.source(base_arc), new_group_nodes_right[rg]);
            }

            // just copy over stuff   
            aic[new_arc] = aic[base_arc];

            if (unsecurityArc[base_arc].has_ref()) {
                std::copy(unsecurityArc[base_arc]->begin(), unsecurityArc[base_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
            }

            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;

                capacity_type capi = group_connections[id][gci->first];
                fsc[new_arc].get_flow(id) = capi;
                fsc[new_arc].get_mean(id) = fsc[base_arc].get_mean(id);
                fsc[new_arc].get_mean(id).reduce(capi/ float(fsc[base_arc].get_flow(id)));

                fsc[base_arc].get_mean(id).reduce((fsc[base_arc].get_flow(id) - capi)/ float(fsc[base_arc].get_flow(id)));
                fsc[base_arc].get_flow(id) -= capi;
            }

            if (rev_in_groups_ev[lg]) {
                kpc[new_arc].bridges.ref() = kpc[base_arc].bridges.ref(); // deep copy
                if (fsc[base_arc].get_flow(guiding_id) == 0) {
                    for (arc_bridge::iterator bi = kpc[base_arc].begin(); bi != kpc[base_arc].end(); ++bi) {
                        kbpc[wc.arcFromId(bi->first)].replace_evidence(wc.id(base_arc), wc.id(new_arc));
                    }
                    kpc[base_arc].clear();
                } else {
                    for (arc_bridge::iterator bi = kpc[base_arc].begin(); bi != kpc[base_arc].end(); ++bi) {
                        kbpc[wc.arcFromId(bi->first)].add_evidence_if(wc.id(base_arc), wc.id(new_arc));
                    }
                }
            } else {
                kbpc[new_arc].bridges.ref() = kbpc[base_arc].bridges.ref(); // deep copy
                if (fsc[base_arc].get_flow(guiding_id) == 0) {
                    for (arc_back_bridge::iterator bi = kbpc[base_arc].begin(); bi != kbpc[base_arc].end(); ++bi) {
                        kpc[wc.arcFromId(*bi)].replace_evidence(wc.id(base_arc), wc.id(new_arc));
                    }
                    kbpc[base_arc].clear();
                } else {
                     for (arc_back_bridge::iterator bi = kbpc[base_arc].begin(); bi != kbpc[base_arc].end(); ++bi) {
                        kpc[wc.arcFromId(*bi)].add_evidence_if(wc.id(base_arc), wc.id(new_arc));
                    }
                }
            }
        } else {
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Single Single \n");
            #endif

            int lga = *rev_in_groups[lg].begin();
            int rga = *rev_out_groups[rg].begin();    
            std::map<int,  capacity_type> caps;
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;

                caps.insert(std::make_pair(id, group_connections[id][gci->first]));
            }

            unravel_ILP(lga, rga, caps, guiding_id, wc, fsc, aic, kpc, kbpc, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);
        }
    }
      
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Add Over and underflow to original edge \n");
    #endif

    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_left.begin(); ngi != new_group_nodes_left.end(); ++ngi) {
        
        ListDigraph::Node n = ngi->second;
        
        struct cinfo {
            float mean = 0;
            float score = 0;
            capacity_type group_flow = 0;
            capacity_type evidenced_flow = 0;
        };
        std::map<int, cinfo> msm;
            
        for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
            int id = *iii;
        
            cinfo& ci = msm[id];
            capacity_type weight = 0;
            
            for (ListDigraph::InArcIt a(wc, n); a!=INVALID; ++a) {
                ci.group_flow += fsc[a].get_flow(id);
                if (!barred[a]) {
                    ci.mean += fsc[a].get_flow(id) * fsc[a].get_mean(id).mean;
                    ci.score += fsc[a].get_flow(id) * fsc[a].get_mean(id).compute_score();
                    weight += fsc[a].get_flow(id);
                }
            }
            ci.mean = ci.mean / (float) weight;
            ci.score = ci.score / (float) weight;
        
            for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
                ci.evidenced_flow += fsc[a].get_flow(id);
            }
        }
        
        if (msm[guiding_id].group_flow > msm[guiding_id].evidenced_flow) {
            ListDigraph::Arc new_arc = wc.addArc(n, node);
            aic[new_arc].edge_type = edge_types::RESOLVE_HELPER;
            aic[new_arc].edge_specifier = exon_edge(size);
            
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;
                cinfo& ci = msm[id];
                
                fsc[new_arc].get_flow(id) = ci.group_flow - ci.evidenced_flow;
                float ratio = (float) ci.evidenced_flow / (float) ci.group_flow;
                fsc[new_arc].get_mean(id).hidden_score = ci.score * ratio;
                fsc[new_arc].get_mean(id).weight = 0;
                fsc[new_arc].get_mean(id).assign_mean(ci.mean * ratio);
            }
        }
    }
  
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_right.begin(); ngi != new_group_nodes_right.end(); ++ngi) {
        
        ListDigraph::Node n = ngi->second;
        
        struct cinfo {
            float mean = 0;
            float score = 0;
            capacity_type group_flow = 0;
            capacity_type evidenced_flow = 0;
        };
        std::map<int, cinfo> msm;
        
        for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
            int id = *iii;
        
            cinfo& ci = msm[id];
            capacity_type weight = 0;
            
            for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
                ci.group_flow += fsc[a].get_flow(id);
                if (!barred[a]) {
                    ci.mean += fsc[a].get_flow(id) * fsc[a].get_mean(id).mean;
                    ci.score += fsc[a].get_flow(id) * fsc[a].get_mean(id).compute_score();
                    weight += fsc[a].get_flow(id);
                }
            }
            ci.mean = ci.mean / (float) weight;
            ci.score = ci.score / (float) weight;
        
            for (ListDigraph::InArcIt a(wc, n); a!=INVALID; ++a) {
                ci.evidenced_flow += fsc[a].get_flow(id);
            }
        }
        
        if (msm[guiding_id].group_flow > msm[guiding_id].evidenced_flow) {
            ListDigraph::Arc new_arc = wc.addArc(node, n);
            aic[new_arc].edge_type = edge_types::RESOLVE_HELPER;
            aic[new_arc].edge_specifier = exon_edge(size);
            
            
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;
                cinfo& ci = msm[id];
                
                fsc[new_arc].get_flow(id) = ci.group_flow - ci.evidenced_flow;
                float ratio = (float) ci.evidenced_flow / (float) ci.group_flow;
                fsc[new_arc].get_mean(id).hidden_score = ci.score * ratio;
                fsc[new_arc].get_mean(id).weight = 0;
                fsc[new_arc].get_mean(id).assign_mean( ci.mean * ratio);
            }
        }
    }
        
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Clean Leftovers New Nodes \n");
    #endif
   

    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        if (rev_in_groups_ev[lg]) {
            for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
                for (arc_bridge::iterator bi = kpc[wc.arcFromId(*ri)].begin(); bi != kpc[wc.arcFromId(*ri)].end(); ++bi) {
                    #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Remove " + std::to_string(*ri) + " AT " + std::to_string(bi->first) + "\n");
                    #endif
                    kbpc[wc.arcFromId(bi->first)].remove_id(*ri);
                }
                kpc[wc.arcFromId(*ri)].clear();
            }
        }
    }
    for (std::set<int>::iterator rg_it = out_groups.begin(); rg_it != out_groups.end(); ++rg_it) {
        int rg = *rg_it;
        if (rev_out_groups_ev[rg]) {
            for (std::set<int>::iterator ri = rev_out_groups[rg].begin(); ri != rev_out_groups[rg].end(); ++ri) {
                for (arc_back_bridge::iterator bi = kbpc[wc.arcFromId(*ri)].begin(); bi != kbpc[wc.arcFromId(*ri)].end(); ++bi) {
                    #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Remove " + std::to_string(*ri) + " AT " + std::to_string(*bi) + "\n");
                    #endif
                    kpc[wc.arcFromId(*bi)].remove_id(*ri);
                }
                kbpc[wc.arcFromId(*ri)].clear();
            }
        }
    }    
        
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_left.begin(); ngi != new_group_nodes_left.end(); ++ngi) {
        clean_barred_leftovers(ngi->second, barred, guiding_id, wc, aic, fsc, kpc, kbpc, unsecurityArc, unsecurityId, input_ids, transcripts);
    }
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_right.begin(); ngi != new_group_nodes_right.end(); ++ngi) {
        clean_barred_leftovers(ngi->second,barred, guiding_id, wc, aic, fsc, kpc, kbpc, unsecurityArc, unsecurityId, input_ids, transcripts);
    }
        
    // we now look at possibly unevidenced leftovers
    // remove all 0 edges here
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fsc[arc].get_flow(guiding_id) == 0) {
            for (arc_bridge::iterator bi = kpc[arc].begin(); bi != kpc[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(bi->first);
                kbpc[a].remove_id(wc.id(arc));
            }
            kpc[arc].clear();
            erase_arc(wc, arc, fsc, input_ids, transcripts);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fsc[arc].get_flow(guiding_id) == 0) {
            for (arc_back_bridge::iterator bi = kbpc[arc].begin(); bi != kbpc[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(*bi);
                kpc[a].remove_id(wc.id(arc));
            }  
            kbpc[arc].clear();
            erase_arc(wc, arc, fsc, input_ids, transcripts);
        }
    }      

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("FINISHED ILP\n"); 
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {

        logger::Instance()->debug("Arc " + std::to_string(wc.id(a))); 
        for (arc_bridge::iterator ab = kpc[a].begin(); ab != kpc[a].end(); ++ab) {
            logger::Instance()->debug(" f " + std::to_string(ab->first)); 
        }
        for (arc_back_bridge::iterator ab = kbpc[a].begin(); ab != kbpc[a].end(); ++ab) {
            logger::Instance()->debug(" b " + std::to_string(*ab)); 
        }
        logger::Instance()->debug("\n" );  
     }
     #endif    
        
    return max_cap;
}

void base_manager::clean_ILP_leftovers(ListDigraph::Node node, int guiding_id,
            ListDigraph::ArcMap<bool> &barred, capacity_type barr_limit,
            ListDigraph &wc,
            ListDigraph::ArcMap<arc_identifier> &ai, ListDigraph::ArcMap<flow_series> &fs,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, std::set<int>& ip_ids, ATC &trans) {
    
    std::set<int> nodes;
    nodes.insert(wc.id(node));
    
    bool has_left = false;
    bool has_right = false;
    
    for (ListDigraph::InArcIt in(wc, node); in != INVALID; ++in) {
        if (fs[in].get_flow(guiding_id) < barr_limit ) { //&& know_paths[in].size() > 0) {
            barred[in] = true;
            unsecurityArc[in]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));  
        }
        nodes.insert(wc.id(wc.source(in))); 
    }
    for (ListDigraph::OutArcIt out(wc, node); out != INVALID; ++out) {
        if (fs[out].get_flow(guiding_id) < barr_limit) { //&& know_back_paths[out].size() > 0) {
            barred[out] = true;
            unsecurityArc[out]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
        }
        nodes.insert(wc.id(wc.target(out)));
    }
   
    if (nodes.size() == 1) { 
        wc.erase(node);
    }
    
    for (std::set<int>::iterator in = nodes.begin(); in != nodes.end(); ++in) { 
        clean_barred_leftovers(wc.nodeFromId(*in), barred, guiding_id, wc, ai, fs, know_paths, know_back_paths, unsecurityArc, unsecurityId, ip_ids, trans);
    }
}

bool base_manager::clean_barred_leftovers(ListDigraph::Node node,
            ListDigraph::ArcMap<bool> &barred, int guiding_id,
            ListDigraph &wc,
            ListDigraph::ArcMap<arc_identifier> &aic,
            ListDigraph::ArcMap<flow_series>& fsc,   
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc, 
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId,
            std::set<int>& ip_ids, ATC &trans) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("clean barred leftovers\n");
    #endif
    
    std::deque<ListDigraph::Arc> inArc, outArc;
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        if (!barred[a]) {
            inArc.push_back(a);   
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        if (!barred[a]) {
            outArc.push_back(a); 
        }
    }
    
    if (inArc.size() > 1 && outArc.size() > 1 || inArc.size() == 0 || outArc.size() == 0) {
        return false;
    }
    
    std::map<int, capacity_type> in_flow, out_flow;
    for (std::set<int>::iterator iii = ip_ids.begin(); iii != ip_ids.end(); ++iii) {
        int id = *iii;
        for (std::deque<ListDigraph::Arc>::iterator a = outArc.begin(); a != outArc.end(); ++a) {
            out_flow[id] += fsc[*a].get_flow(id);
        }
        for (std::deque<ListDigraph::Arc>::iterator a = inArc.begin(); a != inArc.end(); ++a) {
            in_flow[id] += fsc[*a].get_flow(id);
        }
    }

    if (inArc.size() == 1) {
        
        ListDigraph::Arc left = inArc[0];
        std::unordered_map<int, std::map<int, capacity_type> > out_to_caps;
        
        for (std::set<int>::iterator iii = ip_ids.begin(); iii != ip_ids.end(); ++iii) {
            int id = *iii;
        
            float ratio = in_flow[id] / (float) out_flow[id];
            if (ratio > 1) ratio = 1;
        
            std::deque<capacity_type> flows;
            capacity_type total = 0;
            for (std::deque<ListDigraph::Arc>::iterator a = outArc.begin(); a != outArc.end(); ++a) {
                flows.push_back(fsc[*a].get_flow(id) * ratio);
                total += flows.back();
            } 
            capacity_type left_over = fsc[left].get_flow(id) - total;
        
            std::deque<capacity_type>::iterator fi = flows.begin();
            for (std::deque<ListDigraph::Arc>::iterator a = outArc.begin(); a != outArc.end(); ++a, ++fi) {
                capacity_type additional_flow = std::min(left_over, fsc[*a].get_flow(id) - *fi);
                left_over -= additional_flow;
                out_to_caps[wc.id(*a)][id] = *fi + additional_flow;
            }
        }

        for (std::deque<ListDigraph::Arc>::iterator a = outArc.begin(); a != outArc.end(); ++a) {
            unravel_ILP(wc.id(left), wc.id(*a), out_to_caps[wc.id(*a)], guiding_id, wc, fsc, aic, know_paths, know_back_paths, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);
            barred[*a] = true;
            unsecurityArc[*a]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
        } 
    } else {
        
        ListDigraph::Arc right = outArc[0];
        std::unordered_map<int, std::map<int, capacity_type> > in_to_caps;
        
        for (std::set<int>::iterator iii = ip_ids.begin(); iii != ip_ids.end(); ++iii) {
            int id = *iii;
        
            float ratio = out_flow[id] / (float) in_flow[id];
            if (ratio > 1) ratio = 1;
        
            std::deque<capacity_type> flows;
            capacity_type total = 0;
            for (std::deque<ListDigraph::Arc>::iterator a = inArc.begin(); a != inArc.end(); ++a) {
                flows.push_back(fsc[*a].get_flow(id) * ratio);
                total += flows.back();
            } 
            capacity_type left_over = fsc[right].get_flow(id) - total;
        
            std::deque<capacity_type>::iterator fi = flows.begin();
            for (std::deque<ListDigraph::Arc>::iterator a = inArc.begin(); a != inArc.end(); ++a, ++fi) {
                capacity_type additional_flow = std::min(left_over, fsc[*a].get_flow(id) - *fi);
                left_over -= additional_flow;
                in_to_caps[wc.id(*a)][id] = *fi + additional_flow;
            }
        }

        for (std::deque<ListDigraph::Arc>::iterator a = inArc.begin(); a != inArc.end(); ++a) {
            unravel_ILP(wc.id(*a), wc.id(right), in_to_caps[wc.id(*a)], guiding_id, wc, fsc, aic, know_paths, know_back_paths, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);
            barred[*a] = true;
            unsecurityArc[*a]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
        } 
    }
        
     // we now look at possibly unevidenced leftovers
    // remove all 0 edges here
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fsc[arc].get_flow(guiding_id) == 0) {
            for (arc_bridge::iterator bi = know_paths[arc].begin(); bi != know_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(bi->first);
                know_back_paths[a].remove_id(wc.id(arc));
            }
            know_paths[arc].clear();
            erase_arc(wc, arc, fsc, ip_ids, trans);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fsc[arc].get_flow(guiding_id) == 0) {
            for (arc_back_bridge::iterator bi = know_back_paths[arc].begin(); bi != know_back_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(*bi);
                know_paths[a].remove_id(wc.id(arc));
            }  
            know_back_paths[arc].clear();
            erase_arc(wc, arc, fsc, ip_ids, trans);
        }
    }
    
    bool ret = false;
    
    int indeg = 0;
    int outdeg = 0;
    ListDigraph::InArcIt in(wc, node);
    ListDigraph::OutArcIt out(wc, node);
    for (;in != INVALID; ++in) {
        ++indeg;
    }
    for (;out != INVALID; ++out) {
        ++outdeg;
    }
    if (indeg == 1 || outdeg == 1) {
        unravel_unevidenced_leftovers(node, guiding_id, wc, aic, fsc, know_paths, know_back_paths, unsecurityArc, unsecurityId, ip_ids, trans, &barred);
        indeg = 0;
        ret = true;
    }
    
    if (indeg == 0) {
        // we completely resolved this one, remove!
        wc.erase(node);
    }
    return ret;
    
}


void base_manager::unravel_ILP(int left, int right, std::map<int,  capacity_type> &caps, int guiding_id, 
            ListDigraph &wc,
            ListDigraph::ArcMap<flow_series>& fsc,   
            ListDigraph::ArcMap<arc_identifier> &aic,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc, 
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, transcript_unsecurity::evidence evidence) {
    
    
    if (caps[guiding_id] == 0) {
        return;
    }
    
    ListDigraph::Arc left_arc = wc.arcFromId(left);
    ListDigraph::Arc right_arc = wc.arcFromId(right);
    
    // this should never be called on a zero edge, so add in new arcs
    ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));
    aic[new_arc].edge_type = edge_types::EXON;
    
    for (std::map<int, capacity_type>::iterator ci = caps.begin(); ci != caps.end(); ++ci) {
        int id = ci->first;
        capacity_type capacity = ci->second;

        capacity_mean c1 = fsc[left_arc].get_mean(id);
        capacity_mean c2 = fsc[right_arc].get_mean(id);

        capacity_type flow_left = fsc[left_arc].get_flow(id);
        capacity_type flow_right = fsc[right_arc].get_flow(id);

        if (flow_left != 0) {
            c1.reduce(capacity/ float(flow_left));
            fsc[left_arc].get_mean(id).reduce((flow_left - capacity)/ float(flow_left));
        } else {
            c1.reduce(0);
        }
        if (flow_right != 0) {
            c2.reduce(capacity/ float(flow_right));
            fsc[right_arc].get_mean(id).reduce((flow_right - capacity)/ float(flow_right));
        } else {
            c2.reduce(0);
        }
        
        c1.update(c2);

        // we know we can extract the capacity freely
        fsc[left_arc].get_flow(id) -= capacity;
        fsc[right_arc].get_flow(id) -= capacity;
        
        fsc[new_arc].get_flow(id) = capacity;
        fsc[new_arc].get_mean(id) = c1;
    }

    if (aic[left_arc].edge_type != edge_types::EXON && aic[right_arc].edge_type != edge_types::EXON) {
        aic[new_arc].edge_type = edge_types::HELPER;
        if (aic[right_arc].edge_specifier.node_index == -1) {
            aic[new_arc].edge_specifier = aic[left_arc].edge_specifier;
        } else {
            aic[new_arc].edge_specifier = aic[right_arc].edge_specifier;
        }
    } else if (aic[left_arc].edge_type == edge_types::HELPER) {
        // leftarc is a helper
        aic[new_arc].edge_specifier = exon_edge(aic[right_arc].edge_specifier); 
        
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else if (aic[right_arc].edge_type == edge_types::HELPER) {
        // rightarc is a helper
        aic[new_arc].edge_specifier = exon_edge(aic[left_arc].edge_specifier);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else {
        aic[new_arc].edge_specifier = exon_edge(aic[left_arc].edge_specifier);
        aic[new_arc].edge_specifier.join_edge(aic[right_arc].edge_specifier);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->begin()));
        }
        unsecurityArc[new_arc]->insert( transcript_unsecurity( unsecurityId[wc.target(left_arc)].id, evidence));
    }
    
    know_paths[new_arc].bridges.ref() = know_paths[right_arc].bridges.ref(); // deep copy
    know_back_paths[new_arc].bridges.ref() = know_back_paths[left_arc].bridges.ref(); // deep copy
    
    if (fsc[left_arc].get_flow(guiding_id) == 0) {
        
        for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
            know_paths[a].replace_evidence(wc.id(left_arc), wc.id(new_arc));
        }
          
    } else {
        
        for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
            know_paths[a].add_evidence_if(wc.id(left_arc), wc.id(new_arc));
        }   
    }
    
    if (fsc[right_arc].get_flow(guiding_id) == 0) {
        
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].replace_evidence(wc.id(right_arc), wc.id(new_arc));
        }
          
    } else {
        
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].add_evidence_if(wc.id(right_arc), wc.id(new_arc));
        }
    }
     
    aic[new_arc].edge_lengths.first_exon = aic[left_arc].edge_lengths.first_exon;
    aic[new_arc].edge_lengths.middle = aic[left_arc].edge_lengths.middle + aic[left_arc].edge_lengths.last_exon + aic[right_arc].edge_lengths.middle;
    aic[new_arc].edge_lengths.last_exon = aic[right_arc].edge_lengths.last_exon;
    
    aic[new_arc].cycle_id_in = aic[left_arc].cycle_id_in;
    aic[new_arc].cycle_id_out = aic[right_arc].cycle_id_out;
        
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Exit ILP.\n");
    #endif

}

void base_manager::unravel_evidences_groups(ListDigraph::Node node, int guiding_id,
            std::map<int, evidence_group> & left_groups, std::map<int, evidence_group> & right_groups, ListDigraph::ArcMap<bool> &barred,
            ListDigraph &wc,
            ListDigraph::ArcMap<arc_identifier> &aic,
            ListDigraph::ArcMap<flow_series>& fsc, 
            ListDigraph::NodeMap<unsigned int> &nic,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, unsigned int size,
            std::set<int> &input_ids, ATC &transcripts) {
    
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Unravel NODE GROUP " + std::to_string(wc.id(node)) + "\n");
    #endif
   
    struct info {
        
        info() {}
        info(int group, bool is_out, std::map<int, capacity_type> &cap, std::set<int> &s) : group(group), is_out(is_out), cap(cap), know_groups(s) 
            {}
        
        int group;
        bool is_out = false;
        std::map<int, capacity_type> cap;
        std::set<int> know_groups;
    };
    
    std::set<int> in_groups, out_groups;
    std::map<int, std::set<int> > rev_in_groups;
    std::map<int, std::set<int> > rev_out_groups;
    std::map<int, bool > rev_in_groups_ev;
    std::map<int, bool> rev_out_groups_ev;
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (aic[a].edge_type == edge_types::HELPER) {
            // no resolving 
            continue;
        }
        
        int id = wc.id(a);
        in_groups.insert(left_groups[id].id);
        rev_in_groups[left_groups[id].id].insert(id);
        rev_in_groups_ev[left_groups[id].id] = left_groups[id].block;
    }

    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (aic[a].edge_type == edge_types::HELPER) {
            // no resolving 
            continue;
        }
        
        int id = wc.id(a);
        out_groups.insert(right_groups[id].id);
        rev_out_groups[right_groups[id].id].insert(id);
        rev_out_groups_ev[right_groups[id].id] = right_groups[id].block;
    }

    std::deque<info> groups;
    std::map<int, info* > rev_groups_info;
    for (std::set<int>::iterator si = in_groups.begin(); si != in_groups.end(); ++si) {
        std::map<int, capacity_type> group_flow;
        std::set<int> kg;
        for (std::set<int>::iterator fi = rev_in_groups[*si].begin(); fi != rev_in_groups[*si].end(); ++fi) {
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                group_flow[*iii] += fsc[wc.arcFromId(*fi)].get_flow(*iii);
            }
            for (arc_bridge::iterator bi = know_paths[wc.arcFromId(*fi)].begin(); bi != know_paths[wc.arcFromId(*fi)].end(); ++bi) {
                kg.insert(right_groups[bi->first].id);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("FG " + std::to_string(right_groups[bi->first].id) + "\n");
                #endif
            }
        }
        
        #ifdef ALLOW_DEBUG
         logger::Instance()->debug("FGT " + std::to_string(*si) + "\n");
        #endif

        groups.push_back(info(*si, false, group_flow, kg));
        rev_groups_info[*si] = &groups.back();
    }
    for (std::set<int>::iterator si = out_groups.begin(); si != out_groups.end(); ++si) {
        std::map<int, capacity_type> group_flow;
        std::set<int> kg;
        for (std::set<int>::iterator fi = rev_out_groups[*si].begin(); fi != rev_out_groups[*si].end(); ++fi) {
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                group_flow[*iii] += fsc[wc.arcFromId(*fi)].get_flow(*iii);
            }
            for (arc_back_bridge::iterator bi = know_back_paths[wc.arcFromId(*fi)].begin(); bi != know_back_paths[wc.arcFromId(*fi)].end(); ++bi) {
                kg.insert(left_groups[*bi].id);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("FGB " + std::to_string(left_groups[*bi].id) + "\n");
                #endif
            }
        }
        
        #ifdef ALLOW_DEBUG
         logger::Instance()->debug("FGTB " + std::to_string(*si) + "\n");
        #endif
        
        groups.push_back(info(*si, true, group_flow, kg));
        rev_groups_info[*si] = &groups.back();
    }

    std::map<int, ListDigraph::Node > new_group_nodes_left, new_group_nodes_right;
       
    bool change = true;
    while (change) {
        change = false;
        // loop over collected
        
        info* candidate;
        capacity_type max = 0;
        
        for (std::deque< info>::iterator it = groups.begin(); it != groups.end(); ++it) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Group: " + std::to_string(it->group) + "\n");
            for (std::set<int>::iterator it2 = it->know_groups.begin(); it2 != it->know_groups.end(); ++it2) {
                logger::Instance()->debug(std::to_string(*it2) + ", ");
            } 
            logger::Instance()->debug("\n");
            #endif
            
            if (it->know_groups.size() == 1) {
                
                capacity_type min = std::min(it->cap[guiding_id], rev_groups_info[*it->know_groups.begin()]->cap[guiding_id]);
                if (min > max ){
                    candidate = &*it;
                    max = min;
                }
            }
        }
        
        if (max == 0) break;
        
        // we have something to resolve, yay
        change = true;
        int lg, rg, del, del_in;
       
        if (candidate->is_out) {
            rg = candidate->group;
            del = rg;
            lg = *candidate->know_groups.begin();
            del_in = lg;
        } else {
            lg = candidate->group;
            del = lg;
            rg = *candidate->know_groups.begin();
            del_in = rg;
        }
        
        std::map<int,  capacity_type> caps;
        for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
            int id = *iii;
            capacity_type cap = std::min(candidate->cap[id], rev_groups_info[del_in]->cap[id]);
            candidate->cap[id] -= cap;
            rev_groups_info[del_in]->cap[id] -= cap;
            caps[id] = cap;
        }
        
        candidate->know_groups.clear();
        rev_groups_info[del_in]->know_groups.erase(del);
    
        if (rev_in_groups_ev[lg] && new_group_nodes_left.find(lg) == new_group_nodes_left.end()) {
            new_group_nodes_left[lg] = wc.addNode();
            nic[new_group_nodes_left[lg]] = nic[node];
            unsecurityId[new_group_nodes_left[lg]] = unsecurityId[node];
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Migrate Group " + std::to_string(lg) + " to new Group " + std::to_string(wc.id(new_group_nodes_left[lg])) + "\n");
            #endif
            
            //move over the edges
            for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
                wc.changeTarget(wc.arcFromId(*ri), new_group_nodes_left[lg]);
                
                for (arc_bridge::iterator bi = know_paths[wc.arcFromId(*ri)].begin(); bi != know_paths[wc.arcFromId(*ri)].end(); ++bi) {
                    know_back_paths[wc.arcFromId(bi->first)].remove_id(*ri);
                }
                know_paths[wc.arcFromId(*ri)].clear();
            }
        }
        
        if (rev_out_groups_ev[rg] && new_group_nodes_right.find(rg) == new_group_nodes_right.end()) {
            ListDigraph::Node nn = wc.addNode();
            new_group_nodes_right[rg] = nn;
            nic[new_group_nodes_right[rg]] = nic[node];
            unsecurityId[new_group_nodes_right[rg]] = unsecurityId[node];
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Migrate Group " + std::to_string(rg) + " to new Group " + std::to_string(wc.id(new_group_nodes_right[rg])) + " " + std::to_string(wc.id(nn)) + "\n");
            #endif
            
            //move over the edges
            for (std::set<int>::iterator ri = rev_out_groups[rg].begin(); ri != rev_out_groups[rg].end(); ++ri) {
                wc.changeSource(wc.arcFromId(*ri), new_group_nodes_right[rg]);
                for (arc_back_bridge::iterator bi = know_back_paths[wc.arcFromId(*ri)].begin(); bi != know_back_paths[wc.arcFromId(*ri)].end(); ++bi) {
                    know_paths[wc.arcFromId(*bi)].remove_id(*ri);
                }
                know_back_paths[wc.arcFromId(*ri)].clear();
            }
        }
        
        if (rev_in_groups_ev[lg] && rev_out_groups_ev[rg]) {
                
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Double Block \n");
            #endif
            
            ListDigraph::Arc new_arc = wc.addArc(new_group_nodes_left[lg], new_group_nodes_right[rg]);
            aic[new_arc].edge_type = edge_types::RESOLVE_HELPER;
            aic[new_arc].edge_specifier = exon_edge(size);
            
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;

                fsc[new_arc].get_flow(id) = caps[id];

                float mean_left = 0;
                float score_left = 0;
                capacity_type weight_left = 0;
                for (ListDigraph::InArcIt a(wc, new_group_nodes_left[lg]); a!=INVALID; ++a) {
                    if (!barred[a]) {
                        mean_left += fsc[a].get_flow(id) * fsc[a].get_mean(id).mean;
                        score_left += fsc[a].get_flow(id) * fsc[a].get_mean(id).compute_score();
                        weight_left += fsc[a].get_flow(id);
                    }
                }
                mean_left = mean_left / (float) weight_left;
                score_left = score_left / (float) weight_left;

                float mean_right = 0;
                float score_right = 0;
                capacity_type weight_right = 0;
                for (ListDigraph::OutArcIt a(wc, new_group_nodes_right[rg]); a!=INVALID; ++a) {
                    if (!barred[a]) {
                        mean_right += fsc[a].get_flow(id) * fsc[a].get_mean(id).mean;
                        score_right += fsc[a].get_flow(id) * fsc[a].get_mean(id).compute_score();
                        weight_right += fsc[a].get_flow(id);
                    }
                }
                mean_right = mean_right / (float) weight_right;
                score_right = score_right / (float) weight_right;

                float mean;
                float score;
                capacity_type weight;
                if (score_left > score_right) {
                    mean = mean_left;
                    score = score_left;
                    weight = weight_left;
                } else {
                    mean = mean_right;
                    score = score_right;
                    weight = weight_right;
                }

                float ratio = (float) caps[id] / (float) weight;
                fsc[new_arc].get_mean(id).hidden_score = score * ratio;
                fsc[new_arc].get_mean(id).weight = 0;
                fsc[new_arc].get_mean(id).assign_mean(mean * ratio);
            }
            
        } else if (rev_in_groups_ev[lg] || rev_out_groups_ev[rg]) {

            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Single Block \n");
            #endif
            
            ListDigraph::Arc base_arc, new_arc;

            if (rev_in_groups_ev[lg]) {
                base_arc = wc.arcFromId(*rev_out_groups[rg].begin());
                new_arc = wc.addArc(new_group_nodes_left[lg], wc.target(base_arc));
            } else { // rev_out_groups_ev[rg].block
                base_arc = wc.arcFromId(*rev_in_groups[lg].begin());
                new_arc = wc.addArc(wc.source(base_arc), new_group_nodes_right[rg]);
            }

            // just copy over stuff   
            aic[new_arc] = aic[base_arc];
            
            if (unsecurityArc[base_arc].has_ref()) {
                std::copy(unsecurityArc[base_arc]->begin(), unsecurityArc[base_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
            }
            
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;
                
                capacity_type capi = caps[id];
                fsc[new_arc].get_flow(id) = capi;
                fsc[new_arc].get_mean(id) = fsc[base_arc].get_mean(id);
                fsc[new_arc].get_mean(id).reduce(capi/ float(fsc[base_arc].get_flow(id)));

                fsc[base_arc].get_mean(id).reduce((fsc[base_arc].get_flow(id) - capi)/ float(fsc[base_arc].get_flow(id)));
                fsc[base_arc].get_flow(id) -= capi;
            }

            if (rev_in_groups_ev[lg]) {
                know_paths[new_arc].bridges.ref() = know_paths[base_arc].bridges.ref(); // deep copy
                if (fsc[base_arc].get_flow(guiding_id) == 0) {
                    for (arc_bridge::iterator bi = know_paths[base_arc].begin(); bi != know_paths[base_arc].end(); ++bi) {
                        know_back_paths[wc.arcFromId(bi->first)].replace_evidence(wc.id(base_arc), wc.id(new_arc));
                    }
                    know_paths[base_arc].clear();
                } else {
                    for (arc_bridge::iterator bi = know_paths[base_arc].begin(); bi != know_paths[base_arc].end(); ++bi) {
                        know_back_paths[wc.arcFromId(bi->first)].add_evidence_if(wc.id(base_arc), wc.id(new_arc));
                    }
                }
            } else {
                know_back_paths[new_arc].bridges.ref() = know_back_paths[base_arc].bridges.ref(); // deep copy
                if (fsc[base_arc].get_flow(guiding_id) == 0) {
                    for (arc_back_bridge::iterator bi = know_back_paths[base_arc].begin(); bi != know_back_paths[base_arc].end(); ++bi) {
                        know_paths[wc.arcFromId(*bi)].replace_evidence(wc.id(base_arc), wc.id(new_arc));
                    }
                    know_back_paths[base_arc].clear();
                } else {
                     for (arc_back_bridge::iterator bi = know_back_paths[base_arc].begin(); bi != know_back_paths[base_arc].end(); ++bi) {
                        know_paths[wc.arcFromId(*bi)].add_evidence_if(wc.id(base_arc), wc.id(new_arc));
                    }
                }
            }
        } else {
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Single Single \n");
            #endif
            unravel_ILP(*rev_in_groups[lg].begin(), *rev_out_groups[rg].begin(), caps, guiding_id, wc, fsc, aic, know_paths, know_back_paths, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);
        }    
    }
    
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_left.begin(); ngi != new_group_nodes_left.end(); ++ngi) {
        
        ListDigraph::Node n = ngi->second;
        
        struct cinfo {
            float mean = 0;
            float score = 0;
            capacity_type group_flow = 0;
            capacity_type evidenced_flow = 0;
        };
        std::map<int, cinfo> msm;
            
        for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
            int id = *iii;
        
            cinfo& ci = msm[id];
            capacity_type weight = 0;
            
            for (ListDigraph::InArcIt a(wc, n); a!=INVALID; ++a) {
                ci.group_flow += fsc[a].get_flow(id);
                if (!barred[a]) {
                    ci.mean += fsc[a].get_flow(id) * fsc[a].get_mean(id).mean;
                    ci.score += fsc[a].get_flow(id) * fsc[a].get_mean(id).compute_score();
                    weight += fsc[a].get_flow(id);
                }
            }
            ci.mean = ci.mean / (float) weight;
            ci.score = ci.score / (float) weight;
        
            for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
                ci.evidenced_flow += fsc[a].get_flow(id);
            }
        }
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("node L " + std::to_string(wc.id(n)) + "\n");
        #endif
        
        if (msm[guiding_id].group_flow > msm[guiding_id].evidenced_flow) {
            ListDigraph::Arc new_arc = wc.addArc(n, node);
            aic[new_arc].edge_type = edge_types::RESOLVE_HELPER;
            aic[new_arc].edge_specifier = exon_edge(size);
            
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;
                cinfo& ci = msm[id];
                
                fsc[new_arc].get_flow(id) = ci.group_flow - ci.evidenced_flow;
                float ratio = (float) ci.evidenced_flow / (float) ci.group_flow;
                fsc[new_arc].get_mean(id).hidden_score = ci.score * ratio;
                fsc[new_arc].get_mean(id).weight = 0;
                fsc[new_arc].get_mean(id).assign_mean(ci.mean * ratio);
            }
        }
    }
                
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_right.begin(); ngi != new_group_nodes_right.end(); ++ngi) {
        
        ListDigraph::Node n = ngi->second;
        
        struct cinfo {
            float mean = 0;
            float score = 0;
            capacity_type group_flow = 0;
            capacity_type evidenced_flow = 0;
        };
        std::map<int, cinfo> msm;
        
        for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
            int id = *iii;
        
            cinfo& ci = msm[id];
            capacity_type weight = 0;
            
            for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
                ci.group_flow += fsc[a].get_flow(id);
                if (!barred[a]) {
                    ci.mean += fsc[a].get_flow(id) * fsc[a].get_mean(id).mean;
                    ci.score += fsc[a].get_flow(id) * fsc[a].get_mean(id).compute_score();
                    weight += fsc[a].get_flow(id);
                }
            }
            ci.mean = ci.mean / (float) weight;
            ci.score = ci.score / (float) weight;
        
            for (ListDigraph::InArcIt a(wc, n); a!=INVALID; ++a) {
                ci.evidenced_flow += fsc[a].get_flow(id);
            }
        }
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("node R " + std::to_string(wc.id(n))   + "\n");
        #endif
        
        if (msm[guiding_id].group_flow > msm[guiding_id].evidenced_flow) {
            ListDigraph::Arc new_arc = wc.addArc(node, n);
            aic[new_arc].edge_type = edge_types::RESOLVE_HELPER;
            aic[new_arc].edge_specifier = exon_edge(size);
            
            
            for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
                int id = *iii;
                cinfo& ci = msm[id];
                
                fsc[new_arc].get_flow(id) = ci.group_flow - ci.evidenced_flow;
                float ratio = (float) ci.evidenced_flow / (float) ci.group_flow;
                fsc[new_arc].get_mean(id).hidden_score = ci.score * ratio;
                fsc[new_arc].get_mean(id).weight = 0;
                fsc[new_arc].get_mean(id).assign_mean(ci.mean * ratio);
            }
        }
    }
                
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_left.begin(); ngi != new_group_nodes_left.end(); ++ngi) {
        clean_barred_leftovers(ngi->second, barred, guiding_id, wc, aic, fsc, know_paths, know_back_paths, unsecurityArc, unsecurityId, input_ids, transcripts);
    }
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_right.begin(); ngi != new_group_nodes_right.end(); ++ngi) {
        clean_barred_leftovers(ngi->second, barred, guiding_id, wc, aic, fsc, know_paths, know_back_paths, unsecurityArc, unsecurityId, input_ids, transcripts);
    }
    
    // we now look at possibly unevidenced leftovers
    // remove all 0 edges here
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fsc[arc].get_flow(guiding_id) == 0) {
            for (arc_bridge::iterator bi = know_paths[arc].begin(); bi != know_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(bi->first);
                know_back_paths[a].remove_id(wc.id(arc));
                
            }
            know_paths[arc].clear();
            erase_arc(wc, arc, fsc, input_ids, transcripts);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fsc[arc].get_flow(guiding_id) == 0) {
            for (arc_back_bridge::iterator bi = know_back_paths[arc].begin(); bi != know_back_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(*bi);
                know_paths[a].remove_id(wc.id(arc));
            }  
            know_back_paths[arc].clear();
            erase_arc(wc, arc, fsc, input_ids, transcripts);
        }
    }
    
    int indeg = 0;
    int outdeg = 0;
    ListDigraph::InArcIt in(wc, node);
    ListDigraph::OutArcIt out(wc, node);
    for (;in != INVALID; ++in) {
        ++indeg;
    }
    for (;out != INVALID; ++out) {
        ++outdeg;
    }
    
    if (indeg == 1 || outdeg == 1) {
        unravel_unevidenced_leftovers(node, guiding_id, wc, aic, fsc, know_paths, know_back_paths, unsecurityArc, unsecurityId, input_ids, transcripts, &barred);
        indeg = 0;
    }
    if (indeg == 0) {
        // we completely resolved this one, remove!
        wc.erase(node);
    }
}

void base_manager::unravel_evidences(ListDigraph::Node node, int guiding,
            ListDigraph &wc,
            ListDigraph::ArcMap<arc_identifier> &ai, ListDigraph::ArcMap<flow_series> &fs,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId,
            std::set<int>& ip_ids, ATC &trans) {
    
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Unravel NODE " + std::to_string(wc.id(node)) + "\n");
    #endif
    
    struct info {
        
        info() {}
        info(ListDigraph::Arc a, bool is_out, capacity_type cap) : arc(a), is_out(is_out), cap(cap) 
            {}
        
        ListDigraph::Arc arc;
        bool is_out = false;
        capacity_type cap = 0;
    };
    
    std::deque<info> arcs;
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        if (ai[a].edge_type == edge_types::HELPER) {
            // no resolving 
            continue;
        }
        arcs.push_back(info(a, false, fs[a].get_flow(guiding)));
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        if (ai[a].edge_type == edge_types::HELPER) {
            // no resolving 
            continue;
        }
        arcs.push_back(info(a, true, fs[a].get_flow(guiding)));
    }

    bool change = true;
    while (change) {
        change = false;
        // loop over collected
        
        info candidate;
        capacity_type max = 0;
        for (std::deque< info>::iterator it = arcs.begin(); it != arcs.end(); ++it) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("ARC " + std::to_string(wc.id(it->arc))  + " OUT " + std::to_string(it->is_out) + " kp " + std::to_string(know_paths[it->arc].size()) + " kbp " + std::to_string(know_back_paths[it->arc].size()) + "\n");
            #endif
            
            if (!it->is_out && know_paths[it->arc].size() == 1) {
                capacity_type min = std::min(it->cap, fs[ wc.arcFromId(know_paths[it->arc].begin()->first) ].get_flow(guiding));
                if (min > max ){
                    candidate = *it;
                    max = min;
                }
            } else if (it->is_out && know_back_paths[it->arc].size() == 1) {
                capacity_type min = std::min(it->cap, fs[ wc.arcFromId(*know_back_paths[it->arc].begin()) ].get_flow(guiding));
                if (min > max){
                    candidate = *it;
                    max = min;
                }
            }
        }
        
        if (max == 0) break;
        
        // we have something to resolve, yay
        change = true;
        
        ListDigraph::Arc a = candidate.arc;
        if (candidate.is_out) {
            unravel_evidence_path_right(*know_back_paths[a].begin(), wc.id(a), guiding,  wc, ai, fs, know_paths, know_back_paths, unsecurityArc, unsecurityId, transcript_unsecurity::RESOLVED, ip_ids);
        } else {
            unravel_evidence_path_left(wc.id(a), know_paths[a].begin()->first, guiding,  wc, ai, fs, know_paths, know_back_paths, unsecurityArc, unsecurityId, transcript_unsecurity::RESOLVED, ip_ids);
        }   
    }
    
    // we now look at possibly unevidenced leftovers
    // remove all 0 edges here
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fs[arc].get_flow(guiding) == 0) {
            erase_arc(wc, arc, fs, ip_ids, trans);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fs[arc].get_flow(guiding) == 0) {
            erase_arc(wc, arc, fs, ip_ids, trans);
        }
    }
    
    int indeg = 0;
    int outdeg = 0;
    ListDigraph::InArcIt in(wc, node);
    ListDigraph::OutArcIt out(wc, node);
    for (;in != INVALID; ++in) {
        ++indeg;
    }
    for (;out != INVALID; ++out) {
        ++outdeg;
    }
    
    if (indeg == 1 || outdeg == 1) {
        unravel_unevidenced_leftovers(node, guiding, wc, ai, fs, know_paths, know_back_paths, unsecurityArc, unsecurityId, ip_ids, trans, NULL);
        indeg = 0;
    }
    if (indeg == 0) {
        // we completely resolved this one, remove!
        wc.erase(node);
    }
}


void base_manager::unravel_evidence_path_left(int left, int right, int guiding,
            ListDigraph &wc, 
            ListDigraph::ArcMap<arc_identifier> &ai, ListDigraph::ArcMap<flow_series> &fs,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, transcript_unsecurity::evidence evidence,
            std::set<int>& ip_ids) {
    
    ListDigraph::Arc left_arc = wc.arcFromId(left);
    ListDigraph::Arc right_arc = wc.arcFromId(right);
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract Left. " + ai[left_arc].edge_specifier.to_string() + ";" + ai[right_arc].edge_specifier.to_string() + " - " +  std::to_string(left) + " " + std::to_string(right) + "\n");
    #endif
    
    // this should never be called on a zero edge, so add in new arcs
    ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));
    ai[new_arc].edge_type = edge_types::EXON;
    
    for (std::set<int>::iterator idi = ip_ids.begin(); idi != ip_ids.end(); ++idi) {
        int id = *idi;
        
        // we get minimal capacity
        capacity_type cap_evidence = fs[left_arc].get_flow(id);
        if (fs[right_arc].get_flow(id) < cap_evidence) {
            cap_evidence = fs[right_arc].get_flow(id);
        }

        capacity_mean c1 = fs[left_arc].get_mean(id);
        capacity_mean c2 = fs[right_arc].get_mean(id);

        capacity_type flow_left = fs[left_arc].get_flow(id);
        capacity_type flow_right = fs[right_arc].get_flow(id);
        
        if (flow_left != 0) {
            c1.reduce(cap_evidence/ float(flow_left));
            fs[left_arc].get_mean(id).reduce((flow_left - cap_evidence)/ float(flow_left));
        } else {
            c1.reduce(0);
        }
        if (flow_right != 0) {
            c2.reduce(cap_evidence/ float(flow_right));
            fs[right_arc].get_mean(id).reduce((flow_right - cap_evidence)/ float(flow_right));
        } else {
            c2.reduce(0);
        }
        
        c1.update(c2);
    
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Capacity " + std::to_string(cap_evidence) + "\n");
        #endif

        fs[left_arc].get_flow(id) -= cap_evidence;
        fs[right_arc].get_flow(id) -= cap_evidence;
        fs[new_arc].get_flow(id) = cap_evidence;        
        fs[new_arc].get_mean(id) = c1;
    }
    
    if (ai[left_arc].edge_type != edge_types::EXON && ai[right_arc].edge_type != edge_types::EXON) {
        ai[new_arc].edge_type = edge_types::HELPER;
        if (ai[right_arc].edge_specifier.node_index == -1) {
            ai[new_arc].edge_specifier = ai[left_arc].edge_specifier;
        } else {
            ai[new_arc].edge_specifier = ai[right_arc].edge_specifier;
        }
    } else if (ai[left_arc].edge_type == edge_types::HELPER) {
        // leftarc is a helper
        ai[new_arc].edge_specifier = exon_edge(ai[right_arc].edge_specifier); 
        
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else if (ai[right_arc].edge_type == edge_types::HELPER) {
        // rightarc is a helper
        ai[new_arc].edge_specifier = exon_edge(ai[left_arc].edge_specifier);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else {
        ai[new_arc].edge_specifier = exon_edge(ai[left_arc].edge_specifier);
        ai[new_arc].edge_specifier.join_edge(ai[right_arc].edge_specifier);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->begin()));
        }
        unsecurityArc[new_arc]->insert( transcript_unsecurity( unsecurityId[wc.target(left_arc)].id, evidence));
    }
     
    ai[new_arc].edge_lengths.first_exon = ai[left_arc].edge_lengths.first_exon;
    ai[new_arc].edge_lengths.middle = ai[left_arc].edge_lengths.middle + ai[left_arc].edge_lengths.last_exon + ai[right_arc].edge_lengths.middle;
    ai[new_arc].edge_lengths.last_exon = ai[right_arc].edge_lengths.last_exon;
    
    ai[new_arc].cycle_id_in = ai[left_arc].cycle_id_in;
    ai[new_arc].cycle_id_out = ai[right_arc].cycle_id_out;
    
    // we now need to transfer the know paths of the right and left side
    know_paths[new_arc].bridges.ref() = know_paths[right_arc].bridges.ref(); // deep copy
    know_back_paths[new_arc].bridges.ref() = know_back_paths[left_arc].bridges.ref(); // deep copy
    
    // we clear the left arc in any case, as it was resolved
    know_paths[left_arc].clear();
    // but evidence can vary
    if (fs[left_arc].get_flow(guiding) == 0) {
        // so we need a full id switch in evidence
        // we need to transfer id references from all previous edges, somewhat bothersome
        for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
            know_paths[a].replace_evidence(wc.id(left_arc), wc.id(new_arc));
        }
    } else {
        // we have to add in the additional evidence
        for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
            know_paths[a].add_evidence_if(wc.id(left_arc), wc.id(new_arc));
        }
    }
    
    // right is counted down depending on whether it was set zero or not
    if (fs[right_arc].get_flow(guiding) == 0) {
        
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].replace_evidence(wc.id(right_arc), wc.id(new_arc));
        }
        
        for (arc_back_bridge::iterator bi = know_back_paths[right_arc].begin(); bi != know_back_paths[right_arc].end(); ++bi) {
            ListDigraph::Arc a = wc.arcFromId(*bi);
            know_paths[a].remove_id(right);
        }  
        know_back_paths[right_arc].clear();
        
    } else {
        
        // so we have capacity left in right path
        // just count right path down one
        // this should be the default?
        
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].add_evidence_if(wc.id(right_arc), wc.id(new_arc));
        }
        
        know_back_paths[right_arc].remove_id(left);
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Exit Left.\n");
    #endif
}

void base_manager::unravel_evidence_path_right(int left, int right, int guiding,
            ListDigraph &wc, 
            ListDigraph::ArcMap<arc_identifier> &ai, ListDigraph::ArcMap<flow_series> &fs,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, transcript_unsecurity::evidence evidence,
            std::set<int>& ip_ids) {
    
    ListDigraph::Arc left_arc = wc.arcFromId(left);
    ListDigraph::Arc right_arc = wc.arcFromId(right);

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract Right. " + ai[left_arc].edge_specifier.to_string() + ";" + ai[right_arc].edge_specifier.to_string() + " - " +  std::to_string(left) + " " + std::to_string(right) + "\n");
    #endif

    // this should never be called on a zero edge, so add in new arcs
    ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));
    ai[new_arc].edge_type = edge_types::EXON;
    
    for (std::set<int>::iterator idi = ip_ids.begin(); idi != ip_ids.end(); ++idi) {
        int id = *idi;
        
        // we get minimal capacity
        capacity_type cap_evidence = fs[left_arc].get_flow(id);
        if (fs[right_arc].get_flow(id) < cap_evidence) {
            cap_evidence = fs[right_arc].get_flow(id);
        }

        capacity_mean c1 = fs[left_arc].get_mean(id);
        capacity_mean c2 = fs[right_arc].get_mean(id);

        capacity_type flow_left = fs[left_arc].get_flow(id);
        capacity_type flow_right = fs[right_arc].get_flow(id);
        
        if (flow_left != 0) {
            c1.reduce(cap_evidence/ float(flow_left));
            fs[left_arc].get_mean(id).reduce((flow_left - cap_evidence)/ float(flow_left));
        } else {
            c1.reduce(0);
        }
        if (flow_right != 0) {
            c2.reduce(cap_evidence/ float(flow_right));
            fs[right_arc].get_mean(id).reduce((flow_right - cap_evidence)/ float(flow_right));
        } else {
            c2.reduce(0);
        }
        
        c1.update(c2);
    
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Capacity " + std::to_string(cap_evidence) + "\n");
        #endif

        fs[left_arc].get_flow(id) -= cap_evidence;
        fs[right_arc].get_flow(id) -= cap_evidence;
        fs[new_arc].get_flow(id) = cap_evidence;        
        fs[new_arc].get_mean(id) = c1;
    }
    
    if (ai[left_arc].edge_type != edge_types::EXON && ai[right_arc].edge_type != edge_types::EXON) {
        ai[new_arc].edge_type = edge_types::HELPER;
        if (ai[right_arc].edge_specifier.node_index == -1) {
            ai[new_arc].edge_specifier = ai[left_arc].edge_specifier;
        } else {
            ai[new_arc].edge_specifier = ai[right_arc].edge_specifier;
        }
    } else if (ai[left_arc].edge_type == edge_types::HELPER) {
        // leftarc is a helper
        ai[new_arc].edge_specifier = exon_edge(ai[right_arc].edge_specifier);
        
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else if (ai[right_arc].edge_type == edge_types::HELPER) {
        // rightarc is a helper
        ai[new_arc].edge_specifier = exon_edge(ai[left_arc].edge_specifier); 
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else {
        ai[new_arc].edge_specifier = exon_edge(ai[left_arc].edge_specifier);
        ai[new_arc].edge_specifier.join_edge(ai[right_arc].edge_specifier);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->begin()));
        }
        unsecurityArc[new_arc]->insert( transcript_unsecurity(unsecurityId[wc.target(left_arc)].id, evidence));
    }
        
    ai[new_arc].edge_lengths.first_exon = ai[left_arc].edge_lengths.first_exon;
    ai[new_arc].edge_lengths.middle = ai[left_arc].edge_lengths.middle + ai[left_arc].edge_lengths.last_exon + ai[right_arc].edge_lengths.middle;
    ai[new_arc].edge_lengths.last_exon = ai[right_arc].edge_lengths.last_exon;
    
    ai[new_arc].cycle_id_in = ai[left_arc].cycle_id_in;
    ai[new_arc].cycle_id_out = ai[right_arc].cycle_id_out;
    
    // we now need to transfer the know paths of the right and left side
    know_paths[new_arc].bridges.ref() = know_paths[right_arc].bridges.ref(); // deep copy
    know_back_paths[new_arc].bridges.ref() = know_back_paths[left_arc].bridges.ref(); // deep copy
    
    // we clear the right arc in any case, as it was resolved
    know_back_paths[right_arc].clear();
    // but evidence can vary
    if (fs[right_arc].get_flow(guiding) == 0) {
        // so we need a full id switch in evidence
        // we need to transfer id references from all previous edges, somewhat bothersome
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].replace_evidence(wc.id(right_arc), wc.id(new_arc));
        }
    } else {
        // we have to add in the additional evidence
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].add_evidence_if(wc.id(right_arc), wc.id(new_arc));
        }  
    }
    
    // right is counted down depending on whether it was set zero or not
    if (fs[left_arc].get_flow(guiding) == 0) {
        
        for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
            know_paths[a].replace_evidence(wc.id(left_arc), wc.id(new_arc));
        }
        
        for (arc_bridge::iterator bi = know_paths[left_arc].begin(); bi != know_paths[left_arc].end(); ++bi) {
            ListDigraph::Arc a = wc.arcFromId(bi->first);
            know_back_paths[a].remove_id(left);
        }
        know_paths[left_arc].clear();
        
    } else {
        
        // so we have capacity left in left path
        // just count right path down one
        // this should be the default?
        
        for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
            know_paths[a].add_evidence_if(wc.id(left_arc), wc.id(new_arc));
        }
        
        know_paths[left_arc].remove_id(right);
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Exit Right.\n");
    #endif
}

void base_manager::unravel_unevidenced_leftovers(ListDigraph::Node &node, int guiding,
            ListDigraph &wc,
            ListDigraph::ArcMap<arc_identifier> &ai, ListDigraph::ArcMap<flow_series> &fs,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, std::set<int>& ip_ids, ATC &trans, 
            ListDigraph::ArcMap<bool> *barred) {
    
    std::unordered_set<int> left, right; 
    for (ListDigraph::InArcIt left_arc(wc, node); left_arc!=INVALID; ++left_arc) {
        for (ListDigraph::OutArcIt right_arc(wc, node); right_arc!=INVALID; ++right_arc) {

            left.insert(wc.id(left_arc));
            right.insert(wc.id(right_arc));

            ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Join leftover " + std::to_string(wc.id(left_arc)) + " " + std::to_string(wc.id(right_arc)) + " to " + std::to_string(wc.id(new_arc)) + "\n");
            #endif
            
            for (std::set<int>::iterator iii = ip_ids.begin(); iii != ip_ids.end(); ++iii) {
                int id = *iii;
                
                capacity_mean c1 = fs[left_arc].get_mean(id);
                capacity_mean c2 = fs[right_arc].get_mean(id);

                capacity_type flow_left = fs[left_arc].get_flow(id);
                capacity_type flow_right = fs[right_arc].get_flow(id);
                
                capacity_type capacity = std::min(flow_left, flow_right);

                if (flow_left != 0) {
                    c1.reduce(capacity/ float(flow_left));
                    fs[left_arc].get_mean(id).reduce((flow_left - capacity)/ float(flow_left));
                } else {
                    c1.reduce(0);
                }
                if (flow_right != 0) {
                    c2.reduce(capacity/ float(flow_right));
                    fs[right_arc].get_mean(id).reduce((flow_right - capacity)/ float(flow_right));
                } else {
                    c2.reduce(0);
                }

                c1.update(c2);

                fs[new_arc].get_flow(id) = capacity;
                fs[new_arc].get_mean(id) = c1;
            }
            
            ai[new_arc].edge_type = edge_types::EXON;

            if (barred != NULL) {
                (*barred)[new_arc] = (*barred)[left_arc] || (*barred)[right_arc];
            }
            
            if (ai[left_arc].edge_type != edge_types::EXON && ai[right_arc].edge_type != edge_types::EXON) {
                ai[new_arc].edge_type = edge_types::HELPER;
                if (ai[right_arc].edge_specifier.node_index == -1) {
                    ai[new_arc].edge_specifier = ai[left_arc].edge_specifier;
                } else {
                    ai[new_arc].edge_specifier = ai[right_arc].edge_specifier;
                }
            } else if (ai[left_arc].edge_type == edge_types::HELPER) {
                // leftarc is a helper
                if (unsecurityArc[right_arc].has_ref()) {
                    std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
                }
                ai[new_arc].edge_specifier = exon_edge(ai[right_arc].edge_specifier);
                
            } else if (ai[right_arc].edge_type == edge_types::HELPER) {
                // rightarc is a helper
                if (unsecurityArc[left_arc].has_ref()) {
                    std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
                }
                ai[new_arc].edge_specifier = exon_edge(ai[left_arc].edge_specifier);
            } else {
                if (unsecurityArc[left_arc].has_ref()) {
                    std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
                }
                if (unsecurityArc[right_arc].has_ref()) {
                    std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->begin()));
                }
                
                if (unsecurityId[wc.target(left_arc)].resolvable && !know_paths[left_arc].has_path_evidence(wc.id(right_arc)) ) {
                    unsecurityArc[new_arc]->insert( transcript_unsecurity( unsecurityId[wc.target(left_arc)].id, transcript_unsecurity::UNEVIDENCED));
                } else {
                    unsecurityArc[new_arc]->insert( transcript_unsecurity( unsecurityId[wc.target(left_arc)].id, transcript_unsecurity::GUESSED));
                }
                
                ai[new_arc].edge_specifier = exon_edge(ai[left_arc].edge_specifier);
                ai[new_arc].edge_specifier.join_edge(ai[right_arc].edge_specifier);
            }
            
            ai[new_arc].edge_lengths.first_exon = ai[left_arc].edge_lengths.first_exon;
            ai[new_arc].edge_lengths.middle = ai[left_arc].edge_lengths.middle + ai[left_arc].edge_lengths.last_exon + ai[right_arc].edge_lengths.middle;
            ai[new_arc].edge_lengths.last_exon = ai[right_arc].edge_lengths.last_exon;

            ai[new_arc].cycle_id_in = ai[left_arc].cycle_id_in;
            ai[new_arc].cycle_id_out = ai[right_arc].cycle_id_out;

            // we now need to transfer the know paths of the right and left side
            know_paths[new_arc].bridges.ref() = know_paths[right_arc].bridges.ref(); // deep copy
            know_back_paths[new_arc].bridges.ref() = know_back_paths[left_arc].bridges.ref(); // deep copy

            // we to remove the right arc  (del later!)
            for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
                know_back_paths[a].add_evidence_if(wc.id(right_arc), wc.id(new_arc));
            }

            //left arc is killed as well (del later!)
            for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
                know_paths[a].add_evidence_if(wc.id(left_arc), wc.id(new_arc));
            }
            
        }
    }

    for ( std::unordered_set<int>::iterator it = left.begin(); it!= left.end(); ++it) {
       
        for (ListDigraph::InArcIt a(wc, wc.source(wc.arcFromId(*it))); a!=INVALID; ++a) {
            know_paths[a].remove_id(*it);
        }
        
        know_paths[wc.arcFromId(*it)].clear();
        know_back_paths[wc.arcFromId(*it)].clear();
        erase_arc(wc, wc.arcFromId(*it), fs, ip_ids, trans);
    }
    for ( std::unordered_set<int>::iterator it = right.begin(); it!= right.end(); ++it) {
        
        for (ListDigraph::OutArcIt a(wc, wc.target(wc.arcFromId(*it))); a!=INVALID; ++a) {
            know_back_paths[a].remove_id(*it);
        }
        
        know_paths[wc.arcFromId(*it)].clear();
        know_back_paths[wc.arcFromId(*it)].clear();
        erase_arc(wc, wc.arcFromId(*it), fs, ip_ids, trans);
    }
    
}

void base_manager::erase_arc(ListDigraph &wc, ListDigraph::Arc arc, ListDigraph::ArcMap<flow_series>& fsc,  std::set<int>& ip_ids, ATC &transcripts) {
    
    for (std::set<int>::iterator iii = ip_ids.begin(); iii != ip_ids.end(); ++iii) {
        int id = *iii;
        transcripts.total_flow_error[id] += fsc[arc].get_flow(id);
    }
    wc.erase(arc);
}

void base_manager::extend_path_left(std::deque<path>& paths, path* p, unsigned int& border,
        ListDigraph& wc, ListDigraph::ArcMap<arc_identifier> &ai,
	ListDigraph::NodeMap<unsigned int>& cni, ListDigraph::ArcMap<arc_bridge> &kp) {
        
    
    if ( paths.size() <= options::Instance()->get_max_enumerated_paths() && border < p->border_index )  {
        
//        #ifdef ALLOW_DEBUG
//        logger::Instance()->debug("Extend Left.\n");
//        #endif
        
        // extend path at last path node
        
        // we create a parent arc for all current
        ListDigraph::Arc arc_save = INVALID;
	bool first = false;
        for (ListDigraph::InArcIt a(wc, p->last_node) ; a!=INVALID; ++a) {
            
            if (ai[a].edge_type == edge_types::HELPER || ai[a].edge_type == edge_types::RESOLVE_HELPER ||
                    (kp[a].is_evidenced_path() && !kp[a].has_path_evidence(wc.id(p->last_arc))) ) {
                continue;
            }
            
	    if (first) {
                arc_save = a;
                first = false;
                continue;
            }
            
            paths.push_back(path(p));  // makes copy
             
            path *np = &paths.back();
            np->last_node = wc.source(a);
            np->last_arc = a;
            np->border_index = cni[np->last_node];
            np->identifier.join_edge(ai[a].edge_specifier); 
            
            extend_path_left(paths, np, border, wc, ai, cni, kp);
        }
        
        if ( arc_save==INVALID) {
            return;
        }
        
        p->last_node = wc.source(arc_save);
        p->last_arc = arc_save;
        p->border_index = cni[p->last_node];
        p->identifier.join_edge(ai[arc_save].edge_specifier); 
        
        extend_path_left(paths, p, border, wc, ai, cni, kp);
    }
}



void base_manager::extend_path_right(std::deque<path>& paths, path* p, unsigned int& border,
        ListDigraph& wc, ListDigraph::ArcMap<arc_identifier> &ai,
	ListDigraph::NodeMap<unsigned int>& cni, ListDigraph::ArcMap<arc_bridge> &kp) {
    
    if ( paths.size() <=  options::Instance()->get_max_enumerated_paths() && border > p->border_index )  {
        
//        #ifdef ALLOW_DEBUG
//        logger::Instance()->debug("Extend right.\n");
//        #endif
        
        // extend path at last path node
        
        // we use first arc to change path, safe it for later

        ListDigraph::Arc arc_save = INVALID;
        bool first = true;
        for ( ListDigraph::OutArcIt a(wc, p->last_node); a!=INVALID; ++a) {
            
//            logger::Instance()->debug("Test Arc " + std::to_string(wc.id(a))+ "\n");
            
            if (ai[a].edge_type == edge_types::HELPER || ai[a].edge_type == edge_types::RESOLVE_HELPER ||
                  (kp[p->last_arc].is_evidenced_path() && !kp[p->last_arc].has_path_evidence(wc.id(a))) ) {
	       continue;
            }
               
            if (first) {
                arc_save = a;
                first = false;
                continue;
            }

            paths.push_back(path(p));  // makes copy
             
            path *np = &paths.back();
            np->last_node = wc.target(a);
            np->last_arc = a;
            np->border_index = cni[np->last_node];
            np->identifier.join_edge(ai[a].edge_specifier);
            
            extend_path_right(paths, np, border, wc, ai, cni, kp);
        }
        
        if ( arc_save==INVALID) {
            return;
        }
        
        //now change path itself
        p->last_node = wc.target(arc_save);
        p->last_arc = arc_save;
        p->border_index = cni[p->last_node];
        p->identifier.join_edge(ai[arc_save].edge_specifier);

        extend_path_right(paths, p, border, wc, ai, cni, kp);
    }
    
}

void base_manager::print_graph_debug_copy(std::ostream& os,
    ListDigraph &wc,
    ListDigraph::ArcMap<flow_series> &fc,
    ListDigraph::ArcMap<arc_identifier> &ai,   
    ListDigraph::Node &cs, ListDigraph::Node &ct) {

    if (options::Instance()->is_debug()) {
    
        digraphWriter(wc, os)
                .arcMap("ai", ai)
                .arcMap("flow", fc)
                .node("source", cs)
                .node("drain", ct)
                .run();  
    }
}


void base_manager::print_graph_gs(std::ostream& os,
    ListDigraph &wc,
    ListDigraph::ArcMap<capacity_type> &fc,
    ListDigraph::ArcMap<arc_identifier> &ai,  
    ListDigraph::Node &cs, ListDigraph::Node &ct) {

    digraphWriter(wc, os)
            .arcMap("Identifier", ai)
            .arcMap("Flow/Capacity", fc)
            .node("Source", cs)
            .node("Drain", ct)
            .run();  
}
