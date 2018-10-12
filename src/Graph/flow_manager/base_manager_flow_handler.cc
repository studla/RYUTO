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

//void base_manager::extract_transcripts_from_flow(std::ostream &gs) {

void base_manager::extract_transcripts_from_flow() {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract Transcripts.\n");
    #endif
        
    if (meta->size == 1) { // TODO: REAL handling
        return;
    }
   
    #ifdef ALLOW_DEBUG
    print_graph_debug_copy(std::cout, g, capacity, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    
    #ifdef ALLOW_DEBUG
    digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("cap", capacity)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run();  
    #endif  

    // transform cap to mean
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
       if (edge_type[a] == edge_types::NODE) {
           means[a] = capacity_mean(regions[a].get_average(), edge_lengths[a].middle);
       } else if (edge_type[a] == edge_types::EXON) {
           means[a] = capacity_mean(regions[a].get_average(), regions[a].get_length());
       }
    } 
    
    // downsize neighbouring exon junctions, as they are overrepresented often
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (edge_type[a] == edge_types::EXON 
                && edge_specifier[a].left_consecutive && edge_specifier[a].right_consecutive 
                &&  edge_specifier[a].id.count()) {
            
            // downsize this one!
            ListDigraph::InArcIt left_node(g, g.source(a));
            ListDigraph::OutArcIt right_node(g, g.target(a));
            
            float this_ave = regions[a].get_average();
            float left_ave = regions[left_node].get_average();
            float right_ave = regions[right_node].get_average();
            
            float sum_left = 0;
            for (ListDigraph::OutArcIt io(g, g.target(left_node)); io != INVALID; ++io) {
                sum_left += regions[io].get_average();
            }
            
            float sum_right = 0;
            for (ListDigraph::InArcIt io(g, g.source(right_node)); io != INVALID; ++io) {
                sum_right += regions[io].get_average();
            }
            
            if (sum_left == 0 || sum_right == 0) {
                continue;
            }
            
            float min_ave = std::min(left_ave*this_ave/sum_left, right_ave*this_ave/sum_right);
            
            float cap = std::min(min_ave, (float) capacity[a]);
            
            capacity[a] = cap;
            regions[a].subregions.begin()->start_count = cap;
            regions[a].subregions.begin()->end_count = cap;
            regions[a].subregions.begin()->basecount = cap * (regions[a].subregions.begin()->end - regions[a].subregions.begin()->start + 1 );

            means[a] = capacity_mean(cap, regions[a].get_length());
            
        }
    }  
    
 
    #ifdef ALLOW_DEBUG
    digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("cap", capacity)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run();  
    #endif
   
   ListDigraph::ArcMap<bool> guided_saves(g);
   std::deque<std::deque<ListDigraph::Arc> > guided_paths;
   find_guide_sources_and_drains(guided_saves, guided_paths);
    
   // this is a set of exons directly neighbouring other exons
   // we only mark those which are either filled introns or in between sources or drains
   ListDigraph::ArcMap<bool> consecutive_s_t(g); // inits as false
   ListDigraph::ArcMap<bool> marked_source(g);
   ListDigraph::ArcMap<bool> marked_drain(g);
   find_s_t_exons(guided_saves, marked_source, marked_drain, consecutive_s_t);
   
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (edge_type[a] == edge_types::NODE && !marked_source[a] && !marked_drain[a]) {
           means[a].scores.clear();
        }
    }
   
   prune_sources_and_drains_adjacent_starts(guided_saves); // directly adjacent/starts at exactly inner exons 

   #ifdef ALLOW_DEBUG
   print_graph_debug_copy(std::cout, g, capacity, edge_specifier, edge_type, edge_lengths, s, t);
   #endif   
   
   if (!options::Instance()->is_keep_introns()) {
        filter_introns(guided_saves, marked_source, marked_drain);
        filter_broken_introns(guided_saves, marked_source, marked_drain);
   }
//   prune_sources_and_drains4(guided_saves, marked_source, marked_drain, neighbouring);  // inhibit drains/sources on very small inner exons
//   prune_sources_and_drains5(guided_saves, marked_source, marked_drain, neighbouring);   // strong filtering of end
//   prune_sources_and_drains3(guided_saves, marked_source, marked_drain, neighbouring);   // median test
 
   #ifdef ALLOW_DEBUG
   print_graph_debug_copy(std::cout, g, capacity, edge_specifier, edge_type, edge_lengths, s, t);
   #endif
   

   #ifdef ALLOW_DEBUG
   print_graph_debug_copy(std::cout, g, capacity, edge_specifier, edge_type, edge_lengths, s, t);
   #endif
   
   //prune_unevidenced_arcs(guided_saves, marked_source, marked_drain);  // percent left right
   
   
//    std::vector<capacity_type> avgs;
//    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//        if (edge_type[a] == edge_types::EXON) {
//            avgs.push_back(capacity[a]);
//        }
//    }
//    std::sort(avgs.begin(), avgs.end());
////    
//
//    std::vector<capacity_type>::iterator ait = avgs.begin();
//    int j = 0;
//    while ( ait != avgs.end() ){
//        //logger::Instance()->debug(std::to_string(j*10));
//        for (int i = 0 ; i < 10 && ait != avgs.end(); ++i, ++ait) {
//            logger::Instance()->debug("," + std::to_string(*ait));
//        }
//        
//        ++j;
//    }

    float low_mark = options::Instance()->get_low_edge_mark();
    
//    if (!avgs.empty()) {
//
//      //  float median_mark = avgs[avgs.size() * 0.5];
//        
//        //avgs.erase(avgs.begin() + avgs.size() / 2 + 1, avgs.end());
//        //float percent_mark = DKMeans(avgs);
//
//        //if (percent_mark > low_mark) {// && percent_mark < median_mark) {
//        //    low_mark = percent_mark;
//        //}
//        
//        low_mark += avgs[avgs.size() * 0.5] / 100;
//        
//    }
//    avgs.clear();
  
    
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
   
    #ifdef ALLOW_DEBUG
    print_graph_debug_copy(std::cout, g, capacity, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    
    #ifdef ALLOW_DEBUG
    digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("cap", capacity)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run();  
    #endif
    
    //correct_start_end();
    
    if (options::Instance()->is_graph_denoising()) {
        
        if (guided_paths.empty()) {
            ListDigraph::ArcMap<bool> neighbouring2(g); // inits as false, we don't block right now
            denoise_graph(guided_saves, neighbouring2, marked_source, marked_drain);
            //denoise_graph_alt(guided_saves, neighbouring);
        } else {
            denoise_graph_guided(guided_paths);
        }
    }
    
    #ifdef ALLOW_DEBUG
    digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("cap", capacity)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run();  
    #endif
     
    #ifdef ALLOW_DEBUG
    print_graph_debug_copy(std::cout, g, capacity, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    
    finalize_flow_graph();
    compute_flow();
    
    // copy is now a DAG
    
        #ifdef ALLOW_DEBUG
    logger::Instance()->debug("After Flow.\n");
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);

    digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("flow", flow)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run(); 
    #endif
    
    // in some cases we might have lost guides though due to poor performance of flow!
    for(std::deque<std::deque<ListDigraph::Arc> >::iterator it = guided_paths.begin(); it != guided_paths.end(); ++it) {
        for(std::deque<ListDigraph::Arc>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
            ++flow[*it2];
        }
    }
    
    // remove all possible 0 flow edges!
    for (ListDigraph::ArcIt a(g); a != INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        
        if (flow[arc] != 0) { // online erase edges
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
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);

    digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("flow", flow)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run(); 
    #endif
    
    // we remove path node that are one to one, one to many, or many to one
    // this keeps flow rules intact
    while ( contract_composites(in_deg, out_deg) ); // of this runs too often we could make a top. sorting

//    for (ListDigraph::ArcIt a(g); a != INVALID;++a) {
//        means[a].compute_score();
//    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Contracting done.\n");
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
    
    digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("flow", flow)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run();  
    #endif

//    gs << "Compacted graph ---\n";
//    print_graph_gs(gs, g, flow, edge_specifier, edge_type, s, t);
//    
    
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
    logger::Instance()->debug("Guide before.\n");
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    
    path_finder* pf = path_finder::create_path_finder(g, s, t, flow, means, edge_specifier, edge_type, edge_lengths, node_index, know_paths, know_back_paths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, raw->size, meta);
    pf->extract_guided_transcripts(transcripts, raw->guide_transcripts);    // test and extract evidenced transcripts 
    delete pf; 
    

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Guide after.\n");
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    
//
//    if (nid == 3) { // only have s and t! all arcs are parallel transcripts
//        float max_score = 0;
//        for (ListDigraph::ArcIt a(g); a != INVALID;++a) {
//            float score = means[a].compute_score();
//            if (score > max_score) {
//                max_score = score;
//            }
//        }
//        for (ListDigraph::ArcIt a(g); a != INVALID;++a) {
//            float score = means[a].compute_score();
//            if (score < max_score * 0.02) {
//                g.erase(a);
//            }
//        }
//    }
//    
    ListDigraph::ArcMap<bool> barred_high(g);
    while(simplify_ambiguous_ILP(unsecurityArc, unsecurityId, know_paths, know_back_paths, in_deg, out_deg, barred_high, false));
    
//    gs << "Resolved graph ---\n";
//    print_graph_gs(gs, g, flow, edge_specifier, edge_type, s, t);
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Simplify done.\n");
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
    #endif

    // we have a simplified graph that we can now extract transcripts from
    // at this point we get the maximum amount of transcripts possible out of the existing DAG
    pf = path_finder::create_path_finder(g, s, t, flow, means, edge_specifier, edge_type, edge_lengths, node_index, know_paths, know_back_paths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, raw->size, meta);
//    pf->extract_guided_transcripts(transcripts, raw->guide_transcripts);    // test and extract evidenced transcripts 
    pf->extract_transcripts(transcripts); // normal path extraction
    delete pf; 
    
    // set exon positions for unsecurities
    for (graph_list<lazy<transcript> >::iterator t = transcripts.transcripts.begin(); t != transcripts.transcripts.end(); ++t) {
        
        if (!(*t).has_ref()) {
            continue;
        }
        
        for ( std::deque<transcript_unsecurity>::iterator u = (*t)->unsecurity_id.begin(); u!= (*t)->unsecurity_id.end(); ++u) {
            u->position = unsecIdToNodeIndex[u->id].id;
            u->left = unsecIdToNodeIndex[u->id].left;
            u->right = unsecIdToNodeIndex[u->id].right;
        }
    }

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Pathfinder done " + std::to_string(transcripts.size())  + ".\n");
    print_graph_debug_copy(std::cout, g, capacity, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    
}


capacity_type base_manager::DKMeans( std::vector<capacity_type> &in) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Inner Custer Group.\n");
    #endif

    if (in.empty()) {
        return 0;
    }
    if (in.size() == 1) {
         return 0;
    }
    
    if (in.size() == 2) {
         return in[0];
    }
    
    //################ FORWARD ################
    
    unsigned int nUnique = 1;
    for (std::vector<capacity_type>::iterator itr=in.begin()+1; itr!=in.end(); ++itr) {
        if (*itr != *(itr -1)) {
            nUnique++;
        }
    }
    
    if (nUnique == 1) {
        return in[1];
    }
    
    const unsigned int N_size = in.size();
    const unsigned int K_size = std::min(nUnique, 4u);
    
    // 2D dynamic lookup table arrays
    std::vector< std::vector< double > > costs( (K_size), std::vector<double>(N_size) );
    std::vector< std::vector< unsigned int > > back( (K_size), std::vector<unsigned int>(N_size) );
    
    for(unsigned int i = 0; i < K_size; ++i) {
            costs[i][0] = 0.0;
            back[i][0] = 0;
    }
    
    // temporary mean values
    double mean_x1, mean_xj;
    double d;
    
    for(unsigned int k = 0; k < K_size; ++k) {    
        
        mean_x1 = in[0];
        
        for(unsigned int i = std::max(1ul,(unsigned long)k); i < N_size; ++i) {
            
            if (k == 0) {
                costs[0][i] = costs[0][i-1] + i / (double) (i + 1 ) *
                          (in[i] - mean_x1) * (in[i] - mean_x1);
                mean_x1 = (i * mean_x1 + in[i]) / (double) (i + 1);
                back[0][i] = 0;
            } else {
                costs[k][i] = -1.0;
                d = 0.0;
                mean_xj = 0.0;
                
                for(unsigned int j = i; j >= k; --j) {
                    
                    d = d + (i - j) / (double) (i - j + 1) * 
                            (in[j] - mean_xj) * (in[j] - mean_xj);
                    mean_xj = (in[j] + (i - j) * mean_xj) / (double)(i - j + 1);

                    if(costs[k][i] == -1.0) { //initialization of D[k,i]
                        
                        if(j == 1) {
                            costs[k][i] = d;
                            back[k][i] = j;
                        } else {
                            costs[k][i] = d + costs[k - 1][j - 1];
                            back[k][i] = j;
                        }
                    } else {
                        
                        if(j == 1) {
                            if(d <= costs[k][i]) {
                                costs[k][i] = d;
                                back[k][i] = j;
                            }
                        } else {
                            if(d + costs[k - 1][j - 1] < costs[k][i]) {
                                costs[k][i] = d + costs[k - 1][j - 1];
                                back[k][i] = j;
                            }
                        }
                    }

                }
                
            }
        }
    }
    
//    logger::Instance()->debug("Matrix \n");
//    for(unsigned int k = 0; k < K_size; ++k) {    
//        for(unsigned int i = k; i < N_size; ++i) {
//            logger::Instance()->debug(" " + std::to_string(costs[k][i]) + ";" + std::to_string(back[k][i]));
//        }
//        logger::Instance()->debug("\n");
//    }
    
    
    //################ Backwards ################
    
    // find smallest number of clusters k that contains all points
    unsigned int range_start;
    unsigned int seperator;
    unsigned int range_end = N_size - 1;
    for ( int k = K_size - 1; k>=0; --k) {
        range_start = back[k][range_end];
        
        double sum = 0.0;
        for(unsigned int i = range_start; i <= range_end; ++i)
            sum += in[i];

//        logger::Instance()->debug("CENTERX " + std::to_string(sum/(range_end-range_start+1)) + "\n");
//        logger::Instance()->debug("CLUSTERX " + std::to_string(range_start) + " - " + std::to_string(range_end) + "\n");
        
        seperator = range_end;
        
        range_end = range_start - 1;
    }
    
    return in[seperator];
}

void base_manager::find_guide_sources_and_drains(ListDigraph::ArcMap<bool> &guided_saves, std::deque<std::deque<ListDigraph::Arc> > &paths) {
    
    for (graph_list<exon_group *>::iterator eg_it = raw->guide_transcripts.begin(); eg_it != raw->guide_transcripts.end(); ++eg_it) {
        
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
        if (edge_type[a] == edge_types::NODE) { //with opened middle node
            int ni = edge_specifier[a].node_index;
            same_node[ni].push_back(a); 
        }
    }
    
    
    for (std::set<int>::iterator src = sources_to_test.begin(); src != sources_to_test.end(); ++src ) {
        
        ListDigraph::Arc node2 = g.arcFromId(*src);
        
        if (edge_type[node2] != edge_types::NODE) { // cannot happen if we end
            continue;
        }
        
        // do we have a non-source adjourning node?
        bool non_source_adjourning = false;
        for (std::deque<ListDigraph::Arc>::iterator nr = same_node[edge_specifier[node2].node_index].begin(); nr != same_node[edge_specifier[node2].node_index].end(); ++nr) {
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
        for (std::deque<ListDigraph::Arc>::iterator nr = same_node[edge_specifier[node2].node_index].begin(); nr != same_node[edge_specifier[node2].node_index].end(); ++nr) {
            for (ListDigraph::InArcIt i(g, g.source(*nr) ); i != INVALID; ++i) {
                if (marked_source[i] && edge_type[i] == edge_types::EXON) {

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

                    if (edge_specifier[edge1].right_consecutive) {
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
        
        if (edge_type[node2] != edge_types::NODE) { // cannot happen if we end
            continue;
        }
        
        // do we have a non-source adjourning node?
        bool non_drain_adjourning = false;
        for (std::deque<ListDigraph::Arc>::iterator nr = same_node[edge_specifier[node2].node_index].begin(); nr != same_node[edge_specifier[node2].node_index].end(); ++nr) {
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
        for (std::deque<ListDigraph::Arc>::iterator nr = same_node[edge_specifier[node2].node_index].begin(); nr != same_node[edge_specifier[node2].node_index].end(); ++nr) {
            for (ListDigraph::OutArcIt o(g, g.target(*nr) ); o != INVALID; ++o) {
                if (marked_drain[o] && edge_type[o] == edge_types::EXON) {

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

                    if (edge_specifier[edge1].left_consecutive) {               
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
        if (edge_type[a] == edge_types::NODE) { //with opened middle node
            
            int ni = edge_specifier[a].node_index;
            same_node_r[ni].push_back(g.target(a));
            same_node_l[ni].push_back(g.source(a));
            
        }
    }
    
    // test for in between neighbours 
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        if (guided_saves[a]) {
            continue;
        }
        
        if (edge_type[a] == edge_types::NODE) { //with opened middle node

            ListDigraph::InArcIt ni(g, g.source(a));
            ListDigraph::OutArcIt no(g, g.target(a));

            ListDigraph::Arc left_edge(ni);
            ListDigraph::Arc rigth_edge(no);

            ++ni;
            ++no;
            
            if (ni == INVALID && no == INVALID && edge_type[left_edge] == edge_types::EXON && edge_type[rigth_edge] == edge_types::EXON) {
                // first condition: only one connecting edge in each direction, otherwise this can't be an intron

                // each edge can have but one node
                ListDigraph::InArcIt left_node(g, g.source(left_edge));
                ListDigraph::OutArcIt right_node(g, g.target(rigth_edge));

                // check if this is really neighbouring and if bot exons are actually connected!
                if (meta->exons[edge_specifier[left_node].node_index].right + 1 == meta->exons[edge_specifier[a].node_index].left 
                        && meta->exons[edge_specifier[a].node_index].right + 1 == meta->exons[edge_specifier[right_node].node_index].left) {

                    bool breaked = false;
                    // look for direct connection, a single edge has to exist or intron would be further split, which is legally impossible
                    for (std::deque<ListDigraph::Node>::iterator nr = same_node_r[edge_specifier[left_node].node_index].begin(); nr != same_node_r[edge_specifier[left_node].node_index].end() && !breaked; ++nr) {
                        ListDigraph::OutArcIt search(g, *nr); // search is legally an edge and con only have one node output!
                        for (;search != INVALID && !breaked; ++search) {
                            
                            bool hit_target = false;
                            for (std::deque<ListDigraph::Node>::iterator nl = same_node_l[edge_specifier[right_node].node_index].begin(); nl != same_node_l[edge_specifier[right_node].node_index].end(); ++nl) {
                                
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
                                logger::Instance()->debug("Test " + std::to_string(g.id(a)) + " " + std::to_string(regions[a].get_average()) + " " + std::to_string(capacity[search]) + "\n");
                                #endif  

                                if (capacity[search] != 0 && ( (capacity[a] < 45 && regions[a].get_average()*100 / capacity[search] < options::Instance()->get_intron_retention_threshold()) || (capacity[a] <= 2) ) )
                                {

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
        } else if (edge_type[a] == edge_types::EXON) {
            if (edge_specifier[a].id.count() == 3) { // only then can middle node be the intron!
                
                unsigned int left = edge_specifier[a].id.find_first();
                unsigned int middle = edge_specifier[a].id.find_next(left);
                unsigned int right = edge_specifier[a].id.find_next(middle);
                               
                if (meta->exons[left].right + 1 == meta->exons[middle].left 
                        && meta->exons[middle].right + 1 == meta->exons[right].left) {
                    
                    
                    ListDigraph::OutArcIt search(g, g.source(a));
                    for (;search != INVALID; ++search) {
                         
                        if (g.target(search) == g.target(a) && search != a) {
                    
                            ListDigraph::InArcIt left_node(g, g.source(a));
                            ListDigraph::OutArcIt right_node(g, g.target(a));

                            if ( capacity[search] != 0 && capacity[a] < 45 && capacity[a]*100 / capacity[search] < options::Instance()->get_intron_retention_threshold() )
                            {
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
        if (edge_type[a] == edge_types::NODE) { //with opened middle node
            
            int ni = edge_specifier[a].node_index;
            same_node_r[ni].push_back(g.target(a));
            same_node_l[ni].push_back(g.source(a));
            
        }
    }    
        
    ListDigraph::ArcMap<bool> kill_status(g);
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (edge_type[a] == edge_types::EXON && !marked_drain[a] && !marked_source[a]) {
            
            ListDigraph::InArcIt left_node(g, g.source(a));
            ListDigraph::OutArcIt right_node(g, g.target(a));
            
            
            if (left_node != INVALID && right_node != INVALID && edge_type[left_node] == edge_types::NODE && edge_type[right_node] == edge_types::NODE ) {
//                && edge_specifier[right_node].node_index == edge_specifier[left_node].node_index + 3 ) {
            //        && (!marked_drain[left_node] && !marked_drain[right_node]) && (!marked_source[left_node] && !marked_source[right_node]) {
            
                bool has_source = false;
                ListDigraph::Arc source, source_base;

                bool has_drain = false;
                ListDigraph::Arc drain, drain_base;

                for (std::deque<ListDigraph::Node>::iterator nl = same_node_l[edge_specifier[right_node].node_index].begin(); nl != same_node_l[edge_specifier[right_node].node_index].end(); ++nl) {
                    for (ListDigraph::InArcIt oi(g, *nl); oi != INVALID ;++oi) {
                        ListDigraph::InArcIt lr(g, g.source(oi));

                        if (marked_source[oi] && lr != INVALID && meta->exons[edge_specifier[lr].node_index].right + 1 == meta->exons[edge_specifier[right_node].node_index].left) {
                            has_source = true;
                            source = oi;
                            source_base = lr;
                        }
                    }
                }
                for (std::deque<ListDigraph::Node>::iterator nr = same_node_r[edge_specifier[left_node].node_index].begin(); nr != same_node_r[edge_specifier[left_node].node_index].end(); ++nr) {
                    for (ListDigraph::OutArcIt oi(g, *nr); oi != INVALID ;++oi) {
                        ListDigraph::OutArcIt rl(g, g.target(oi));

                        if (marked_drain[oi] && rl != INVALID && meta->exons[edge_specifier[left_node].node_index].right + 1 == meta->exons[edge_specifier[rl].node_index].left) {
                            has_drain = true;
                            drain = oi;
                            drain_base = rl;
                        }
                    }
                }
                
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Broken Intron Test " + std::to_string(g.id(a)) +" s " + std::to_string(g.id(source)) +" d "+ std::to_string(g.id(drain)) + "\n");
                #endif

                if (has_drain && has_source &&  capacity[source] < capacity[a] && capacity[drain] < capacity[a]) {
                // likely we can kill both

                    rcount source_cap = capacity[source];
                    rcount drain_cap = capacity[drain];
                    
                    rcount source_cap_base = means[source_base].mean;
                    rcount drain_cap_base = means[drain_base].mean;
                    rcount min = std::min(source_cap_base, drain_cap_base);
                    rcount max = std::max(source_cap_base, drain_cap_base);
                    
                    if (source_cap_base > 30 && drain_cap_base > 30) {
                        continue;
                    }
                    
                    if ( (source_cap_base == drain_cap_base || (max - min)*100/(max) < 30 || (max < 15 && min > 5)) && capacity[a] !=0 && source_cap_base *100 / capacity[a] < options::Instance()->get_broken_intron_retention_threshold() 
                            && drain_cap_base*100 / capacity[a] < options::Instance()->get_broken_intron_retention_threshold()) {
                    
                            kill_status[source] = true;
                            kill_status[drain] = true;

                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Erase Source and Drain Broken Intron " + std::to_string(g.id(source)) + " " + std::to_string(g.id(drain)) + "\n");
                            #endif
                        
                    } else if (source_cap_base > drain_cap_base && drain_cap_base*100 / capacity[a] < options::Instance()->get_broken_intron_retention_threshold()) {
                        
//                        capacity[source] = source_cap - drain_cap;
//                        capacity[source_base] = capacity[source];  
                        
                        kill_status[drain] = true;
                        
                        #ifdef ALLOW_DEBUG
                         logger::Instance()->debug("Equalize source\n");
                        #endif
                        
                    } else if (source_cap_base < drain_cap_base && source_cap_base *100 / capacity[a] < options::Instance()->get_broken_intron_retention_threshold()) {
                        
//                        capacity[drain] = drain_cap - source_cap;
//                        capacity[drain_base] = capacity[drain];  
                        
                        kill_status[source] = true;
                        #ifdef ALLOW_DEBUG
                         logger::Instance()->debug("Equalize drain\n");
                        #endif
                    }
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


void base_manager::prune_sources_and_drains5(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &push_block) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Prune 5 \n");
    #endif

    struct sd_mark {
        
        sd_mark() {}
        sd_mark(ListDigraph::Arc a, ListDigraph::Arc n, bool ne) : arc(a), node(n), neighbouring(ne) 
            {}

        ListDigraph::Arc arc;
        ListDigraph::Arc node;
        bool neighbouring;
    };
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        if (edge_type[a] != edge_types::NODE) {
            continue;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Test " + std::to_string(g.id(a)) + "\n");
        #endif
        
        /// ---------- SOURCE ----------
    
        if(!marked_source[a]) {
            std::deque<sd_mark> sources;
            rcount arc_count = 0;
            rcount node_count = 0;
            bool has_neighbouring = false;
            for (ListDigraph::InArcIt oi(g, g.source(a)); oi != INVALID ;++oi) {

                if (edge_type[oi] == edge_types::HELPER) {
                    // helper can't possibly be disproven/proven
                    continue;
                }

                ListDigraph::InArcIt node(g, g.source(oi));
                if (node == INVALID) {
                    continue;
                }

                if (marked_source[oi]) {
                    bool neighbouring = meta->exons[edge_specifier[a].node_index].left == meta->exons[edge_specifier[node].node_index].right + 1;
                    sources.push_back(sd_mark(oi, node, neighbouring));
                    has_neighbouring = has_neighbouring || neighbouring;
                } else {
                    arc_count += capacity[oi];
                    node_count += capacity[node];
                }
            }

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Sloop done " + std::to_string(node_count) + " " + std::to_string(arc_count) + "\n");
            for (std::deque<sd_mark>::iterator a_it = sources.begin(); a_it != sources.end(); ++a_it) {
                logger::Instance()->debug("Source Adjoint " + std::to_string(g.id(a_it->arc)) + " " + std::to_string(g.id(a_it->node)) + " " + std::to_string(a_it->neighbouring) + "\n");
            }
            #endif

            if (sources.size() > 1) {
                // multiple_source_case
                 rcount max_node_b = 0;
                rcount max_arc_b = 0;
                rcount max_node_n = 0;
                rcount max_arc_n = 0;
                std::deque<sd_mark>::iterator n_it;
                for (std::deque<sd_mark>::iterator a_it = sources.begin(); a_it != sources.end(); ++a_it) {
                    if (a_it->neighbouring) {
                        n_it = a_it;

                        if (max_node_n < capacity[a_it->node]) {
                            max_node_n = capacity[a_it->node];
                        }
                        if (max_arc_n < capacity[a_it->arc]) {
                            max_arc_n = capacity[a_it->arc];
                        }
                    } else {
                        if (max_node_b < capacity[a_it->node]) {
                            max_node_b = capacity[a_it->node];
                        }
                        if (max_arc_b < capacity[a_it->arc]) {
                            max_arc_b = capacity[a_it->arc];
                        }
                    }    
                }
                rcount max_node = std::max(max_node_n, max_node_b);
                if (max_node < node_count) {
                    max_node = node_count;
                }
                rcount max_arc = std::max(max_arc_n, max_arc_b);
                if (max_arc < arc_count) {
                    max_arc = arc_count;
                }

                for (std::deque<sd_mark>::iterator a_it = sources.begin(); a_it != sources.end(); ++a_it) {
                    
                    if (guided_saves[a_it->arc] || guided_saves[a_it->node]) {
                        continue;
                    }
                    
                    if ( (max_node >=2 && capacity[a_it->node] < 2) || (max_node >=4 && capacity[a_it->node] < 3) || (max_node >=8 && capacity[a_it->node] < 4)
                            || (max_arc >=3 && capacity[a_it->arc] < 2) || (max_arc >=6 && capacity[a_it->arc] < 3) || (max_arc >=10 && capacity[a_it->arc] < 4) 
                            || capacity[a_it->node] * 100 / float(max_node) < options::Instance()->get_arc_filter()
                            || (a_it == n_it && max_arc <= 3 && (max_arc_b >3 || max_arc_b >= max_arc_n) && (max_node_b > 3 || max_node_b >= max_node_n)) ) {

                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Erase Source dn5 " + std::to_string(g.id(a_it->arc)) +"\n");
                        #endif
                        g.erase(a_it->arc);
                        g.erase(a_it->node);

                    }
                } 
            } else if (sources.size() == 1 && !guided_saves[sources.begin()->arc] && !guided_saves[sources.begin()->node]) {
                if ((capacity[sources.begin()->node] <= 3 && node_count > 3 )
                        || capacity[sources.begin()->node] * 100 / float(node_count+capacity[sources.begin()->node]) < options::Instance()->get_arc_filter() ) {

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Erase Source2 dn5 " + std::to_string(g.id(sources.begin()->arc)) +"\n");
                    #endif
                    g.erase(sources.begin()->arc);
                    g.erase(sources.begin()->node);
                }
            } 
        }
        
        /// ---------- DRAIN ----------
        if (!marked_drain[a]) {
            std::deque<sd_mark> drains;
            rcount arc_count = 0;
            rcount node_count = 0;
            bool has_neighbouring = false;
            for (ListDigraph::OutArcIt oi(g, g.target(a)); oi != INVALID ;++oi) {

                if (edge_type[oi] == edge_types::HELPER) {
                    // helper can't possibly be disproven/proven
                    continue;
                }

                ListDigraph::OutArcIt node(g, g.target(oi));
                if (node == INVALID) {
                    continue;
                }

                if (marked_drain[oi]) {
                    bool neighbouring = meta->exons[edge_specifier[node].node_index].left == meta->exons[edge_specifier[a].node_index].right + 1;
                    drains.push_back(sd_mark(oi, node, neighbouring));
                    has_neighbouring = has_neighbouring || neighbouring;
                } else {
                    arc_count += capacity[oi];
                    node_count += capacity[node];
                }
            }

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Dloop done " + std::to_string(node_count) + " " + std::to_string(arc_count) + "\n");
            for (std::deque<sd_mark>::iterator a_it = drains.begin(); a_it != drains.end(); ++a_it) {
                logger::Instance()->debug("Drain Adjoint " + std::to_string(g.id(a_it->arc)) + " " + std::to_string(g.id(a_it->node)) + " " + std::to_string(a_it->neighbouring) + "\n");
            }
            #endif

            if (drains.size() > 1) {
                // multiple_source_case
                rcount max_node_b = 0;
                rcount max_arc_b = 0;
                rcount max_node_n = 0;
                rcount max_arc_n = 0;
                std::deque<sd_mark>::iterator n_it;
                for (std::deque<sd_mark>::iterator a_it = drains.begin(); a_it != drains.end(); ++a_it) {
                    if (a_it->neighbouring) {
                        n_it = a_it;

                        if (max_node_n < capacity[a_it->node]) {
                            max_node_n = capacity[a_it->node];
                        }
                        if (max_arc_n < capacity[a_it->arc]) {
                            max_arc_n = capacity[a_it->arc];
                        }
                    } else {
                        if (max_node_b < capacity[a_it->node]) {
                            max_node_b = capacity[a_it->node];
                        }
                        if (max_arc_b < capacity[a_it->arc]) {
                            max_arc_b = capacity[a_it->arc];
                        }
                    }    
                }
                rcount max_node = std::max(max_node_n, max_node_b);
                if (max_node < node_count) {
                    max_node = node_count;
                }
                rcount max_arc = std::max(max_arc_n, max_arc_b);
                if (max_arc < arc_count) {
                    max_arc = arc_count;
                }

                for (std::deque<sd_mark>::iterator a_it = drains.begin(); a_it != drains.end(); ++a_it) {
                    
                    if (guided_saves[a_it->arc] || guided_saves[a_it->node]) {
                        continue;
                    }
                    
                    if ( (max_node >=2 && capacity[a_it->node] < 2) || (max_node >=4 && capacity[a_it->node] < 3) || (max_node >=8 && capacity[a_it->node] < 4)
                            || (max_arc >=3 && capacity[a_it->arc] < 2) || (max_arc >=6 && capacity[a_it->arc] < 3) || (max_arc >=10 && capacity[a_it->arc] < 4) 
                            || capacity[a_it->node] * 100 / float(max_node) < options::Instance()->get_arc_filter()
                            || (a_it == n_it && max_arc <= 3 && (max_arc_b >3 || max_arc_b >= max_arc_n) && (max_node_b > 3 || max_node_b >= max_node_n)) ) {

                        #ifdef ALLOW_DEBUG
                        logger::Instance()->info("Erase Drain dn5 " + std::to_string(g.id(a_it->arc)) +"\n");
                        #endif
                        g.erase(a_it->arc);
                        g.erase(a_it->node);
                    }
                } 
            } else if (drains.size() == 1 && !guided_saves[drains.begin()->arc] && !guided_saves[drains.begin()->node]) {
                if ( (capacity[drains.begin()->node] <= 3 && node_count > 3)
                        || capacity[drains.begin()->node] * 100 / float(node_count+capacity[drains.begin()->node]) < options::Instance()->get_arc_filter()) {
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->info("Erase Drain2 dn5 " + std::to_string(g.id(drains.begin()->arc)) +"\n");
                    #endif
                    g.erase(drains.begin()->arc);
                    g.erase(drains.begin()->node);
                }
            } 
        }
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Prune 5 Out\n");
    #endif
}


void base_manager::prune_sources_and_drains4(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &push_block) {
    
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        if (edge_type[a] != edge_types::NODE || marked_source[a] || marked_drain[a]) {
            continue;
        }
        
        bool has_dubious_source = false;
        ListDigraph::Arc dubious_source;
        for (ListDigraph::InArcIt oi(g, g.source(a)); oi != INVALID ;++oi) {
            if (!guided_saves[oi] && marked_source[oi]) {
                    has_dubious_source = true;
                    dubious_source = oi;
            }
        }
        
        bool has_dubious_drain = false;
        ListDigraph::Arc dubious_drain;
        
        for (ListDigraph::OutArcIt oi(g, g.target(a)); oi != INVALID ;++oi) {
            if (!guided_saves[oi] && marked_drain[oi]) {
                    has_dubious_drain = true;
                    dubious_drain = oi;
            }
        }
        
        if (capacity[a] < 11) {
            
             if (has_dubious_source) {
                g.erase(dubious_source);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Erase Source d4l " + std::to_string(g.id(dubious_source)) + "\n");
                #endif
            }
            if (has_dubious_drain) {
                g.erase(dubious_drain);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Erase Drain d4l " + std::to_string(g.id(dubious_drain)) + "\n");
                #endif
            }
        } 

//        if (has_dubious_source && has_dubious_drain) {
//            
//            rcount source_cap = capacity[dubious_source];
//            rcount drain_cap = capacity[dubious_drain];
//            rcount min = std::min(source_cap, drain_cap);
//            rcount max = std::max(source_cap, drain_cap);
//            
//            
//            if ( (source_cap == drain_cap || (max - min)*100/(max) < 30 ) &&  regions[a].get_left()!=0 && regions[a].get_right()!=0 && source_cap*100 / regions[a].get_left() < options::Instance()->get_intron_retention_threshold() 
//                            && drain_cap*100 / regions[a].get_right() < options::Instance()->get_intron_retention_threshold()) {
//                    
//                g.erase(dubious_source);
//                g.erase(dubious_drain);
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Erase Both d4 " + std::to_string(g.id(dubious_source)) + "\n");
//                #endif
//                continue;
//            } else if (source_cap > drain_cap) {
//                ListDigraph::InArcIt left_node(g, g.source(dubious_source));
//                capacity[dubious_source] = source_cap - drain_cap;
//                capacity[left_node] = capacity[dubious_source];    
//                g.erase(dubious_drain);
//                has_dubious_drain = false;
//                #ifdef ALLOW_DEBUG
//                 logger::Instance()->debug("Equalize source d4\n");
//                #endif
//                
//            } else {
//                ListDigraph::InArcIt right_node(g, g.target(dubious_drain));
//                capacity[dubious_drain] = drain_cap - source_cap;
//                capacity[right_node] = capacity[dubious_drain];    
//                g.erase(dubious_source);
//                has_dubious_source = false;
//                #ifdef ALLOW_DEBUG
//                 logger::Instance()->debug("Equalize drain d4\n");
//                #endif
//            }
//        }
//        
//        if (has_dubious_source) {
//            rcount left_count = 0;
//            rcount left_capacity = 0;
//            
//            for (ListDigraph::InArcIt edge(g, g.source(a)); edge != INVALID ;++edge) {
//                ListDigraph::InArcIt node(g, g.source(edge));
//                if ( node != INVALID && edge != dubious_source ) {
//                    left_count += regions[node].get_right();
//                    left_capacity += capacity[node];
//                }
//            }
//                    
//            rcount right_count = regions[a].get_left();
//            rcount right_capacity = capacity[a];    
//
//            if (left_count > right_count || left_capacity > right_capacity) {
//                g.erase(dubious_source);
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Erase Source d4 " + std::to_string(g.id(dubious_source)) + "\n");
//                #endif
//            }   
//        }
//        
//        if (has_dubious_drain) {
//            rcount right_count = 0;
//            rcount right_capacity = 0;
//            
//            for (ListDigraph::OutArcIt edge(g, g.target(a)); edge != INVALID ;++edge) {
//                ListDigraph::OutArcIt node(g, g.target(edge));
//                if ( node != INVALID && edge != dubious_drain ) {
//                    right_count += regions[node].get_left();
//                    right_capacity += capacity[node];
//                }
//            }
//            
//            rcount left_count = regions[a].get_right();
//            rcount left_capacity = capacity[a];    
//            if (right_count > left_count || right_capacity > left_capacity) {
//                g.erase(dubious_drain);
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Erase Drain d4 " + std::to_string(g.id(dubious_drain)) + "\n");
//                #endif
//            }   
//        }
    }
}

void base_manager::prune_sources_and_drains3(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &push_block) {
    
    // we want to get the median left and median right capacity of each node-arc
    std::deque<ListDigraph::Node> top_order;
    pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > visitor(top_order);
    DfsVisit<ListDigraph, pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > > dfs(g, visitor);
    dfs.init();
    dfs.addSource(s);
    dfs.start();
    
    // now do the actual checking with median
    for (std::deque<ListDigraph::Node>::iterator n = top_order.begin(); n != top_order.end(); ++n) {
        
        
        ListDigraph::OutArcIt onnt(g, *n);
        bool node_right = (onnt != INVALID && edge_type[onnt] == edge_types::NODE);
        
        if (!node_right || (push_block[onnt] && !guided_saves[onnt]) ) {
            continue;
        }
        
        bool has_dubious_source = false;
        ListDigraph::Arc dubious_source;
        
        bool has_dubious_drain = false;
        ListDigraph::Arc dubious_drain;
        
        for (ListDigraph::InArcIt oi(g, g.source(onnt)); oi != INVALID ;++oi) {
            if (!guided_saves[oi] && push_block[oi] && marked_source[oi]) {
                    has_dubious_source = true;
                    dubious_source = oi;
            }
        }
        
        for (ListDigraph::OutArcIt oi(g, g.target(onnt)); oi != INVALID ;++oi) {
            if (!guided_saves[oi] && push_block[oi] && marked_drain[oi]) {
                    has_dubious_drain = true;
                    dubious_drain = oi;
            }
        }
             
//        if (has_dubious_source && capacity[dubious_source] <= 2) {
//            logger::Instance()->debug("Erase dubious source by capacity: " + std::to_string(g.id(dubious_source)) + "\n");
//            g.erase(dubious_source);
//            has_dubious_source = false;
//        }
//        if (has_dubious_drain && capacity[dubious_drain] <= 2) {
//            logger::Instance()->debug("Erase dubious drain by capacity: " + std::to_string(g.id(dubious_drain)) + "\n");
//            g.erase(dubious_drain);
//            has_dubious_drain = false;
//        }
//        
        
        if (!has_dubious_source && !has_dubious_drain) {
            continue;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Test Arc Median: " + std::to_string(g.id(onnt)) + "\n");
        #endif
        
        // functional node, get valid iterator in Top Sorting!
        std::deque<ListDigraph::Node>::reverse_iterator left(n);
        std::deque<ListDigraph::Node>::reverse_iterator left_end = top_order.rend();
        std::deque<ListDigraph::Node>::iterator right = n;
        std::deque<ListDigraph::Node>::iterator right_end = top_order.end();
        ++right;
        --left;
        
        rcount median_right = get_forward_median(right, right_end, guided_saves, push_block, marked_drain);
        rcount median_left = get_backward_median(left, left_end, guided_saves, push_block, marked_source);
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Medians: " + std::to_string(median_left) + " " + std::to_string(median_right) + "\n");
        #endif
        
	bool erase_drain = false;
	bool erase_source = false;        

	if (median_left != 0 && median_right != 0) {
		// we actually have a valid position here
		if (median_left < median_right || (median_left - median_right) * 100 / float(median_right) < options::Instance()->get_coverage_change_limit() ) {
		    //this does NOT support a drain
		   erase_drain = true;
		} 
		
		if (median_left > median_right || (median_right - median_left) * 100 / float(median_left) < options::Instance()->get_coverage_change_limit() ) {
		    //this does NOT support a SOURCE
		   erase_source = true;
		} 
	}

        if (has_dubious_drain) {
            if (erase_drain) {
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Erase Drain " + std::to_string(g.id(dubious_drain)) + "\n");
                    #endif
                    g.erase(dubious_drain);
            } else {
                rcount other_side = 0;
                for (ListDigraph::OutArcIt edge(g, *right); edge != INVALID ;++edge) {
                    ListDigraph::OutArcIt node(g, g.target(edge));
                    if ( node != INVALID && edge != dubious_drain ) {
                        other_side += regions[node].get_left();
                        }
                    }
                rcount right = regions[onnt].get_right();
                if (other_side > right) {
                    g.erase(dubious_drain);
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Erase Drain2 " + std::to_string(g.id(dubious_drain)) + "\n");
                    #endif
                } else if (other_side != 0) {
                    capacity[dubious_drain] = std::min(right-other_side, capacity[dubious_drain]);
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Lower Capacity Drain " + std::to_string(g.id(dubious_drain)) + "\n");
                    #endif
                }
            }
        }

        if (has_dubious_source) {
            if (erase_source) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Erase Source " + std::to_string(g.id(dubious_source)) + "\n");
                #endif
                g.erase(dubious_source);
            } else {
                rcount other_side = 0;
                for (ListDigraph::InArcIt edge(g, *left); edge != INVALID ;++edge) {
                    ListDigraph::InArcIt node(g, g.source(edge));
                    if ( node != INVALID && edge != dubious_source ) {
                        other_side += regions[node].get_right();
                    }
                }
                rcount left = regions[onnt].get_left();
                if (other_side > left) {
                    g.erase(dubious_source);
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Erase Source 2 " + std::to_string(g.id(dubious_source)) + "\n");
                    #endif
                } else if (other_side != 0) {
                    capacity[dubious_source] = std::min(left-other_side, capacity[dubious_source]);
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Lower Capacity Source " + std::to_string(g.id(dubious_source)) + "\n");
                    #endif
                } 
            }
        }
    } 
}

rcount base_manager::get_forward_median(std::deque<ListDigraph::Node>::iterator &start, std::deque<ListDigraph::Node>::iterator &end,
        ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &push_block, ListDigraph::ArcMap<bool> &marked_drain) {
    
    std::set<int> contributing_nodes;
    ListDigraph::ArcMap<float> contribution_factor(g); // init as 0
    ListDigraph::ArcMap<rcount> bias(g); // init as 0
    ListDigraph::ArcMap<bool> visited(g);

    ListDigraph::InArcIt sn1(g, *start);
    contribution_factor[sn1] = 1;    
    visited[sn1] = true;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Entry: " + std::to_string(g.id(sn1)) + "\n");
    #endif
    
    for (std::deque<ListDigraph::Node>::iterator n = start; n!= end; ++n) { 
        // first test if we can propagate
        ListDigraph::InArcIt innt(g, *n);
        bool node_left= (innt != INVALID && edge_type[innt] == edge_types::NODE);
        if (node_left && visited[innt]) {
              for (ListDigraph::OutArcIt oi(g, *n); oi != INVALID ;++oi) {
                  ListDigraph::OutArcIt node(g, g.target(oi));
                  if ( node != INVALID && capacity[node] != 0 && !marked_drain[node] && (guided_saves[node] || !push_block[node]) ) {
                      visited[node] = true;
                  }
              }
	}
    }
    
    for (std::deque<ListDigraph::Node>::iterator n = start; n!= end; ++n) {

	ListDigraph::InArcIt innt(g, *n);
	bool node_left= (innt != INVALID && edge_type[innt] == edge_types::NODE);
        if (node_left && visited[innt]) {
            
            if (innt != sn1) {
                
                rcount b = 0;		
                for (ListDigraph::InArcIt in_edge(g, g.source(innt)); in_edge != INVALID ;++in_edge) {
                    ListDigraph::InArcIt in_node(g, g.source(in_edge));
                    if (in_node != INVALID && !visited[in_node]) { //&& (guided_saves[in_node] || !push_block[in_node]) ) {
                            b += capacity[in_edge];
                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Bias From: " + std::to_string(g.id(in_edge)) + "\n");
                            #endif
                    }
                }

                bias[innt] += b; 
                bias[innt] = std::min(bias[innt], capacity[innt]);
            }
            
            std::vector<ListDigraph::Arc> adjourning;
            rcount total_edge_count = 0;
            for (ListDigraph::OutArcIt oi(g, *n); oi != INVALID ;++oi) {
                // all edges have by definition a node handy
                ListDigraph::OutArcIt node(g, g.target(oi));
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Node1: " + std::to_string(g.id(innt))  + " Edge1: " + std::to_string(g.id(oi)) + " Node2: " + std::to_string(g.id(node)) + " " + std::to_string(capacity[oi])+ "\n");
                #endif
                
                if ( node != INVALID && capacity[node] != 0 && visited[node]) { //&& (guided_saves[node] || !push_block[node]) ) {
                    total_edge_count+= capacity[oi];
                    adjourning.push_back(oi);
                }
            }
            
            if (total_edge_count == 0) {
                continue;
            }
            
            for (std::vector<ListDigraph::Arc>::iterator edge = adjourning.begin(); edge != adjourning.end(); ++edge) {
                
                ListDigraph::OutArcIt target(g, g.target(*edge));
                                
                contribution_factor[target] += contribution_factor[innt] * capacity[*edge]/total_edge_count;
                bias[target] += bias[innt] * capacity[*edge]/total_edge_count;
                
                contributing_nodes.insert(g.id(target));
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("fwd: " + std::to_string(g.id(target)) + " " + std::to_string(contribution_factor[target]) + " " + std::to_string(bias[target]) + "\n"); 
                #endif
            }
        }
    }
    
    if (contributing_nodes.empty()) {
        return 0;
    }
    
    std::vector<rcount> median;
    for (std::set<int>::iterator ni = contributing_nodes.begin(); ni != contributing_nodes.end(); ++ni) {
	
       ListDigraph::Arc node = g.arcFromId(*ni);
    
        if (!visited[node]) { // truly skip ends
            continue;
        }
       
       rcount value = 0;
       if (capacity[node] > bias[node]) {
       	   value = (capacity[node] - bias[node]) / contribution_factor[node];
       }
       median.push_back(value);
       #ifdef ALLOW_DEBUG
       logger::Instance()->debug("Add to Median List FWD: " + std::to_string(g.id(node)) + " " + std::to_string(value) + " of " + std::to_string(capacity[node]) + " at factor " + std::to_string(contribution_factor[node]) + " with bias " + std::to_string(bias[node])+ "\n");
       #endif
    }
    
    if (median.empty()) {
        return 0;
    }
    
    std::sort (median.begin(), median.end());
    return median[ median.size() / 2 ];
}

rcount base_manager::get_backward_median(std::deque<ListDigraph::Node>::reverse_iterator &start, std::deque<ListDigraph::Node>::reverse_iterator &end,
        ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &push_block, ListDigraph::ArcMap<bool> &marked_source) {
    
    std::set<int> contributing_nodes;
    ListDigraph::ArcMap<float> contribution_factor(g); // init as 0
    ListDigraph::ArcMap<rcount> bias(g); // init as 0
    ListDigraph::ArcMap<bool> visited(g);

    ListDigraph::OutArcIt sn1(g, *start);
    contribution_factor[sn1] = 1;    
    visited[sn1] = true;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Entry: " + std::to_string(g.id(sn1)) + "\n");
    #endif
    
    for (std::deque<ListDigraph::Node>::reverse_iterator n = start; n!= end; ++n) { 
        // first test if we can propagate
        ListDigraph::OutArcIt onnt(g, *n);
        bool node_right= (onnt != INVALID && edge_type[onnt] == edge_types::NODE);
        if (node_right && visited[onnt]) {
              for (ListDigraph::InArcIt oi(g, *n); oi != INVALID ;++oi) {
                  ListDigraph::InArcIt node(g, g.source(oi));
                  if ( node != INVALID && capacity[node] != 0 && !marked_source[node] && (guided_saves[node] || !push_block[node]) ) {
                      visited[node] = true;
                  }
              }
	}
    }
    
    for (std::deque<ListDigraph::Node>::reverse_iterator n = start; n!= end; ++n) {

	ListDigraph::OutArcIt onnt(g, *n);
        bool node_right= (onnt != INVALID && edge_type[onnt] == edge_types::NODE);
        if (node_right && visited[onnt]) {
        
            if (onnt != sn1) {
                
                rcount b = 0;		
                for (ListDigraph::OutArcIt out_edge(g, g.target(onnt)); out_edge != INVALID ;++out_edge) {
                    ListDigraph::OutArcIt out_node(g, g.target(out_edge));
                    if (out_node != INVALID && !visited[out_node]) { //&& (guided_saves[out_node] || !push_block[out_node]) ) {
                            b += capacity[out_edge];
                            
                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Bias From: " + std::to_string(g.id(out_edge)) + "\n");
                            #endif
                    }
                }

                bias[onnt] += b; 
                bias[onnt] = std::min(bias[onnt], capacity[onnt]);
            }
            
            
            std::vector<ListDigraph::Arc> adjourning;
            rcount total_edge_count = 0;
            for (ListDigraph::InArcIt oi(g, *n); oi != INVALID ;++oi) {
                // all edges have by definition a node handy
                ListDigraph::InArcIt node(g, g.source(oi));
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Node1: " + std::to_string(g.id(onnt))  + " Edge1: " + std::to_string(g.id(oi)) + " Node2: " + std::to_string(g.id(node)) + "\n");
                #endif
                
                if (node != INVALID &&  capacity[node] != 0 && visited[node]) { // && (guided_saves[node] || !push_block[node]) ) {
		
                    total_edge_count+= capacity[oi];
                    adjourning.push_back(oi);
                }
            }
            
            if (total_edge_count == 0) {
                continue;
            }
            
            for (std::vector<ListDigraph::Arc>::iterator edge = adjourning.begin(); edge != adjourning.end(); ++edge) {
             
                ListDigraph::InArcIt target(g, g.source(*edge));
                
                contribution_factor[target] += contribution_factor[onnt] * (capacity[*edge])/total_edge_count;
		bias[target] += bias[onnt] * (capacity[*edge])/total_edge_count;
                
                contributing_nodes.insert(g.id(target));
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("fwd: " + std::to_string(g.id(target)) + " " + std::to_string(contribution_factor[target]) + + " " + std::to_string(bias[target]) + "\n"); 
                #endif
            }       
        }
    }
    
    if (contributing_nodes.empty()) {
        return 0;
    }
    
    std::vector<rcount> median;
    for (std::set<int>::iterator ni = contributing_nodes.begin(); ni != contributing_nodes.end(); ++ni) {
	
       ListDigraph::Arc node = g.arcFromId(*ni);
    
       if (!visited[node]) { // truly skip ends
            continue;
        }
       
       rcount value = 0;
       if (capacity[node] > bias[node]) {
       	   value = (capacity[node] - bias[node]) / contribution_factor[node];
       }
       median.push_back(value);
       
       #ifdef ALLOW_DEBUG
       logger::Instance()->debug("Add to Median List BWD: " + std::to_string(g.id(node)) + " " + std::to_string(value) + " of " + std::to_string(capacity[node]) + " at factor " + std::to_string(contribution_factor[node]) + " with bias " + std::to_string(bias[node])+ "\n");
       #endif

    }
    
    if (median.empty()) {
        return 0;
    }
    
    std::sort (median.begin(), median.end());
    return median[ median.size() / 2 ];
}


bool base_manager::test_consistent_sloping_source(ListDigraph::Node p, ListDigraph::ArcMap<rcount> &maximum, ListDigraph::Arc last, rpos length) {
    
   #ifdef ALLOW_DEBUG 
   logger::Instance()->debug("Follow Sloping Source " + std::to_string(g.id(p)) + " " + std::to_string(g.id(last)) + " Len " + std::to_string(length) +"\n"); 
   #endif 
   
   // test incoming
   ListDigraph::OutArcIt onnt(g, p);
   bool node_right= (onnt != INVALID && edge_type[onnt] == edge_types::NODE);
   if (node_right) {

	bool fail = false;
        capacity_type incoming = 0;
   	for (ListDigraph::InArcIt i(g, p); i!=INVALID; ++i) {
           
            incoming += maximum[i];
            
   	    if (last == i || edge_type[i] == edge_types::HELPER) {
      	   	continue;
   	    }	
	
	    if ( regions[i].is_increasing_region(length) 
                && maximum[i] == capacity[i]
                && maximum[last] <  maximum[i]) {
		fail = true;
		break;
	    }
        } 
	if (fail || incoming > (regions[onnt].get_left() + capacity[onnt]) / 2) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Failed Incoming\n"); 
            #endif
	    return false;
   	}
    }

    // test forward
    bool found = false;
    for ( ListDigraph::OutArcIt a(g, p); a!=INVALID; ++a) {
        
//        logger::Instance()->debug("Test Inner Source " + std::to_string(g.id(a)) +"\n");
        
        // ignoring circular references and drains, as well as normally corrected edges
        if (edge_type[a] == edge_types::BACKLINK 
            || (edge_type[a] == edge_types::HELPER && g.target(a) == t) 
            || capacity[a] == 0 
            || !regions[a].is_increasing_region(length)  
                ) { 
            continue;
        }
        
        if (length > regions[a].total_length) {
            found = found || test_consistent_sloping_source(g.target(a), maximum, a, length -regions[a].total_length);
        } else {
            found = true;
            break;
        }
    }
    
    return found;
}


bool base_manager::test_consistent_sloping_drain(ListDigraph::Node p, ListDigraph::ArcMap<rcount> &maximum, ListDigraph::Arc last, rpos length) {
    
    #ifdef ALLOW_DEBUG
   logger::Instance()->debug("Follow Sloping Drain " + std::to_string(g.id(p)) + " " + std::to_string(g.id(last)) +"\n");  
    #endif
    
   // test incoming
   ListDigraph::InArcIt innt(g, p);
   
   bool node_left= (innt != INVALID && edge_type[innt] == edge_types::NODE);
   if (node_left) {

	bool fail = false;
        capacity_type outgoing = 0;
   	for (ListDigraph::OutArcIt i(g, p); i!=INVALID; ++i) {
            
           outgoing += maximum[i];  
            
   	   if (last == i || edge_type[i] == edge_types::HELPER) {
      	   	continue;
   	   }	
	
	   if ( regions[i].is_decreasing_region(length) 
                && maximum[i] == capacity[i]
                && maximum[last] <  maximum[i]) {
		fail = true;
		break;
	   }
         }
//        logger::Instance()->debug("Outgoing " + std::to_string(fail) + " "+ std::to_string(outgoing) + " "+ std::to_string(regions[innt].get_right()) + " "+ std::to_string(capacity[innt])  + "\n");
	 if (fail || outgoing > (regions[innt].get_right() + capacity[innt]) / 2 ) {
             #ifdef ALLOW_DEBUG
             logger::Instance()->debug("Failed Outgoing\n");
             #endif
	     return false;
   	 }
    }

    // test backward
    bool found = false;
    for ( ListDigraph::InArcIt a(g, p); a!=INVALID; ++a) {
        
//        logger::Instance()->debug("Test Inner Drain " + std::to_string(g.id(a)) + " " + std::to_string(maximum[a]) + " " + std::to_string(capacity[a]) + "\n");

        // ignoring circular references and drains, as well as normally corrected edges
        if (edge_type[a] == edge_types::BACKLINK 
            || (edge_type[a] == edge_types::HELPER && g.source(a) == s) 
            || capacity[a] == 0 
            || !regions[a].is_decreasing_region(length)      
                ) { 
//            logger::Instance()->debug("Continue\n");
            continue;
        }
        
//        logger::Instance()->debug("Next\n");
        
        if (length > regions[a].total_length) {
            found = found ||  test_consistent_sloping_drain(g.source(a), maximum, a, length -regions[a].total_length);
        } else {
            found = true;
            break;
        }
    }
    
//    logger::Instance()->debug("Return " + std::to_string(found) +"\n");
    return found;
}


void base_manager::prune_unevidenced_arcs(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Normal Arc Rest.\n");
    #endif
    
    ListDigraph::ArcMap<bool> kill_status(g); // only valid for EXON edges
    ListDigraph::ArcMap<bool> kill_fix(g); 
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        if (edge_type[a] == edge_types::EXON) {
            if (regions[a].get_average() < options::Instance()->get_minimal_edge_capacity()) {
                kill_status[a] = true;
                kill_fix[a] = true;
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Kill Exon by filter " + std::to_string(g.id(a)) + ".\n");
                #endif
                continue;
            }
        }
        
        if (edge_type[a] != edge_types::NODE) {
            continue;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Prune NodeArc " + std::to_string(g.id(a)) + ".\n");
        #endif
        
        // we go over every node an check for possible removal of arcs
        
        rcount left = regions[a].get_left();
        rcount right = regions[a].get_right();
        rcount total_max = std::max(left, right);
        
        rcount in_count_norm = 0;
        rcount in_count_source = 0;
        rcount max_value_in_norm = 0;
        rcount max_value_in_source = 0;
        
        rcount out_count_norm = 0;
        rcount out_count_drain = 0;
        rcount max_value_out_norm = 0;
        rcount max_value_out_drain = 0;
        for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
            if (edge_type[i] == edge_types::HELPER) {
                // helper can't possibly be disproven/proven
                continue;
            }
            
            if (marked_source[i]) {
                if (capacity[i] > max_value_in_source) {
                    max_value_in_source = capacity[i];
                }
                in_count_source += capacity[i];
            } else {
                if (capacity[i] > max_value_in_norm) {
                    max_value_in_norm = capacity[i];
                }
                in_count_norm += capacity[i];
            }
        }
        rcount max_in = std::max(max_value_in_source, max_value_in_norm);
        
        for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
            if (edge_type[o] == edge_types::HELPER) {
                // helper can't possibly be disproven/proven
                continue;
            }
            
            if (marked_drain[o]) {
                if (capacity[o] > max_value_out_drain) {
                    max_value_out_drain = capacity[o];
                }
                out_count_drain += capacity[o];
            } else {
                if (capacity[o] > max_value_out_norm) {
                    max_value_out_norm = capacity[o];
                }
                out_count_norm += capacity[o];
            }
        }
        rcount max_out = std::max(max_value_out_drain, max_value_out_norm);
        
        for (int i=2; i < 4; ++i) {
            
            rcount reference = 1;
            rcount max_single = 1;
            switch(i) {
                case 0:
                    reference = total_max;
                    max_single = total_max;
                    break;
                case 1:
                    reference = std::max(out_count_drain+out_count_norm, in_count_source+in_count_norm);
                    max_single = std::max(max_in, max_out);
                    break;
                case 2:
                    reference = in_count_source+in_count_norm;
                    max_single = max_in;
                    break;
                case 3:
                    reference = in_count_norm;
                    max_single = max_value_in_norm; 
                    break;
            }
            
            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                
                if (edge_type[i] == edge_types::HELPER) {
                    // helper can't possibly be disproven/proven
                    continue;
                }      
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("CapTest L " + std::to_string(reference) + ": " + std::to_string(capacity[i]) + " " + std::to_string(capacity[i] * 100 / float(reference)) + " " + std::to_string(max_single)  + ".\n");
                #endif
                // || (max_single >= 2 && capacity[i] < 2 ) || (max_single > 5 && capacity[i] < 3 )  // || (max_single >= 20 && capacity[i] < 4 //  )
                if (capacity[i] * 100 / float(reference) < options::Instance()->get_arc_filter() || (max_single >= 4 && capacity[i] < 2 ) || (max_single >= 15 && capacity[i] < 3 ) || (max_single >= 25 && capacity[i] < 4) ) {
                    kill_status[i] = true;
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Kill in " + std::to_string(g.id(i)) + ".\n");
                    #endif
                }
            }
            
            break;
            
            bool source_open = false;
            bool norm_open = false;
            bool norm_arc_adjoining = false;
            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                if (edge_type[i] == edge_types::HELPER) {
                    continue;
                }

                if (marked_source[i] && !kill_status[i]) {
                    source_open = true;
                } else if (!marked_source[i]) {
                    if (!kill_status[i]) {
                        norm_open = true;
                    }
                    norm_arc_adjoining = true;
                } 
            }
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Status L " + std::to_string(source_open) + " " + std::to_string(norm_open) + " " + std::to_string(norm_arc_adjoining)  + ".\n");
            #endif
            
            if ( norm_arc_adjoining && norm_open || !norm_arc_adjoining && source_open) {
                for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                    if (edge_type[i] == edge_types::HELPER ) {
                        continue;
                    }

                    if (kill_status[i]) kill_fix[i] = true;
                }  
                break;
            }

            // reset
            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                if (edge_type[i] == edge_types::HELPER ) {
                    continue;
                }

                if (!kill_fix[i]) kill_status[i] = false;
            }  

        }
        
        for (int i=2; i < 4; ++i) {
            
            rcount reference = 1;
            rcount max_single = 1;
            switch(i) {
                case 0:
                    reference = total_max;
                    max_single = total_max;
                    break;
                case 1:
                    reference = std::max(out_count_drain+out_count_norm, in_count_source+in_count_norm);
                    max_single = std::max(max_in, max_out);
                    break;
                case 2:
                    reference = out_count_drain+out_count_norm;
                    max_single = max_out;
                    break;
                case 3:
                    reference = out_count_norm;
                    max_single = max_value_out_norm;      
                    break;
            }
            
            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                
                if (edge_type[o] == edge_types::HELPER) {
                    // helper can't possibly be disproven/proven
                    continue;
                }      
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("CapTest R " + std::to_string(reference) + ": " + std::to_string(capacity[o]) + " " + std::to_string(capacity[o] * 100 / float(reference)) + " " + std::to_string(max_single)  + ".\n");
                #endif
                // || (max_single >= 2 && capacity[o] < 2 ) || (max_single > 5 && capacity[o] < 3 ) // || (max_single >= 10 && capacity[o] < 3 ) // 
                if (capacity[o] * 100 / float(reference) < options::Instance()->get_arc_filter() || (max_single >= 4 && capacity[o] < 2 ) || (max_single >= 15 && capacity[o] < 3 ) || (max_single >= 25 && capacity[o] < 4 ) ) {
                    kill_status[o] = true;
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Kill out " + std::to_string(g.id(o)) + ".\n");
                    #endif
                }
            }

            break;
            
            bool drain_open = false;
            bool norm_open = false;
            bool norm_arc_adjoining = false;
            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                if (edge_type[o] == edge_types::HELPER) {
                    continue;
                }

                if (marked_drain[o] && !kill_status[o]) {
                    drain_open = true;
                } else if (!marked_drain[o]) {
                    if (!kill_status[o]) {
                        norm_open = true;
                    }
                    norm_arc_adjoining = true;
                } 
            }
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Status R " + std::to_string(drain_open) + " " + std::to_string(norm_open) + " " + std::to_string(norm_arc_adjoining)  + ".\n");
            #endif
            
            if ( norm_arc_adjoining && norm_open || !norm_arc_adjoining && drain_open) {
                for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                    if (edge_type[o] == edge_types::HELPER ) {
                        continue;
                    }

                    if (kill_status[o]) kill_fix[o] = true;
                }  
                break;
            }

            // reset
            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                if (edge_type[o] == edge_types::HELPER ) {
                    continue;
                }

                if (!kill_fix[o]) kill_status[o] = false;
            }  
        }  
    }
    
    // test if this separates, if yes we can't do it
    ListDigraph wc;
    ListDigraph::ArcMap<bool> killc(wc);
    ListDigraph::ArcMap<bool> guided_saves_c(wc);
    ListDigraph::Node cs, ct;
    DigraphCopy<ListDigraph, ListDigraph> copy(g,wc);
    copy.arcMap(kill_status, killc).arcMap(guided_saves, guided_saves_c)
     .node(s,cs).node(t,ct).run();
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ) {
        if (killc[a] == true && guided_saves_c[a] == false) {
            ListDigraph::Arc c(a);
            ++a;
            wc.erase(c);
        } else {
            ++a;
        }
    }
    
    Dfs<ListDigraph> dfs_b(wc);
    if (!dfs_b.run(cs,ct)) {
        return;
    }
     
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("ALL Marked\n");
    #endif
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ) {
        if (kill_status[a] == true && guided_saves[a] == false) {
            ListDigraph::Arc c(a);
            ++a;
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Erase " + std::to_string(g.id(c)) +"\n");
            #endif
            g.erase(c);
        } else {
            ++a;
        }
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Initial Erased Done \n");
    #endif
    
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
        logger::Instance()->debug("Test " + std::to_string(g.id(no)) + " " + std::to_string(regions[no].get_length_to_first_zero_from_right()) + ".\n");
        #endif
        
        if (regions[no].get_length_to_first_zero_from_right() <= 20 || regions[no].get_average() < 2) {
            //erase
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
        logger::Instance()->debug("Test " + std::to_string(g.id(ni)) + " " + std::to_string(regions[ni].get_length_to_first_zero_from_left()) + ".\n");
        #endif
        
        if (regions[ni].get_length_to_first_zero_from_left() <= 20 || regions[ni].get_average() < 2) {
            //erase
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
        
        if (edge_type[arc] == edge_types::NODE) {
           rpos length = regions[arc].total_length;
           float average = regions[arc].get_average();
           
            bool left_neighbour = false;
            bool left_source = false;
            bool left_other = false;
            for(ListDigraph::InArcIt ei(g, g.source(arc)); ei != INVALID; ++ei) {
               if (edge_specifier[ei].right_consecutive) {
                   left_neighbour = true;
               } else if (edge_type[ei] == edge_types::HELPER) {
                   left_source = true;
               }
            }
           
            bool right_neighbour = false;
            bool right_drain = false;
            for(ListDigraph::OutArcIt eo(g, g.target(arc)); eo != INVALID; ++eo) {
                if (edge_specifier[eo].left_consecutive ) {
                    right_neighbour = true;
                } else if (edge_type[eo] == edge_types::HELPER) {
                   right_drain = true;
                }
            }
           
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test Arc " +  std::to_string(g.id(arc))  + " " + std::to_string(left_source) + " " + std::to_string(right_drain) + " " + std::to_string(length) + " " + std::to_string(average) + "\n");
            #endif

            if ( (left_source || right_drain) && (!left_neighbour || !right_neighbour) && (length <= 15 || average <= 1.5) && !guided_saves[arc]) {
               
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
        if (edge_type[a] == edge_types::NODE) { //with opened middle node
            
            int ni = edge_specifier[a].node_index;
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
        
        if (edge_type[arc] == edge_types::NODE ) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Check Dev"  + std::to_string(g.id(arc)) + " " + std::to_string(regions[arc].get_deviation())  + "\n");
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
            
            if (edge_type[arc_left] != edge_types::HELPER && edge_type[arc_right] != edge_types::HELPER) {
                continue;
            }
            
            if (regions[arc].get_deviation() > 0.01) {
                continue;
            }
            
            if (edge_type[arc_left] == edge_types::HELPER) {
                // we test if this is unique!
                unsigned int counter = 0;
                for (ListDigraph::InArcIt oi(g, g.target(arc_right)); oi != INVALID ;++oi) {
                    ++counter;
                }
                if (counter == 1) {
                    continue;
                }
            }
            
            if (edge_type[arc_right] == edge_types::HELPER) {
                // we test if this is unique!
                unsigned int counter = 0;
                for (ListDigraph::OutArcIt oi(g, g.source(arc_left)); oi != INVALID ;++oi) {
                    ++counter;
                }
                if (counter == 1) {
                    continue;
                }
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Erase by Dev " + std::to_string(g.id(arc)) + ".\n");
            #endif
            
            res = true;
            
            g.erase(arc);
            g.erase(arc_left);
            g.erase(arc_right);
        }
    }
    return res;
}
    

bool base_manager::erase_low_boundary_st(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain, ListDigraph::ArcMap<bool> &consecutive_s_t) {

    std::vector<std::deque<ListDigraph::Node> > same_node_r(meta->size);
    std::vector<std::deque<ListDigraph::Node> > same_node_l(meta->size);
    std::vector<std::deque<ListDigraph::Arc> > same_node(meta->size);
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (edge_type[a] == edge_types::NODE) { //with opened middle node
            
            int ni = edge_specifier[a].node_index;
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
        
        if (edge_type[arc] == edge_types::EXON && g.source(arc) != s && g.target(arc) != t) {
            
            // we go over every node an check for possible removal of arcs

            ListDigraph::InArcIt node_left(g, g.source(arc));
            ListDigraph::OutArcIt node_right(g, g.target(arc));
            
            if (node_left == INVALID || node_right == INVALID) {
                    continue;
            }
            
            unsigned int out_degree = 0;
            bool drain_exists = false;
            for (std::deque<ListDigraph::Node>::iterator nr = same_node_r[edge_specifier[node_left].node_index].begin(); nr != same_node_r[edge_specifier[node_left].node_index].end(); ++nr) {
                for (ListDigraph::OutArcIt oi(g, *nr); oi != INVALID ;++oi) {
                    if (edge_type[oi] == edge_types::HELPER) {
                        drain_exists = true;
                        continue;
                    }
                    ++out_degree;
                }
            }
            
            unsigned int in_degree = 0;
            bool source_exists = false;
            for (std::deque<ListDigraph::Node>::iterator nl = same_node_l[edge_specifier[node_right].node_index].begin(); nl != same_node_l[edge_specifier[node_right].node_index].end(); ++nl) {
                for (ListDigraph::InArcIt oi(g, *nl); oi != INVALID ;++oi) {
                    if (edge_type[oi] == edge_types::HELPER) {
                        source_exists = true;
                        continue;
                    }
                    ++in_degree;
                }
            }


            if (meta->exons[edge_specifier[node_left].node_index].right + 1 == meta->exons[edge_specifier[node_right].node_index].left) {
                continue; // disregard neighbouring here!
            }
            
            rcount we = capacity[arc];
            rcount left = 0;
            for (std::deque<ListDigraph::Arc>::iterator ni = same_node[edge_specifier[node_left].node_index].begin(); ni != same_node[edge_specifier[node_left].node_index].end(); ++ni) {
                left += regions[*ni].get_average();
            }
            
            rcount right = 0;
            for (std::deque<ListDigraph::Arc>::iterator ni = same_node[edge_specifier[node_right].node_index].begin(); ni != same_node[edge_specifier[node_right].node_index].end(); ++ni) {
                right += regions[*ni].get_average();
            }

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test N Scallop " + std::to_string(g.id(arc)) + " : " + std::to_string(we) + " " + std::to_string(left) + " " + std::to_string(right) + " : " + std::to_string(meta->exons[edge_specifier[node_left].node_index].right) + " " + std::to_string(meta->exons[edge_specifier[node_right].node_index].left)  + ".\n");
            #endif

            if ((out_degree == 1 && left>= 10.0 * we * we + 10.0) || (in_degree == 1 && right>= 10.0 * we * we + 10.0)) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Erase by boundary " + std::to_string(g.id(arc)) + ".\n");
                #endif

                res = true;
                
                g.erase(arc);
                if (!drain_exists && out_degree == 1) {
                    ListDigraph::Arc narc = g.addArc(g.target(node_left), t);
                    edge_type[narc] = edge_types::HELPER;
                    initialize_source_drain_arc(narc);
                }
               if (!source_exists && in_degree == 1)  {
                    ListDigraph::Arc narc = g.addArc(s, g.source(node_right));
                    edge_type[narc] = edge_types::HELPER;
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
        if (edge_type[a] == edge_types::NODE) { //with opened middle node
            int ni = edge_specifier[a].node_index;
            same_node[ni].push_back(a);  
        }
    }
    
    bool res = false;
    for (ListDigraph::ArcIt a(g); a != INVALID;) {
        
        ListDigraph::Arc arc(a);
        ++a;
        
        capacity_type current = regions[arc].get_average();
        
        if (edge_type[arc] == edge_types::NODE) {
            
            capacity_type max_in = 0;
            bool erased = false;
            for (std::deque<ListDigraph::Arc>::iterator ni = same_node[edge_specifier[arc].node_index].begin(); ni != same_node[edge_specifier[arc].node_index].end(); ++ni) {
                for (ListDigraph::InArcIt oi(g, g.source(*ni)); oi != INVALID ;++oi) {

                    ListDigraph::InArcIt node_left(g, g.source(oi));
                    if (node_left == INVALID) continue;

                    if (!edge_specifier[oi].right_consecutive) continue;

                    capacity_type nla = regions[node_left].get_average();
                    capacity_type wl = std::max(nla , capacity[oi]);

                    if (wl > max_in) {
                        max_in = wl;
                    }
                }
            }
            for (std::deque<ListDigraph::Arc>::iterator ni = same_node[edge_specifier[arc].node_index].begin(); ni != same_node[edge_specifier[arc].node_index].end(); ++ni) {
                for (ListDigraph::InArcIt oi(g, g.source(*ni)); oi != INVALID ;) {

                    ListDigraph::Arc oik(oi);
                    ++oi;

                    if (edge_specifier[oik].right_consecutive) continue;

                    capacity_type w = capacity[oik];

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Test scallop L " + std::to_string(g.id(oik)) + " to " + std::to_string(g.id(arc)) + " : " + std::to_string(w) + " " + std::to_string(max_in) + " " + std::to_string(current) + ".\n");
                    #endif

                     if ( max_in >= 2.0 * w * w + 18.0 && current >= 2.0 * w * w + 18.0 && !guided_saves[oik]) {
                        g.erase(oik);
                        
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Kill scallop L " + std::to_string(g.id(arc)) + ".\n");
                        #endif
                        res = true;
                    }

                }
            }
            
            capacity_type max_out = 0;
            for (std::deque<ListDigraph::Arc>::iterator ni = same_node[edge_specifier[arc].node_index].begin(); ni != same_node[edge_specifier[arc].node_index].end(); ++ni) {
                for (ListDigraph::OutArcIt oi(g, g.target(*ni)); oi != INVALID ;++oi) {

                    ListDigraph::OutArcIt node_right(g, g.target(oi));
                    if (node_right == INVALID) continue;
                    if (!edge_specifier[oi].left_consecutive) continue;

                    capacity_type nra = regions[node_right].get_average();
                    capacity_type wr = std::max(nra , capacity[oi]);
                    if (wr > max_out) {
                        max_out = wr;
                    }
                }
            }
            
            for (std::deque<ListDigraph::Arc>::iterator ni = same_node[edge_specifier[arc].node_index].begin(); ni != same_node[edge_specifier[arc].node_index].end(); ++ni) {
                for (ListDigraph::OutArcIt oi(g, g.target(*ni)); oi != INVALID ;) {

                    ListDigraph::Arc oik(oi);
                    ++oi;

                    if (edge_specifier[oik].left_consecutive) continue;

                    capacity_type w = capacity[oik];

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Test scallop R " + std::to_string(g.id(oik)) + " to " + std::to_string(g.id(arc)) + " : " + std::to_string(w) + " " + std::to_string(max_out) + " " + std::to_string(current) + ".\n");
                    #endif
                    
                    if ( max_out >= 2.0 * w * w + 18.0 && current >= 2.0 * w * w + 18.0 && !guided_saves[oik]) {
                        g.erase(oik);
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Kill scallop R " + std::to_string(g.id(arc)) + ".\n");
                        #endif
                        res = true;
                    }

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
        if (edge_type[a] == edge_types::NODE) { //with opened middle node
            int ni = edge_specifier[a].node_index;
            same_node[ni].push_back(a);  
        }
    }
    
    ListDigraph::ArcMap<bool> filtered_sd(g);
    ListDigraph::ArcMap<bool> kill_status(g); // potentially killed edges!
//    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//        if (edge_type[a] == edge_types::NODE) { 
//            
//            
//            rcount max_in = 0;
//            rcount total_in = 0;
//            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
//                
//                if (edge_type[i] == edge_types::HELPER) {
//                    // helper can't possibly be disproven/proven
//                    continue;
//                }
//
//                if (capacity[i] > max_in) {
//                    max_in = capacity[i];
//                }
//                total_in += capacity[i];
//            }
//            
//            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
//                
//                if (edge_type[i] == edge_types::HELPER || guided_saves[i]) {
//                    // helper can't possibly be disproven/proven
//                    continue;
//                }      
//                
//                //if ( 1.75 * capacity[i] * capacity[i] < max_in || capacity[i] * capacity[i] * 100 / float(total_in) < options::Instance()->get_arc_filter()) {
//                if (marked_source[i] && edge_specifier[i].right_consecutive && capacity[i] * 100 / float(max_in) < 5 && capacity[i] < 20) {
//                    kill_status[i] = true;
//                    logger::Instance()->debug("Mark Left " + std::to_string(g.id(i)) + ".\n");
//
//                    ListDigraph::InArcIt is(g, g.source(i)); // get the last node
//                    if (is != INVALID)  kill_status[is] = true;
//
//                }
//            }
//            
//            rcount max_out = 0;
//            rcount total_out = 0;
//            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
//                
//                if (edge_type[o] == edge_types::HELPER) {
//                    // helper can't possibly be disproven/proven
//                    continue;
//                }
//
//                if (capacity[o] > max_out) {
//                    max_out = capacity[o];
//                }
//                total_out += capacity[o];
//            }
//            
//            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
//                
//                if (edge_type[o] == edge_types::HELPER || guided_saves[o] ) {
//                    // helper can't possibly be disproven/proven
//                    continue;
//                }      
//                //if ( 1.75 * capacity[o] * capacity[o] < max_out || capacity[o] * capacity[o] * 100 / float(total_out) < options::Instance()->get_arc_filter()  ) {
//                if (marked_drain[o] && edge_specifier[o].left_consecutive && capacity[o] * 100 / float(max_out) < 5 && capacity[o] < 20) {
//                    kill_status[o] = true;
//                    logger::Instance()->debug("Mark Right " + std::to_string(g.id(o)) + ".\n");
//                    
//                    ListDigraph::OutArcIt ot(g, g.target(o)); // get the last node
//                    if (ot != INVALID)  kill_status[ot] = true;
//
//                }
//            }
//        } 
//    }
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
	logger::Instance()->debug("Test low Source " + std::to_string(g.id(no)) + " " + std::to_string(regions[no].get_average_to_first_zero_from_right()) + ".\n");
        #endif
             
        if (regions[no].get_average_to_first_zero_from_right() <= low_mark || regions[no].total_length <= 20 ) {
            //erase
            #ifdef ALLOW_DEBUG
	    logger::Instance()->debug("Mark Short Source " + std::to_string(g.id(arc)) + " " + std::to_string(regions[no].get_average()) + ".\n");
            #endif
            kill_status[no] = true;
            filtered_sd[no] = true;
//            kill_status[no2] = true;
            for (ListDigraph::OutArcIt no3(g, g.target(no)); no3 != INVALID; ++no3) {
                kill_status[no3] = true;
                filtered_sd[no3] = true;
            }
            
//            std::deque<ListDigraph::Arc> open;
//            for (ListDigraph::OutArcIt no3(g, g.target(no)); no3 != INVALID; ++no3) {
//                open.push_back(no3);
//            }
//            while (true) {
//                if (open.empty()) {
//                    break;
//                }
//                std::deque<ListDigraph::Arc> next;
//                for(std::deque<ListDigraph::Arc>::iterator it = open.begin(); it != open.end(); ++it) {
//                    kill_status[*it] = true;
//                    ListDigraph::InArcIt ni(g, g.source(*it));
//                    if (ni == INVALID) continue;
//                    ++ni;
//                    if (ni != INVALID) continue;
//                    
//                    for (ListDigraph::OutArcIt o(g,g.target(*it)); o != INVALID; ++o) {
//                        next.push_back(o);
//                    }
//                }
//                open = next;
//            }
            
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
	logger::Instance()->debug("Test low Drain " + std::to_string(g.id(ni)) + " " + std::to_string(regions[ni].get_average_to_first_zero_from_left()) + ".\n");
        #endif
        
        if (regions[ni].get_average_to_first_zero_from_left() <= low_mark  || regions[ni].total_length <= 20 ) {
            //erase
            #ifdef ALLOW_DEBUG
	    logger::Instance()->debug("Mark Short Drain " + std::to_string(g.id(arc)) + " " + std::to_string(regions[ni].get_average()) + ".\n");
            #endif
            kill_status[ni] = true;
            filtered_sd[ni] = true;
//            kill_status[ni2] = true;
            for (ListDigraph::InArcIt ni3(g, g.source(ni)); ni3 != INVALID; ++ni3) {
                kill_status[ni3] = true;
                filtered_sd[ni] = true;
            }
            
//            std::deque<ListDigraph::Arc> open;
//            for (ListDigraph::InArcIt ni3(g, g.source(ni)); ni3 != INVALID; ++ni3) {
//                open.push_back(ni3);
//            }
//            while (true) {
//                if (open.empty()) {
//                    break;
//                }
//                std::deque<ListDigraph::Arc> next;
//                for(std::deque<ListDigraph::Arc>::iterator it = open.begin(); it != open.end(); ++it) {
//                    kill_status[*it] = true;
//                    ListDigraph::OutArcIt no(g, g.target(*it));
//                    if (no == INVALID) continue;
//                    ++no;
//                    if (no != INVALID) continue;
//                    
//                    for (ListDigraph::InArcIt o(g,g.source(*it)); o != INVALID; ++o) {
//                        next.push_back(o);
//                    }
//                }
//                open = next;
//            }
            
        }
    }    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        if (edge_type[a] == edge_types::EXON && capacity[a] <= low_mark) {
            kill_status[a] = true;
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Mark Exon by threshold " + std::to_string(g.id(a)) + ".\n");
             #endif
        } // else if (edge_type[a] == edge_types::NODE && (regions[a].get_length_to_first_zero_from_left() <= 20 ||  regions[a].get_length_to_first_zero_from_right() <= 20)) {
//            kill_status[a] = true;
//            logger::Instance()->debug("Mark Possible bad intron " + std::to_string(g.id(a)) + ".\n");
//        } 
    }

    ListDigraph::Arc max_node;
    capacity_type max_cap = 0;
    bool node_leftover = false;
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        if (edge_type[a] == edge_types::NODE) { 

            if (capacity[a] > max_cap) {
                max_cap = capacity[a];
                max_node = a;
            }
            
            bool has_in_arc = false;
            bool has_in_helper = false;
            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                if (edge_type[i] == edge_types::HELPER) {
                    has_in_helper = true;
                } else if (kill_status[i] != true) {
                    has_in_arc = true;
                }
            }

            bool has_out_arc = false;
            bool has_out_helper = false;
            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                if (edge_type[o] == edge_types::HELPER) {
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
        kill_status[max_node] = false;
    }
    
//    digraphWriter(g, std::cout)
//                .arcMap("edge_specifier", edge_specifier)
//                .arcMap("edge_type", edge_type)
//                .arcMap("cap", capacity)
//                .arcMap("kill", kill_status)
//                .run();  
    
    // rescue mission!
    
    while (true) {
        bool changed = false;
        ListDigraph::ArcMap<bool> next_rescue(g); // mark all rescues here!
                
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
            
            if (kill_status[a] || edge_type[a] == edge_types::HELPER) {
                continue;
            }
            
            bool has_in_arc = false;
            bool has_in_helper = false;
            for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {
                if (edge_type[i] == edge_types::HELPER) {
                    has_in_helper = true;
                } else if (kill_status[i] != true) {
                    has_in_arc = true;
                }
            }

            bool has_out_arc = false;
            bool has_out_helper = false;
            for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                if (edge_type[o] == edge_types::HELPER) {
                    has_out_helper = true;
                } else if(kill_status[o] != true) {
                    has_out_arc = true;
                }
            }
            bool has_in = has_in_arc || has_in_helper;
            bool has_out = has_out_arc || has_out_helper; 
            
            bool has_non_sd_filtered_in = false;
            float max_count_in = 0;
            bool has_non_sd_filtered_out = false;
            float max_count_out = 0;
            if (edge_type[a] == edge_types::NODE) {
                for (std::deque<ListDigraph::Arc>::iterator ni = same_node[edge_specifier[a].node_index].begin(); ni != same_node[edge_specifier[a].node_index].end(); ++ni) {
                    for (ListDigraph::InArcIt i(g, g.source(*ni)); i!=INVALID; ++i) {
                        
                        if (!filtered_sd[i]) {
                            has_non_sd_filtered_in = true;
                        }
                        
                        float cap = means[i].mean;

                        if (max_count_in < cap) {
                            max_count_in = cap;
                        }
                    }
                    for (ListDigraph::OutArcIt o(g, g.target(*ni)); o!=INVALID; ++o) {
                        if (!filtered_sd[o]) {
                            has_non_sd_filtered_out = true;
                        }
                        
                        float cap = means[o].mean;
                        
                        if (max_count_out < cap) {
                            max_count_out = cap;
                        }
                    }
                }
            }
            
            if (!has_in) {               
                bool has_max = false;
                float max_count = 0;                
                std::deque<ListDigraph::Arc> max_arc_all;

                for (ListDigraph::InArcIt i(g, g.source(a)); i!=INVALID; ++i) {

                    float cap = means[i].mean;

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
                    
                    if (!filtered_sd[*it] || !has_non_sd_filtered_in || max_count_in == max_count) { // arc is not low drain or source or only filtered ones exist
                    
                        next_rescue[*it] = true;
                        changed = true;       
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Rescue In " + std::to_string(g.id(a)) + " " + std::to_string(g.id(*it)) +".\n");
                        #endif
                    }
                }
            }
            

            if (!has_out) {                
                bool has_max = false;
                float max_count = 0;                
                std::deque<ListDigraph::Arc> max_arc_all;
                
                for (ListDigraph::OutArcIt o(g, g.target(a)); o!=INVALID; ++o) {
                    
                    float cap = means[o].mean;

                    logger::Instance()->debug("Cap " + std::to_string(g.id(o)) + " " + std::to_string(cap) +".\n");
                    
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
                     if (!filtered_sd[*it] || !has_non_sd_filtered_out || max_count_out == max_count) { // arc is not low drain or source or only filtered ones exist
                        next_rescue[*it] = true;
                        changed = true;
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Rescue Out " + std::to_string(g.id(a)) + " " + std::to_string(g.id(*it)) +".\n");
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


void base_manager::denoise_graph(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &push_block, ListDigraph::ArcMap<bool> &marked_source, ListDigraph::ArcMap<bool> &marked_drain) {
    
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
        
        average_push[a]= means[a].mean;
         
        if (marked_source[a]) {
            average_push[a] = rec_score_backward(a);
        }
        if (marked_drain[a]) {
            average_push[a] = rec_score_forward(a);
        }
         
//        average_push_fwd[a]= means[a].compute_score();
//        average_push_bwd[a]= means[a].compute_score();
//        if (edge_type[a] == edge_types::EXON) {
//            average_push_fwd[a] = rec_score_forward(a);
//            average_push_bwd[a] = rec_score_backward(a);
//        }
//        
//        #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("fwd/bwd " + std::to_string(g.id(a)) + " " + std::to_string( average_push_fwd[a]) + " / " + std::to_string( average_push_bwd[a]) + "\n");
//        #endif
//        
//        if (edge_type[a] == edge_types::NODE) {
//            average_push_fwd[a]= regions[a].get_average();
//            average_push_bwd[a]= regions[a].get_average();
//            
//            logger::Instance()->debug("Origin " + std::to_string(g.id(a)) + "\n");
//            
//            for (ListDigraph::InArcIt i(g, g.source(a)); i != INVALID ;++i) {
//                capacity_mean mean = means[i];
//                ListDigraph::Arc c(i); 
//                while (true) {
//                    ListDigraph::InArcIt next(g, g.source(c));
//                    if (next == INVALID) break;
//                    
//                    logger::Instance()->debug("LN " + std::to_string(g.id(next)) + "\n");
//                    
//                    ListDigraph::InArcIt ti(g, g.source(next));
//                    if (ti == INVALID) break;
//                    ++ti;
//                    if (ti != INVALID) break;  // we need to stop if non or more than two arcs adjourn!
//                    
//                    ListDigraph::OutArcIt to(g, g.target(next));
//                    if (to == INVALID) break;
//                    ++to;
//                    if (to != INVALID) break;  // we need to stop if non or more than two arcs adjourn!
//                      
//                    if (edge_type[next] == edge_types::EXON) {
//                        mean.update(means[next]);
//                    }
//                    
//                    c = next;     
//                }
//                #ifdef ALLOW_DEBUG
//                    logger::Instance()->debug("fwd " + std::to_string(g.id(i)) + " " + std::to_string(mean.mean) + "\n");
//                #endif
//                average_push_bwd[i] = mean.mean;
//            }
//            
//            for (ListDigraph::OutArcIt i(g, g.target(a)); i != INVALID ;++i) {
//                capacity_mean mean = means[i];
//                ListDigraph::Arc c(i); 
//                while (true) {
//                    ListDigraph::OutArcIt next(g, g.target(c));
//                    if (next == INVALID) break;
//                                        
//                    ListDigraph::InArcIt ti(g, g.source(next));
//                    if (ti == INVALID) break;
//                    ++ti;
//                    if (ti != INVALID) break;  // we need to stop if non or more than two arcs adjourn!
//                    
//                    ListDigraph::OutArcIt to(g, g.target(next));
//                    if (to == INVALID) break;
//                    ++to;
//                    if (to != INVALID) break;  // we need to stop if non or more than two arcs adjourn!
//                    
//                    if (edge_type[next] == edge_types::EXON) {
//                        mean.update(means[next]);
//                    }
//                    
//                    c = next;     
//                }
//                #ifdef ALLOW_DEBUG
//                    logger::Instance()->debug("bwd " + std::to_string(g.id(i)) + " " + std::to_string(mean.mean) + "\n");
//                #endif
//                average_push_fwd[i] = mean.mean;
//            }
//            
//        }
    }
//    
//    for (std::deque<ListDigraph::Node>::iterator n = top_order.begin(); n!= top_order.end(); ++n) {
//        
//        ListDigraph::InArcIt in(g, *n);
//        ListDigraph::OutArcIt out(g, *n);
//        
//        if (in == INVALID || out == INVALID) {
//            continue;
//        }
//        
//        ListDigraph::InArcIt in_next(in);
//        ListDigraph::OutArcIt out_next(out);
//        ++in_next;
//        ++out_next;
//        
//        if (in_next == INVALID && out_next == INVALID && edge_type[in] != edge_types::HELPER ) { 
//            // push IN to OUT
//            rcount in_cout = average_push[in];
//            rcount out_cout = average_push[out];
//            if (in_cout < out_cout) {
//                average_push[in] = in_cout;
//                average_push[out] = in_cout;
//                
//                logger::Instance()->debug("PushAverage R " + std::to_string(g.id(in)) + " " + std::to_string(g.id(out)) + " : " + std::to_string(in_cout) + "\n");
//                
//            }
//        }
//    }
//    
//    for (std::deque<ListDigraph::Node>::reverse_iterator n = top_order.rbegin(); n!= top_order.rend(); ++n) {
//        ListDigraph::InArcIt in(g, *n);
//        ListDigraph::OutArcIt out(g, *n);
//        
//        if (in == INVALID || out == INVALID) {
//            continue;
//        }
//        
//        ListDigraph::InArcIt in_next(in);
//        ListDigraph::OutArcIt out_next(out);
//        ++in_next;
//        ++out_next;
//        
//        if (in_next == INVALID && out_next == INVALID && edge_type[out] != edge_types::HELPER ) { 
//            // push OUT to IN
//            rcount in_cout = average_push[in]; 
//            rcount out_cout = average_push[out];
//            if (in_cout > out_cout) {
//                average_push[in] = out_cout;
//                average_push[out] = out_cout;
//                
//                logger::Instance()->debug("PushAverage L " + std::to_string(g.id(in)) + " " + std::to_string(g.id(out)) + " : " + std::to_string(out_cout) + "\n");
//            }
//        }
//    }
    
    // forward potential!
    ListDigraph::ArcMap<rcount> cp_fwd(g); // init as 0
    push_potential_forward(cp_fwd, top_order, average_push, push_block);
//    std::set<ListDigraph::Node> additional_drain_candidates;
//    find_over_capacitated_end(t, cp_fwd, additional_drain_candidates);
    
    // backward potential!
    ListDigraph::ArcMap<rcount> cp_bw(g); // init as 0
    push_potential_backward(cp_bw, top_order, average_push, push_block);
//    std::set<ListDigraph::Node> additional_source_candidates;
//    find_over_capacitated_start(s, cp_bw, additional_source_candidates);
    
    // we add in an abort if increases are too high!
    // we do this ad hoc!
//    rcount total_max = 0;
//    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//        total_max = std::max(regions[a].get_max(), total_max);
//    }
//    total_max = total_max * 1.5;
//    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//        if ( std::max(cp_fwd[a], cp_bw[a]) >  total_max) {
//            logger::Instance()->debug("Abort Denoise!");
//            return;
//        }
//    }
    
    
    
    
//    rpos testlength = raw->average_fragment_length/2;
//    if (testlength < 100) {
//        testlength = 100;
//    }
//    
//    bool source_change = false;
//    for (std::set<ListDigraph::Node>::iterator n = additional_source_candidates.begin(); n != additional_source_candidates.end(); ++n) {  
//        
//        ListDigraph::OutArcIt onnt(g, *n);
//        bool node_right= (onnt != INVALID && edge_type[onnt] == edge_types::NODE);
//        if(!node_right) {
//            continue;
//        }
//        
//        rcount incoming = 0;
//        for ( ListDigraph::InArcIt a(g, *n); a!=INVALID; ++a) {
//            incoming += capacity[a] + cp_fwd[a];  
//        }
//        
//        
//        if ( cp_fwd[onnt] > 10 
//                || incoming < capacity[onnt] || !(incoming - capacity[onnt]) * 100 / float(capacity[onnt]) >= options::Instance()->get_coverage_change_limit() 
//                || !test_consistent_sloping_source(*n, testlength)
//                || test_no_pre_sloping_source(*n, testlength, cp_fwd)) {
//            continue;
//        }
//        
//        
//        logger::Instance()->debug("Add Source " + std::to_string(g.id(*n)) + "\n");
//
//        ListDigraph::Arc na = g.addArc(s, *n);
//        capacity[na] = 0;
//        edge_type[na] = edge_types::HELPER;
//        source_change = true;
//    }
//    
//    
//    bool drain_change = false;
//    for (std::set<ListDigraph::Node>::iterator n = additional_drain_candidates.begin(); n != additional_drain_candidates.end(); ++n) {  
//        
//        ListDigraph::InArcIt innt(g, *n);
//        bool node_left= (innt != INVALID && edge_type[innt] == edge_types::NODE);
//        if(!node_left) {
//            continue;
//        }
//        
//        rcount outgoing = 0;
//        for ( ListDigraph::InArcIt a(g, *n); a!=INVALID; ++a) {
//            outgoing += capacity[a] + cp_bw[a];  
//        }
//        
//        
//        if ( cp_bw[innt] > 10 
//                || outgoing < capacity[innt] || !(outgoing - capacity[innt]) * 100 / float(capacity[innt]) >= options::Instance()->get_coverage_change_limit() 
//                || !test_consistent_sloping_drain(*n, testlength)
//                || test_no_pre_sloping_drain(*n, testlength, cp_bw)) {
//            continue;
//        }
//        
//        
//        logger::Instance()->debug("Add Drain " + std::to_string(g.id(*n)) + "\n");
//
//        ListDigraph::Arc na = g.addArc(*n, t);
//        capacity[na] = 0;
//        edge_type[na] = edge_types::HELPER;
//        drain_change = true;
//    }
//       
//    if (drain_change) push_potential_forward(cp_fwd, top_order);
//    if (source_change) push_potential_backward(cp_bw, top_order);
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("add " + std::to_string(g.id(a)) + " " + std::to_string(cp_fwd[a]) + " " + std::to_string(cp_bw[a]) + "\n");
        #endif
        
        capacity_type max = std::max(cp_fwd[a], cp_bw[a]);
        if (edge_type[a] != edge_types::HELPER) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Correct Capacity " + std::to_string(g.id(a)) + " from " + std::to_string(capacity[a]) + " to " + std::to_string(capacity[a]+max) + ".\n");
            #endif
            capacity[a] += max;
        }
    }
}



bool base_manager::test_no_pre_sloping_source(ListDigraph::Node p, rpos length, ListDigraph::ArcMap<rcount> &cp_fwd) {
    
    for ( ListDigraph::InArcIt a(g, p); a!=INVALID; ++a) {
        
        if (cp_fwd[a] < 10 && regions[a].is_increasing_region(length)) {
            
            if (length < regions[a].total_length) {
                return true;
            }
            
            // test also adjourning
            for ( ListDigraph::InArcIt b(g, g.source(a)); b!=INVALID; ++b) {
                if (cp_fwd[b] < 10 && regions[b].is_increasing_region(length - regions[a].total_length)) {
                    return true;
                }
            }
        } 
    }
    
    return false;
}


bool base_manager::test_consistent_sloping_source(ListDigraph::Node p, rpos length) {
    
    // test forward
    bool found = false;
    for ( ListDigraph::OutArcIt a(g, p); a!=INVALID; ++a) {
        
//        logger::Instance()->debug("Test Inner Source " + std::to_string(g.id(a)) +"\n");
        
        // ignoring circular references and drains, as well as normally corrected edges
        if (edge_type[a] == edge_types::BACKLINK 
            || (edge_type[a] == edge_types::HELPER && g.target(a) == t) 
            || capacity[a] == 0 
            || !regions[a].is_increasing_region(length)  
                ) { 
            continue;
        }
        
        if (length > regions[a].total_length) {
            found = found || test_consistent_sloping_source(g.target(a), length -regions[a].total_length);
        } else {
            found = true;
            break;
        }
    }
    
    return found;
}

bool base_manager::test_no_pre_sloping_drain(ListDigraph::Node p, rpos length, ListDigraph::ArcMap<rcount> &cp_bwd) {
    
    for ( ListDigraph::OutArcIt a(g, p); a!=INVALID; ++a) {
        
        if (cp_bwd[a] < 10 && regions[a].is_decreasing_region(length)) {
            
            if (length < regions[a].total_length) {
                return true;
            }
            
            // test also adjourning
            for ( ListDigraph::OutArcIt b(g, g.target(a)); b!=INVALID; ++b) {
                if (cp_bwd[b] < 10 && regions[b].is_decreasing_region(length - regions[a].total_length)) {
                    return true;
                }
            }
        } 
    }
    
    return false;
}

bool base_manager::test_consistent_sloping_drain(ListDigraph::Node p, rpos length) {
    
    // test forward
    bool found = false;
    for ( ListDigraph::InArcIt a(g, p); a!=INVALID; ++a) {
        
//        logger::Instance()->debug("Test Inner Source " + std::to_string(g.id(a)) +"\n");
        
        // ignoring circular references and drains, as well as normally corrected edges
        if (edge_type[a] == edge_types::BACKLINK 
            || (edge_type[a] == edge_types::HELPER && g.source(a) == s) 
            || capacity[a] == 0 
            || !regions[a].is_decreasing_region(length)  
                ) { 
            continue;
        }
        
        if (length > regions[a].total_length) {
            found = found || test_consistent_sloping_drain(g.source(a), length -regions[a].total_length);
        } else {
            found = true;
            break;
        }
    }
    
    return found;
}

float base_manager::rec_score_forward(ListDigraph::Arc a) {
    
    unsigned int count = 0;
    for (ListDigraph::InArcIt i(g, g.source(a)); i != INVALID ;++i) {
         ++count;
    }
    if (count > 1) {
        return 0;
    }
    
    float score_left = means[a].compute_score();
    
    float score_right = 0;
    for (ListDigraph::OutArcIt i(g, g.target(a)); i != INVALID ;++i) {
           score_right+=rec_score_forward(i);
    }
    if (score_right == 0) {
        return score_left;
    }
    
    return sqrt(score_left * score_right);
}


float base_manager::rec_score_backward(ListDigraph::Arc a) {
    
    unsigned int count = 0;
    for (ListDigraph::OutArcIt i(g, g.target(a)); i != INVALID ;++i) {
         ++count;
    }
    if (count > 1) {
        return 0;
    }
    
    float score_right = means[a].compute_score();
    
    float score_left = 0;
    for (ListDigraph::InArcIt i(g, g.source(a)); i != INVALID ;++i) {
           score_left+=rec_score_forward(i);
    }
    if (score_left == 0) {
        return score_right;
    }
    
    return sqrt(score_left * score_right);
}

void base_manager::push_potential_forward(ListDigraph::ArcMap<rcount> &cp_fwd, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, ListDigraph::ArcMap<bool> &push_block) {
        
     ListDigraph::ArcMap<rcount> cs(g);  

     for (std::deque<ListDigraph::Node>::iterator n = top_order.begin(); n!= top_order.end(); ++n) {
        
//          logger::Instance()->debug("Push_Node " + std::to_string(g.id(*n)) + "\n");
                  
        // first test if we can propagate
        ListDigraph::InArcIt innt(g, *n);
        bool node_left= (innt != INVALID && edge_type[innt] == edge_types::NODE);
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
            rcount max = regions[node].get_max();
            rcount right = regions[node].get_right();
            left = std::min(regions[node].get_left(), left);
            
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
                if (!push_block[oi]) {
                    total_out += average_push[oi];
                }
                
                if (edge_type[oi] == edge_types::HELPER) {
                    helper = true;
                }
            }
            
//            if (helper) {
//                // find how many sinks we share!
//                unsigned int count = 0;
//                unsigned int calls = 0;
//                count_drains_in_distance(*n, raw->average_fragment_length, count, calls);
//                
////                logger::Instance()->debug("Helper Count " + std::to_string(count) + "\n");
//                
//                right_potential = right_potential / float(count);
//            }
            
//            logger::Instance()->debug("Right Potential " + std::to_string(g.id(*n)) + " " + std::to_string(right_potential) + "\n");
            
            for (ListDigraph::OutArcIt oi(g, *n); oi != INVALID ;++oi) {
                if (push_block[oi]) {
                    cp_fwd[oi] = 0;
                    cs[oi] = 0;
                } else {
                    if (total_out != 0) {
                        cp_fwd[oi] = right_potential * average_push[oi]/float(total_out); // upscaling according to own size
                        cs[oi] = right * average_push[oi]/float(total_out);
                    } else {
                        cp_fwd[oi] = right_potential; // upscaling according to own size
                        cs[oi] = right;
                    }
                }
                
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Set " + std::to_string(g.id(oi)) + " " + std::to_string(cp_fwd[oi]) + "\n");
//                #endif
            }
        }
    }
}


void base_manager::push_potential_backward(ListDigraph::ArcMap<rcount> &cp_bw, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, ListDigraph::ArcMap<bool> &push_block) {
   
     ListDigraph::ArcMap<rcount> cs(g);  

     for (std::deque<ListDigraph::Node>::reverse_iterator n = top_order.rbegin(); n!= top_order.rend(); ++n) {
         // first test if we can propagate
        ListDigraph::OutArcIt onnt(g, *n);
        bool node_right= (onnt != INVALID && edge_type[onnt] == edge_types::NODE);
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
            rcount max = regions[node].get_max();
            right = std::min(regions[node].get_right(), right);
            rcount left = regions[node].get_left();
            
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
                
                if (!push_block[oi]) {
                    total_out += average_push[oi];
                }
                if (edge_type[oi] == edge_types::HELPER) {
                    helper = true;
                }
            }
            
//            if (helper) {
//                // find how many sinks we share!
//                unsigned int count = 0;
//                unsigned int calls = 0;
//                count_sources_in_distance(*n, raw->average_fragment_length, count, calls);
//                
////                logger::Instance()->debug("Helper " + std::to_string(count) + "\n");
//                
//                left_potential = left_potential / float(count);
//            }
//            
//            logger::Instance()->debug("Left Potential " + std::to_string(g.id(*n)) + " " + std::to_string(left_potential) + "\n");
            
            for (ListDigraph::InArcIt oi(g, *n); oi != INVALID ;++oi) {
                if ( push_block[oi]) {
                    cp_bw[oi] = 0;
                    cs[oi] = 0;
                } else {
                    if (total_out != 0) {
                        cp_bw[oi] = left_potential * average_push[oi]/float(total_out); // upscaling according to own size
                        cs[oi] = left * average_push[oi]/float(total_out);
                    } else {
                        cp_bw[oi] = left_potential; // upscaling according to own size
                        cs[oi] = left;                        
                    }
                }
//                
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Set " + std::to_string(g.id(oi)) + " " + std::to_string(cp_bw[oi]) + "\n");
//                #endif
            }
        }
    }
}

void base_manager::count_drains_in_distance(ListDigraph::Node p, rpos average_fragment_length, unsigned int &count, unsigned int &calls) {
    

    if (calls > options::Instance()->get_max_enumerated_paths()) {
        count = 1;
        return;
    }
    
    ListDigraph::OutArcIt a(g, p);
    for ( ; a!=INVALID; ++a) {
        
        if (edge_type[a] == edge_types::BACKLINK) { // ignoring circular references
            continue;
        }
        
        if (edge_type[a] == edge_types::HELPER) {
            ++count;
        }
        
        if (average_fragment_length > regions[a].total_length) {
            count_drains_in_distance(g.target(a), average_fragment_length - regions[a].total_length, count, ++calls);
        }
    }
}


void base_manager::count_sources_in_distance(ListDigraph::Node p, rpos average_fragment_length, unsigned int &count, unsigned int &calls) {
    
    if (calls > options::Instance()->get_max_enumerated_paths()) {
        count = 1;
        return;
    }
    
    ListDigraph::InArcIt a(g, p);
    for ( ; a!=INVALID; ++a) {
        
        if (edge_type[a] == edge_types::BACKLINK) { // ignoring circular references
            continue;
        }
        
        if (edge_type[a] == edge_types::HELPER) {
            ++count;
        }
        
        if (average_fragment_length > regions[a].total_length) {
            count_sources_in_distance(g.source(a), average_fragment_length - regions[a].total_length, count, ++calls);
        }
    }
}


void base_manager::find_over_capacitated_end(ListDigraph::Node p, ListDigraph::ArcMap<rcount> &cp_fwd, std::set<ListDigraph::Node> &nodes) {

    nodes.insert(p);
    ListDigraph::InArcIt a(g, p);
    for ( ; a!=INVALID; ++a) {
        
        if (edge_type[a] == edge_types::BACKLINK) { // ignoring circular references
            continue;
        }
        
        if (edge_type[a] == edge_types::HELPER || capacity[a] == 0 || cp_fwd[a] * 100 / capacity[a] > options::Instance()->get_coverage_bias_limit()) {
            if (nodes.find(g.source(a)) != nodes.end()) { // no double explore to safe time!
                continue;
            }
            find_over_capacitated_end(g.source(a), cp_fwd, nodes);
        }
    }
}

void base_manager::find_over_capacitated_start(ListDigraph::Node p, ListDigraph::ArcMap<rcount> &cp_fwd, std::set<ListDigraph::Node> &nodes) {

    nodes.insert(p);
    ListDigraph::OutArcIt a(g, p);
    for ( ; a!=INVALID; ++a) {
        
        if (edge_type[a] == edge_types::BACKLINK) { // ignoring circular references
            continue;
        }
        
//        logger::Instance()->debug("Find START " + std::to_string(g.id(a)) + " " + std::to_string(cp_fwd[a])+ " " + std::to_string(capacity[a])+";\n");
        
        if (edge_type[a] == edge_types::HELPER || capacity[a] == 0 || cp_fwd[a] * 100 / capacity[a] > options::Instance()->get_coverage_bias_limit()) {
            if (nodes.find(g.target(a)) != nodes.end()) { // no double explore to safe time!
                continue;
            }
            find_over_capacitated_start(g.target(a), cp_fwd, nodes);
        }
    }
}


void base_manager::denoise_graph_alt(ListDigraph::ArcMap<bool> &guided_saves, ListDigraph::ArcMap<bool> &push_block) {
    
    // first we get a topological sorting of all nodes
    std::deque<ListDigraph::Node> top_order;
    pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > visitor(top_order);
    DfsVisit<ListDigraph, pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > > dfs(g, visitor);
    dfs.init();
    dfs.addSource(s);
    dfs.start();

    ListDigraph::ArcMap<rcount> average_push(g); // init as 0, we set a value
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) { // init as the average!
        average_push[a]= regions[a].get_average();
    }
    
    // forward potential!
    ListDigraph::ArcMap<rcount> cp_fwd(g); // init as 0
    push_potential_forward_alt(cp_fwd, top_order, average_push, push_block);
    // backward potential!
    ListDigraph::ArcMap<rcount> cp_bw(g); // init as 0
    push_potential_backward_alt(cp_bw, top_order, average_push, push_block);
 
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("add " + std::to_string(g.id(a)) + " " + std::to_string(cp_fwd[a]) + " " + std::to_string(cp_bw[a]) + "\n");
        #endif
        
        capacity_type max = std::max(cp_fwd[a], cp_bw[a]);
        if (edge_type[a] != edge_types::HELPER) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Correct Capacity " + std::to_string(g.id(a)) + " from " + std::to_string(capacity[a]) + " to " + std::to_string(capacity[a]+max) + ".\n");
            #endif
            capacity[a] = max;
        }
    }
}

void base_manager::push_potential_forward_alt(ListDigraph::ArcMap<rcount> &cp_fwd, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, ListDigraph::ArcMap<bool> &push_block) {
        
     for (std::deque<ListDigraph::Node>::iterator n = top_order.begin(); n!= top_order.end(); ++n) {
        
         #ifdef ALLOW_DEBUG 
        logger::Instance()->debug("Push_Node " + std::to_string(g.id(*n)) + "\n");
        #endif
                  
        // first test if we can propagate
        ListDigraph::InArcIt innt(g, *n);
        bool node_left= (innt != INVALID && edge_type[innt] == edge_types::NODE);
        if (node_left) {
            
            // this means we can only have one valid left, to of course multiple outs edges
            ListDigraph::Arc node(innt);
            
            rcount in_max = 0;
            for (ListDigraph::InArcIt inpi(g, g.source(node)); inpi != INVALID ;++inpi) {
                in_max += cp_fwd[inpi];
            }
            
            rcount node_max = regions[node].get_max();
            rcount max = std::max(node_max, in_max);
           
            // we set the push
            cp_fwd[node] = max; 
            
            // we push the potential to the next adjourning edges
            rcount total_out = 0;
            bool helper = false;
            for (ListDigraph::OutArcIt oi(g, *n); oi != INVALID ;++oi) {
                if (!push_block[oi]) {
                    total_out += average_push[oi];
                }
                
                if (edge_type[oi] == edge_types::HELPER) {
                    helper = true;
                }
            }
            
            if (helper) {
                // find how many sinks we share!
                if (regions[node].get_right() > total_out) {
                    total_out += regions[node].get_right() - total_out;
                }
            }
                        
            for (ListDigraph::OutArcIt oi(g, *n); oi != INVALID ;++oi) {
                if (push_block[oi]) {
                    cp_fwd[oi] = 0;
                } else {
                    if (total_out != 0) {
                        capacity_type nmax = max * average_push[oi]/float(total_out);
                        cp_fwd[oi] = std::max( nmax , capacity[oi]); // upscaling according to own size
                    } else {
                        cp_fwd[oi] = std::max( max, capacity[oi]); 
                    }
                }
                
                #ifdef ALLOW_DEBUG 
                logger::Instance()->debug("Set " + std::to_string(g.id(oi)) + " " + std::to_string(cp_fwd[oi]) + "\n");
                #endif
            }
        }
    }
}

void base_manager::push_potential_backward_alt(ListDigraph::ArcMap<rcount> &cp_bwd, std::deque<ListDigraph::Node> &top_order, ListDigraph::ArcMap<rcount> &average_push, ListDigraph::ArcMap<bool> &push_block) {
        
     for (std::deque<ListDigraph::Node>::reverse_iterator n = top_order.rbegin(); n!= top_order.rend(); ++n) {
         // first test if we can propagate
        ListDigraph::OutArcIt onnt(g, *n);
        bool node_right= (onnt != INVALID && edge_type[onnt] == edge_types::NODE);
        if (node_right) {
            
             #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Push_Node " + std::to_string(g.id(*n)) + "\n");
            #endif

            // this means we can only have one valid left, to of course multiple outs edges
            ListDigraph::Arc node(onnt);
            
            rcount out_max = 0;
            for (ListDigraph::OutArcIt inpi(g, g.target(node)); inpi != INVALID ;++inpi) {
                out_max += cp_bwd[inpi];
            }
            
            rcount node_max = regions[node].get_max();
            rcount max = std::max(node_max, out_max);
           
            // we set the push
            cp_bwd[node] = max; 
            
            // we push the potential to the next adjourning edges
            rcount total_in = 0;
            bool helper = false;
            for (ListDigraph::InArcIt oi(g, *n); oi != INVALID ;++oi) {
                
                if (!push_block[oi]) {
                    total_in += average_push[oi];
                }
                if (edge_type[oi] == edge_types::HELPER) {
                    helper = true;
                }
            }
            
            if (helper) {
                // find how many sinks we share!
                if (regions[node].get_right() > total_in) {
                    total_in += regions[node].get_left() - total_in;
                }
            }
                        
            for (ListDigraph::InArcIt oi(g, *n); oi != INVALID ;++oi) {
                if (push_block[oi]) {
                    cp_bwd[oi] = 0;
                } else {
                    if (total_in != 0) {
                        capacity_type nmax = max * average_push[oi]/float(total_in);
                        cp_bwd[oi] = std::max( nmax, capacity[oi]); // upscaling according to own size
                    } else {
                        cp_bwd[oi] = std::max( max, capacity[oi]); 
                    }
                }
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Set " + std::to_string(g.id(oi)) + " " + std::to_string(cp_bwd[oi]) + "\n");
                #endif
            }
        }
    }
}

void base_manager::denoise_graph_guided(std::deque<std::deque<ListDigraph::Arc> > &guides) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Graph Denoise Guided \n");
    #endif
    
    ListDigraph::NodeMap<float> ratios(g);
    ListDigraph::ArcMap<rcount> guided_capacity(g);
    ListDigraph::ArcMap<float> guided_percentage(g);
    
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
        
        bool helper = false;
        
        rcount in_push = 0;
        for (ListDigraph::InArcIt inpi(g, n); inpi != INVALID ;++inpi) {
            if (edge_type[inpi] == edge_types::HELPER) {
                helper = true;
            }
            in_push += capacity[inpi];
        }
        
        rcount out_push = 0;
        for (ListDigraph::OutArcIt onpi(g, n); onpi != INVALID ;++onpi) {
            if (edge_type[onpi] == edge_types::HELPER) {
                helper = true;
            }
            out_push += capacity[onpi];
        }
        if (helper || out_push == 0 || in_push == 0) {
            ratios[n] = 1;
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
            
            if (edge_type[*a_it] == edge_types::HELPER) {
                continue;
            }
            
            float fac = ratios[g.source(*a_it)];
            way_factor = way_factor * fac;
            rcount local_limit = std::ceil(capacity[*a_it] * 1/way_factor);
            
//            logger::Instance()->debug("Overshoot Test " + std::to_string(g.id(*a_it)) + " " + std::to_string(local_limit) + " " + std::to_string(fac) + " " + std::to_string(way_factor) + "\n");
            if (local_limit < guide_cap || guide_cap == 0) {
                guide_cap = local_limit ;
            }
        }
        
//        logger::Instance()->debug("Guide Cap " + std::to_string(guide_cap) + "\n");
        
        // set guide_increase in second and third run!
        rcount mean = 0;
        rcount guide_current = guide_cap;
        a_it = g_it->begin();
        mean += guide_current;
        
        if (edge_type[*a_it] == edge_types::HELPER) {
            guided_percentage[*a_it] = 1;
        } else {
            if (capacity[*a_it] > 0) {
                guided_percentage[*a_it] += guide_current/float(capacity[*a_it]); 
            } else {
                guided_percentage[*a_it] = 1;
            }
        }
        
        ++a_it;
        for (; a_it != g_it->end(); ++a_it) {
            
            float fac = ratios[g.source(*a_it)];
            guide_current = std::ceil(guide_current * fac);
            mean += guide_current;
            
//            logger::Instance()->debug("Current " + std::to_string(g.id(*a_it)) + " " + std::to_string(guide_current) + " " + std::to_string(fac) + "\n");
            
            if (edge_type[*a_it] == edge_types::HELPER) {
                guided_percentage[*a_it] = 1;
            } else {
                if (capacity[*a_it] > 0) {
                    guided_percentage[*a_it] += guide_current/float(capacity[*a_it]); 
                } else {
                    guided_percentage[*a_it] = 1;
                }
            }
        }
        mean = mean / (g_it->size() - 1);
        
//        logger::Instance()->debug("Guide Mean " + std::to_string(mean) + "\n");
        
        for (std::deque<ListDigraph::Arc>::iterator f_it = g_it->begin(); f_it != g_it->end(); ++f_it) {
            guided_capacity[*f_it] += mean;
        }
    }
    
    float truth_limit = 1 - options::Instance()->get_guide_trust_level();
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Set " + std::to_string(g.id(a)) + " " + std::to_string(guided_percentage[a]) + " " + std::to_string(guided_capacity[a]) + "\n");
        #endif
        
        if (guided_percentage[a] >= truth_limit) {
            capacity[a] = guided_capacity[a];
        } else if (guided_percentage[a] > 0) {
            capacity[a] = guided_capacity[a] * 1/guided_percentage[a];
        }
    } 
}

void base_manager::correct_start_end() {

    bool correctable = true;
    
    // as a first step, build up ALL paths that need correction
    // STARTS
    std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > start_region;
    std::deque<float> initial_start;
    for(ListDigraph::OutArcIt a(g, s); a != INVALID; ++a) {
        start_region.push_back(std::make_pair(std::deque<ListDigraph::Arc>(), 0) );
        correctable &= find_starts(g.target(a), start_region, raw->average_fragment_length, start_region.end()-1);
    }
    for (std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator si = start_region.begin(); si != start_region.end(); ++si) {
        initial_start.push_back(si->second);
    }
    
    // ENDS
    std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > end_region; 
    std::deque<float> initial_end;
    for(ListDigraph::InArcIt a(g, t); a != INVALID; ++a) {
        end_region.push_back(std::make_pair(std::deque<ListDigraph::Arc>(), 0) );
        correctable &= find_ends(g.source(a), end_region, raw->average_fragment_length, end_region.end()-1);
    }
    for (std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator si = end_region.begin(); si != end_region.end(); ++si) {
        initial_end.push_back(si->second);
    }
    
    // next we get the estimated corrections for the first run!
    
    ListDigraph::ArcMap<float > capacity_correction_start(g);
    ListDigraph::ArcMap<float > capacity_correction_end(g);
    if (correctable) {
        recursive_start_correction(start_region, capacity_correction_start);
        recursive_end_correction(end_region, capacity_correction_end);
    }
    
    bool changes = false;
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) { // just for now!    //TODODODODODODDO
        
        capacity_type in = 0;
        capacity_type out = 0;
        bool skip = false;
        for(ListDigraph::InArcIt i(g,n);i!=INVALID;++i) {
            
            if (edge_type[i] == edge_types::HELPER) {
                skip = true;
                break;
            }
            capacity_type max_alt = std::round(std::max(capacity_correction_start[i], capacity_correction_end[i]));
            in += std::max(capacity[i], max_alt);
        }
        if (skip) {
            continue;
        }
        for(ListDigraph::OutArcIt o(g,n);o!=INVALID;++o) {
            if (edge_type[o] == edge_types::HELPER) {
                skip = true;
                break;
            }
            capacity_type max_alt = std::round(std::max(capacity_correction_start[o], capacity_correction_end[o]));
            out += std::max(capacity[o], max_alt);
        }
        if (skip) {
            continue;
        }
        
        rcount n1 = std::max(in, out);
        rcount n2 = std::min(in, out);

        ListDigraph::InArcIt innt(g, n);
        bool node_left= (innt != INVALID && edge_type[innt] == edge_types::NODE);
        if ( (n1 - n2) * 100 / float(n2) >= options::Instance()->get_coverage_change_limit() ) {

//            logger::Instance()->debug("In " + std::to_string(in) + " Out " + std::to_string(out) + " at " + std::to_string(g.id(n)) + "\n" );
            
            if (in > out) {
                if (node_left) {
//                    logger::Instance()->debug("Insert new Drain: " + std::to_string(g.id(n)) + "\n" );

                    // down trend
                    ListDigraph::Arc na = g.addArc(n, t);
                    capacity[na] = 0;
                    edge_type[na] = edge_types::HELPER;
                    end_region.push_back(std::make_pair(std::deque<ListDigraph::Arc>(), 0) );
                    correctable &= find_ends(n, end_region, raw->average_fragment_length, end_region.end()-1);
                    changes = true;
                }
            } else {
                if (!node_left) {
//                    logger::Instance()->debug("Insert new Source: " + std::to_string(g.id(n)) + "\n" );

                    ListDigraph::Arc na = g.addArc(s, n);
                    capacity[na] = 0;
                    edge_type[na] = edge_types::HELPER;
                    start_region.push_back(std::make_pair(std::deque<ListDigraph::Arc>(), 0) );
                    correctable &= find_starts(n, start_region, raw->average_fragment_length, start_region.end()-1);
                    changes = true;
                }
            }
        }
    }
    
    if (changes && correctable) {
        
        for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
            capacity_correction_start[a] = 0;
            capacity_correction_end[a] = 0;
        }
        
        std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator sri = start_region.begin();
        std::deque<float>::iterator sci = initial_start.begin();
        for (;sci != initial_start.end(); ++sri, ++sci) {
            sri->second = *sci;
        }
        
        std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator eri = end_region.begin();
        std::deque<float>::iterator eci = initial_end.begin();
        for (;eci != initial_end.end(); ++eri, ++eci) {
            eri->second = *eci;
        }
        
        // so we need to do the
        recursive_start_correction(start_region, capacity_correction_start);
        recursive_end_correction(end_region, capacity_correction_end);
    }
    
    // finally upcorrect those that are too small!
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        capacity_type max = std::max(std::round(capacity_correction_start[a]), std::round(capacity_correction_end[a]));
        if (edge_type[a] != edge_types::HELPER && capacity[a] < max) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Correct Capacity " + std::to_string(g.id(a)) + " from " + std::to_string(capacity[a]) + " to " + std::to_string(max) + ".\n");
            #endif
            capacity[a] = max;
        }
    }
}

void base_manager::recursive_start_correction(std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > &start_region, ListDigraph::ArcMap<float > &capacity_correction) {

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Start Correction " + std::to_string(start_region.size()) + " " + std::to_string(raw->average_fragment_length) + ".\n");
    #endif
    
    ListDigraph::ArcMap<std::vector<unsigned int> > startcount(g);
    ListDigraph::NodeMap<unsigned int > node_available(g); // a "resolved by" relation
    ListDigraph::NodeMap<unsigned int > node_used(g); // a "resolved by" relation
    std::set<int> order_tmp;
    
    // we have all possible lines
    unsigned int i = 0;
    for ( std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator it = start_region.begin(); it!= start_region.end(); ++it, ++i) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("S-Path " + std::to_string(it->second) + " a:");
        #endif
        
        std::deque<ListDigraph::Arc>::iterator in_it = it->first.begin();
        ++node_available[g.source(*in_it)];
        for (; in_it != it->first.end(); ++in_it) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug(" " + std::to_string(g.id(*in_it)));
            #endif
            
            startcount[*in_it].push_back(i);
            ++node_available[g.target(*in_it)];
            if (in_it+1 == it->first.end()) {
                // so this is the last one
                ++node_used[g.target(*in_it)];
                order_tmp.insert(g.id(g.target(*in_it)));
                
                 #ifdef ALLOW_DEBUG
                logger::Instance()->debug(" I(" + std::to_string(g.id(g.target(*in_it)))+")");
                #endif
                
            }
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("\n");
        #endif
    }   
    
    std::deque<ListDigraph::Node> order;
    while (!order_tmp.empty()) {
        
        int next_id;
        ListDigraph::Node n;
        for(std::set<int>::iterator n_it = order_tmp.begin(); n_it != order_tmp.end(); ++n_it) {
            
            next_id = *n_it;
            n = g.nodeFromId(next_id);
            if (node_used[n] == node_available[n]) {
                order_tmp.erase(n_it);
                break;
            }
        }
        
        order.push_back(n);
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Order Node " + std::to_string(g.id(n))+ "\n");
        #endif
        
        for (ListDigraph::InArcIt a(g, n); a != INVALID; ++a) {
            
            if (startcount[a].size() == 0) {
                continue;
            }
            
            node_used[g.source(a)] += startcount[a].size();
//            logger::Instance()->debug("Test Add " + std::to_string(g.id(g.source(a))) + " " + std::to_string(node_used[g.source(a)]) + " " + std::to_string(node_available[g.source(a)]) +"\n");
            if (node_used[g.source(a)] == node_available[g.source(a)]) {
                order_tmp.insert(g.id(g.source(a)));
            }
        }
    }
    
    for (std::deque<ListDigraph::Node >::iterator ni = order.begin(); ni != order.end(); ++ni) {
        
        rcount outcov = 0;
        for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
            outcov += regions[a].get_average();
        }
        rcount incov_unused = 0;
        rcount incov_used = 0;
        bool has_source = false;
        bool only_source = true;
        for(ListDigraph::InArcIt a(g, *ni); a!=INVALID; ++a) {
            if (edge_type[a] == edge_types::HELPER) {
                has_source = true;
            } else if (startcount[a].empty()) {
                incov_unused += regions[a].get_average();
                only_source = false;
            } else {
                incov_used += regions[a].get_average();
                only_source = false;
            }
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("NODE " + std::to_string(g.id(*ni)) + " " + std::to_string(incov_used) + " " + std::to_string(incov_unused) + " " + std::to_string(outcov) + ".\n");
        #endif
        
        if (only_source || incov_used == 0 || outcov == 0) {
            continue;
        }
        
        rcount offset = 0; // flow that needs to be removed in total!
        if(incov_unused < outcov) {
            offset = incov_unused;
        } else {
            // in this case disable ALL
            for(ListDigraph::InArcIt a(g, *ni); a!=INVALID; ++a) {
                if (!startcount[a].empty() && edge_type[a] != edge_types::HELPER) {
                    for (std::vector<unsigned int>::iterator ids = startcount[a].begin(); ids != startcount[a].end(); ++ids) {
                        start_region[*ids].second = 0;
                    }
                }
            }  
            continue;
        }
        
        if (has_source) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Has Source .\n");
            #endif
            
            rcount max_ending = 0;
            rcount total_ending = 0;
            rcount max_overlap = 0;
            for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
                for (std::vector<unsigned int>::iterator ids = startcount[a].begin(); ids != startcount[a].end(); ++ids) {
                    if (start_region[*ids].first.front() == a) { // this one ends here!
                        if (start_region[*ids].second > max_ending) {
                            max_ending = start_region[*ids].second;
                        }
                        total_ending += start_region[*ids].second;
                    } else {
                        if (start_region[*ids].second > max_overlap) {
                            max_overlap = start_region[*ids].second;
                        }
                    }
                }
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Found " + std::to_string(max_overlap) + " " + std::to_string(max_ending) + " " + std::to_string(total_ending) +"\n");
            #endif
            if (total_ending > max_overlap && total_ending != 0) { // else something is really wrong in mapping, we can't fix that...
                rcount dissipation = total_ending - max_overlap;
                for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
                    for (std::vector<unsigned int>::iterator ids = startcount[a].begin(); ids != startcount[a].end(); ++ids) {
                        if (start_region[*ids].first.front() == a) {
                            start_region[*ids].second = dissipation * start_region[*ids].second/total_ending;
                        }
                    }
                }
                offset += dissipation;
            } else {
                // we kill the stopping arcs??
                for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
                    for (std::vector<unsigned int>::iterator ids = startcount[a].begin(); ids != startcount[a].end(); ++ids) {
                        if (start_region[*ids].first.front() == a) {
                            start_region[*ids].second = 0;
                        }
                    }
                }
            }
        }
        
        // now go over the individual splits to left
        
        float unreduced_offset = 0;
        float global_leftovers = 0; 
        for(ListDigraph::InArcIt a(g, *ni); a!=INVALID; ++a) {
            if (!startcount[a].empty() && edge_type[a] != edge_types::HELPER) {
                float percentage_modifier = regions[a].get_average()/float(incov_used);
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Modifier " + std::to_string(regions[a].get_average()) + "/" + std::to_string(incov_used) + "=" + std::to_string(percentage_modifier) + " o" + std::to_string(offset) + " cov" + std::to_string(outcov) + " s " + std::to_string(startcount[a].size())  + ".\n");
                #endif
               
                for (std::vector<unsigned int>::iterator ids = startcount[a].begin(); ids != startcount[a].end(); ++ids) {
                    float dedicated_offset = offset/float(startcount[a].size()) * percentage_modifier;
                    float new_value = start_region[*ids].second * percentage_modifier;
                    
                    if (dedicated_offset > new_value) {
                        start_region[*ids].second = 0;
                        unreduced_offset += dedicated_offset - new_value;
                    } else {
                        start_region[*ids].second = new_value - dedicated_offset;
                    }
                    global_leftovers += start_region[*ids].second;
                }
            }
        }  
        
        if (unreduced_offset >= global_leftovers) {
            for(ListDigraph::InArcIt a(g, *ni); a!=INVALID; ++a) {
                if (!startcount[a].empty() && edge_type[a] != edge_types::HELPER) {
                    for (std::vector<unsigned int>::iterator ids = startcount[a].begin(); ids != startcount[a].end(); ++ids) {
                        start_region[*ids].second = 0;
                    }
                }
            }  
        } else {
            for(ListDigraph::InArcIt a(g, *ni); a!=INVALID; ++a) {
                if (!startcount[a].empty() && edge_type[a] != edge_types::HELPER) {
                    for (std::vector<unsigned int>::iterator ids = startcount[a].begin(); ids != startcount[a].end(); ++ids) {
                        start_region[*ids].second -= unreduced_offset * start_region[*ids].second/global_leftovers;
                    }
                }
            }  
        }
        
        
        #ifdef ALLOW_DEBUG
        int j = 0;
        for ( std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator it = start_region.begin(); it!= start_region.end(); ++it, ++j) {
            logger::Instance()->debug("Intermediate Path Value " + std::to_string(j) + " " + std::to_string(it->second) + ".\n");
        }
        #endif
    }
    
    // now all needs to be added up
    for ( std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator it = start_region.begin(); it!= start_region.end(); ++it, ++i) {

        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Finalized Path Value " + std::to_string(it->second) + ".\n");
        #endif
        
        for (std::deque<ListDigraph::Arc>::iterator in_it = it->first.begin(); in_it != it->first.end(); ++in_it) {
            capacity_correction[*in_it] += it->second;
        }
    } 
    
}

bool base_manager::find_starts(ListDigraph::Node p, std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > &start_region, rpos average_fragment_length, std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator sri) {
    
    if ( start_region.size() > options::Instance()->get_max_enumerated_paths() ) {
        return false;
    }
    
    ListDigraph::OutArcIt a(g, p);
    
    while (a != INVALID && (edge_type[a] == edge_types::BACKLINK || (edge_type[a] == edge_types::HELPER && g.target(a) == t)) ) {
        ++a;
    }
    
    if (a == INVALID) {
        return true;
    }
    ListDigraph::Arc arc_save(a);
    ++a;
    for ( ; a!=INVALID; ++a) {
        
        if (edge_type[a] == edge_types::BACKLINK || (edge_type[a] == edge_types::HELPER && g.target(a) == t)) { // ignoring circular references and drains, otherwise we have a BAD time
            continue;
        }
        
        start_region.push_back(*sri);
        start_region.rbegin()->first.push_back(a);
        start_region.rbegin()->second = regions[a].get_pos_max(average_fragment_length);
        if (average_fragment_length > regions[a].total_length) {
            find_starts(g.target(a), start_region, average_fragment_length -regions[a].total_length, start_region.end()-1);
        }
    }
    
    sri->first.push_back(arc_save);
    sri->second = regions[arc_save].get_pos_max(average_fragment_length);
    if (average_fragment_length > regions[arc_save].total_length) {
        find_starts(g.target(arc_save), start_region, average_fragment_length -regions[arc_save].total_length, sri);
    }
    
    return true;
}

void base_manager::recursive_end_correction( std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > &end_region, ListDigraph::ArcMap<float > &capacity_correction) {

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("End Correction " + std::to_string(end_region.size())  + " " + std::to_string(raw->average_fragment_length) + ".\n");
    #endif
    
    ListDigraph::ArcMap<std::vector<unsigned int> > endcount(g);
    ListDigraph::NodeMap<unsigned int > node_available(g); // a "resolved by" relation
    ListDigraph::NodeMap<unsigned int > node_used(g); // a "resolved by" relation
    std::set<int> order_tmp;
    
    // we have all possible lines
    unsigned int i = 0;
    for ( std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator it = end_region.begin(); it!= end_region.end(); ++it, ++i) {
        
         #ifdef ALLOW_DEBUG
        logger::Instance()->debug("E-Path " + std::to_string(it->second) + " a:");
        #endif
        
        std::deque<ListDigraph::Arc>::iterator in_it = it->first.begin();
        ++node_available[g.target(*in_it)];
        for (; in_it != it->first.end(); ++in_it) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug(" " + std::to_string(g.id(*in_it)));
            #endif
            
            endcount[*in_it].push_back(i);
            ++node_available[g.source(*in_it)];
            if (in_it == it->first.end()-1) {
                // so this is the last one
                ++node_used[g.source(*in_it)];
                order_tmp.insert(g.id(g.source(*in_it)));
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("I(" + std::to_string(g.id(g.source(*in_it)))+ ")");
                #endif
            }
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("\n");
        #endif 
    }
    
    std::deque<ListDigraph::Node> order;
    while (!order_tmp.empty()) {
        
        
        int next_id;
        ListDigraph::Node n;
        for(std::set<int>::iterator n_it = order_tmp.begin(); n_it != order_tmp.end(); ++n_it) {
            
            next_id = *n_it;
            n = g.nodeFromId(next_id);
            if (node_used[n] == node_available[n]) {
                order_tmp.erase(n_it);
                break;
            }
        }

        order.push_back(n);
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Order Node " + std::to_string(g.id(n))+ "\n");
        #endif
        
        for (ListDigraph::OutArcIt a(g, n); a != INVALID; ++a) {
            
            if (endcount[a].size() == 0) {
                continue;
            }
            
            node_used[g.target(a)] += endcount[a].size();
//            logger::Instance()->debug("Test Add " + std::to_string(g.id(g.target(a))) + " " + std::to_string(node_used[g.target(a)]) + " " + std::to_string(node_available[g.target(a)]) +"\n");
            if (node_used[g.target(a)] == node_available[g.target(a)]) {
                order_tmp.insert(g.id(g.target(a)));
            }
        }
    }
    
    for (std::deque<ListDigraph::Node >::iterator ni = order.begin(); ni != order.end(); ++ni) {
        
        rcount incov = 0;
        for(ListDigraph::InArcIt a(g, *ni); a!=INVALID; ++a) {
            incov += regions[a].get_average();
        }
        rcount outcov_unused = 0;
        rcount outcov_used = 0;
        bool has_drain = false;
        bool only_drain = true;
        for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
            if (edge_type[a] == edge_types::HELPER) {
                has_drain = true;
            } else if (endcount[a].empty()) {
                outcov_unused += regions[a].get_average();
                only_drain = false;
            } else {
                outcov_used += regions[a].get_average();
                only_drain = false;
            }
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("NODE " + std::to_string(g.id(*ni)) + " " + std::to_string(outcov_used) + " " + std::to_string(outcov_unused) + " " + std::to_string(incov) + ".\n");
        #endif
        
        if (only_drain || outcov_used == 0 || incov == 0) {
            continue;
        }
        
        rcount offset = 0;
        if(outcov_unused < incov) {
            offset = outcov_unused;
        } else {
            for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
                if (!endcount[a].empty() && edge_type[a] != edge_types::HELPER) {
                    for (std::vector<unsigned int>::iterator ids = endcount[a].begin(); ids != endcount[a].end(); ++ids) {
                        end_region[*ids].second = 0;
                    }
                }
            }  
            continue;
        }
        
         if (has_drain) {
             
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Has drain .\n");
            #endif
             
            rcount max_ending = 0;
            rcount total_ending = 0;
            rcount max_overlap = 0;
            for(ListDigraph::InArcIt a(g, *ni); a!=INVALID; ++a) {
                for (std::vector<unsigned int>::iterator ids = endcount[a].begin(); ids != endcount[a].end(); ++ids) {
                    if (end_region[*ids].first.front() == a) { // this one ends here!
                        if (end_region[*ids].second > max_ending) {
                            max_ending = end_region[*ids].second;
                        }
                        total_ending += end_region[*ids].second;
                    } else {
                        if (end_region[*ids].second > max_overlap) {
                            max_overlap = end_region[*ids].second;
                        }
                    }
                }
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Found " + std::to_string(max_overlap) + " " + std::to_string(max_ending) + " " + std::to_string(total_ending) + "\n");
            #endif
            if (total_ending > max_overlap && total_ending != 0) { // something is wrong here otherwise
                rcount dissipation = total_ending - max_overlap;
                for(ListDigraph::InArcIt a(g, *ni); a!=INVALID; ++a) {
                    for (std::vector<unsigned int>::iterator ids = endcount[a].begin(); ids != endcount[a].end(); ++ids) {
                        if (end_region[*ids].first.front() == a) {
                             end_region[*ids].second = dissipation * end_region[*ids].second/total_ending;
                        }
                    }
                }
                offset += dissipation;
            } else {
                for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
                    for (std::vector<unsigned int>::iterator ids = endcount[a].begin(); ids != endcount[a].end(); ++ids) {
                        if (end_region[*ids].first.front() == a) {
                            end_region[*ids].second = 0;
                        }
                    }
                }
            }
        }
        
        // now go over the individual splits to left
        float unreduced_offset = 0;
        float global_leftovers = 0; 
        for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
            if (!endcount[a].empty() && edge_type[a] != edge_types::HELPER) {
                float percentage_modifier = regions[a].get_average()/float(outcov_used);
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Modifier " + std::to_string(regions[a].get_average()) + "/" + std::to_string(outcov_used) + "=" + std::to_string(percentage_modifier) + " o" + std::to_string(offset) + " cov" + std::to_string(incov) + " s " + std::to_string(endcount[a].size())  + ".\n");
                #endif
               
                for (std::vector<unsigned int>::iterator ids = endcount[a].begin(); ids != endcount[a].end(); ++ids) {
                    float dedicated_offset = offset/float(endcount[a].size()) * percentage_modifier;
                    float new_value = end_region[*ids].second * percentage_modifier;
                    
                    if (dedicated_offset > new_value) {
                        end_region[*ids].second = 0;
                        unreduced_offset += dedicated_offset - new_value;
                    } else {
                        end_region[*ids].second = new_value - dedicated_offset;
                    }
                    global_leftovers += end_region[*ids].second;
                }
            }
        }  
        
         #ifdef ALLOW_DEBUG
         logger::Instance()->debug("UnRevStats " + std::to_string(unreduced_offset) + " " + std::to_string(global_leftovers) + "\n");
         #endif
        
        if (unreduced_offset >= global_leftovers) {
            for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
                if (!endcount[a].empty() && edge_type[a] != edge_types::HELPER) {
                    for (std::vector<unsigned int>::iterator ids = endcount[a].begin(); ids != endcount[a].end(); ++ids) {
                        end_region[*ids].second = 0;
                    }
                }
            }  
        } else {
            for(ListDigraph::OutArcIt a(g, *ni); a!=INVALID; ++a) {
                if (!endcount[a].empty() && edge_type[a] != edge_types::HELPER) {
                    for (std::vector<unsigned int>::iterator ids = endcount[a].begin(); ids != endcount[a].end(); ++ids) {
                        end_region[*ids].second -= unreduced_offset * end_region[*ids].second/global_leftovers;
                    }
                }
            }  
        }
        
        #ifdef ALLOW_DEBUG
        for ( std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator it = end_region.begin(); it!= end_region.end(); ++it, ++i) {
            logger::Instance()->debug("Intermediate Path Value " + std::to_string(it->second) + ".\n");
        }
        #endif
        
    }
    
    // now all needs to be added up
    for ( std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator it = end_region.begin(); it!= end_region.end(); ++it, ++i) {
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Finalized Path Value " + std::to_string(it->second) + ".\n");
        #endif
        
        for (std::deque<ListDigraph::Arc>::iterator in_it = it->first.begin(); in_it != it->first.end(); ++in_it) {
            capacity_correction[*in_it] += it->second;
        }
    } 
    
}

bool base_manager::find_ends(ListDigraph::Node p, std::deque< std::pair<std::deque<ListDigraph::Arc>, float> > &end_region, rpos average_fragment_length, std::deque< std::pair<std::deque<ListDigraph::Arc>, float> >::iterator sri) {
    
    
    if ( end_region.size() > options::Instance()->get_max_enumerated_paths() ) {
        return false;
    }
    
    ListDigraph::InArcIt a(g, p);
    
    while (a != INVALID && (edge_type[a] == edge_types::BACKLINK || (edge_type[a] == edge_types::HELPER && g.source(a) == s)) ) {
        ++a;
    }
    
    if (a == INVALID) {
        return true;
    }
    ListDigraph::Arc arc_save(a);
    ++a;
    for ( ; a!=INVALID; ++a) {
        
        if (edge_type[a] == edge_types::BACKLINK || (edge_type[a] == edge_types::HELPER && g.source(a) == s)) { // ignoring circular references and drains, otherwise we have a BAD time
            continue;
        }
        
        end_region.push_back(*sri);
        end_region.rbegin()->first.push_back(a);
        end_region.rbegin()->second = regions[a].get_pos_max_reverse(average_fragment_length);
        if (average_fragment_length > regions[a].total_length) {
            find_ends(g.source(a), end_region, average_fragment_length -regions[a].total_length, end_region.end()-1);
        }
    }
    
    sri->first.push_back(arc_save);
    sri->second = regions[arc_save].get_pos_max_reverse(average_fragment_length);
    if (average_fragment_length > regions[arc_save].total_length) {

        find_ends(g.source(arc_save), end_region, average_fragment_length -regions[arc_save].total_length, sri);
    }
    
    return true;
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
        
        if (edge_type[a] == edge_types::BACKLINK) { // ignoring circular references, otherwise we have a BAD time
            continue;
        }
        
        bool correct_edge = false;
        if (edge_type[a] != edge_types::EXON) { // all edges except for circular and labeled exon edges
            // just jump add
            
            bool correct = exon;
            if (edge_specifier[a].node_index == next_start && edge_specifier[a].node_index == goal_end_index) {
                correct = true;
            }
            
            correct_edge = recursive_arc_backtracing(goal, next_start, goal_end_index, g.target(a), a, path, correct);
        } else {
            
            unsigned int start_index = edge_specifier[a].id.find_first();
            
            // this is an EXON here
            if (start_index != next_start) {
                continue;
            }
            
            // for guides path evidences do not matter!
            
            unsigned int end_index; // this is slow but meh, doesn't happen too often
            boost::dynamic_bitset<>::size_type index = start_index;
            while(index != boost::dynamic_bitset<>::npos) {
                end_index = index;
                index = edge_specifier[a].id.find_next(index);
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test Cont " + edge_specifier[a].to_string() + " vs " + goal.to_string() + " ; " + std::to_string(start_index) + "-"+  std::to_string(end_index) +".\n");
            #endif
            // exon edge, we only follow this if it is compliant to the goal
            if ( edge_specifier[a].is_contained_in(goal, start_index, end_index) ) {
                
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



bool base_manager::contract_composites(InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg) {
    
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("contract_composites.\n");
    #endif
    
    bool changed = false;
    
     // loop over all node that can be a potential start
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
        // we can erase while looping, but be careful
       #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Before\n");
        digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("flow", flow)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run();  
        #endif
        
        ListDigraph::Node node(n);
                
        if ( out_deg[node] != 1 || in_deg[node] != 1) { // not an in between node from in perspective
            // test backwards arcs
            
            for (ListDigraph::OutArcIt a(g, node); a!=INVALID; ) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Iterator " + std::to_string(g.id(a)) + "\n");
                #endif
    
                ListDigraph::Node next = g.target(a);

                if ( (in_deg[next] == 1 || out_deg[next] == 1) && out_deg[next] != 0) { // not an in between node from out perspective

                    // this is the perfect seed to get a path

                    changed = true;

                    // we do the merge
                    // this means we build a composite edge from all encountered

                    capacity_mean path_mean = means[a];
                    capacity_type path_flow = flow[a];
                    exon_edge path_ident = edge_specifier[a];
                    edge_length path_length = edge_lengths[a];
                    unsigned int cidin = cycle_id_in[a] ;
                    
                    if (edge_type[a] == edge_types::NODE) {
                        path_length.first_exon = path_length.middle;
                        path_length.middle = 0;
                    }
                    
                    // if in-degree is > 1 we can not yet delete nodes and arcs, only update them
                    bool del = in_deg[next] <= 1;
                    
                    ListDigraph::Arc arc_tmp(a);
                    ++a;
                    g.erase(arc_tmp);
                                       
                    // follow contraction
                    follow_contraction(g, flow, edge_specifier, edge_type, edge_lengths, means, in_deg, out_deg, cycle_id_in, cycle_id_out, path_flow, path_ident, path_length, path_mean, cidin, node, next, del);
                            
                } else {
                    ++a;
                }
            }
            
        }
    }
    
    return changed;
}


void base_manager::follow_contraction(ListDigraph &wc,
        ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::ArcMap<capacity_mean> &mc,
        InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg,
        ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
        capacity_type path_flow, exon_edge path_ident, edge_length path_length, capacity_mean path_mean,
        unsigned int cidin,
        ListDigraph::Node first, ListDigraph::Node current, bool del) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Follow Contraction " + path_ident.to_string() + " " + std::to_string(del) + "; "+  std::to_string(wc.id(first)) + " " + std::to_string(wc.id(current)) +  "\n");
    #endif
    
    for (ListDigraph::OutArcIt a(wc, current); a!=INVALID; ) {
        
        ListDigraph::Node next = wc.target(a);
    
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("iterator " + ces[a].to_string() + " " +  std::to_string(wc.id(a)) +  " ; "+ std::to_string(wc.id(next)) + "\n");
            logger::Instance()->debug("means " + std::to_string(path_mean.mean) + " " +  std::to_string(path_mean.weight) +  " ; "+ std::to_string(fc[a]) + " " + std::to_string(path_flow)+ "\n");
        #endif
        
        capacity_mean next_mean = path_mean;
        if (fc[a] >= path_flow) { // this means we first enter with a second incoming edge !
            capacity_mean new_mean = mc[a];
            new_mean.reduce(path_flow/(float) fc[a]);
            next_mean.update(new_mean); 
            
        } else { //one to many
            next_mean.reduce(fc[a]/(float) path_flow);            
            next_mean.update(mc[a]); 
        }

        capacity_type up_path_flow = path_flow;
        if (fc[a] < path_flow) {
            up_path_flow = fc[a];
        }
          
        exon_edge path_ident_c = path_ident;
        if (!ces[a].id.empty()) {
            if (path_ident.id.empty()) {
                path_ident_c = ces[a];
            } else {
                path_ident_c.join_edge(ces[a]);
            }
        } else {
            if (path_ident.id.empty() && path_ident.node_index == -1) {
                path_ident_c = ces[a];
            }
        }
        
        
        edge_length up_path_length;
        if (cet[a] == edge_types::EXON) {
            up_path_length = cel[a];
            up_path_length.middle += path_length.middle + path_length.last_exon;
            up_path_length.first_exon = path_length.first_exon;
        } else {
            up_path_length = path_length;
            if (up_path_length.first_exon == 0) {
                up_path_length.first_exon = cel[a].middle;
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
            wc.erase(arc_tmp);
            if (out_deg[current] == 0 ) {
                wc.erase(current);
            }
            
            follow_contraction(wc, fc, ces, cet, cel, mc, in_deg, out_deg, cycle_id_in, cycle_id_out, up_path_flow, path_ident_c, up_path_length, next_mean, cidin, first, next, up_del);
            
        } else if ( !del &&  out_deg[next] == 1 ) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Next 2 \n");
            #endif
            
            // can continue to merge this one :)
            // however, another contraction blocks this one
            
            ++a;
            follow_contraction(wc, fc, ces, cet, cel, mc, in_deg, out_deg, cycle_id_in, cycle_id_out, up_path_flow, path_ident_c, up_path_length, next_mean, cidin, first, next, del);
            
        } else {
 
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Next 3 Arc: " + path_ident_c.to_string() + " " + std::to_string(wc.id(first)) + " " + std::to_string(wc.id(next)) + "\n");
            #endif
            
            // I won’t give up, no I won’t give in
            // Till I reach the end
            // And then I’ll start again
            
            // we have found the end, create new edges if we deleted to up here and have fun :)
            
            ListDigraph::Arc arc_tmp(a);
                ++a;
            
            ListDigraph::Arc narc = wc.addArc(first, next);
            fc[narc] = up_path_flow;
            ces[narc] = path_ident_c;
            if (path_ident_c.id.empty()) {
                cet[narc] = edge_types::HELPER; // in rare cases helper are only contracted into the node 
            } else {
                cet[narc] = edge_types::EXON;  
            }
            cel[narc] = up_path_length;
            cycle_id_in[narc] = cidin;
            cycle_id_out[narc] = cycle_id_out[arc_tmp];
            mc[narc] = next_mean;
            
            if (del) {    
                 wc.erase(arc_tmp);
                 
                 if (a==INVALID) {
                     wc.erase(current);
                 }
            }
            
        }    
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Exit \n");

        digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("flow", flow)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .node("source", s)
                .node("drain", t)
                .run();  
    #endif
    
}

bool base_manager::simplify_ambiguous(
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap<unsecurity_id> &unsecurityId,
        ListDigraph::NodeMap<bool> &resolved,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg) {
    
    return false;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Simplify ambiguous \n");
    #endif
    
//    bool changed = false;
//    
//    for (ListDigraph::NodeIt n(g); n != INVALID;) {
//        // we can erase while looping, but be careful
//        
//        
//        #ifdef ALLOW_DEBUG
//        logger::Instance()->debug("Process Node " + std::to_string(g.id(n)) + "\n");
//        #endif
//        
//        ListDigraph::Node node(n);
//                
//        if ( out_deg[node] > 1 && in_deg[node] > 1 && !resolved[node]) {
//            // this is an unresolvable node
//            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Unresolved Node " + std::to_string(g.id(n)) + "\n");
//            #endif
//            
//            changed = true;
//            
//            std::deque<lpath> pre_path;
//            std::deque<rpath> post_path;
//            
//            bool helper = false;
//            
//            // we init pre and post with directly adjacent edges only
//            for (ListDigraph::InArcIt a(g, node); a!=INVALID; ++a) {
//                if (edge_type[a] == edge_types::HELPER) {
//                    // helper can't possibly be disproven/proven
//                    helper = true;
//                    continue;
//                }
//                pre_path.push_back( lpath(&edge_specifier[a], node_index[g.source(a)], a, g.source(a)) );
//            }
//            
//            for (ListDigraph::OutArcIt a(g, node); a!=INVALID; ++a) {
//                if (edge_type[a] == edge_types::HELPER) {
//                    // helper can't possibly be disproven/proven
//                    helper = true;
//                    continue;
//                }
//                post_path.push_back( rpath(&edge_specifier[a], node_index[g.target(a)], a, g.target(a), edge_lengths[a]) );
//            }
//            
//            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("l-Path: ");
//            for ( std::deque<lpath>::iterator it = pre_path.begin(); it!= pre_path.end(); ++it) {
//                logger::Instance()->debug("(" + it->identifier.to_string()+")" );
//            }
//            logger::Instance()->debug("\n");
//            logger::Instance()->debug("r-Path: ");
//            for ( std::deque<rpath>::iterator it = post_path.begin(); it!= post_path.end(); ++it) {
//                logger::Instance()->debug("(" + it->identifier.to_string()+")" );
//            }
//            logger::Instance()->debug("\n");
//            #endif
//            
//            // we need to get the id of the exon we bridge here
//            ListDigraph::OutArcIt a(g, node);
//            while (edge_type[a] == edge_types::HELPER) {
//                ++a;
//            }
//            unsigned int exon_index = edge_specifier[a].id.find_first();
//            
//                        exon_edge last;
//            rpos min_right_hit = 0;
//            rcount evidence_count = 0;
//            rcount threshold_count = post_path.size();// pre_path.size() + post_path.size();
//            bool matches_all_left = true;
//            
//            // we are only interested in the direct arcs
//            // all paths going further cannot be resolved, because more path through these nodes might exist
//            std::deque<lpath* > pre_path_match;
//            path_evidence_set<int> post_path_match;
//            path_evidence_map<int, rpos> post_path_min_length;
//            path_evidence_map<int, rcount> post_path_count;
//            
//            // +++++++++++++++++++++  Overarching  +++++++++++++++++++++
//            // we start with fragments overarching this node before dealing with paired data
//            for (graph_list<std::pair<exon_group *, rcount>  >::iterator si = raw->single_for_exon[exon_index].begin(); si != raw->single_for_exon[exon_index].end() ;++si) {
//                
//                exon_edge left, right;
//                si->first->bin_mask.left_split(exon_index, left);
//                si->first->bin_mask.right_split(exon_index, right);
//                
//                rcount count = si->first->frag_count + si->second; 
//                
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Test Overarch: " + left.to_string() + " " + right.to_string() + " " + std::to_string(count) + "\n");
//                #endif
//                
//                unsigned int leftmost = si->first->range_start;
//                unsigned int rightmost = si->first->range_end;
//                
//                if (last != left) {
//                    // finalize last set :)
//                    
//                    if (!last.id.empty()) {
//
//                        // so now we have all matches collected
//                        // we need to put in evidences
//                        
//                        if (evidence_count > threshold_count && !matches_all_left) resolve_overarching(pre_path_match, post_path_match, post_path_min_length, post_path_count, min_right_hit, g, edge_type, node);
//
//                        pre_path_match.clear();
//                        post_path_match.clear();
//                        post_path_count.clear();
//                        min_right_hit = 0;
//                        evidence_count = 0;
//                        matches_all_left = true;
//                        // we keep paths, and therefore the min_length_list
//                    }
//                
//                    // go over all left paths first
//                    for ( unsigned i = 0; i < pre_path.size(); ++i) {
//
//                        lpath *pp = &pre_path[i];
//                        if (pp->is_parent) { // only work on maximal paths
//                            continue;
//                        }
//                        
//                        // can this match? path long enough?
//                        if (pp->border_index > leftmost) {
//                            // extendpath
//                            extend_path_left(pre_path, pp, leftmost, g, edge_specifier, edge_type, node_index, know_paths);
//                            continue; // so this was made a parent or can't hit when we ran out of paths
//                        }
//                        
//                        if (left.is_contained_in(pp->identifier, leftmost, exon_index)) {
//                            pre_path_match.push_back( pp );
//                        } else {
//                            matches_all_left = false;
//                        }
//                    }
//                    last = left;
//                }
//                
//                // we have no distance in the middle here, so we can directly look at exons that should match on right
//                rpos len = 0;
//                boost::dynamic_bitset<>::size_type i = right.id.find_first();
//                while(i < meta->size ) {
//                    len += meta->exons[i].exon_length;
//                    i = right.id.find_next(i);
//                }
//                if (min_right_hit == 0 || min_right_hit > len ) {
//                    min_right_hit = len;
//                }
//                
//                // we keep all counts for later variables this time!
//                resolve_count_set next;
//                next.count = count;
//                for (std::deque<lpath* >::iterator li = pre_path_match.begin(); li != pre_path_match.end(); ++li) {
//                    next.left.insert( g.id((*li)->starting_arc) );
//                }
//                
//                // then we go over right paths
//                bool matches_all_right = true;
//                std::deque<int> right_ids;
//                for ( rpos i = 0; i < post_path.size(); ++i ) {
//                    
//                    rpath *pp = &post_path[i];
//                    
//                    // can this match? path long enough?
//                    if (pp->border_index < rightmost) {
//                        // extendpath
//                        extend_path_right(post_path, pp, rightmost, g, edge_specifier, node_index, edge_lengths, edge_type, post_path_min_length, know_paths);
//                    }
//                    
//                    if (right.is_contained_in(pp->identifier, exon_index, rightmost)) {
//                        int id = g.id(pp->starting_arc);
//                        right_ids.push_back(id);
//                    } else {
//                        matches_all_right = false;
//                    }
//                }
//                
//                if (!matches_all_right) {
//                    for (std::deque<int>::iterator id_it = right_ids.begin(); id_it != right_ids.end(); ++id_it) {
//                        int id = *id_it;
//                        post_path_match.insert(id);
//                        post_path_count[id] += count;
//                        next.right.insert(id);  
//                        evidence_count += count;
//                    }
//                }
//                     
//                // now we have everything we need to extend evidences,
//                // try if next has same left arcs !
//            }
//            
//            // fit in last list and clean evidence lists;
//            if (evidence_count > threshold_count && !matches_all_left) resolve_overarching(pre_path_match, post_path_match, post_path_min_length, post_path_count, min_right_hit, g, edge_type, node);
//            pre_path_match.clear();
//            post_path_match.clear();
//            post_path_count.clear();
//            matches_all_left = true;
//            evidence_count = 0;
//            
//
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Overarching done \n");
//            for ( std::deque<lpath>::iterator  pp = pre_path.begin(); pp!= pre_path.end(); ++pp) {
//                logger::Instance()->debug("(" + pp->identifier.to_string() +")" + std::to_string(pp->is_parent) + "\n" );
//                for (path_evidence_map<int, rcount>::iterator ei = pp->evidence.begin(); ei!= pp->evidence.end(); ++ei) {
//                    logger::Instance()->debug("Linked " + std::to_string(ei->first) +" " + std::to_string(pp->evidence.is_blocked(ei->first)) + "\n" );
//                }
//            }
//            #endif
//            
//            
//            
//            // +++++++++++++++++++++  Paired  +++++++++++++++++++++
//            // in the second step we apply paired end data
//                                                       //min read    , max
//            path_evidence_map<int, std::pair<std::pair< rpos, rpos > ,rpos> > post_path_range;
//            std::deque<rpos> hit_distance;
//            exon_group* lastp = NULL;
//            for (graph_list<paired_exon_group *>::iterator pi = raw->pairs_for_exon[exon_index].begin(); pi != raw->pairs_for_exon[exon_index].end() ; ++pi) {
//                
//                unsigned int leftmost = (*pi)->left_read->range_start;
//                unsigned int rightmost = (*pi)->right_read->range_end;
//                                
//                rcount count = (*pi)->count; 
//                
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Test Paired\n");
//                logger::Instance()->debug("Test Paired: " + (*pi)->left_read->bin_mask.to_string() + " " + (*pi)->right_read->bin_mask.to_string() + " " + std::to_string(count) + "\n");
//                #endif
//                
//                if (lastp == NULL || lastp != (*pi)->left_read ) {
//                    // finalize last set :)
//                    
//                    if (lastp != NULL) {
//             
//                        // so now we have all matches collected
//                        // we need to put in evidences
//                        if (evidence_count > threshold_count && !matches_all_left) resolve_pair(pre_path_match, hit_distance, post_path_match, post_path_min_length, post_path_range, post_path_count, g, edge_type, node);
//
//                        pre_path_match.clear();
//                        post_path_match.clear();
//                        post_path_range.clear();
//                        post_path_count.clear();
//                        min_right_hit = 0;
//                        evidence_count = 0;
//                        matches_all_left = true;
//                        // we keep paths, and therefore the min_length_list
//                    }
//                
//                    // go over all left paths first
//                    for ( unsigned i = 0; i < pre_path.size(); ++i) {
//
//                        lpath *pp = &pre_path[i];
//                        
//                        if (pp->is_parent) { // only work on maximal paths
//                            continue;
//                        }
//                        
//                        // can this match? path long enough?
//                        if (pp->border_index > leftmost) {
//                            // extendpath
//                            extend_path_left(pre_path, pp, leftmost, g, edge_specifier, edge_type, node_index, know_paths);
//                            continue; // so this was made a parent or can't hit when we ran out of paths
//                        }
//                        
//                        if ((*pi)->left_read->bin_mask.is_contained_in(pp->identifier, leftmost, (*pi)->left_read->range_end )) {
//                            pre_path_match.push_back( pp );
//                            
//                            // we have a winner, get hit distance for filtering
//                            rpos dist = 0;
//                            boost::dynamic_bitset<>::size_type  i =  pp->identifier.id.find_next((*pi)->left_read->range_end);
//                            while(i < exon_index) { // EXcluding middle node
//                                dist += meta->exons[i].exon_length;
//                                i = pp->identifier.id.find_next(i);
//                            }
//                            
//                            // the distance to the middle!
//                            hit_distance.push_back(dist);
//                            
//                        } else {
//                            matches_all_left = false;
//                        }
//                    }
//                    lastp = (*pi)->left_read;
//                }
//                
//                
//                // the length of the right hit itself
//                rpos len = 0;
//                boost::dynamic_bitset<>::size_type  i =  (*pi)->right_read->bin_mask.id.find_first();
//                while(i < meta->size ) {
//                    len += meta->exons[i].exon_length;
//                    i = (*pi)->right_read->bin_mask.id.find_next(i);
//                }
//                
//                resolve_count_set next;
//                next.count = count;
//                for (std::deque<lpath* >::iterator li = pre_path_match.begin(); li != pre_path_match.end(); ++li) {
//                    next.left.insert( g.id((*li)->starting_arc) );
//                }
//                
//                // then we go over right paths
//                bool matches_all_right = true;
//                std::deque<int> right_ids;
//                std::deque<int> js;
//                for ( unsigned j = 0; j < post_path.size(); ++j ) {
//                    
//                    rpath *pp = &post_path[j];
//                    
//                    // can this match? path long enough?
//                    if (pp->border_index < rightmost) {
//                        // extendpath
//                        extend_path_right(post_path, pp, rightmost, g, edge_specifier, node_index, edge_lengths, edge_type, post_path_min_length, know_paths);
//                    }
//                    
//                    if ((*pi)->right_read->bin_mask.is_contained_in(pp->identifier, (*pi)->right_read->range_start, rightmost)) {
//                        int id = g.id(pp->starting_arc);
//                        right_ids.push_back(id);
//                        js.push_back(j);
//                    } else {
//                        matches_all_right = false;
//                    }
//                }
//                
//                if (!matches_all_right) {
//                    std::deque<int>::iterator js_it = js.begin();
//                    for (std::deque<int>::iterator id_it = right_ids.begin(); id_it != right_ids.end(); ++id_it, ++js_it) {
//                        
//                        int id = *id_it;
//                        rpath *pp = &post_path[*js_it];
//                         
//                        post_path_match.insert(id);
//                        next.right.insert(id);
//                        
//                        // so now we need to mark the hit distances
//                        boost::dynamic_bitset<>::size_type  i = std::max(exon_index, (*pi)->left_read->range_end+1);
//                        rpos dist = 0;      // INcluding middle node
//                        while( i < (*pi)->right_read->range_start ) { 
//                            dist += meta->exons[i].exon_length;
//                            i = pp->identifier.id.find_next(i);
//                        }
//                        
//                        path_evidence_map<int, std::pair< std::pair< rpos, rpos > ,rpos> >::iterator mark = post_path_range.find(id);
//                        if (mark == post_path_range.end()) {
//                            post_path_range[id] = std::make_pair(std::make_pair(dist, len),dist);
//                        } else {
//                            if (mark->second.first.first > dist) {
//                                mark->second.first.first = dist;
//                                mark->second.first.second = len;
//                            } else if (mark->second.second < dist) {
//                                mark->second.second = dist;
//                            }
//                        }
//                        
//                        post_path_count[id] += count;
//                        
//                        evidence_count += count;
//                    }
//                }
//                
//                // now we have everything we need to extend evidences,
//                // try if next has same left arcs !
//            }
//            // fit in last list and clean evidence lists;
//            if (evidence_count > threshold_count && !matches_all_left) resolve_pair(pre_path_match, hit_distance, post_path_match, post_path_min_length, post_path_range, post_path_count, g, edge_type, node);
//            pre_path_match.clear();
//            post_path_range.clear();
//            post_path_match.clear();
//            post_path_count.clear();
//            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Paired done \n");
//            for ( std::deque<lpath>::iterator  pp = pre_path.begin(); pp!= pre_path.end(); ++pp) {
//                logger::Instance()->debug("(" + pp->identifier.to_string() +")" + std::to_string(pp->is_parent) + "\n" );
//                for (path_evidence_map<int, rcount>::iterator ei = pp->evidence.begin(); ei!= pp->evidence.end(); ++ei) {
//                    logger::Instance()->debug("Linked " + std::to_string(ei->first) +" " + std::to_string(pp->evidence.is_blocked(ei->first)) + "\n" );
//                }
//            }
//            #endif
//            
//            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("l-Path-extended: ");
//            for ( std::deque<lpath>::iterator it = pre_path.begin(); it!= pre_path.end(); ++it) {
//                logger::Instance()->debug("(" + it->identifier.to_string() + " " + std::to_string(it->is_parent) +")\n" );
//            }
//            logger::Instance()->debug("r-Path-extended: ");
//            for ( std::deque<rpath>::iterator it = post_path.begin(); it!= post_path.end(); ++it) {
//                logger::Instance()->debug("(" + it->identifier.to_string() +")\n" );
//            }
//            #endif
//            
//            // +++++++++++++++++++++  Add in Info  +++++++++++++++++++++
//            
//            // now resolve what can be resolved
//            // we use know_paths and know_back_paths
//            
//            // we iterate over all possible paths and now mark the found on direct arcs
//            bool has_evidence = false;
//            for ( std::deque<lpath>::iterator pp = pre_path.begin(); pp != pre_path.end(); ++pp) {
//
//                if (pp->is_parent) { // we only look at the extended paths
//                    continue;
//                }
//                
//                for (path_evidence_map<int, rcount>::iterator ei = pp->evidence.begin(); ei!= pp->evidence.end(); ++ei) {
//                    if (!pp->evidence.is_blocked(ei->first)) {
//                        know_paths[pp->starting_arc][ei->first] = ei->second;
//                        know_back_paths[g.arcFromId(ei->first)].bridges->insert(g.id(pp->starting_arc));
//                        has_evidence = true;
//                    }
//                }
//            }
//            
//            unsecurityId[n].resolvable = has_evidence;
//            if (!has_evidence) {
//                ++n;
//                resolved[node] = true;
//                continue;
//            }
//            
//            
//            // +++++++++++++++++++++  Cycle  +++++++++++++++++++++
//            
//            // this can handle newly added edges that resolve cycles in the graph
//            // as they are added as helpers, they are not in the enumerated Paths
//            
//            for (ListDigraph::InArcIt a(g, node); a!=INVALID; ++a) {
//                 
//                if (edge_type[a] == edge_types::HELPER && cycle_id_in[a] != 0) { // this is a helper added by
//                    
//                    for (graph_list<paired_exon_group>::iterator circle = raw->chim_circle_bin_list.begin(); circle != raw->chim_circle_bin_list.end(); ++circle) {
//                        
//                        // unintuitively, the in helper edge connects to the OUT 
//                        
//                        // bounds
//                        unsigned int out_left = circle->right_read->range_start; // jump
//                        unsigned int out_right = circle->right_read->range_end;
//
//                        if ( exon_index != out_left) {
//                            continue;
//                        }  
//                     
//                        // if we are here there should exist at least one match to out_path
//                        // as all fragments are as well in single bin, paths are guaranteed long enough!
//                        // we go over right paths to see which matches the evidence
//                        for ( std::deque<rpath>::iterator pp = post_path.begin(); pp != post_path.end(); ++pp ) {
//                            
//                            if (out_right > pp->border_index) { // paths NOT containing the fragment might be shorter, but we don't care for those
//                                continue;
//                            }
//                            
//                            if (circle->right_read->bin_mask.is_contained_in(pp->identifier, out_left, out_right) ) {
//                                // we have a correct hit, add the known info directly
//                                int idb = g.id(pp->starting_arc);
//                                know_paths[a][idb] += 1; // evidence of one means just found, arbitrary number, could also be zero
//                                know_back_paths[pp->starting_arc].bridges->insert(g.id(a));
//                            }
//                            
//                        }
//                    }
//                }
//            }
//            // now we do the same for outward arcs
//            for (ListDigraph::OutArcIt a(g, node); a!=INVALID; ++a) {
//                 
//                if (edge_type[a] == edge_types::HELPER && cycle_id_out[a] != 0) { // this is a helper added by
//                    
//                    for (graph_list<paired_exon_group>::iterator circle = raw->chim_circle_bin_list.begin(); circle != raw->chim_circle_bin_list.end(); ++circle) {
//                        
//                        // unintuitively, the out helper edge connects to the IN 
//                        
//                        // bounds
//                        unsigned int in_right = circle->left_read->range_end; // jump
//                        unsigned int in_left = circle->left_read->range_start;
//
//                        if ( exon_index != in_right) {
//                            continue;
//                        }  
//                     
//                        // if we are here there should exist at least one match to out_path
//                        // as all fragments are as well in single bin, paths are guaranteed long enough!
//                        // we go over left paths to see which matches the evidence
//                        for ( std::deque<lpath>::iterator pp = pre_path.begin(); pp != pre_path.end(); ++pp ) {
//                            
//                            if (in_left < pp->border_index) { // paths NOT containing the fragment might be shorter, but we don't care for those
//                                continue;
//                            }
//                            
//                            if (circle->left_read->bin_mask.is_contained_in(pp->identifier, in_left, in_right) ) {
//                                // we have a correct hit, add the known info directly
//                                int idb = g.id(a);
//                                know_paths[pp->starting_arc][idb] += 1; // evidence of one means just found, arbitrary number, could also be zero
//                                know_back_paths[a].bridges->insert(g.id(pp->starting_arc));
//                            }
//                            
//                        }
//                    }
//                }
//            }
//            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Cycle done \n");
//            #endif
//            
//            // +++++++++++++++++++++  Empty Evidences and Resolve  +++++++++++++++++++++
//            
////            
////              #ifdef ALLOW_DEBUG
////            logger::Instance()->debug("Pre empty " + std::to_string(g.id(n)) + "\n");
////            for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
////
////               logger::Instance()->debug("Arc " + std::to_string(g.id(a))); 
////               for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
////                   logger::Instance()->debug(" f " + std::to_string(ab->first)); 
////               }
////               for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
////                   logger::Instance()->debug(" b " + std::to_string(*ab)); 
////               }
////               logger::Instance()->debug("\n" );  
////            }
////            #endif
////            
////            // TODO: this should not be needed
//////            // if we have arcs with no info, we don't know what happened
////            for (ListDigraph::InArcIt a(g, node); a!=INVALID; ++a) {
////                int id = g.id(a);
////                logger::Instance()->debug("A " + std::to_string(id) + "\n");
////                if (edge_type[a] != edge_types::HELPER && !know_paths[a].is_evidenced_path()) { 
////                    
////                    ListDigraph::Arc max;
////                    rcount max_value = 0;
////                    for (ListDigraph::OutArcIt b(g, node); b!=INVALID; ++b) {
////                        
////                        if (edge_type[b] == edge_types::HELPER) {
////                            // no direct join of source and sinks
////                            continue;
////                        }
////                         logger::Instance()->debug("B " + std::to_string(g.id(b)) + " " + std::to_string(flow[b]) + "\n");
////                        if (flow[b] > max_value) {
////                            logger::Instance()->debug("setmax\n");
////                            max_value = capacity[b];
////                            max = b;
////                        }
////                        
////                    }
////                    
////                    int idb = g.id(max);
////                    
////                    logger::Instance()->debug("IDB1 " + std::to_string(idb) + " " + std::to_string(max_value) + "\n");
////                    
////                    know_paths[a][idb] += 1; // evidence of one means just found, arbitrary number, could also be zero
////                    know_back_paths[max].bridges->insert(id);
////                }
////            }    
////            for (ListDigraph::OutArcIt a(g, node); a!=INVALID; ++a) {  
////                int id = g.id(a);
////                if (edge_type[a] != edge_types::HELPER && !know_back_paths[a].is_evidenced_path()) { 
////                    
////                    ListDigraph::Arc max;
////                    rcount max_value = 0;
////                    for (ListDigraph::InArcIt b(g, node); b!=INVALID; ++b) {
////                        
////                        if (edge_type[b] == edge_types::HELPER) {
////                            // no direct join of source and sinks
////                            continue;
////                        }
////                         
////                        if (flow[b] > max_value) {
////                            max_value = capacity[b];
////                            max = b;
////                        }
////                        
////                    }
////
////                    int idb = g.id(max);
////                    
////                    logger::Instance()->debug("IDB2 " + std::to_string(idb) + " " + std::to_string(max_value) + "\n");
////                    
////                    know_paths[max][id] += 1; // evidence of one means just found, arbitrary number, could also be zero
////                    know_back_paths[a].bridges->insert(idb);
////                }
////                
////            }
////            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Pre empty " + std::to_string(g.id(n)) + "\n");
//            for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//
//               logger::Instance()->debug("Arc " + std::to_string(g.id(a))); 
//               for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
//                   logger::Instance()->debug(" f " + std::to_string(ab->first)); 
//               }
//               for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
//                   logger::Instance()->debug(" b " + std::to_string(*ab)); 
//               }
//               logger::Instance()->debug("\n" );  
//            }
//            #endif
//            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Empty Evidences Added \n");
//            #endif
//
////            if (helper) {
////                resolved[node] = true;
////                ++n;
////                continue;
////            }
//            
//            // now we can try and extract nodes that have unambiguous info
//            ++n;
//            unravel_evidences(node, know_paths, know_back_paths, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, &resolved);
//            resolved[node] = true;
//                        
//        } else {
//            // resolved node (source, drain)
//             ++n;
//        }
//        
//        
//            // logging stuff to remove
////        logger::Instance()->debug("Node Done " + std::to_string(g.id(n)) + "\n");
////        for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
////            
////           logger::Instance()->debug("Arc " + std::to_string(g.id(a))); 
////           for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
////               logger::Instance()->debug(" f " + std::to_string(ab->first)); 
////           }
////           for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
////               logger::Instance()->debug(" b " + std::to_string(*ab)); 
////           }
////           logger::Instance()->debug("\n" );  
////        }
//        
//        // end logging stuff
//        
//        #ifdef ALLOW_DEBUG
//        logger::Instance()->debug("Step done\n");
//        if (n != INVALID) print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, n, n);
//        #endif
//        
//    }
//    
//    return changed;
}


bool base_manager::simplify_ambiguous_ILP(
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap<unsecurity_id> &unsecurityId,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        InDegMap<ListDigraph> &in_deg, OutDegMap<ListDigraph> &out_deg, 
        
        ListDigraph::ArcMap<bool> &barred, bool final_resolve) {
        
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Simplify ambiguous ILP \n");
    #endif
     
    for (ListDigraph::ArcIt a(g); a != INVALID;++a) {
        means[a].reduce_score();
    }
    
    ListDigraph::ArcMap<arc_bridge> temp_know_paths(g);
    ListDigraph::ArcMap<arc_back_bridge> temp_know_back_paths(g);
    ListDigraph::ArcMap<arc_bridge> temp_know_paths_add_max(g);
    ListDigraph::ArcMap<arc_back_bridge> temp_know_back_paths_add_max(g);
    ListDigraph::NodeMap<classify_component> classify(g);
    ListDigraph::NodeMap<std::unordered_set<int> > single_to_single(g);
    ListDigraph::NodeMap<std::unordered_set<int> > block_delete_high_res(g);
    ListDigraph::NodeMap<std::unordered_set<int> > block_delete(g);
    
    ListDigraph::NodeMap<std::map<int, evidence_group> > left_groups_ids(g);
    ListDigraph::NodeMap<std::map<int, evidence_group> > right_groups_ids(g);
    
    // ######### Classify Each Node ######### //
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Classify ##################### \n");
    #endif
  
    for (ListDigraph::NodeIt n(g); n != INVALID;++n) {
        // we can erase while looping, but be careful

        ListDigraph::Node node(n);       
        
        // out_deg and in_deg unfunction because lemon are fucking stupip fucks
        
        unsigned int pre_size = 0;
        unsigned int post_size = 0;
        
        for (ListDigraph::InArcIt a(g, node); a!=INVALID; ++a) {
                if (barred[a] && !final_resolve) { // do NOT include barred edges unless for final resolutions!
                    continue;
                }
                ++pre_size;
            }
            
            for (ListDigraph::OutArcIt a(g, node); a!=INVALID; ++a) {
                if (barred[a] && !final_resolve) { // do NOT include barred edges unless for final resolutions!
                    continue;
                }
                ++post_size;
            }
        
        if ( pre_size > 1 && post_size > 1) {
            // this is an unresolvable node
                        
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Classify Node " + std::to_string(g.id(n)) + "\n");
            #endif
                        
            std::deque<path> pre_path;
            std::deque<path> post_path;
            
            bool helper = false;
            
            // we init pre and post with directly adjacent edges only
            for (ListDigraph::InArcIt a(g, node); a!=INVALID; ++a) {
                if (edge_type[a] == edge_types::HELPER) {
                    // helper can't possibly be disproven/proven
                    helper = true;
                    continue;
                } else if (barred[a] && !final_resolve  || edge_type[a] == edge_types::RESOLVE_HELPER) { // do NOT include barred edges unless for final resolutions!
                    continue;
                }
                pre_path.push_back( path(&edge_specifier[a], node_index[g.source(a)], a, g.source(a)) );
            }
            
            for (ListDigraph::OutArcIt a(g, node); a!=INVALID; ++a) {
                if (edge_type[a] == edge_types::HELPER) {
                    // helper can't possibly be disproven/proven
                    helper = true;
                    continue;
                } else if (barred[a] && !final_resolve || edge_type[a] == edge_types::RESOLVE_HELPER) { // do NOT include barred edges unless for final resolutions!
                    continue;
                }
                post_path.push_back( path(&edge_specifier[a], node_index[g.target(a)], a, g.target(a)) );
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
            ListDigraph::OutArcIt a(g, node);
            while (a != INVALID && (edge_type[a] == edge_types::HELPER || edge_type[a] == edge_types::RESOLVE_HELPER)) {
                ++a;
            }
            ei = a;
            if (a==INVALID) {
                ListDigraph::InArcIt a(g, node);
                while (edge_type[a] == edge_types::HELPER || edge_type[a] == edge_types::RESOLVE_HELPER) {
                    ++a;
                }
                ei = a;
            }
            
            unsigned int exon_index = edge_specifier[ei].id.find_first();
            
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
            for (graph_list<std::pair<exon_group *, rcount>  >::iterator si = raw->single_for_exon[exon_index].begin(); si != raw->single_for_exon[exon_index].end() ;++si) {
                
                exon_edge left, right;
                si->first->bin_mask.left_split(exon_index, left);
                si->first->bin_mask.right_split(exon_index, right);
                
                rcount count = si->first->frag_count + si->second; 
                
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
                            extend_path_left(pre_path, pp, leftmost, g, edge_specifier, edge_type, node_index, know_paths);
                        }
                        
                        if (left.is_contained_in(pp->identifier, leftmost, exon_index)) {
                            pre_path_match.insert( g.id(pp->starting_arc) );
//                            #ifdef ALLOW_DEBUG
//                            logger::Instance()->debug("Match L: " + std::to_string(g.id(pp->starting_arc)) + "\n");
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
                        extend_path_right(post_path, pp, rightmost, g, edge_specifier, edge_type, node_index, know_paths);
                    }
                    
                    if (right.is_contained_in(pp->identifier, exon_index, rightmost)) {
                        post_path_match.insert(g.id(pp->starting_arc));
//                        #ifdef ALLOW_DEBUG
//                        logger::Instance()->debug("Match R: " + std::to_string(g.id(pp->starting_arc)) + "\n");
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
                
                rcount count = (*pi)->count; 
                
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
                            extend_path_left(pre_path, pp, leftmost, g, edge_specifier, edge_type, node_index, know_paths);
                        }
                        
                        if (left.is_contained_in(pp->identifier, leftmost, exon_index)) {
                            pre_path_match.insert( g.id(pp->starting_arc) );
//                            #ifdef ALLOW_DEBUG
//                            logger::Instance()->debug("Match L: " + std::to_string(g.id(pp->starting_arc)) + "\n");
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
                        extend_path_right(post_path, pp, rightmost, g, edge_specifier, edge_type, node_index, know_paths);
                    }
                    
                    if (right.is_contained_in(pp->identifier, exon_index, rightmost)) {
                        post_path_match.insert(g.id(pp->starting_arc));
//                        #ifdef ALLOW_DEBUG
//                        logger::Instance()->debug("Match R: " + std::to_string(g.id(pp->starting_arc)) + "\n");
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
                                
                rcount count = (*pi)->count; 
                
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
                            extend_path_left(pre_path, pp, leftmost, g, edge_specifier, edge_type, node_index, know_paths);
                        }
                        
                        if ( left_a.is_contained_in(pp->identifier, la1, la2 ) && (ra1 == 0 || right_a.is_contained_in(pp->identifier, ra1, ra2)) ) {
                            pre_path_match.insert( g.id(pp->starting_arc) );
//                            #ifdef ALLOW_DEBUG
//                            logger::Instance()->debug("Match L: " + std::to_string(g.id(pp->starting_arc)) + "\n");
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
                        extend_path_right(post_path, pp, rightmost, g, edge_specifier, edge_type, node_index, know_paths);
                    }
                    
                    if ( right_b.is_contained_in(pp->identifier, rb1, rb2) && (lb1 == 0 || left_b.is_contained_in(pp->identifier, lb1, lb2 )) ) {
                        post_path_match.insert(g.id(pp->starting_arc));
//                        #ifdef ALLOW_DEBUG
//                        logger::Instance()->debug("Match R " + std::to_string(g.id(pp->starting_arc)) + "\n");
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
                        rpos pos = std::max(node_index[g.source(g.arcFromId(*first))], node_index[g.source(g.arcFromId(*pi))]);
                        if (!edge_specifier[g.arcFromId(*first)].is_contained_in(edge_specifier[g.arcFromId(*pi)], std::max(la1, pos), exon_index)) {
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
                        rpos pos = std::max(node_index[g.target(g.arcFromId(*first))], node_index[g.target(g.arcFromId(*pi))]);
                        if (!edge_specifier[g.arcFromId(*first)].is_contained_in(edge_specifier[g.arcFromId(*pi)], std::min(pos, rb2), exon_index)) {
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
                    for (ListDigraph::OutArcIt a(g, n); a!=INVALID; ++a) {
                        if (edge_type[a] == edge_types::HELPER) {
                            continue;  // helper can't possibly be disproven/proven
                        } 
                        int id = g.id(a);
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
//            for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
//                int pp = g.id(a);
//                for (path_evidence_map<int, rcount>::iterator ei = evidences[pp].begin(); ei!= evidences[pp].end(); ++ei) {
//                    if (!evidences[pp].is_blocked(ei->first)) {
//                        count_left[pp] += ei->second;
//                        count_right[ei->first] += ei->second;
//                    }
//                }
//            }
//            
            
            for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
                int pp = g.id(a);
                for (path_evidence_map<int, rcount>::iterator ei = evidences[pp].begin(); ei!= evidences[pp].end(); ++ei) {
                    if (!evidences[pp].is_blocked(ei->first) ) {//&& (ei->second * 100 / count_left[pp] >= 5 || ei->second * 100 / count_right[ei->first] >= 5) ) {
                        temp_know_paths[g.arcFromId(pp)][ei->first] = ei->second;
                        temp_know_back_paths[g.arcFromId(ei->first)].bridges->insert(pp);
                        
                        block_delete[n].insert(pp);
                        block_delete[n].insert(ei->first);
                    }
                }
            }
                 
            unsecurityId[n].resolvable = has_evidence;
            
            // we compute the groups first 
            std::deque<std::set<int> > left_groups;
            std::deque<std::set<int> > right_groups;
            compute_edge_groups(left_raw_groups, left_groups);
            compute_edge_groups(right_raw_groups, right_groups);
            
            for (std::deque<std::set<int> >::iterator lgi = left_groups.begin(); lgi != left_groups.end(); ++lgi) {
                std::set<int>::iterator s = lgi->begin();
                std::set<int>::iterator s2 = s;
                int max = *s;
                float max_score = means[g.arcFromId(*s)].compute_score();
                ++s;
                for(; s != lgi->end(); ++s) {
                    float score = means[g.arcFromId(*s)].compute_score();
                    if (score > max_score) {
                        max_score = score;
                        max = *s;
                    }
                }
//                for(; s2 != lgi->end(); ++s2) {
//                    float score = means[g.arcFromId(*s2)].compute_score();
//                    if ( max_score * 0.5 <= score && block_delete[n].find(*s2) != block_delete[n].end()) block_delete_high_res[n].insert(*s2);
//                }
                if (block_delete[n].find(max) != block_delete[n].end()) block_delete_high_res[n].insert(max);
            }
            
            for (std::deque<std::set<int> >::iterator lgi = right_groups.begin(); lgi != right_groups.end(); ++lgi) {
                std::set<int>::iterator s = lgi->begin();
                std::set<int>::iterator s2 = s;
                int max = *s;
                float max_score = means[g.arcFromId(*s)].compute_score();
                ++s;
                for(; s != lgi->end(); ++s) {
                    float score = means[g.arcFromId(*s)].compute_score();
                    if (score > max_score) {
                        max_score = score;
                        max = *s;
                    }
                }
//                for(; s2 != lgi->end(); ++s2) {
//                    float score = means[g.arcFromId(*s2)].compute_score();
//                    if ( max_score * 0.5 <= score && block_delete[n].find(*s2) != block_delete[n].end()) block_delete_high_res[n].insert(*s2);
//                }
                if (block_delete[n].find(max) != block_delete[n].end()) block_delete_high_res[n].insert(max);
            }
            
            int id_iterator = 0;
            for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
                left_groups_ids[n][g.id(a)] = evidence_group(id_iterator, false);
                ++id_iterator;
            }
            for (ListDigraph::OutArcIt a(g, n); a!=INVALID; ++a) {
                right_groups_ids[n][g.id(a)] = evidence_group(id_iterator, false);
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
//                        for(std::set<int>::iterator s2 = lgi->begin(); s2 != lgi->end(); ++s2) {
//                            for (arc_bridge::iterator ab = temp_know_paths[g.arcFromId(*s)].begin(); ab != temp_know_paths[g.arcFromId(*s)].end(); ++ab) {
//                                 temp_know_back_paths[g.arcFromId(ab->first)].bridges->insert(*s2);
//                                 temp_know_paths[g.arcFromId(*s2)][ab->first] += 1;
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
//                        for(std::set<int>::iterator s2 = rgi->begin(); s2 != rgi->end(); ++s2) {
//                            for (arc_back_bridge::iterator ab = temp_know_back_paths[g.arcFromId(*s)].begin(); ab != temp_know_back_paths[g.arcFromId(*s)].end(); ++ab) {
//                                 temp_know_back_paths[g.arcFromId(*s2)].bridges->insert(*ab);
//                                 temp_know_paths[g.arcFromId(*ab)][*s2] += 1;
//                            }
//                        }
                        //if (block_delete[n].find(*s) != block_delete[n].end()) block_delete_high_res[n].insert(*s);
                    }
                    ++id_iterator;
                    ++rgi;
                }
            }
            
            // identify the components
            for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
                
                if(barred[a]) {
                    continue;
                }
                
                bool matched = false;
                for (std::deque<component>::iterator ci = classify[n].components.begin(); ci !=  classify[n].components.end(); ++ci) {
                    if ( ci->nodes.find(g.id(a)) != ci->nodes.end()) {
                        matched = true;
                        break;
                    }
                }
                
                if (!matched) {
                    classify[n].components.push_back(component());
                    find_component_left( classify[n].components.back().nodes, g.id(a), g, temp_know_paths, temp_know_back_paths);
                }
            }
            // now set up the right numbers!
            for (std::deque<component >::iterator ci = classify[n].components.begin(); ci !=  classify[n].components.end(); ) {
                        
                if (ci->nodes.size() == 1) {
//                     #ifdef ALLOW_DEBUG
//                    logger::Instance()->debug("CI Single"+ std::to_string(g.id(n)) + " as " + std::to_string(*ci->nodes.begin())+"\n");
//                    #endif
                    ci = classify[n].components.erase(ci);
                    continue;
                }
                
                unsigned int in = 0;
                unsigned int out = 0;
                unsigned int edges = 0;
                for ( std::set<int>::iterator it = ci->nodes.begin(); it != ci->nodes.end(); ++it) {                   
                    ListDigraph::Arc a = g.arcFromId(*it);
                    if (g.target(a) == node ) {
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
                
                classify[n].assigned_in += in;
                classify[n].assigned_out += out;
                classify[n].assigned_edges += edges;
//                
//                #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("CI "+ std::to_string(g.id(n)) + " as " + std::to_string(in)+"-"+std::to_string(out) + " " + std::to_string(edges) +"\n");
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
                            ListDigraph::Arc a = g.arcFromId(*it);
                            if (g.target(a) == node ) {
                               temp_know_paths[a].clear();
                            } else {
                               temp_know_back_paths[a].clear();
                            }
                            block_delete_high_res[n].erase(*it);
                            block_delete[n].erase(*it);
                        }     
                        // rewind
                        classify[n].assigned_in -= in;
                        classify[n].assigned_out -= out;
                        classify[n].assigned_edges -= edges;
                        
                        ci->in_nodes = 0;
                        ci->out_nodes = 0;
                        ci->edges = 0;
                        
                    } else {    
                        // we are fully connected!
                        classify[n].fully_connected_in += in;
                        classify[n].fully_connected_out += out;
                        classify[n].fully_connected_edges += edges;
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
            logger::Instance()->debug("CLASSIFY SET On "+ std::to_string(g.id(n)) + " with " + std::to_string(pre_size)  + "-" + std::to_string(classify[n].assigned_in) + "  " + std::to_string(post_size)+ "-" + std::to_string(classify[n].assigned_out)+"\n");
            #endif
            
            classify[n].unassigned_in = pre_size - classify[n].assigned_in;
            classify[n].unassigned_out = post_size - classify[n].assigned_out;   
         
            if (classify[n].components.size() == 1) {
                
                std::deque<component >::iterator ci = classify[n].components.begin();
                if (classify[n].unassigned_in == 0) {
                    // join out to max in!
                    capacity_type max = 0;
                    ListDigraph::Arc marc;
                    for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
                        if (ci->nodes.find(g.id(a)) != ci->nodes.end()) {
                            if (flow[a] > max) {
                                max = flow[a];
                                marc = a;
                            }
                        }
                    }
                    for (ListDigraph::OutArcIt a(g, n); a!=INVALID; ++a) {
                        if (ci->nodes.find(g.id(a)) == ci->nodes.end() && !barred[a]) {

                            temp_know_paths_add_max[marc][g.id(a)] += 1;
                            temp_know_back_paths_add_max[a].bridges->insert(g.id(marc));
                            ci->fragments.push_back(resolve_count_set());
                            ci->fragments.back().left.insert(g.id(marc));
                            ci->fragments.back().right.insert(g.id(a));
                            ci->fragments.back().count = 1;
                            ci->nodes.insert(g.id(a));
                        }
                    }
                } else if (classify[n].unassigned_out == 0) {
                     // join in to max out!
                    capacity_type max = 0;
                    ListDigraph::Arc marc;
                    for (ListDigraph::OutArcIt a(g, n); a!=INVALID; ++a) {
                        if (ci->nodes.find(g.id(a)) != ci->nodes.end()) {
                            if (flow[a] > max) {
                                max = flow[a];
                                marc = a;
                            }
                        }
                    }
                    for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
                        if (ci->nodes.find(g.id(a)) == ci->nodes.end() && !barred[a]) {

                            temp_know_paths_add_max[a][g.id(marc)] += 1;
                            temp_know_back_paths_add_max[marc].bridges->insert(g.id(a));
                            ci->fragments.push_back(resolve_count_set());
                            ci->fragments.back().left.insert(g.id(a));
                            ci->fragments.back().right.insert(g.id(marc));
                            ci->fragments.back().count = 1;
                            ci->nodes.insert(g.id(a));
                        }
                    }
                }
            }

            
            unsigned int edges_in_groups = 0;
            unsigned int nodes_offset = 0;
            std::set<int> lgs;
            std::set<int> rgs;
            std::set<std::pair<int, int> > combinations;
            for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
                ListDigraph::Arc l = a;
                for (arc_bridge::iterator bi = know_paths[l].begin(); bi != know_paths[l].end(); ++bi) {
                    ListDigraph::Arc r = g.arcFromId(bi->first);

                    // get the original group combos!
                    int i1 = 0;
                    for ( ;i1 < left_groups.size(); ++i1) {
                        if(left_groups[i1].find(g.id(l)) != left_groups[i1].end()) {
                            lgs.insert(i1);
                            break;
                        }
                    }
                    int i2 = 0;
                    for ( ;i2 < right_groups.size(); ++i2) {
                        if(right_groups[i2].find(g.id(r)) != right_groups[i2].end()) {
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
                        combinations.insert(std::make_pair(i1 + g.id(l), i2));
                    } else if (i1 != left_groups.size() && i2 == right_groups.size()) {
                        combinations.insert(std::make_pair(i1, i2 + g.id(r)));
                    }
                    
                }
            }
            for ( std::set<int>::iterator ig = lgs.begin(); ig != lgs.end(); ++ig) {
                nodes_offset += left_groups[*ig].size() - 1;
            }
            for ( std::set<int>::iterator ig = rgs.begin(); ig != rgs.end(); ++ig) {
                nodes_offset += right_groups[*ig].size() - 1;
            }
            
            classify[n].edge_group_offset = edges_in_groups - combinations.size();
            classify[n].node_group_offset = nodes_offset;
            
//                         save the maximal member of each group!
//            for(std::set<std::set<int> >::iterator group = left_raw_groups.begin(); group != left_raw_groups.end(); ++group) {
//                float max = 0;
//                ListDigraph::Arc arc; 
//                
//                for (std::set<int>::iterator gi = group->begin(); gi != group->end(); ++gi) {    
//                    ListDigraph::Arc a = g.arcFromId(*gi);
//                    
//                    float score = means[a].compute_score();
//                    if (score > max) {
//                        max = score;
//                        arc = a;
//                    }
//                }
//                
//                if (block_delete[n].find(g.id(arc)) != block_delete[n].end()) block_delete_high_res[n].insert(g.id(arc));
//            }
//                    
//            for(std::set<std::set<int> >::iterator group = right_raw_groups.begin(); group != right_raw_groups.end(); ++group) {
//                float max = 0;
//                ListDigraph::Arc arc; 
//                
//                for (std::set<int>::iterator gi = group->begin(); gi != group->end(); ++gi) {    
//                    ListDigraph::Arc a = g.arcFromId(*gi);
//                    
//                    float score = means[a].compute_score();
//                    if (score > max) {
//                        max = score;
//                        arc = a;
//                    }
//                }
//                
//                if (block_delete[n].find(g.id(arc)) != block_delete[n].end()) block_delete_high_res[n].insert(g.id(arc));
//            }
            
        }
    }
    
    
    // ######### Test Actions! ######### //
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Actions ##################### \n");
    digraphWriter(g, std::cout)
                .arcMap("edge_specifier", edge_specifier)
                .arcMap("edge_type", edge_type)
                .arcMap("cap", flow)
                .arcMap("length", edge_lengths)
                .arcMap("means", means)
                .arcMap("barred", barred)
                .run();  

//    logger::Instance()->debug("Has Block Mark \n");
//    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//        logger::Instance()->debug("Arc " + std::to_string(g.id(a)) + "\n"); 
//        bool barred = false;
//        for(std::set<transcript_unsecurity>::iterator u_it = unsecurityArc[a]->begin(); u_it != unsecurityArc[a]->end() ; ++u_it) {
//
//            logger::Instance()->debug(std::to_string(u_it->id) + "-" + std::to_string(u_it->evidenced) + ", " );
//        }
//        logger::Instance()->debug("\n"); 
//    }
//    logger::Instance()->debug("Nodes \n");
//    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
//        logger::Instance()->debug("Node " + std::to_string(g.id(n)) + " : " + std::to_string(unsecurityId[n].id) + "\n"); 
//    }

    logger::Instance()->debug("Classify done \n");
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
       logger::Instance()->debug("Arc " + std::to_string(g.id(a))); 
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
//    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//       logger::Instance()->debug("Arc " + std::to_string(g.id(a))); 
//       for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
//           logger::Instance()->debug(" cf " + std::to_string(ab->first)); 
//       }
//       for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
//           logger::Instance()->debug(" cb " + std::to_string(*ab)); 
//       }
//       logger::Instance()->debug("\n" );  
//    }
    
    for (ListDigraph::NodeIt n(g); n != INVALID;++n) {
        logger::Instance()->debug("For n " + std::to_string(g.id(n)) +": \n");
        logger::Instance()->debug("Left Groups\n");
        for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
            logger::Instance()->debug("Arc " + std::to_string(g.id(a)) + " G " + std::to_string(left_groups_ids[n][g.id(a)].id)  + "\n"); 
        }
        logger::Instance()->debug("Right Groups\n");
        for (ListDigraph::OutArcIt a(g, n); a!=INVALID; ++a) {
            logger::Instance()->debug("Arc " + std::to_string(g.id(a)) + " G " + std::to_string(right_groups_ids[n][g.id(a)].id)  + "\n"); 
        }
    }
    #endif

    
//    std::unordered_set<int> delete_block;
//    for (ListDigraph::ArcIt a(g); a != INVALID;++a) {  
//        
//        ListDigraph::Node ns = g.source(a);
//        ListDigraph::Node nt = g.target(a);
//        
//        bool lblock = block_delete_high_res[ns].find(g.id(a)) !=  block_delete_high_res[ns].end();
//        bool rblock = block_delete_high_res[nt].find(g.id(a)) !=  block_delete_high_res[nt].end();
//       
//        if ( rblock || lblock ) {
//            delete_block.insert(g.id(a));
//        }
//    }
    
    
    std::deque<ListDigraph::Node> top_order;

    pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > visitor(top_order);
    DfsVisit<ListDigraph, pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > > dfs(g, visitor);
    dfs.init();
    dfs.addSource(s);
    dfs.start();
    
    std::deque<std::pair<ListDigraph::Node, float> > ratio_order;
    for (ListDigraph::NodeIt n(g); n != INVALID;++n) {  
        float error = 0;
        capacity_type f_count = 0;
        for ( std::deque<component>::iterator ci = classify[n].components.begin(); ci != classify[n].components.end(); ++ci) {
            compute_ratio(n, temp_know_paths, temp_know_back_paths,  ci->nodes, g, flow, error, f_count);
        }
        ratio_order.push_back(std::make_pair(n, error / (float) f_count));
    }
    std::sort(ratio_order.begin(), ratio_order.end(), [](const std::pair<ListDigraph::Node, float> &a1, const std::pair<ListDigraph::Node, float> &a2){return a1.second < a2.second;});
    
    
    const float low_barr_ratio = 0.40;
        
    // we have a perfect match
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
        
        if ( classify[n].unassigned_in + classify[n].assigned_in < 2 || classify[n].unassigned_out + classify[n].assigned_out < 2 ) {
            continue;
        }
        
        unsigned int max = classify[n].assigned_in;
        if (classify[n].assigned_out > max) {
            max = classify[n].assigned_out;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Testing "+ std::to_string(g.id(n)) + " un " + std::to_string(classify[n].unassigned_in)+"-"+std::to_string(classify[n].unassigned_out) + " as " + std::to_string(classify[n].assigned_in)+"-"+std::to_string(classify[n].assigned_out) + " " + std::to_string(classify[n].assigned_edges) +"\n");
        #endif
        
        if (classify[n].unassigned_in == 0 && classify[n].unassigned_out == 0 && classify[n].assigned_edges / 2 == max && classify[n].components.size() == 1) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option A "+ std::to_string(g.id(n)) +"\n");
            #endif
            //this is perfectly resolved! use ALL
            insert_local_known(n, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
            add_local_known(n, g, temp_know_paths_add_max, temp_know_back_paths_add_max, know_paths, know_back_paths);

            capacity_type max = 0;
            for ( std::deque<component>::iterator ci = classify[n].components.begin(); ci != classify[n].components.end(); ++ci) {
                capacity_type count = unravel_evidences_ILP(n, know_paths, know_back_paths, left_groups_ids[n], right_groups_ids[n], barred, ci->fragments, ci->nodes, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, node_index, meta->size);
                if (count > max) {
                    max = count;
                }
            }
            clean_ILP_leftovers(n, know_paths, know_back_paths, barred, max * low_barr_ratio, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
            return true;
        }  
    }

    // we have a near perfect match
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
    
         if ( classify[n].unassigned_in + classify[n].assigned_in < 2 || classify[n].unassigned_out + classify[n].assigned_out < 2 ) {
            continue;
        }
        
        unsigned int edge_count = classify[n].assigned_edges / 2 - classify[n].edge_group_offset;
        unsigned int nodes = classify[n].assigned_in + classify[n].assigned_out - classify[n].node_group_offset;
        
        if (edge_count + 2 <= 1 + nodes
                && classify[n].unassigned_in == 0 && classify[n].unassigned_out == 0 && classify[n].components.size() <= 2) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option B "+ std::to_string(g.id(n)) +"\n");
            #endif
            // one component and easy degree, also use ILP!
            insert_local_known(n, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
            add_local_known(n, g, temp_know_paths_add_max, temp_know_back_paths_add_max, know_paths, know_back_paths);
            
            capacity_type max = 0;
            for ( std::deque<component>::iterator ci = classify[n].components.begin(); ci != classify[n].components.end(); ++ci) {
                capacity_type count = unravel_evidences_ILP(n, know_paths, know_back_paths, left_groups_ids[n], right_groups_ids[n], barred, ci->fragments, ci->nodes, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, node_index, meta->size);
                if (count > max) {
                    max = count;
                }
            }
            clean_ILP_leftovers(n, know_paths, know_back_paths,  barred, max * low_barr_ratio, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
            return true;
        }
        
    }
      
    if (!final_resolve) {
        bool barr_edge = false;
        ListDigraph::Node node;
        ListDigraph::Arc barc;
        float dr = 0;
        for (std::deque<ListDigraph::Node>::iterator tn = top_order.begin(); tn != top_order.end(); ++tn) {
            ListDigraph::Node n(*tn);
            
            
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("CHECK BARR NODE"+ std::to_string(g.id(n)) +"\n");
                #endif
            
            if ( classify[n].unassigned_in + classify[n].assigned_in < 2 || classify[n].unassigned_out + classify[n].assigned_out < 2 ) {
                continue;
            }
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("IN "+ std::to_string(g.id(n)) +"\n");
            #endif
                
            ListDigraph::Arc nbarc;
            float ndr = 0;
            if (bar_smallest_edge(n, g, flow, means, barred, 0.33, nbarc, ndr, block_delete_high_res[n])) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Report "+ std::to_string(g.id(nbarc)) + " " + std::to_string(ndr) +"\n");
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
                if (g.target(barc) == node) {
                    node2 = g.source(barc);
                } else {
                    node2 = g.target(barc);
                }
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Option C.1 "+ std::to_string(g.id(node)) + " Arc " + std::to_string(g.id(barc)) +"\n");
                #endif
                barred[barc] = true;
                unsecurityArc[barc]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
                insert_local_known(node, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
                clean_barred_leftovers(node, know_paths, know_back_paths, barred, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
                clean_barred_leftovers(node2, know_paths, know_back_paths, barred, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);

                return true;
        }
    }
             
    if (!final_resolve) {
        bool barr_edge = false;
        ListDigraph::Node node;
        ListDigraph::Arc barc;
        float dr = 0;
        for (std::deque<ListDigraph::Node>::iterator tn = top_order.begin(); tn != top_order.end(); ++tn) {
            ListDigraph::Node n(*tn);
            
            if ( classify[n].unassigned_in + classify[n].assigned_in < 2 || classify[n].unassigned_out + classify[n].assigned_out < 2 ) {
                continue;
            }
            
            ListDigraph::Arc nbarc;
            float ndr = 0;
            if (bar_negligible_edge(n, g, flow, means, barred, 0.05, nbarc, ndr, block_delete_high_res[n])) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Report "+ std::to_string(g.id(nbarc)) + " " + std::to_string(ndr) +"\n");
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
                if (g.target(barc) == node) {
                    node2 = g.source(barc);
                } else {
                    node2 = g.target(barc);
                }
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Option N.1 "+ std::to_string(g.id(node)) + " Arc " + std::to_string(g.id(barc)) +"\n");
                #endif
                barred[barc] = true;
                unsecurityArc[barc]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
                insert_local_known(node, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
                clean_barred_leftovers(node, know_paths, know_back_paths, barred, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
                clean_barred_leftovers(node2, know_paths, know_back_paths, barred, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);

                return true;
        }
    }
    
    
        // include near perfect multi component matches
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
        
        if ( classify[n].unassigned_in + classify[n].assigned_in < 2 || classify[n].unassigned_out + classify[n].assigned_out < 2 ) {
            continue;
        }
        
        unsigned int components = classify[n].unassigned_in + classify[n].unassigned_out + classify[n].components.size();
        
        unsigned int edge_count = classify[n].assigned_edges / 2 - classify[n].edge_group_offset;
        unsigned int nodes = classify[n].assigned_in + classify[n].assigned_out - classify[n].node_group_offset;

        if (edge_count + 2 * components <= 3 + classify[n].unassigned_in + classify[n].unassigned_out + nodes
            && components > 1 && classify[n].assigned_edges > 0) {
            // multi component and easy degree, also use ILP!
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option D "+ std::to_string(g.id(n)) +"\n");
            #endif
            insert_local_known(n, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
            add_local_known(n, g, temp_know_paths_add_max, temp_know_back_paths_add_max, know_paths, know_back_paths);
            
            capacity_type max = 0;
            for ( std::deque<component>::iterator ci = classify[n].components.begin(); ci != classify[n].components.end(); ++ci) {
                capacity_type count = unravel_evidences_ILP(n, know_paths, know_back_paths, left_groups_ids[n], right_groups_ids[n], barred, ci->fragments, ci->nodes, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, node_index, meta->size);
                if (count > max) {
                    max = count;
                }
            }
            clean_ILP_leftovers(n, know_paths, know_back_paths,  barred, max * low_barr_ratio, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
            return true;
        }
    }

        // clean up bad multi unique
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
        
        if ( classify[n].unassigned_in + classify[n].assigned_in < 2 || classify[n].unassigned_out + classify[n].assigned_out < 2 ) {
            continue;
        }
        
        if (classify[n].assigned_in > classify[n].fully_connected_in && classify[n].assigned_out > classify[n].fully_connected_out && classify[n].assigned_edges > classify[n].fully_connected_edges ) {
            
            bool has_loose_ends = false;
            
            std::set<int> in_groups, out_groups;
            std::map<int, std::set<int> > rev_in_groups;
            std::map<int, std::set<int> > rev_out_groups;
            for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
                if (edge_type[a] == edge_types::HELPER) continue;
                int id = g.id(a);
                in_groups.insert(left_groups_ids[n][id].id);
                rev_in_groups[left_groups_ids[n][id].id].insert(id);
            }
            for (ListDigraph::OutArcIt a(g, n); a!=INVALID; ++a) {
                if (edge_type[a] == edge_types::HELPER) continue;
                int id = g.id(a);
                out_groups.insert(right_groups_ids[n][id].id);
                rev_out_groups[right_groups_ids[n][id].id].insert(id);
            }
            for (std::set<int>::iterator si = in_groups.begin(); si != in_groups.end() && !has_loose_ends; ++si) {
                std::set<int> kg;
                for (std::set<int>::iterator fi = rev_in_groups[*si].begin(); fi != rev_in_groups[*si].end(); ++fi) {
                    for (arc_bridge::iterator bi = know_paths[g.arcFromId(*fi)].begin(); bi != know_paths[g.arcFromId(*fi)].end(); ++bi) {
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
                    for (arc_back_bridge::iterator bi = know_back_paths[g.arcFromId(*fi)].begin(); bi != know_back_paths[g.arcFromId(*fi)].end(); ++bi) {
                        kg.insert(left_groups_ids[n][*bi].id);
                    }
                }
                if (kg.size() == 1) {
                    has_loose_ends = true;
                }
            }
            
            if (has_loose_ends) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Option F "+ std::to_string(g.id(n)) +"\n");
                #endif
                insert_local_known(n, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
                
                unravel_evidences_groups(n, know_paths, know_back_paths, left_groups_ids[n], right_groups_ids[n], barred, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, node_index, meta->size);

//                unravel_single(n, know_paths, know_back_paths, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
                clean_barred_leftovers(n, know_paths, know_back_paths, barred, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);

                return true;
            }
        }
    }
        
    // we have a match although bad
    for (std::deque<std::pair<ListDigraph::Node, float> >::iterator tn = ratio_order.begin(); tn != ratio_order.end(); ++tn) {
        ListDigraph::Node n(tn->first);
        
         if ( classify[n].unassigned_in + classify[n].assigned_in < 2 || classify[n].unassigned_out + classify[n].assigned_out < 2 ) {
            continue;
        }
        
         #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option G TEST "+  std::to_string(g.id(n)) + ": "+ std::to_string(classify[n].unassigned_in) + " " + std::to_string(classify[n].unassigned_out) + " " + std::to_string(classify[n].components.size())  +"\n");
            #endif
         
        if (classify[n].unassigned_in + classify[n].unassigned_out + classify[n].components.size() == 1) { // one component, although bad edgecount
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Option G "+ std::to_string(g.id(n)) +"\n");
            #endif
            // one component and easy degree, also use ILP!
            insert_local_known(n, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
            add_local_known(n, g, temp_know_paths_add_max, temp_know_back_paths_add_max, know_paths, know_back_paths);
            
            capacity_type max = 0;
            for ( std::deque<component>::iterator ci = classify[n].components.begin(); ci != classify[n].components.end(); ++ci) {
                capacity_type count = unravel_evidences_ILP(n, know_paths, know_back_paths, left_groups_ids[n], right_groups_ids[n], barred, ci->fragments, ci->nodes, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, node_index, meta->size);
                if (count > max) {
                    max = count;
                }
            }
            clean_ILP_leftovers(n, know_paths, know_back_paths, barred, max * low_barr_ratio, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
            return true;
        }
    }
        
    if (!final_resolve) {
        bool barr_edge = false;
        ListDigraph::Node node;
        ListDigraph::Arc barc;
        float dr = 0;
        for (std::deque<ListDigraph::Node>::iterator tn = top_order.begin(); tn != top_order.end(); ++tn) {
            ListDigraph::Node n(*tn);
        
            if ( classify[n].unassigned_in + classify[n].assigned_in < 2 || classify[n].unassigned_out + classify[n].assigned_out < 2 ) {
                continue;
            }
            
            ListDigraph::Arc nbarc;
            float ndr = 0;
            if (bar_negligible_edge(n, g, flow, means, barred, 0.75, nbarc, ndr, block_delete_high_res[n])) {
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
                if (g.target(barc) == node) {
                    node2 = g.source(barc);
                } else {
                    node2 = g.target(barc);
                }
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Option H "+ std::to_string(g.id(node)) + " Arc " + std::to_string(g.id(barc)) +"\n");
                #endif
                barred[barc] = true;
                unsecurityArc[barc]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
                insert_local_known(node, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
                clean_barred_leftovers(node, know_paths, know_back_paths, barred, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
                clean_barred_leftovers(node2, know_paths, know_back_paths, barred, g, flow, means, edge_specifier, edge_type, edge_lengths, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
                return true;
        }
    } 
    
    for (ListDigraph::NodeIt n(g); n != INVALID;++n) {
         insert_local_known(n, g, temp_know_paths, temp_know_back_paths, know_paths, know_back_paths);
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

bool base_manager::bar_negligible_edge(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<bool> &barred, float ratio, ListDigraph::Arc &barc, float &value, std::unordered_set<int> &block_delete) {
        
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
        
        float score = mc[a].compute_score();
        
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

        float score = mc[a].compute_score();
        
        logger::Instance()->debug("RA " + std::to_string(wc.id(a)) + " " + std::to_string(score) + " " + std::to_string(block_delete.find(wc.id(a)) == block_delete.end()) + "\n");

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


bool base_manager::bar_smallest_edge(ListDigraph::Node node, ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<bool> &barred, float ratio, ListDigraph::Arc &barc, float &value, std::unordered_set<int> &block_delete) {
        
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
        
        float score = mc[a].compute_score();
        
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

        float score = mc[a].compute_score();
        
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

bool base_manager::bar_smallest_edge_by_group(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, std::set<int> &group, float ratio, ListDigraph::Arc &barc, float &value) {
        
    float max = 0;
    float sum = 0;
    float min = -1;
    ListDigraph::Arc arc;
    
    
    for (std::set<int>::iterator gi = group.begin(); gi != group.end(); ++gi) {    
        ListDigraph::Arc a = wc.arcFromId(*gi);
        
        float score = mc[a].mean;
        if (score > max) {
            max = score;
        }
        if (score < min || min < 0) {
            min = score;
            arc = a;
        }
        sum += score;
    }   
    
    bool bar = max != 0 && min > 0 && (float) min <= ratio * (float) max;
    float r = 1;
    if (sum != 0 && min > 0) r = min / sum;
   
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Barr " + std::to_string(wc.id(arc)) + " "  + std::to_string(bar) + "\n");
    #endif
    
    if (bar) {
        barc = arc;
        value = r;
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
    
ListDigraph::Arc base_manager::extract_path(ListDigraph::Arc left_arc, ListDigraph::Arc right_arc, capacity_type capacity,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc,
            ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId) {
    
    // the idea is to extract the flow until it is a full separate egde that can be removed without any further problems
    // this is one call on two edges
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract Path. " + std::to_string(capacity) + " " + std::to_string(wc.id(left_arc)) + " " + std::to_string(wc.id(right_arc))  + " \n");
    #endif

    
    capacity_mean c1 = mc[left_arc];
    capacity_mean c2 = mc[right_arc];
    
    c1.reduce(capacity/ float(fc[left_arc]));
    c2.reduce(capacity/ float(fc[right_arc]));
    mc[left_arc].reduce((fc[left_arc] - capacity)/ float(fc[left_arc]));
    mc[right_arc].reduce((fc[right_arc] - capacity)/ float(fc[right_arc]));
    c1.update(c2);
    
    
    // we know we can extract the capacity freely
    fc[left_arc] -= capacity;
    fc[right_arc] -= capacity;
    
    // new we create a new edge out of the two (this is the extraction part, in the end we have a full s-t arc)
    ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));
    
    fc[new_arc] = capacity;
    cet[new_arc] = edge_types::EXON;
    mc[new_arc] = c1;

//    TODO: this should NOT be needed
//    if (cet[left_arc] == edge_types::HELPER) {
//        // leftarc is a helper
//        ces[new_arc] = exon_edge(ces[right_arc]); 
//    } else if (cet[right_arc] == edge_types::HELPER) {
//        // rightarc is a helper
//        ces[new_arc] = exon_edge(ces[left_arc]); 
//    } else {
//        ces[new_arc] = exon_edge(ces[left_arc]);
//        ces[new_arc].join_edge(ces[right_arc]);
//    }
//    
//    cel[new_arc].first_exon = cel[left_arc].first_exon;
//    cel[new_arc].middle = cel[left_arc].middle + cel[left_arc].last_exon + cel[right_arc].middle;
//    cel[new_arc].last_exon = cel[right_arc].last_exon;
//    
//    cycle_id_in[new_arc] = cycle_id_in[left_arc];
//    cycle_id_out[new_arc] = cycle_id_out[right_arc];
    
    // we now need to transfer the know paths of the right and left side
    know_paths[new_arc].bridges.ref() = know_paths[right_arc].bridges.ref(); // deep copy 
   // know_back_paths[new_arc] = know_back_paths[left_arc]; // deep copy, if called right this is always empty!
    
    if (fc[left_arc] != 0 && fc[right_arc] != 0) {
        // we are done already, because nothing changed in the evidence of the edges!
        return new_arc;
    }
    
    // if we are here, one of the edges has been removed
    // first we remove all invalidated
    
    if (fc[left_arc] == 0) {
        for (arc_bridge::iterator bi = know_paths[left_arc].begin(); bi != know_paths[left_arc].end(); ++bi) {
            ListDigraph::Arc a = wc.arcFromId(bi->first);
            know_back_paths[a].remove_id(wc.id(left_arc));
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("REMOVE BACK " + std::to_string(wc.id(a)) + " " + std::to_string(wc.id(left_arc))  + " \n");
            #endif
        }
    } 
    if (fc[right_arc] == 0) {
        
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
   
    unravel_evidences(wc.target(left_arc), know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, NULL);

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
            ListDigraph::ArcMap<capacity_type> &fc, float &error, capacity_type &flow) {
    
    std::deque<ListDigraph::Arc> inArc, outArc;

    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        if (component.find(wc.id(a)) != component.end()) {
            flow += fc[a];
             inArc.push_back(a);
        }
    }

    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        if (component.find(wc.id(a)) != component.end()) {
            flow += fc[a];
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
        abs1 += fc[left];
        for (arc_bridge::iterator bi = know_paths[left].begin(); bi != know_paths[left].end(); ++bi) {
            ListDigraph::Arc right = wc.arcFromId(bi->first);
            
            abs1 -= bi_edges[std::make_pair(wc.id(left), wc.id(right))];
            abs2 += bi_edges[std::make_pair(wc.id(left), wc.id(right))]; 
        }
        
        left_z.push_back(lp.addCol());
        lp.addRow(left_z.back() >= abs1);
        lp.addRow( abs2 <= fc[left]); 
        optimize += left_z.back();// * fc[left];
    }

    for (std::deque<ListDigraph::Arc>::iterator a_it = outArc.begin(); a_it != outArc.end(); ++a_it) {
        ListDigraph::Arc right = *a_it;
        
        Lp::Expr abs1, abs2;
        abs1 += fc[right];
        for (arc_back_bridge::iterator bi = know_back_paths[right].begin(); bi != know_back_paths[right].end(); ++bi) {
            ListDigraph::Arc left = wc.arcFromId(*bi);
            
            abs1 -= bi_edges[std::make_pair(wc.id(left), wc.id(right))];
            abs2 += bi_edges[std::make_pair(wc.id(left), wc.id(right))]; 
        }
        
        right_z.push_back(lp.addCol());
        lp.addRow(right_z.back() >= abs1);
        lp.addRow(abs2 <= fc[right]); 
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
        error += std::sqrt(fc[left] - xee_sum);
    }
    
    for (std::deque<ListDigraph::Arc>::iterator a_it = outArc.begin(); a_it != outArc.end(); ++a_it) {
        ListDigraph::Arc right = *a_it;
        
        int xee_sum = 0; 
        for (arc_back_bridge::iterator bi = know_back_paths[right].begin(); bi != know_back_paths[right].end(); ++bi) {
            ListDigraph::Arc left = wc.arcFromId(*bi);
            
            xee_sum += std::round(lp.primal(bi_edges[std::make_pair(wc.id(left), wc.id(right))]));      
        }
        error += std::sqrt(fc[right] - xee_sum);
    }  
     
}


capacity_type base_manager::unravel_evidences_ILP(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            std::map<int, evidence_group> & left_groups, std::map<int, evidence_group> & right_groups, ListDigraph::ArcMap<bool> &barred,
            std::deque< resolve_count_set > &hit_counter_all,
            std::set<int> &component,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, ListDigraph::NodeMap<unsigned int> &ni, unsigned int size) {
    
    using namespace lemon;
    
    capacity_type flow_reference_in = 0;
    capacity_type flow_reference_out = 0;
    capacity_type flow_evidence = 0;
    
    std::deque<ListDigraph::Arc> inArc;
    
    std::set<int> in_groups, out_groups;
    std::map<int, std::set<int> > rev_in_groups;
    std::map<int, std::set<int> > rev_out_groups;
    
    std::map<int, evidence_group* > rev_in_groups_ev;
    std::map<int, evidence_group* > rev_out_groups_ev;
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        int id = wc.id(a);
        if (component.find(id) != component.end()) {
            inArc.push_back(a);
            flow_reference_in += fc[a];
            in_groups.insert(left_groups[id].id);
            rev_in_groups[left_groups[id].id].insert(id);
            rev_in_groups_ev[left_groups[id].id] = &left_groups[id];
        }
    }

    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        int id = wc.id(a);
        if (component.find(id) != component.end()) { 
            flow_reference_out += fc[a];
            out_groups.insert(right_groups[id].id);
            rev_out_groups[right_groups[id].id].insert(id);
            rev_out_groups_ev[right_groups[id].id] = &right_groups[id];
        }
    }

    capacity_type flow_reference = std::min(flow_reference_in, flow_reference_out);

    for (std::deque< resolve_count_set >::iterator ha = hit_counter_all.begin(); ha != hit_counter_all.end(); ++ha) {
        flow_evidence += ha->count;
    }
    
    // now vamp up the evidences to match!
    for (std::deque< resolve_count_set >::iterator ha = hit_counter_all.begin(); ha != hit_counter_all.end(); ++ha) {
        
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
        ha->count = ha->count * flow_reference/ (float) flow_evidence;
        
        logger::Instance()->debug("Correct to " + std::to_string(ha->count) + "\n");
    }
    
    Lp lp;
    std::map<std::pair< int, int > , Lp::Col> bi_edges;
    std::deque<Lp::Col> left_z;
    std::deque<Lp::Col> right_z;
    Lp::Expr optimize;
    
    for (std::deque<ListDigraph::Arc>::iterator a_it = inArc.begin(); a_it != inArc.end(); ++a_it) {
        ListDigraph::Arc left = *a_it;
        
        for (arc_bridge::iterator bi = know_paths[left].begin(); bi != know_paths[left].end(); ++bi) {
            ListDigraph::Arc right = wc.arcFromId(bi->first);
            
            int lg = left_groups[wc.id(left)].id;
            int rg = right_groups[wc.id(right)].id;
            
            logger::Instance()->debug("Bi " + std::to_string(wc.id(left)) + " " + std::to_string(wc.id(right)) + " : " + std::to_string(lg) + " " + std::to_string(rg) + "\n");
            
            if ( bi_edges.find(std::make_pair(lg, rg)) != bi_edges.end() ) {
                continue; // only set up this variable once
            }
            
            bi_edges.insert(std::make_pair(std::make_pair(lg, rg), lp.addCol()));
            
            lp.colLowerBound(bi_edges[std::make_pair(lg, rg)], 0);
            lp.colUpperBound(bi_edges[std::make_pair(lg, rg)], Lp::INF);
        }
    }

    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        
        capacity_type group_cap = 0;
        std::set<int> know_groups;
        for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
            ListDigraph::Arc left = wc.arcFromId(*ri);
            group_cap += fc[left];
            for (arc_bridge::iterator bi = know_paths[left].begin(); bi != know_paths[left].end(); ++bi) {
                know_groups.insert(right_groups[bi->first].id);
            }
        }
        
        Lp::Expr abs1, abs2;
        abs1 += group_cap;
        for (std::set<int>::iterator ki = know_groups.begin(); ki != know_groups.end(); ++ki) {
            int rg = *ki;
            
            abs1 -= bi_edges[std::make_pair(lg, rg)];
            abs2 += bi_edges[std::make_pair(lg, rg)]; 
        }
        
        left_z.push_back(lp.addCol());
        lp.addRow(left_z.back() >= abs1);
        lp.addRow( abs2 <= group_cap); 
        optimize += left_z.back();
    }
    
    for (std::set<int>::iterator rg_it = out_groups.begin(); rg_it != out_groups.end(); ++rg_it) {
        int rg = *rg_it;
        
        capacity_type group_cap = 0;
        std::set<int> know_groups;
        for (std::set<int>::iterator ri = rev_out_groups[rg].begin(); ri != rev_out_groups[rg].end(); ++ri) {
            ListDigraph::Arc right = wc.arcFromId(*ri);
            group_cap += fc[right];
            for (arc_back_bridge::iterator bi = know_back_paths[right].begin(); bi != know_back_paths[right].end(); ++bi) {
                know_groups.insert(left_groups[*bi].id);
            }
        }
        
        Lp::Expr abs1, abs2;
        abs1 += group_cap;
        for (std::set<int>::iterator ki = know_groups.begin(); ki != know_groups.end(); ++ki) {
            int lg = *ki;
            
            abs1 -= bi_edges[std::make_pair(lg, rg)];
            abs2 += bi_edges[std::make_pair(lg, rg)]; 
        }
        
        left_z.push_back(lp.addCol());
        lp.addRow(left_z.back() >= abs1);
        lp.addRow( abs2 <= group_cap); 
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
                lp.colUpperBound(constraints.back(), ha->count);
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
            lp.addRow(all_combos == ha->count);
        }
    }
    
    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;

        std::set<int> know_groups;
        for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
            for (arc_bridge::iterator bi = know_paths[wc.arcFromId(*ri)].begin(); bi != know_paths[wc.arcFromId(*ri)].end(); ++bi) {
                ListDigraph::Arc right = wc.arcFromId(bi->first);
                know_groups.insert(right_groups[wc.id(right)].id);
            }
        }
        
        for (std::set<int>::iterator ki = know_groups.begin(); ki != know_groups.end(); ++ki) {
            int rg = *ki;
            
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
        return 0;
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
    
    for (std::deque<ListDigraph::Arc>::iterator a_it = inArc.begin(); a_it != inArc.end(); ++a_it) {
        ListDigraph::Arc left = *a_it;
        
        for (arc_bridge::iterator bi = know_paths[left].begin(); bi != know_paths[left].end(); ++bi) {
            ListDigraph::Arc right = wc.arcFromId(bi->first);
            
            int lg = left_groups[wc.id(left)].id;
            int rg = right_groups[wc.id(right)].id;
            
            if ( bi_edges2.find(std::make_pair(lg, rg)) != bi_edges2.end() ) {
                continue; // only set up this variable once
            }
            
            bi_edges2.insert(std::make_pair(std::make_pair(lg, rg), lp2.addCol()));
            
            lp2.colLowerBound(bi_edges2[std::make_pair(lg, rg)], lp.primal(bi_edges[std::make_pair(lg, rg)]));
            lp2.colUpperBound(bi_edges2[std::make_pair(lg, rg)], Lp::INF);
        }
    }

    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        
        capacity_type group_cap = 0;
        std::set<int> know_groups;
        for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
            ListDigraph::Arc left = wc.arcFromId(*ri);
            group_cap += fc[left];
            for (arc_bridge::iterator bi = know_paths[left].begin(); bi != know_paths[left].end(); ++bi) {
                know_groups.insert(right_groups[bi->first].id);
            }
        }
        
        Lp::Expr abs1, abs2;
        abs1 += group_cap;
        for (std::set<int>::iterator ki = know_groups.begin(); ki != know_groups.end(); ++ki) {
            int rg = *ki;
            
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
        std::set<int> know_groups;
        for (std::set<int>::iterator ri = rev_out_groups[rg].begin(); ri != rev_out_groups[rg].end(); ++ri) {
            ListDigraph::Arc right = wc.arcFromId(*ri);
            group_cap += fc[right];
            for (arc_back_bridge::iterator bi = know_back_paths[right].begin(); bi != know_back_paths[right].end(); ++bi) {
                know_groups.insert(left_groups[*bi].id);
            }
        }
        
        Lp::Expr abs1, abs2;
        abs1 += group_cap;
        for (std::set<int>::iterator ki = know_groups.begin(); ki != know_groups.end(); ++ki) {
            int lg = *ki;
            
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
        return 0;
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Costs " + std::to_string(lp2.primal()) + " as " + std::to_string(lp2.primal(optimize)) + "\n");
    for (std::map<std::pair< int, int > , Lp::Col>::iterator it = bi_edges2.begin(); it != bi_edges2.end(); ++it) {
         logger::Instance()->debug("XX2 " + std::to_string(it->first.first) + "-" + std::to_string(it->first.second) + " " + std::to_string(lp2.primal(it->second)) + "\n");
    }
    logger::Instance()->debug("LP B Done\n");
    #endif
    
    // primal solutions of Xees are now the interesting characteristic!
    // we extract exactly with them!
    
    // move around blocked nodes
    std::map<int, ListDigraph::Node > new_group_nodes_left;
    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        if (rev_in_groups_ev[lg]->block) {
            new_group_nodes_left[lg] = wc.addNode();
            ni[new_group_nodes_left[lg]] = ni[node];
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
        if (rev_out_groups_ev[rg]->block) {
            new_group_nodes_right[rg] = wc.addNode();
            ni[new_group_nodes_right[rg]] = ni[node];
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
    logger::Instance()->debug("Extract Groups\n");
    #endif
    
    capacity_type max_cap = 0;
    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        
        std::set<int> know_groups;
        for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
            ListDigraph::Arc left = wc.arcFromId(*ri);
            for (arc_bridge::iterator bi = know_paths[left].begin(); bi != know_paths[left].end(); ++bi) {
                ListDigraph::Arc right = wc.arcFromId(bi->first);
                if (right_groups[wc.id(right)].id != -1) know_groups.insert(right_groups[wc.id(right)].id);
            }
        }
        
        for (std::set<int>::iterator ki = know_groups.begin(); ki != know_groups.end(); ++ki) {
            int rg = *ki;
            
            capacity_type cap = lp2.primal(bi_edges2[std::make_pair(lg, rg)]) + 0.001;
            
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
            if (rev_in_groups_ev[lg]->block && rev_out_groups_ev[rg]->block) {
                
                #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Double Block \n");
                #endif

                ListDigraph::Arc new_arc = wc.addArc(new_group_nodes_left[lg], new_group_nodes_right[rg]);
                cet[new_arc] = edge_types::RESOLVE_HELPER;
                fc[new_arc] = cap;
                ces[new_arc] = exon_edge(size);
                
                float mean_left = 0;
                float score_left = 0;
                capacity_type weight_left = 0;
                for (ListDigraph::InArcIt a(wc, new_group_nodes_left[lg]); a!=INVALID; ++a) {
                    if (!barred[a]) {
                        mean_left += fc[a] * mc[a].mean;
                        score_left += fc[a] * mc[a].compute_score();
                        weight_left += fc[a];
                    }
                }
                mean_left = mean_left / (float) weight_left;
                score_left = score_left / (float) weight_left;
                
                float mean_right = 0;
                float score_right = 0;
                capacity_type weight_right = 0;
                for (ListDigraph::OutArcIt a(wc, new_group_nodes_right[rg]); a!=INVALID; ++a) {
                    if (!barred[a]) {
                        mean_right += fc[a] * mc[a].mean;
                        score_right += fc[a] * mc[a].compute_score();
                        weight_right += fc[a];
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

                float ratio = (float) cap / (float) weight;
                mc[new_arc].hidden_score = score * ratio;
                mc[new_arc].weight = 0;
                mc[new_arc].mean = mean * ratio;
                
            } else if (rev_in_groups_ev[lg]->block || rev_out_groups_ev[rg]->block) {
                
                #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Single Block \n");
                #endif
                
                ListDigraph::Arc base_arc, new_arc;
                
                if (rev_in_groups_ev[lg]->block) {
                    base_arc = wc.arcFromId(*rev_out_groups[rg].begin());
                    new_arc = wc.addArc(new_group_nodes_left[lg], wc.target(base_arc));
                } else { // rev_out_groups_ev[rg].block
                    base_arc = wc.arcFromId(*rev_in_groups[lg].begin());
                    new_arc = wc.addArc(wc.source(base_arc), new_group_nodes_right[rg]);
                }
            
                // just copy over stuff   
                ces[new_arc] = ces[base_arc];
                cet[new_arc] = cet[base_arc];
                cel[new_arc] = cel[base_arc];
                fc[new_arc] = cap;
                mc[new_arc] = mc[base_arc];
                mc[new_arc].reduce(cap/ float(fc[base_arc]));
                
                mc[base_arc].reduce((fc[base_arc] - cap)/ float(fc[base_arc]));
                fc[base_arc] -= cap;
                
                
                if (unsecurityArc[base_arc].has_ref()) {
                    std::copy(unsecurityArc[base_arc]->begin(), unsecurityArc[base_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
                }
                
                cycle_id_in[new_arc] = cycle_id_in[base_arc];
                cycle_id_out[new_arc] = cycle_id_out[base_arc];
                
                if (rev_in_groups_ev[lg]->block) {
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("In Know \n");
                    #endif
                    know_paths[new_arc].bridges.ref() = know_paths[base_arc].bridges.ref(); // deep copy
                    if (fc[base_arc] == 0) {
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
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Out Know \n");
                    #endif
                    know_back_paths[new_arc].bridges.ref() = know_back_paths[base_arc].bridges.ref(); // deep copy
                    if (fc[base_arc] == 0) {
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
                unravel_ILP(*rev_in_groups[lg].begin(), *rev_out_groups[rg].begin(), cap, know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);
            }
             
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("STEP\n"); 
            digraphWriter(wc, std::cout)
                        .arcMap("edge_specifier", ces)
                        .arcMap("edge_type", cet)
                        .arcMap("fc", fc)
                        .arcMap("length", cel)
                        .arcMap("means", mc)
                        .arcMap("barred", barred)
                        .run();   
            for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {

                logger::Instance()->debug("Arc " + std::to_string(wc.id(a))); 
                for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
                    logger::Instance()->debug(" f " + std::to_string(ab->first)); 
                }
                for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
                    logger::Instance()->debug(" b " + std::to_string(*ab)); 
                }
                logger::Instance()->debug("\n" );  
             }
             #endif
             
        }
    }
    
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Add Over and underflow to original edge \n");
    #endif
    
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_left.begin(); ngi != new_group_nodes_left.end(); ++ngi) {
        
        ListDigraph::Node n = ngi->second;
        
        float mean = 0;
        float score = 0;
        capacity_type weight = 0;
        capacity_type group_flow = 0;
        for (ListDigraph::InArcIt a(wc, n); a!=INVALID; ++a) {
            group_flow += fc[a];
            if (!barred[a]) {
                mean += fc[a] * mc[a].mean;
                score += fc[a] * mc[a].compute_score();
                weight += fc[a];
            }
        }
        mean = mean / (float) weight;
        score = score / (float) weight;
        
        capacity_type evidenced_flow = 0;
        for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
            evidenced_flow += fc[a];
        }
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("node " + std::to_string(wc.id(n)) + " GroupFlow " +  std::to_string(group_flow) + " EvFlow " +  std::to_string(evidenced_flow)  + "\n");
        #endif
        
        if (group_flow > evidenced_flow) {
            ListDigraph::Arc new_arc = wc.addArc(n, node);
            cet[new_arc] = edge_types::RESOLVE_HELPER;
            fc[new_arc] = group_flow - evidenced_flow;
            ces[new_arc] = exon_edge(size);
            
            float ratio = (float) evidenced_flow / (float) group_flow;
            mc[new_arc].hidden_score = score * ratio;
            mc[new_arc].weight = 0;
            mc[new_arc].mean = mean * ratio;
        }
    }
  
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_right.begin(); ngi != new_group_nodes_right.end(); ++ngi) {
        
        ListDigraph::Node n = ngi->second;
        
        float mean = 0;
        float score = 0;
        capacity_type weight = 0;
        capacity_type group_flow = 0;
        for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
            group_flow += fc[a];
            if (!barred[a]) {
                mean += fc[a] * mc[a].mean;
                score += fc[a] * mc[a].compute_score();
                weight += fc[a];
            }
        }
        mean = mean / (float) weight;
        score = score / (float) weight;
        
        capacity_type evidenced_flow = 0;
        for (ListDigraph::InArcIt a(wc, n); a!=INVALID; ++a) {
            evidenced_flow += fc[a];
        }
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("node " + std::to_string(wc.id(n)) + " GroupFlow " +  std::to_string(group_flow) + " EvFlow " +  std::to_string(evidenced_flow)  + "\n");
        #endif
        
        if (group_flow > evidenced_flow) {
            ListDigraph::Arc new_arc = wc.addArc(node, n);
            cet[new_arc] = edge_types::RESOLVE_HELPER;
            fc[new_arc] = group_flow - evidenced_flow;
            ces[new_arc] = exon_edge(size);
            
            float ratio = (float) evidenced_flow / (float) group_flow;
            mc[new_arc].hidden_score = score * ratio;
            mc[new_arc].weight = 0;
            mc[new_arc].mean = mean * ratio;
        }
    }
        
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Clean Leftovers New Nodes \n");
    #endif
   

    for (std::set<int>::iterator lg_it = in_groups.begin(); lg_it != in_groups.end(); ++lg_it) {
        int lg = *lg_it;
        if (rev_in_groups_ev[lg]->block) {
            for (std::set<int>::iterator ri = rev_in_groups[lg].begin(); ri != rev_in_groups[lg].end(); ++ri) {
                for (arc_bridge::iterator bi = know_paths[wc.arcFromId(*ri)].begin(); bi != know_paths[wc.arcFromId(*ri)].end(); ++bi) {
                    #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Remove " + std::to_string(*ri) + " AT " + std::to_string(bi->first) + "\n");
                    #endif
                    know_back_paths[wc.arcFromId(bi->first)].remove_id(*ri);
                }
                know_paths[wc.arcFromId(*ri)].clear();
            }
        }
    }
    for (std::set<int>::iterator rg_it = out_groups.begin(); rg_it != out_groups.end(); ++rg_it) {
        int rg = *rg_it;
        if (rev_out_groups_ev[rg]->block) {
            for (std::set<int>::iterator ri = rev_out_groups[rg].begin(); ri != rev_out_groups[rg].end(); ++ri) {
                for (arc_back_bridge::iterator bi = know_back_paths[wc.arcFromId(*ri)].begin(); bi != know_back_paths[wc.arcFromId(*ri)].end(); ++bi) {
                    #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Remove " + std::to_string(*ri) + " AT " + std::to_string(*bi) + "\n");
                    #endif
                    know_paths[wc.arcFromId(*bi)].remove_id(*ri);
                }
                know_back_paths[wc.arcFromId(*ri)].clear();
            }
        }
    }    
        
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_left.begin(); ngi != new_group_nodes_left.end(); ++ngi) {
        clean_barred_leftovers(ngi->second, know_paths, know_back_paths, barred, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
    }
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_right.begin(); ngi != new_group_nodes_right.end(); ++ngi) {
        clean_barred_leftovers(ngi->second, know_paths, know_back_paths, barred, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
    }
    
        
    // we now look at possibly unevidenced leftovers
    // remove all 0 edges here
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            for (arc_bridge::iterator bi = know_paths[arc].begin(); bi != know_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(bi->first);
                know_back_paths[a].remove_id(wc.id(arc));
            }
            know_paths[arc].clear();
            wc.erase(arc);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            for (arc_back_bridge::iterator bi = know_back_paths[arc].begin(); bi != know_back_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(*bi);
                know_paths[a].remove_id(wc.id(arc));
            }  
            know_back_paths[arc].clear();
            wc.erase(arc);
        }
    }      

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("FINISHED ILP\n"); 
    digraphWriter(wc, std::cout)
                .arcMap("edge_specifier", ces)
                .arcMap("edge_type", cet)
                .arcMap("fc", fc)
                .arcMap("length", cel)
                .arcMap("means", mc)
                .arcMap("barred", barred)
                .run();   
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {

        logger::Instance()->debug("Arc " + std::to_string(wc.id(a))); 
        for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
            logger::Instance()->debug(" f " + std::to_string(ab->first)); 
        }
        for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
            logger::Instance()->debug(" b " + std::to_string(*ab)); 
        }
        logger::Instance()->debug("\n" );  
     }
     #endif    
        
    return max_cap;
}

void base_manager::clean_ILP_leftovers(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap<bool> &barred, capacity_type barr_limit,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId) {
    
    std::set<int> nodes;
    nodes.insert(wc.id(node));
    
    bool has_left = false;
    bool has_right = false;
    
    for (ListDigraph::InArcIt in(wc, node); in != INVALID; ++in) {
        if (fc[in] < barr_limit ) { //&& know_paths[in].size() > 0) {
            barred[in] = true;
            unsecurityArc[in]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));  
        }
        nodes.insert(wc.id(wc.source(in))); 
    }
    for (ListDigraph::OutArcIt out(wc, node); out != INVALID; ++out) {
        if (fc[out] < barr_limit) { //&& know_back_paths[out].size() > 0) {
            barred[out] = true;
            unsecurityArc[out]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
        }
        nodes.insert(wc.id(wc.target(out)));
    }
   
    if (nodes.size() == 1) { 
        wc.erase(node);
    }
    
    for (std::set<int>::iterator in = nodes.begin(); in != nodes.end(); ++in) { 
        clean_barred_leftovers(wc.nodeFromId(*in), know_paths, know_back_paths, barred, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
    }
}

bool base_manager::clean_barred_leftovers(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph::ArcMap<bool> &barred,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId) {
    
    logger::Instance()->debug("clean barred leftovers\n");
    
    std::deque<ListDigraph::Arc> inArc, outArc;
    capacity_type in_flow = 0;
    capacity_type out_flow = 0;
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        if (!barred[a]) {
            inArc.push_back(a);
            in_flow += fc[a];
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        if (!barred[a]) {
            outArc.push_back(a);
            out_flow += fc[a];
        }
    }
    
    if (inArc.size() > 1 && outArc.size() > 1 || inArc.size() == 0 || outArc.size() == 0) {
        return false;
    }
    
    if (inArc.size() == 1) {
        
        ListDigraph::Arc left = inArc[0];
        
        float ratio = in_flow / (float) out_flow;
        if (ratio > 1) ratio = 1;
        
//        std::sort(outArc.begin(), outArc.end(), [&fc](const ListDigraph::Arc a1, const ListDigraph::Arc a2){return fc[a1] < fc[a2];});
       
        std::deque<capacity_type> flows;
        capacity_type total = 0;
        for (std::deque<ListDigraph::Arc>::iterator a = outArc.begin(); a != outArc.end(); ++a) {
            flows.push_back(fc[*a] * ratio);
            total += flows.back();
        } 
        capacity_type left_over = fc[left] - total;
        
        std::deque<capacity_type>::iterator fi = flows.begin();
        for (std::deque<ListDigraph::Arc>::iterator a = outArc.begin(); a != outArc.end(); ++a, ++fi) {
            
            
            capacity_type additional_flow = std::min(left_over, fc[*a] - *fi);
            left_over -= additional_flow;
            
            logger::Instance()->debug("left " + std::to_string(*fi) + " " + std::to_string(additional_flow) + " on " + std::to_string(fc[*a]) + " " + std::to_string(ratio) +"\n");
            
            unravel_ILP(wc.id(left), wc.id(*a), *fi + additional_flow, know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);
            barred[*a] = true;
            unsecurityArc[*a]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
            if (fc[left] == 0) {
                break;
            }
        } 
    } else {
        
        ListDigraph::Arc right = outArc[0];
        
        float ratio = out_flow / (float) in_flow;
        if (ratio > 1) ratio = 1;
                
//        std::sort(inArc.begin(), inArc.end(), [&fc](const ListDigraph::Arc a1, const ListDigraph::Arc a2){return fc[a1] < fc[a2];});
        
        std::deque<capacity_type> flows;
        capacity_type total = 0;
        for (std::deque<ListDigraph::Arc>::iterator a = inArc.begin(); a != inArc.end(); ++a) {
            flows.push_back(fc[*a] * ratio);
            total += flows.back();
        }
        capacity_type left_over = fc[right] - total;
        
        std::deque<capacity_type>::iterator fi = flows.begin(); 
        for (std::deque<ListDigraph::Arc>::iterator a = inArc.begin(); a != inArc.end(); ++a, ++fi) {
            
            capacity_type additional_flow = std::min(left_over, fc[*a] - *fi);
            left_over -= additional_flow;
            
            logger::Instance()->debug("right " + std::to_string(*fi) + " " + std::to_string(additional_flow) + " on " + std::to_string(fc[*a]) + " " + std::to_string(ratio) +"\n");

            unravel_ILP(wc.id(*a), wc.id(right), *fi + additional_flow, know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);
            barred[*a] = true;
            unsecurityArc[*a]->insert( transcript_unsecurity( 0, transcript_unsecurity::BARRED));
            if (fc[right] == 0) {
                break;
            }
        }
    }
    
        #ifdef ALLOW_DEBUG
    digraphWriter(wc, std::cout)
                .arcMap("edge_specifier", ces)
                .arcMap("edge_type", cet)
                .arcMap("fc", fc)
                .arcMap("length", cel)
                .arcMap("means", mc)
                .arcMap("barred", barred)
                .run();   
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {

       logger::Instance()->debug("Arc " + std::to_string(wc.id(a))); 
       for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
           logger::Instance()->debug(" f " + std::to_string(ab->first)); 
       }
       for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
           logger::Instance()->debug(" b " + std::to_string(*ab)); 
       }
       logger::Instance()->debug("\n" );  
    }
     #endif
    
     // we now look at possibly unevidenced leftovers
    // remove all 0 edges here
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            for (arc_bridge::iterator bi = know_paths[arc].begin(); bi != know_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(bi->first);
                know_back_paths[a].remove_id(wc.id(arc));
            }
            know_paths[arc].clear();
            wc.erase(arc);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            for (arc_back_bridge::iterator bi = know_back_paths[arc].begin(); bi != know_back_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(*bi);
                know_paths[a].remove_id(wc.id(arc));
            }  
            know_back_paths[arc].clear();
            wc.erase(arc);
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
        unravel_unevidenced_leftovers(node, know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, &barred);
        indeg = 0;
        ret = true;
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("After NEAR OUT \n");
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {

       logger::Instance()->debug("Arc " + std::to_string(wc.id(a))); 
       for (arc_bridge::iterator ab = know_paths[a].begin(); ab != know_paths[a].end(); ++ab) {
           logger::Instance()->debug(" f " + std::to_string(ab->first)); 
       }
       for (arc_back_bridge::iterator ab = know_back_paths[a].begin(); ab != know_back_paths[a].end(); ++ab) {
           logger::Instance()->debug(" b " + std::to_string(*ab)); 
       }
       logger::Instance()->debug("\n" );  
    }
     #endif
    
    
    if (indeg == 0) {
        // we completely resolved this one, remove!
        wc.erase(node);
    }
    return ret;
    
}


void base_manager::unravel_ILP(int left, int right, capacity_type cap_evidence,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths, 
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, transcript_unsecurity::evidence evidence) {
    
    
    if (cap_evidence == 0) {
        return;
    }
    
    ListDigraph::Arc left_arc = wc.arcFromId(left);
    ListDigraph::Arc right_arc = wc.arcFromId(right);
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract ILP. " + ces[left_arc].to_string() + ";" + ces[right_arc].to_string() + " - " +  std::to_string(left) + " " + std::to_string(right) + "\n");
    logger::Instance()->debug("Capacity " + std::to_string(cap_evidence) + "\n");
    #endif
    
    capacity_mean c1 = mc[left_arc];
    capacity_mean c2 = mc[right_arc];
    
    logger::Instance()->debug(c1.to_string()+ "\n");
    logger::Instance()->debug(c2.to_string()+ "\n");
    
    c1.reduce(cap_evidence/ float(fc[left_arc]));
    c2.reduce(cap_evidence/ float(fc[right_arc]));
    
    logger::Instance()->debug(c1.to_string()+ "\n");
    logger::Instance()->debug(c2.to_string()+ "\n");
    
    mc[left_arc].reduce((fc[left_arc] - cap_evidence)/ float(fc[left_arc]));
    mc[right_arc].reduce((fc[right_arc] - cap_evidence)/ float(fc[right_arc]));
    c1.update(c2);
    
    logger::Instance()->debug(c1.to_string()+ "\n");
    
    fc[left_arc] -= cap_evidence;
    fc[right_arc] -= cap_evidence;
    
    // this should never be called on a zero edge, so add in new arcs
    ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));
    
    fc[new_arc] = cap_evidence;
    cet[new_arc] = edge_types::EXON;
    mc[new_arc] = c1;
    
    if (cet[left_arc] != edge_types::EXON && cet[right_arc] != edge_types::EXON) {
        cet[new_arc] = edge_types::HELPER;
        if (ces[right_arc].node_index == -1) {
            ces[new_arc] = ces[left_arc];
        } else {
            ces[new_arc] = ces[right_arc];
        }
    } else if (cet[left_arc] == edge_types::HELPER) {
        // leftarc is a helper
        ces[new_arc] = exon_edge(ces[right_arc]); 
        
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else if (cet[right_arc] == edge_types::HELPER) {
        // rightarc is a helper
        ces[new_arc] = exon_edge(ces[left_arc]);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else {
        ces[new_arc] = exon_edge(ces[left_arc]);
        ces[new_arc].join_edge(ces[right_arc]);
        
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
    
    if (fc[left_arc] == 0) {
        
        for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
            know_paths[a].replace_evidence(wc.id(left_arc), wc.id(new_arc));
        }
          
    } else {
        
        for (ListDigraph::InArcIt a(wc, wc.source(left_arc)); a!=INVALID; ++a) {
            know_paths[a].add_evidence_if(wc.id(left_arc), wc.id(new_arc));
        }   
    }
    
    if (fc[right_arc] == 0) {
        
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].replace_evidence(wc.id(right_arc), wc.id(new_arc));
        }
          
    } else {
        
        for (ListDigraph::OutArcIt a(wc, wc.target(right_arc)); a!=INVALID; ++a) {
            know_back_paths[a].add_evidence_if(wc.id(right_arc), wc.id(new_arc));
        }
    }
     
    cel[new_arc].first_exon = cel[left_arc].first_exon;
    cel[new_arc].middle = cel[left_arc].middle + cel[left_arc].last_exon + cel[right_arc].middle;
    cel[new_arc].last_exon = cel[right_arc].last_exon;
    
    cycle_id_in[new_arc] = cycle_id_in[left_arc];
    cycle_id_out[new_arc] = cycle_id_out[right_arc];
        
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Exit ILP.\n");
    #endif

}

void base_manager::unravel_evidences_groups(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            std::map<int, evidence_group> & left_groups, std::map<int, evidence_group> & right_groups, ListDigraph::ArcMap<bool> &barred,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc,  ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId, ListDigraph::NodeMap<unsigned int> &ni, unsigned int size) {
    
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Unravel NODE GROUP" + std::to_string(wc.id(node)) + "\n");
    #endif
   
    struct info {
        
        info() {}
        info(int group, bool is_out, capacity_type cap, std::set<int> &s) : group(group), is_out(is_out), cap(cap), know_groups(s) 
            {}
        
        int group;
        bool is_out = false;
        capacity_type cap = 0;
        std::set<int> know_groups;
    };
    
    std::set<int> in_groups, out_groups;
    std::map<int, std::set<int> > rev_in_groups;
    std::map<int, std::set<int> > rev_out_groups;
    std::map<int, evidence_group* > rev_in_groups_ev;
    std::map<int, evidence_group* > rev_out_groups_ev;
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (cet[a] == edge_types::HELPER) {
            // no resolving 
            continue;
        }
        
        int id = wc.id(a);
        in_groups.insert(left_groups[id].id);
        rev_in_groups[left_groups[id].id].insert(id);
        rev_in_groups_ev[left_groups[id].id] = &left_groups[id];
    }

    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        
        if (cet[a] == edge_types::HELPER) {
            // no resolving 
            continue;
        }
        
        int id = wc.id(a);
        out_groups.insert(right_groups[id].id);
        rev_out_groups[right_groups[id].id].insert(id);
        rev_out_groups_ev[right_groups[id].id] = &right_groups[id];
    }

    std::deque<info> groups;
    std::map<int, info* > rev_groups_info;
    for (std::set<int>::iterator si = in_groups.begin(); si != in_groups.end(); ++si) {
        capacity_type group_flow = 0;
        std::set<int> kg;
        for (std::set<int>::iterator fi = rev_in_groups[*si].begin(); fi != rev_in_groups[*si].end(); ++fi) {
            group_flow += fc[wc.arcFromId(*fi)];
            for (arc_bridge::iterator bi = know_paths[wc.arcFromId(*fi)].begin(); bi != know_paths[wc.arcFromId(*fi)].end(); ++bi) {
                kg.insert(right_groups[bi->first].id);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("FG " + std::to_string(right_groups[bi->first].id) + "\n");
                #endif
            }
        }
        
        #ifdef ALLOW_DEBUG
         logger::Instance()->debug("FGT " + std::to_string(*si) + " " + std::to_string(group_flow) + "\n");
        #endif

        groups.push_back(info(*si, false, group_flow, kg));
        rev_groups_info[*si] = &groups.back();
    }
    for (std::set<int>::iterator si = out_groups.begin(); si != out_groups.end(); ++si) {
        capacity_type group_flow = 0;
        std::set<int> kg;
        for (std::set<int>::iterator fi = rev_out_groups[*si].begin(); fi != rev_out_groups[*si].end(); ++fi) {
            group_flow += fc[wc.arcFromId(*fi)];
            for (arc_back_bridge::iterator bi = know_back_paths[wc.arcFromId(*fi)].begin(); bi != know_back_paths[wc.arcFromId(*fi)].end(); ++bi) {
                kg.insert(left_groups[*bi].id);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("FGB " + std::to_string(left_groups[*bi].id) + "\n");
                #endif
            }
        }
        
        #ifdef ALLOW_DEBUG
         logger::Instance()->debug("FGTB " + std::to_string(*si) + " " + std::to_string(group_flow) + "\n");
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
                
                capacity_type min = std::min(it->cap, rev_groups_info[*it->know_groups.begin()]->cap);
                if (min > max ){
                    candidate = &*it;
                    max = min;
                }
            }
        }
        
        if (max == 0) break;
        
        // we have something to resolve, yay
        capacity_type cap = max;
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
        candidate->know_groups.clear();
        candidate->cap -= cap;
        rev_groups_info[del_in]->know_groups.erase(del);
        rev_groups_info[del_in]->cap -= cap;
        
        if (rev_in_groups_ev[lg]->block && new_group_nodes_left.find(lg) == new_group_nodes_left.end()) {
            new_group_nodes_left[lg] = wc.addNode();
            ni[new_group_nodes_left[lg]] = ni[node];
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
        
        if (rev_out_groups_ev[rg]->block && new_group_nodes_right.find(rg) == new_group_nodes_right.end()) {
            new_group_nodes_right[rg] = wc.addNode();
            ni[new_group_nodes_right[rg]] = ni[node];
            unsecurityId[new_group_nodes_right[rg]] = unsecurityId[node];
            
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Migrate Group " + std::to_string(rg) + " to new Group " + std::to_string(wc.id(new_group_nodes_left[rg])) + "\n");
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
 
        if (rev_in_groups_ev[lg]->block && rev_out_groups_ev[rg]->block) {
                
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Double Block \n");
            #endif
            
            ListDigraph::Arc new_arc = wc.addArc(new_group_nodes_left[lg], new_group_nodes_right[rg]);
            cet[new_arc] = edge_types::RESOLVE_HELPER;
            fc[new_arc] = cap;
            ces[new_arc] = exon_edge(size);
            
            float mean_left = 0;
            float score_left = 0;
            capacity_type weight_left = 0;
            for (ListDigraph::InArcIt a(wc, new_group_nodes_left[lg]); a!=INVALID; ++a) {
                if (!barred[a]) {
                    mean_left += fc[a] * mc[a].mean;
                    score_left += fc[a] * mc[a].compute_score();
                    weight_left += fc[a];
                }
            }
            mean_left = mean_left / (float) weight_left;
            score_left = score_left / (float) weight_left;

            float mean_right = 0;
            float score_right = 0;
            capacity_type weight_right = 0;
            for (ListDigraph::OutArcIt a(wc, new_group_nodes_right[rg]); a!=INVALID; ++a) {
                if (!barred[a]) {
                    mean_right += fc[a] * mc[a].mean;
                    score_right += fc[a] * mc[a].compute_score();
                    weight_right += fc[a];
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

            float ratio = (float) cap / (float) weight;
            mc[new_arc].hidden_score = score * ratio;
            mc[new_arc].weight = 0;
            mc[new_arc].mean = mean * ratio;

        } else if (rev_in_groups_ev[lg]->block || rev_out_groups_ev[rg]->block) {

            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("LG " + std::to_string(lg) + " RG " + std::to_string(rg) + "Single Block \n");
            #endif
            
            ListDigraph::Arc base_arc, new_arc;

            if (rev_in_groups_ev[lg]->block) {
                base_arc = wc.arcFromId(*rev_out_groups[rg].begin());
                new_arc = wc.addArc(new_group_nodes_left[lg], wc.target(base_arc));
            } else { // rev_out_groups_ev[rg].block
                base_arc = wc.arcFromId(*rev_in_groups[lg].begin());
                new_arc = wc.addArc(wc.source(base_arc), new_group_nodes_right[rg]);
            }

            // just copy over stuff   
            ces[new_arc] = ces[base_arc];
            cet[new_arc] = cet[base_arc];
            cel[new_arc] = cel[base_arc];
            fc[new_arc] = cap;
            mc[new_arc] = mc[base_arc];
            mc[new_arc].reduce(cap/ float(fc[base_arc]));

            fc[base_arc] -= cap;
            mc[base_arc].reduce((fc[base_arc] - cap)/ float(fc[base_arc]));

            if (unsecurityArc[base_arc].has_ref()) {
                std::copy(unsecurityArc[base_arc]->begin(), unsecurityArc[base_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
            }

            cycle_id_in[new_arc] = cycle_id_in[base_arc];
            cycle_id_out[new_arc] = cycle_id_out[base_arc];

            if (rev_in_groups_ev[lg]->block) {
                know_paths[new_arc].bridges.ref() = know_paths[base_arc].bridges.ref(); // deep copy
                if (fc[base_arc] == 0) {
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
                if (fc[base_arc] == 0) {
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
            unravel_ILP(*rev_in_groups[lg].begin(), *rev_out_groups[rg].begin(), cap, know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);
        }
    }
    

    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_left.begin(); ngi != new_group_nodes_left.end(); ++ngi) {
        
        ListDigraph::Node n = ngi->second;
        
        float mean = 0;
        float score = 0;
        capacity_type weight = 0;
        capacity_type group_flow = 0;
        for (ListDigraph::InArcIt a(wc, n); a!=INVALID; ++a) {
            group_flow += fc[a];
            if (!barred[a]) {
                mean += fc[a] * mc[a].mean;
                score += fc[a] * mc[a].compute_score();
                weight += fc[a];
            }
        }
        mean = mean / (float) weight;
        score = score / (float) weight;
        
        capacity_type evidenced_flow = 0;
        for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
            evidenced_flow += fc[a];
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Group Node Left " + std::to_string(wc.id(n)) + " GF " + std::to_string(group_flow) + " EF " + std::to_string(evidenced_flow) + "\n");
        #endif

        if (group_flow > evidenced_flow) {
            ListDigraph::Arc new_arc = wc.addArc(n, node);
            cet[new_arc] = edge_types::RESOLVE_HELPER;
            fc[new_arc] = group_flow - evidenced_flow;
            ces[new_arc] = exon_edge(size);
            
            float ratio = (float) evidenced_flow / (float) group_flow;
            mc[new_arc].hidden_score = score * ratio;
            mc[new_arc].weight = 0;
            mc[new_arc].mean = mean * ratio;
        }
    }
  
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_right.begin(); ngi != new_group_nodes_right.end(); ++ngi) {
        
        ListDigraph::Node n = ngi->second;
        
        float mean = 0;
        float score = 0;
        capacity_type weight = 0;
        capacity_type group_flow = 0;
        for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
            group_flow += fc[a];
            if (!barred[a]) {
                mean += fc[a] * mc[a].mean;
                score += fc[a] * mc[a].compute_score();
                weight += fc[a];
            }
        }
        mean = mean / (float) weight;
        score = score / (float) weight;
        
        capacity_type evidenced_flow = 0;
        for (ListDigraph::InArcIt a(wc, n); a!=INVALID; ++a) {
            evidenced_flow += fc[a];
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Group Node Left " + std::to_string(wc.id(n)) + " GF " + std::to_string(group_flow) + " EF " + std::to_string(evidenced_flow) + "\n");
        #endif

        if (group_flow > evidenced_flow) {
            ListDigraph::Arc new_arc = wc.addArc(node, n);
            cet[new_arc] = edge_types::RESOLVE_HELPER;
            fc[new_arc] = group_flow - evidenced_flow;
            ces[new_arc] = exon_edge(size);
            
            float ratio = (float) evidenced_flow / (float) group_flow;
            mc[new_arc].hidden_score = score * ratio;
            mc[new_arc].weight = 0;
            mc[new_arc].mean = mean * ratio;
        }
    }
    
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_left.begin(); ngi != new_group_nodes_left.end(); ++ngi) {
        clean_barred_leftovers(ngi->second, know_paths, know_back_paths, barred, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
    }
    for( std::map<int, ListDigraph::Node >::iterator ngi = new_group_nodes_right.begin(); ngi != new_group_nodes_right.end(); ++ngi) {
        clean_barred_leftovers(ngi->second, know_paths, know_back_paths, barred, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
    }
    
    // we now look at possibly unevidenced leftovers
    // remove all 0 edges here
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            for (arc_bridge::iterator bi = know_paths[arc].begin(); bi != know_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(bi->first);
                know_back_paths[a].remove_id(wc.id(arc));
                
            }
            know_paths[arc].clear();
            wc.erase(arc);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            for (arc_back_bridge::iterator bi = know_back_paths[arc].begin(); bi != know_back_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(*bi);
                know_paths[a].remove_id(wc.id(arc));
            }  
            know_back_paths[arc].clear();
            wc.erase(arc);
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
        unravel_unevidenced_leftovers(node, know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, NULL);
        indeg = 0;
    }
    if (indeg == 0) {
        // we completely resolved this one, remove!
        wc.erase(node);
    }
}

void base_manager::unravel_single(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId) {
    
    int left, right;
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ++a) {
        if (know_paths[a].size() == 1 && know_back_paths[wc.arcFromId(know_paths[a].begin()->first)].size() == 1 ) {
        
            left = wc.id(a);
            right = know_paths[a].begin()->first; 

        }
    }

    unravel_ILP(left, right, std::min(fc[wc.arcFromId(left)], fc[wc.arcFromId(right)]), know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, transcript_unsecurity::EVIDENCED);    
    
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            for (arc_bridge::iterator bi = know_paths[arc].begin(); bi != know_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(bi->first);
                know_back_paths[a].remove_id(wc.id(arc));
            }
            know_paths[arc].clear();
            wc.erase(arc);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            for (arc_back_bridge::iterator bi = know_back_paths[arc].begin(); bi != know_back_paths[arc].end(); ++bi) {
                ListDigraph::Arc a = wc.arcFromId(*bi);
                know_paths[a].remove_id(wc.id(arc));
            }  
            know_back_paths[arc].clear();
            wc.erase(arc);
        }
    }
}

void base_manager::unravel_evidences(ListDigraph::Node node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc,  ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id> &unsecurityId,
            ListDigraph::NodeMap<bool>* resolved) {
    
    
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
        if (cet[a] == edge_types::HELPER) {
            // no resolving 
            continue;
        }
        arcs.push_back(info(a, false, fc[a]));
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID; ++a) {
        if (cet[a] == edge_types::HELPER) {
            // no resolving 
            continue;
        }
        arcs.push_back(info(a, true, fc[a]));
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
                capacity_type min = std::min(it->cap, fc[ wc.arcFromId(know_paths[it->arc].begin()->first) ]);
                if (min > max ){
                    candidate = *it;
                    max = min;
                }
            } else if (it->is_out && know_back_paths[it->arc].size() == 1) {
                capacity_type min = std::min(it->cap, fc[ wc.arcFromId(*know_back_paths[it->arc].begin()) ]);
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
            unravel_evidence_path_right(*know_back_paths[a].begin(), wc.id(a), know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, transcript_unsecurity::RESOLVED, resolved);
        } else {
            unravel_evidence_path_left(wc.id(a), know_paths[a].begin()->first, know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, transcript_unsecurity::RESOLVED, resolved);
        }   
    }
    
    // we now look at possibly unevidenced leftovers
    // remove all 0 edges here
    for (ListDigraph::InArcIt a(wc, node); a!=INVALID; ) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            wc.erase(arc);
        }
    }
    for (ListDigraph::OutArcIt a(wc, node); a!=INVALID;) {
        ListDigraph::Arc arc(a);
        ++a;
        if (fc[arc] == 0) {
            wc.erase(arc);
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
        unravel_unevidenced_leftovers(node, know_paths, know_back_paths, wc, fc, mc, ces, cet, cel, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, NULL);
        indeg = 0;
    }
    if (indeg == 0) {
        // we completely resolved this one, remove!
        wc.erase(node);
    }
}


void base_manager::unravel_evidence_path_left(int left, int right,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, transcript_unsecurity::evidence evidence,
            ListDigraph::NodeMap<bool>* resolved) {
    
    ListDigraph::Arc left_arc = wc.arcFromId(left);
    ListDigraph::Arc right_arc = wc.arcFromId(right);
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Extract Left. " + ces[left_arc].to_string() + ";" + ces[right_arc].to_string() + " - " +  std::to_string(left) + " " + std::to_string(right) + "\n");
    #endif
    
    // we get minimal capacity
    capacity_type cap_evidence = fc[left_arc];
    if (fc[right_arc] < cap_evidence) {
        cap_evidence = fc[right_arc];
    }
    
    capacity_mean c1 = mc[left_arc];
    capacity_mean c2 = mc[right_arc];
    
    c1.reduce(cap_evidence/ float(fc[left_arc]));
    c2.reduce(cap_evidence/ float(fc[right_arc]));
    mc[left_arc].reduce((fc[left_arc] - cap_evidence)/ float(fc[left_arc]));
    mc[right_arc].reduce((fc[right_arc] - cap_evidence)/ float(fc[right_arc]));
    c1.update(c2);
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Capacity " + std::to_string(cap_evidence) + "\n");
    #endif
    
    fc[left_arc] -= cap_evidence;
    fc[right_arc] -= cap_evidence;
    
    // this should never be called on a zero edge, so add in new arcs
    ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));
    if (resolved != NULL) {
        (*resolved)[wc.source(left_arc)] = false;
        (*resolved)[wc.target(right_arc)] = false;
    }
    
    fc[new_arc] = cap_evidence;
    cet[new_arc] = edge_types::EXON;
    mc[new_arc] = c1;
    
    if (cet[left_arc] != edge_types::EXON && cet[right_arc] != edge_types::EXON) {
        cet[new_arc] = edge_types::HELPER;
        if (ces[right_arc].node_index == -1) {
            ces[new_arc] = ces[left_arc];
        } else {
            ces[new_arc] = ces[right_arc];
        }
    } else if (cet[left_arc] == edge_types::HELPER) {
        // leftarc is a helper
        ces[new_arc] = exon_edge(ces[right_arc]); 
        
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else if (cet[right_arc] == edge_types::HELPER) {
        // rightarc is a helper
        ces[new_arc] = exon_edge(ces[left_arc]);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else {
        ces[new_arc] = exon_edge(ces[left_arc]);
        ces[new_arc].join_edge(ces[right_arc]);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->begin()));
        }
        unsecurityArc[new_arc]->insert( transcript_unsecurity( unsecurityId[wc.target(left_arc)].id, evidence));
    }
     
    cel[new_arc].first_exon = cel[left_arc].first_exon;
    cel[new_arc].middle = cel[left_arc].middle + cel[left_arc].last_exon + cel[right_arc].middle;
    cel[new_arc].last_exon = cel[right_arc].last_exon;
    
    cycle_id_in[new_arc] = cycle_id_in[left_arc];
    cycle_id_out[new_arc] = cycle_id_out[right_arc];
    
    // we now need to transfer the know paths of the right and left side
    know_paths[new_arc].bridges.ref() = know_paths[right_arc].bridges.ref(); // deep copy
    know_back_paths[new_arc].bridges.ref() = know_back_paths[left_arc].bridges.ref(); // deep copy
    
    // we clear the left arc in any case, as it was resolved
        know_paths[left_arc].clear();
    // but evidence can vary
    if (fc[left_arc] == 0) {
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
    if (fc[right_arc] == 0) {
        
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
    
//                    digraphWriter(wc, std::cout)
//                .arcMap("edge_specifier", ces)
//                .arcMap("edge_type", cet)
//                .arcMap("flow", fc)
//                .arcMap("length", cel)
//                .arcMap("means", mc)
//                .run();   
                    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Exit Left.\n");
    #endif
}

void base_manager::unravel_evidence_path_right(int left, int right,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths, 
            ListDigraph &wc, 
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, transcript_unsecurity::evidence evidence,
            ListDigraph::NodeMap<bool>* resolved) {
    
        ListDigraph::Arc left_arc = wc.arcFromId(left);
        ListDigraph::Arc right_arc = wc.arcFromId(right);

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Extract Right. " + ces[left_arc].to_string() + ";" + ces[right_arc].to_string() + " - " +  std::to_string(left) + " " + std::to_string(right) + "\n");
        #endif

        // we get minimal capacity
        capacity_type cap_evidence = fc[left_arc];
        if (fc[right_arc] < cap_evidence) {
            cap_evidence = fc[right_arc];
        }
    
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Capacity " + std::to_string(cap_evidence) + "\n");
    #endif

    capacity_mean c1 = mc[left_arc];
    capacity_mean c2 = mc[right_arc];
    
    c1.reduce(cap_evidence/ float(fc[left_arc]));
    c2.reduce(cap_evidence/ float(fc[right_arc]));
    mc[left_arc].reduce((fc[left_arc] - cap_evidence)/ float(fc[left_arc]));
    mc[right_arc].reduce((fc[right_arc] - cap_evidence)/ float(fc[right_arc]));
    c1.update(c2);
    
    fc[left_arc] -= cap_evidence;
    fc[right_arc] -= cap_evidence;
    
    // this should never be called on a zero edge, so add in new arcs
    ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));
    if (resolved != NULL) {
        (*resolved)[wc.source(left_arc)] = false;
        (*resolved)[wc.target(right_arc)] = false;
    }
    
    fc[new_arc] = cap_evidence;
    cet[new_arc] = edge_types::EXON;
    mc[new_arc] = c1;
    
    if (cet[left_arc] != edge_types::EXON && cet[right_arc] != edge_types::EXON) {
        cet[new_arc] = edge_types::HELPER;
        if (ces[right_arc].node_index == -1) {
            ces[new_arc] = ces[left_arc];
        } else {
            ces[new_arc] = ces[right_arc];
        }
    } else if (cet[left_arc] == edge_types::HELPER) {
        // leftarc is a helper
        ces[new_arc] = exon_edge(ces[right_arc]);
        
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else if (cet[right_arc] == edge_types::HELPER) {
        // rightarc is a helper
        ces[new_arc] = exon_edge(ces[left_arc]); 
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        
    } else {
        ces[new_arc] = exon_edge(ces[left_arc]);
        ces[new_arc].join_edge(ces[right_arc]);
        
        if (unsecurityArc[left_arc].has_ref()) {
            std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
        }
        if (unsecurityArc[right_arc].has_ref()) {
            std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->begin()));
        }
        unsecurityArc[new_arc]->insert( transcript_unsecurity(unsecurityId[wc.target(left_arc)].id, evidence));
    }
        
    cel[new_arc].first_exon = cel[left_arc].first_exon;
    cel[new_arc].middle = cel[left_arc].middle + cel[left_arc].last_exon + cel[right_arc].middle;
    cel[new_arc].last_exon = cel[right_arc].last_exon;
    
    cycle_id_in[new_arc] = cycle_id_in[left_arc];
    cycle_id_out[new_arc] = cycle_id_out[right_arc];
    
    // we now need to transfer the know paths of the right and left side
    know_paths[new_arc].bridges.ref() = know_paths[right_arc].bridges.ref(); // deep copy
    know_back_paths[new_arc].bridges.ref() = know_back_paths[left_arc].bridges.ref(); // deep copy
    
    // we clear the right arc in any case, as it was resolved
    know_back_paths[right_arc].clear();
    // but evidence can vary
    if (fc[right_arc] == 0) {
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
    if (fc[left_arc] == 0) {
        
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
    
//                    digraphWriter(wc, std::cout)
//                .arcMap("edge_specifier", ces)
//                .arcMap("edge_type", cet)
//                .arcMap("flow", fc)
//                .arcMap("length", cel)
//                .arcMap("means", mc)
//                .run();  
                    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Exit Right.\n");
    #endif
}

void base_manager::unravel_unevidenced_leftovers(ListDigraph::Node &node,
            ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
            ListDigraph &wc,
            ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<capacity_mean> &mc, ListDigraph::ArcMap<exon_edge> &ces,
            ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
            ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
            ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
            ListDigraph::NodeMap< unsecurity_id > &unsecurityId, ListDigraph::ArcMap<bool> *barred) {
    
    std::unordered_set<int> left, right; 
    for (ListDigraph::InArcIt left_arc(wc, node); left_arc!=INVALID; ++left_arc) {
        for (ListDigraph::OutArcIt right_arc(wc, node); right_arc!=INVALID; ++right_arc) {

            left.insert(wc.id(left_arc));
            right.insert(wc.id(right_arc));

            ListDigraph::Arc new_arc = wc.addArc(wc.source(left_arc), wc.target(right_arc));

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Join leftover " + std::to_string(wc.id(left_arc)) + " " + std::to_string(wc.id(right_arc)) + " to " + std::to_string(wc.id(new_arc)) + "\n");
            #endif
            
            capacity_type capacity = std::min(fc[left_arc], fc[right_arc]);
            
            capacity_mean c1 = mc[left_arc];
            capacity_mean c2 = mc[right_arc];

            c1.reduce(capacity/ float(fc[left_arc]));
            c2.reduce(capacity/ float(fc[right_arc]));
            mc[left_arc].reduce((fc[left_arc] - capacity)/ float(fc[left_arc]));
            mc[right_arc].reduce((fc[right_arc] - capacity)/ float(fc[right_arc]));
            c1.update(c2);
            
            
            fc[new_arc] =  capacity;
            cet[new_arc] = edge_types::EXON;
            mc[new_arc] = c1;
            if (barred != NULL) {
                (*barred)[new_arc] = (*barred)[left_arc] || (*barred)[right_arc];
            }
                    
            if (cet[left_arc] != edge_types::EXON && cet[right_arc] != edge_types::EXON) {
                cet[new_arc] = edge_types::HELPER;
                if (ces[right_arc].node_index == -1) {
                    ces[new_arc] = ces[left_arc];
                } else {
                    ces[new_arc] = ces[right_arc];
                }
            } else if (cet[left_arc] == edge_types::HELPER) {
                // leftarc is a helper
                if (unsecurityArc[right_arc].has_ref()) {
                    std::copy(unsecurityArc[right_arc]->begin(), unsecurityArc[right_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
                }
                ces[new_arc] = exon_edge(ces[right_arc]);
                
            } else if (cet[right_arc] == edge_types::HELPER) {
                // rightarc is a helper
                if (unsecurityArc[left_arc].has_ref()) {
                    std::copy(unsecurityArc[left_arc]->begin(), unsecurityArc[left_arc]->end(), std::inserter(unsecurityArc[new_arc].ref(), unsecurityArc[new_arc]->end()));
                }
                ces[new_arc] = exon_edge(ces[left_arc]);
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
                
                ces[new_arc] = exon_edge(ces[left_arc]);
                ces[new_arc].join_edge(ces[right_arc]);
            }
            
            cel[new_arc].first_exon = cel[left_arc].first_exon;
            cel[new_arc].middle = cel[left_arc].middle + cel[left_arc].last_exon + cel[right_arc].middle;
            cel[new_arc].last_exon = cel[right_arc].last_exon;

            cycle_id_in[new_arc] = cycle_id_in[left_arc];
            cycle_id_out[new_arc] = cycle_id_out[right_arc];

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
            
            
//                    digraphWriter(wc, std::cout)
//                .arcMap("edge_specifier", ces)
//                .arcMap("edge_type", cet)
//                .arcMap("flow", fc)
//                .arcMap("length", cel)
//                .arcMap("means", mc)
//                .run();  
        }
    }

    for ( std::unordered_set<int>::iterator it = left.begin(); it!= left.end(); ++it) {
       
        for (ListDigraph::InArcIt a(wc, wc.source(wc.arcFromId(*it))); a!=INVALID; ++a) {
            know_paths[a].remove_id(*it);
        }
        
        know_paths[wc.arcFromId(*it)].clear();
        know_back_paths[wc.arcFromId(*it)].clear();
        wc.erase(wc.arcFromId(*it));
    }
    for ( std::unordered_set<int>::iterator it = right.begin(); it!= right.end(); ++it) {
        
        for (ListDigraph::OutArcIt a(wc, wc.target(wc.arcFromId(*it))); a!=INVALID; ++a) {
            know_back_paths[a].remove_id(*it);
        }
        
        know_paths[wc.arcFromId(*it)].clear();
        know_back_paths[wc.arcFromId(*it)].clear();
        wc.erase(wc.arcFromId(*it));
    }
    
}

void base_manager::resolve_overarching(path_evidence_map<int, path_evidence > &evidences,
            path_evidence_set<int> &pre_path_match, path_evidence_set<int> &post_path_match,
            path_evidence_map<int, rpos> &post_path_min_length, path_evidence_map<int, rcount> &post_path_count, rpos min_right_hit,
            ListDigraph &wc, ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::Node n, bool allow_block) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Resolve Overarch \n");
    #endif

    // first we create an evidence set to avoid recomputation in inner loop;
    // we already know the supported arcs, so we need to collect blocked only
    
    path_evidence_set<int> post_path_blocked;
    for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
       if (cet[a] == edge_types::HELPER) {
            // helper can't possibly be disproven/proven
            continue;
       } 
       int id = wc.id(a);
       if ( post_path_match.find(id) == post_path_match.end()) {
           // meaning this arc was not detected in any way
           // we test length if we can reject it
           
           path_evidence_map<int, rpos>::iterator ml = post_path_min_length.find(id);
           if (ml == post_path_min_length.end() || min_right_hit < ml->second) {
               // we can reject this, as we SHOULD have seen a hit otherwise
               post_path_blocked.insert(id);
           }
       }
    }
    
    // we have a blocked list and a match list, that we now need to incorporate to ALL found left arcs
    for (path_evidence_set<int>::iterator li = pre_path_match.begin(); li != pre_path_match.end(); ++li) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Pre Id " + std::to_string(*li)  + "\n");
        #endif
        for (path_evidence_set<int>::iterator mi = post_path_match.begin(); mi != post_path_match.end(); ++mi ) {
            evidences[*li].add_evidence(*mi, post_path_count[*mi]);
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Evidence Id " + std::to_string(*mi) + "\n");
            #endif
        }
        if (allow_block) {
            for (path_evidence_set<int>::iterator bi = post_path_blocked.begin(); bi != post_path_blocked.end(); ++bi ) {
                evidences[*li].add_blocked(*bi);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Block Id " + std::to_string(*bi) + "\n");
                #endif
            }
        }
    }
}


void base_manager::resolve_pair(path_evidence_map<int, path_evidence > &evidences,
            path_evidence_set<int> &pre_path_match,  path_evidence_map<int, rpos> &hit_distance, 
            path_evidence_set<int> &post_path_match, path_evidence_map<int, rpos> &post_path_min_length,
            path_evidence_map<int, std::pair< std::pair< rpos, rpos > ,rpos> > &post_path_range,
            path_evidence_map<int, rcount> &post_path_count,
            ListDigraph &wc, ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::Node n, bool allow_block) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Resolve Pair \n");
    #endif
    
    // evidence depends on pair distance, so we need to look all right arcs for all left individually
    // avoid constant hashing at least
    std::vector<std::tuple< int, bool, rpos, rpos, rpos, rpos, rcount> > arcs;
    for (ListDigraph::OutArcIt a(wc, n); a!=INVALID; ++a) {
       if (cet[a] == edge_types::HELPER) {
            // helper can't possibly be disproven/proven
            continue;
       } 
       int id = wc.id(a);
       rpos lowesthit = 0;
       rpos lowest_hit_read_length = 0;
       rpos highesthit = 0;
       rpos min_found_range = 0;
       rcount count = 0;
       bool hit = false;
       
       if ( post_path_match.find(id) != post_path_match.end()) {
            // hit, so we have the info
           lowesthit = post_path_range[id].first.first;
           lowest_hit_read_length = post_path_range[id].first.second;
           highesthit = post_path_range[id].second; 
           count = post_path_count[id];
           hit = true;
       }
       
       path_evidence_map<int, rpos>::iterator pi = post_path_min_length.find(id);
       if (pi == post_path_min_length.end()) {
           min_found_range = std::numeric_limits<rpos>::max(); // I know, I hate myself for this, but branch prediction rejoices...
       } else {
           min_found_range = pi->second;
       }
       
       arcs.push_back(std::make_tuple(id, hit, lowesthit, lowest_hit_read_length, highesthit, min_found_range, count));
    }

    for (path_evidence_set<int>::iterator li = pre_path_match.begin(); li != pre_path_match.end(); ++li) {
        
        rpos hd = hit_distance[*li];
        
//       logger::Instance()->debug("Pre Id " + (*li)->identifier.to_string()  + "\n");     
        rpos critical_length = 0;
        
        // first run, we need the max min hit length that supports the arcs
        for (std::vector<std::tuple< int, bool, rpos, rpos, rpos, rpos, rcount> >::iterator ai = arcs.begin(); ai!= arcs.end(); ++ai ) {
//             logger::Instance()->debug("In\n");
            if (std::get<1>(*ai)) {
//                logger::Instance()->debug("hit\n");
                // we have a hit, we have a range!
                rpos low = std::get<2>(*ai);
                rpos high = std::get<4>(*ai);
                
                if (low + hd  <= options::Instance()->get_max_paired_distance() && options::Instance()->get_min_paired_distance() <= high + hd) {
                    
                    // this means we have an actual valid hit on this arc
                    if (low + hd <  options::Instance()->get_min_paired_distance()) {
                        low = options::Instance()->get_min_paired_distance() - hd;
                    }
                    
                    if (critical_length < low + std::get<3>(*ai)) {
                        critical_length = low + std::get<3>(*ai);
                    }
                }
            }
        }
        
//        logger::Instance()->debug("Crit Len " + std::to_string(critical_length)  + "\n");
        
        // second run with know critical length for this left path
        for (std::vector<std::tuple< int, bool, rpos, rpos, rpos, rpos, rcount> >::iterator ai = arcs.begin(); ai!= arcs.end(); ++ai ) {
            if (std::get<1>(*ai)) { // potential hit
                // we have a hit, we have a range!
                rpos low = std::get<2>(*ai);
                rpos high = std::get<4>(*ai);
                
//                logger::Instance()->debug("Len Found a) " + std::to_string(low) + " " + std::to_string(high) + "\n");
                
                if (low + hd  <= options::Instance()->get_max_paired_distance() && options::Instance()->get_min_paired_distance() <= high + hd) {
                    // this means we have an actual valid hit on this arc, double computation but not expensive
                    evidences[*li].add_evidence(std::get<0>(*ai), std::get<6>(*ai));
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Evidence Id " + std::to_string(std::get<0>(*ai)) + "\n");
                    #endif
                } else if (allow_block){
                    // invalid hit that was long enough though, block!
                    evidences[*li].add_blocked(std::get<0>(*ai));
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Block Id " + std::to_string(std::get<0>(*ai)) + "\n");
                    #endif
                }
            } else {
                  #ifdef ALLOW_DEBUG
//                logger::Instance()->debug("Len Found b) " + std::to_string(std::get<5>(*ai))  + "\n");
                  #endif
                if (critical_length <= std::get<5>(*ai) && allow_block) {
                    // this path is long enough to block it!
                     evidences[*li].add_blocked(std::get<0>(*ai));
                     #ifdef ALLOW_DEBUG
                     logger::Instance()->debug("Block Id " + std::to_string(std::get<0>(*ai)) + "\n");
                     #endif
                }
            }
        }
    }
}

void base_manager::extend_path_left(std::deque<path>& paths, path* p, unsigned int& border,
        ListDigraph& wc, ListDigraph::ArcMap<exon_edge>& ces, ListDigraph::ArcMap<edge_types::edge_type> &cet,
	ListDigraph::NodeMap<unsigned int>& cni, ListDigraph::ArcMap<arc_bridge> &know_paths) {
        
    
    if ( paths.size() <= options::Instance()->get_max_enumerated_paths() && border < p->border_index )  {
        
//        #ifdef ALLOW_DEBUG
//        logger::Instance()->debug("Extend Left.\n");
//        #endif
        
        // extend path at last path node
        
        // we create a parent arc for all current
        ListDigraph::Arc arc_save = INVALID;
	bool first = false;
        for (ListDigraph::InArcIt a(wc, p->last_node) ; a!=INVALID; ++a) {
            
            if (cet[a] == edge_types::HELPER || edge_type[a] == edge_types::RESOLVE_HELPER ||
                    (know_paths[a].is_evidenced_path() && !know_paths[a].has_path_evidence(wc.id(p->last_arc))) ) {
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
            np->identifier.join_edge(ces[a]); 
            
            extend_path_left(paths, np, border, wc, ces, cet, cni, know_paths);
        }
        
        if ( arc_save==INVALID) {
            return;
        }
        
        p->last_node = wc.source(arc_save);
        p->last_arc = arc_save;
        p->border_index = cni[p->last_node];
        p->identifier.join_edge(ces[arc_save]); 
        
        extend_path_left(paths, p, border, wc, ces, cet, cni, know_paths);
    }
}



void base_manager::extend_path_right(std::deque<path>& paths, path* p, unsigned int& border,
        ListDigraph& wc, ListDigraph::ArcMap<exon_edge>& ces, ListDigraph::ArcMap<edge_types::edge_type> &cet,
	ListDigraph::NodeMap<unsigned int>& cni, ListDigraph::ArcMap<arc_bridge> &know_paths) {
    
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
            
            if (cet[a] == edge_types::HELPER || edge_type[a] == edge_types::RESOLVE_HELPER ||
                  (know_paths[p->last_arc].is_evidenced_path() && !know_paths[p->last_arc].has_path_evidence(wc.id(a))) ) {
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
            np->identifier.join_edge(ces[a]);
            
            extend_path_right(paths, np, border, wc, ces, cet, cni, know_paths);
        }
        
        if ( arc_save==INVALID) {
            return;
        }
        
        //now change path itself
        p->last_node = wc.target(arc_save);
        p->last_arc = arc_save;
        p->border_index = cni[p->last_node];
        p->identifier.join_edge(ces[arc_save]);

        extend_path_right(paths, p, border, wc, ces, cet, cni, know_paths);
    }
    
}

void base_manager::print_graph_debug_copy(std::ostream& os,
    ListDigraph &wc,
    ListDigraph::ArcMap<capacity_type> &fc,
    ListDigraph::ArcMap<exon_edge> &ces,
    ListDigraph::ArcMap<edge_types::edge_type> &cet,
    ListDigraph::ArcMap<edge_length> &cel,    
    ListDigraph::Node &cs, ListDigraph::Node &ct) {

    if (options::Instance()->is_debug()) {
    
        digraphWriter(wc, os)
                .arcMap("edge_specifier", ces)
                .arcMap("edge_type", cet)
                .arcMap("flow", fc)
                .arcMap("length", cel)
                .node("source", cs)
                .node("drain", ct)
                .run();  
    }
}


void base_manager::print_graph_gs(std::ostream& os,
    ListDigraph &wc,
    ListDigraph::ArcMap<capacity_type> &fc,
    ListDigraph::ArcMap<exon_edge> &ces,
    ListDigraph::ArcMap<edge_types::edge_type> &cet,   
    ListDigraph::Node &cs, ListDigraph::Node &ct) {

    digraphWriter(wc, os)
            .arcMap("Edge", ces)
            .arcMap("Type", cet)
            .arcMap("Flow/Capacity", fc)
            .node("source", cs)
            .node("drain", ct)
            .run();  
}
