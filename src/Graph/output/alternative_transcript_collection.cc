/* 
 * File:   alternative_transcript_collection.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 13, 2016, 9:51 AM
 */

#include "alternative_transcript_collection.h"
#include "../../Datatype_Templates/move.h"
#include "../../Options/options.h"
#include "../../Datatype_Templates/reader_list.h"
#include "../../Datatype_Templates/graph_list.h"

alternative_transcript_collection::alternative_transcript_collection() {
}


alternative_transcript_collection::~alternative_transcript_collection() {
}

bool alternative_transcript_collection::empty() {
    return transcripts.empty();
}

int alternative_transcript_collection::size() {
    return transcripts.size();
}

void alternative_transcript_collection::join(alternative_transcript_collection &t) {
     _MOVE_RANGE(t.transcripts.begin(), t.transcripts.end(), std::back_inserter(transcripts));
     for (std::map<int, capacity_type>::iterator tfi = t.total_flow_error.begin(); tfi != t.total_flow_error.end(); ++tfi) {
         total_flow_error[tfi->first] += tfi->second;
     }
     _MOVE_RANGE(t.flow_error.begin(), t.flow_error.end(), std::back_inserter(flow_error));
}

void alternative_transcript_collection::print(std::ostream& os) {
    
    // header
//    os << "Cycle In\tCycle Out\tEdge\tFlow\tUnsecurities\n";
//    for (graph_list<lazy<transcript> >::iterator it = transcripts.begin(); it!= transcripts.end(); ++it) {
//        if (it->ref().flow > options::Instance()->get_capacity_filter()) { 
//            it->ref().print(os);
//        }
//    } 
    
}

void alternative_transcript_collection::print_gtf(std::ostream &os, std::string &gene_id) {
    
    unsigned int transcript_id = 1;
    for (graph_list<lazy<transcript> >::iterator it = transcripts.begin(); it!= transcripts.end(); ++it) {
        (*it)->print_gtf_entry(os, gene_id, transcript_id, input_main_id);
        ++transcript_id;
    }
}

void alternative_transcript_collection::finalize_borders(exon_meta* me) {

    for (graph_list<lazy<transcript> >::iterator it = transcripts.begin(); it!= transcripts.end(); ++it) {
        (*it)->finalize_borders(me);
    } 
    
    for (std::map<unsigned int, std::map<int, capacity_type> >::iterator fi = exon_flow_error.begin(); fi!= exon_flow_error.end(); ++fi) {
        unsigned int pos = fi->first;
        
        flow_error.push_back(lazy<flow_error_t>());
        flow_error.back()->is_exon = true;
        flow_error.back()->start = me->exons[pos].left;
        flow_error.back()->end = me->exons[pos].right;
        flow_error.back()->errors.insert(fi->second.begin(), fi->second.end());
        flow_error.back()->chromosome = me->chromosome;
        flow_error.back()->strand = me->strand;
        
    }
    
    for (std::map<std::pair<unsigned int, unsigned int>, std::map<int, capacity_type> >::iterator fi = junction_flow_error.begin(); fi!= junction_flow_error.end(); ++fi) {
        unsigned int p1 = fi->first.first;
        unsigned int p2 = fi->first.second;
        
        flow_error.push_back(lazy<flow_error_t>());
        flow_error.back()->is_exon = false;
        flow_error.back()->start = me->exons[p1].right;
        flow_error.back()->end = me->exons[p2].left;
        flow_error.back()->errors.insert(fi->second.begin(), fi->second.end());
        flow_error.back()->chromosome = me->chromosome;
        flow_error.back()->strand = me->strand;
    } 
     
}

void alternative_transcript_collection::print_count_matrix(std::ostream &os, std::string &gene_id, std::set<int> &ids) {
    
    unsigned int transcript_id = 1;
    for (graph_list<lazy<transcript> >::iterator it = transcripts.begin(); it!= transcripts.end(); ++it) {
        (*it)->print_count_matrix_entry(os, gene_id, transcript_id, input_main_id, ids);
        ++transcript_id;
    }
}

void alternative_transcript_collection::print_error_matrix(std::ostream &os, std::string &gene_id, std::set<int> &ids) {
    
    for (graph_list<lazy<flow_error_t> >::iterator it = flow_error.begin(); it!= flow_error.end(); ++it) {
        os << (*it)->chromosome << "\t" << (*it)->strand << "\t" << std::to_string((*it)->start) << "\t" << std::to_string((*it)->end) << "\t";
        if ((*it)->is_exon) {
            os << "Exon";
        } else {
            os << "Junction";
        }
        for(std::set<int>::iterator ii = ids.begin(); ii!= ids.end(); ++ii) {
            if (*ii == -1) {
                continue;
            }
            
            os << "\t" << std::to_string((*it)->errors[*ii]);
        }
        os << "\n";
    }
    
}

void alternative_transcript_collection::compute_region(std::list<std::pair<rpos, rpos> > &regions) {

    for (graph_list<lazy<transcript> >::iterator it = transcripts.begin(); it!= transcripts.end(); ++it) {
                
        rpos left =  (*it)->exons[0].first;
        rpos right = (*it)->exons[(*it)->exons.size() - 1].second;
                
        if (regions.empty()) {                         
            // we add a new on
            regions.push_back(std::make_pair(left, right)); 
        } else {

            greader_list<std::pair<rpos, rpos>* > matched_regions;
            for (std::list<std::pair<rpos, rpos> >::iterator reg_it = regions.begin(); reg_it != regions.end(); ++reg_it) {
                if (reg_it->first <= right && reg_it->second >= left) {
                    // this is an overlap!
                    matched_regions.push_back(&*reg_it);
                }
            } 

            if (matched_regions.size() == 0) {
                 regions.push_back(std::make_pair(left, right));
            } else if (matched_regions.size() == 1) {
                 if (left < matched_regions.back()->first) matched_regions.back()->first = left;
                 if (right > matched_regions.back()->second) matched_regions.back()->second = right;
            } else { // matched_regions.size() > 1

                greader_list<std::pair<rpos, rpos>*  >::iterator fm = matched_regions.begin();
                greader_list<std::pair<rpos, rpos>*  >::iterator match_it = fm;
                ++match_it;

                if (left < (*fm)->first) (*fm)->first = left;
                if (right > (*fm)->second) (*fm)->second = right;

                for (; match_it != matched_regions.end(); ++match_it) {

                    if ( (*match_it)->first < (*fm)->first) (*fm)->first = (*match_it)->first;
                    if ( (*match_it)->second > (*fm)->second) (*fm)->second = (*match_it)->second;

                    for (std::list<std::pair<rpos, rpos>>::iterator reg_it = regions.begin(); reg_it != regions.end(); ++reg_it) {
                        if (&*reg_it == *match_it) {
                            regions.erase(reg_it);
                            break;
                        }
                    }
                } 
            }    
        }
    }
}


void alternative_transcript_collection::vote(int id, graph_list<lazy<transcript> > &keep, std::list<std::pair<rpos, rpos> > &regions) {

    unsigned int reg_count = 0;
    for (std::list<std::pair<rpos, rpos>>::iterator reg_it = regions.begin(); reg_it != regions.end(); ++reg_it, ++reg_count) {
    
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("GTF Region " + std::to_string(reg_it->first) + " " + std::to_string(reg_it->second) + ".\n");
        #endif
        
        graph_list<lazy<transcript> > regional;

        for (graph_list<lazy<transcript> >::iterator it = transcripts.begin(); it!= transcripts.end(); ++it) {
            
            rpos left =  (*it)->exons[0].first;
            rpos right = (*it)->exons[(*it)->exons.size() - 1].second;
            
            if (left <= reg_it->second && right >= reg_it->first) {
                regional.push_back(*it);
            }
        }
        
        capacity_type uniform_max = 0;
        float mean_max = 0;
        float score_max = 0;
//        std::deque<float> all_scores;

        for (graph_list<lazy<transcript> >::iterator it = regional.begin(); it!= regional.end(); ++it) {
            if (it->ref().series[id].flow > uniform_max) {
                uniform_max = it->ref().series[id].flow;
            }
            if (it->ref().series[id].mean > mean_max) {
                mean_max = it->ref().series[id].mean;
            }
            float sco = it->ref().series[id].score;
            if (sco > score_max) {
                score_max = sco;
            }       
        }
        
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Base_values " + std::to_string(uniform_max) + " " + std::to_string(mean_max) + ".\n");
        #endif
        
        
        for (graph_list<lazy<transcript> >::iterator it = regional.begin(); it!= regional.end(); ++it) {

            (*it)->post_filter_regional_group = reg_count;
            
            rpos length = (*it)->length;
            unsigned int exon_count = (*it)->exons.size();
   
            bool unevidenced = false;
            bool barred = false;
            for(std::deque<transcript_unsecurity>::iterator u_it = it->ref().unsecurity_id.begin(); u_it != it->ref().unsecurity_id.end() ; ++u_it) {
                if (u_it->evidenced == transcript_unsecurity::BARRED) {
                    barred = true;
                }
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Transcript: "+it->ref().found_edge.to_string() + " ec " + std::to_string(exon_count) + " l " + std::to_string(it->ref().length) + " f " + std::to_string(it->ref().series[id].flow) + " m " + std::to_string(it->ref().series[id].mean) + " s " + std::to_string(it->ref().series[id].score) + " bar " + std::to_string(barred) +"\n");
            #endif 
             
            if ( it->ref().guided 
                    ||  
                    (
//                    uniform_max > options::Instance()->get_group_filter_base()    
//                   ( (highlander_mode && it->ref().mean >= options::Instance()->get_group_filter_base())
//                      || 
//                      (!highlander_mode) )
                    it->ref().series[id].mean >= mean_max * options::Instance()->get_percentage_filter() / 100.0
     //               && ( it->ref().score > 8 || it->ref().score > score_median * 0.1 )
//                    && it->ref().flow > options::Instance()->get_capacity_filter()
                    && (it->ref().series[id].mean > options::Instance()->get_mean_filter() || it->ref().series[id].score > options::Instance()->get_mean_filter())
                    && it->ref().series[id].score > options::Instance()->get_minimal_score()
                    && mean_max > options::Instance()->get_group_mean_min()  
                    //&& length >= options::Instance()->get_min_transcript_length()
                    && length >= (options::Instance()->get_min_transcript_length_base() + options::Instance()->get_min_transcript_length_extension() * exon_count)    //
                    && !unevidenced
                    && !barred
                    && (exon_count > 1 || ( it->ref().series[id].mean > options::Instance()->get_min_single_coverage() && length >= options::Instance()->get_min_single_length()) )
                    )) {
                         #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("keep\n");
                        #endif 
                        keep.push_back(*it);
            }
        } 
    }
}


void alternative_transcript_collection::multi_vote(std::set<int> &ids, graph_list<lazy<transcript> > &keep, std::list<std::pair<rpos, rpos> > &regions) {

    unsigned int reg_count = 0;
    for (std::list<std::pair<rpos, rpos>>::iterator reg_it = regions.begin(); reg_it != regions.end(); ++reg_it, ++reg_count) {
            
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("GTF Region " + std::to_string(reg_it->first) + " " + std::to_string(reg_it->second) + ".\n");
        #endif

        graph_list<lazy<transcript> > regional;

        for (graph_list<lazy<transcript> >::iterator it = transcripts.begin(); it!= transcripts.end(); ++it) {

            rpos left =  (*it)->exons[0].first;
            rpos right = (*it)->exons[(*it)->exons.size() - 1].second;

            if (left <= reg_it->second && right >= reg_it->first) {
                regional.push_back(*it);
            }
        }

        std::map<int, capacity_type> uniform_max;
        std::map<int, float> mean_max;
        std::map<int, float> score_max;

        for(std::set<int>::iterator ii = ids.begin(); ii != ids.end(); ++ii) {
            int id = *ii;
        
            uniform_max[id] = 0;
            mean_max[id] = 0;
            score_max[id] = 0;
    //        std::deque<float> all_scores;

            for (graph_list<lazy<transcript> >::iterator it = regional.begin(); it!= regional.end(); ++it) {
                if (it->ref().series[id].flow > uniform_max[id]) {
                    uniform_max[id] = it->ref().series[id].flow;
                }
                if (it->ref().series[id].mean > mean_max[id]) {
                    mean_max[id] = it->ref().series[id].mean;
                }
                float sco = it->ref().series[id].score;
                if (sco > score_max[id]) {
                    score_max[id] = sco;
                }       
            }
        }

        std::set<std::string> ref_genes;
        bool has_unguided = false;
        for (graph_list<lazy<transcript> >::iterator it = regional.begin(); it!= regional.end(); ++it) {
             if (it->ref().guided) {
                 ref_genes.insert(it->ref().guide_gene);
             } else {
                 has_unguided = true;
             }
        }

        if (has_unguided) { // nothing to do if all are refs perfectly
            if (ref_genes.size() == 1) { //
                std::string gene_name = *ref_genes.begin();
                for (graph_list<lazy<transcript> >::iterator it = regional.begin(); it!= regional.end(); ++it) {
                   if (!it->ref().guided) {
                       it->ref().guide_grouped = true;
                       it->ref().guide_gene = gene_name;
                   } 
                } 
            } else {
                for (graph_list<lazy<transcript> >::iterator it = regional.begin(); it!= regional.end(); ++it) {
                   if (it->ref().guided) {
                       it->ref().ignore_guide_gene = true;
                   }
                }
            }
        }

        for (graph_list<lazy<transcript> >::iterator it = regional.begin(); it!= regional.end(); ++it) {
            
            (*it)->post_filter_regional_group = reg_count;
            
              #ifdef ALLOW_DEBUG
               logger::Instance()->debug("Transcript: "+it->ref().found_edge.to_string() + "\n");
               #endif 
            
            unsigned int vote_count = 0;
            for(std::set<int>::iterator ii = ids.begin(); ii != ids.end(); ++ii) {
                int id = *ii;

              #ifdef ALLOW_DEBUG
               logger::Instance()->debug("ID: "+std::to_string(id) + "\n");
               #endif 

                rpos length = (*it)->length;
                unsigned int exon_count = (*it)->exons.size();

                bool unevidenced = false;
                bool barred = false;
                for(std::deque<transcript_unsecurity>::iterator u_it = it->ref().unsecurity_id.begin(); u_it != it->ref().unsecurity_id.end() ; ++u_it) {
                    if (u_it->evidenced == transcript_unsecurity::BARRED) {
                        barred = true;
                    }
                }

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug( " ec " + std::to_string(exon_count) + " l " + std::to_string(it->ref().length) + " f " + std::to_string(it->ref().series[id].flow) + " m " + std::to_string(it->ref().series[id].mean) + " s " + std::to_string(it->ref().series[id].score) + " bar " + std::to_string(barred) +"\n");
                #endif 

                if ( it->ref().guided 
                        ||  
                        (
    //                    uniform_max > options::Instance()->get_group_filter_base()    
    //                   ( (highlander_mode && it->ref().mean >= options::Instance()->get_group_filter_base())
    //                      || 
    //                      (!highlander_mode) )
                        it->ref().series[id].mean >= mean_max[id] * options::Instance()->get_percentage_filter() / 100.0
         //               && ( it->ref().score > 8 || it->ref().score > score_median * 0.1 )
    //                    && it->ref().flow > options::Instance()->get_capacity_filter()
                        && (it->ref().series[id].mean > options::Instance()->get_mean_filter() || it->ref().series[id].score > options::Instance()->get_mean_filter())
                        && it->ref().series[id].score > options::Instance()->get_minimal_score()
                        && mean_max[id] > options::Instance()->get_group_mean_min()  
                        //&& length >= options::Instance()->get_min_transcript_length()
                        && length >= (options::Instance()->get_min_transcript_length_base() + options::Instance()->get_min_transcript_length_extension() * exon_count)    //
                        && !unevidenced
                        && !barred
                        && (exon_count > 1 || ( it->ref().series[id].mean > options::Instance()->get_min_single_coverage() && length >= options::Instance()->get_min_single_length()) )
                        )) {
                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Vote\n");
                            #endif 
                            vote_count += 1;
                }
            }
            if (vote_count * 100 / ids.size() > options::Instance()->get_vote_percentage_low()) {
                keep.push_back(*it);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Keep\n");
                #endif 
            }
        }

        
 
        
    }
    
    
    
}

void alternative_transcript_collection::filter_nested(int id, graph_list<lazy<transcript> > &keep) {

    for (graph_list<lazy<transcript> >::iterator i = keep.begin(); i!= keep.end();) {  // all transcripts
        
        if ((*i)->exons.size() <= 1 || i->ref().guided) { // keep single and guided!
            ++i;
            continue;
        }
        
        bool kill = false;
        for(int k = 1; k < (*i)->exons.size(); k++) {   // look over all introns
            rpos is = (*i)->exons[k-1].second;
            rpos ie = (*i)->exons[k].first;
            
            for (graph_list<lazy<transcript> >::iterator j = keep.begin(); j != keep.end(); ++j) {
                
                if ((*j)->exons.size() <= 1) { // keep single
                    continue;
                }
                
                rpos left = (*j)->exons[0].first;
                rpos right = (*j)->exons[(*j)->exons.size()-1].second;
                if ((*j)->series[id].flow >= (*i)->series[id].flow && (*j)->series[id].score >= (*i)->series[id].score && left > is && right < ie) {
                    kill = true;
                    break;
                }
            }
            if(kill) break;
        }
        if(kill) {
            i = keep.erase(i);
        } else {
            ++i;
        }
    }
}


void alternative_transcript_collection::filter_transcripts(std::set<int> &ids) {

     #ifdef ALLOW_DEBUG
      logger::Instance()->debug("Filter Set All\n");
     #endif
    
    // this requires FINALIZED transcripts 
    
    std::list<std::pair<rpos, rpos> > regions;
    compute_region(regions);
    
    graph_list<lazy<transcript> > keep;
    multi_vote(ids, keep, regions);
    
    // now the nested filter on the rest!
    filter_nested(-1, keep);
    
    std::sort(keep.begin(), keep.end(), []( lazy<transcript> a, lazy<transcript> b) -> bool {return (a->exons.begin()->first == b->exons.begin()->first && a->exons.back().second < b->exons.back().second) || a->exons.begin()->first < b->exons.begin()->first;});

    transcripts = keep; 
}


void alternative_transcript_collection::filter_transcripts(int id) {
    
    
     #ifdef ALLOW_DEBUG
      logger::Instance()->debug("Filter Set " + std::to_string(id) + "\n");
     #endif
    
    // this requires FINALIZED transcripts 
    
    std::list<std::pair<rpos, rpos> > regions;
    compute_region(regions);
    
    graph_list<lazy<transcript> > keep;
    vote(id, keep, regions);
    
    // now the nested filter on the rest!
    filter_nested(id, keep);
    
    std::sort(keep.begin(), keep.end(), []( lazy<transcript> a, lazy<transcript> b) -> bool {return (a->exons.begin()->first == b->exons.begin()->first && a->exons.back().second < b->exons.back().second) || a->exons.begin()->first < b->exons.begin()->first;});
//
//    for (graph_list<lazy<transcript> >::iterator i = keep.begin(); i!= keep.end();) {  // all transcripts
//        
//        if ((*i)->guided) {
//            continue;
//        }
//        
//        unsigned int min_dist = 50;
//        graph_list<lazy<transcript> >::iterator min = keep.end();
//        for (graph_list<lazy<transcript> >::iterator j = keep.begin(); j != keep.end() && (*j)->exons.front().first < (*i)->exons.front().first; ++j) {
//           
//            if ((*j)->exons.back().second > (*i)->exons.front().first) {
//                continue;
//            }
//            
//            unsigned int d = (*i)->exons.front().first - (*j)->exons.back().second < min_dist;
//            if (d < min_dist) {
//                min_dist = d;
//                min = j;
//            }
//        }
//        
//        if (min != keep.end() && (*min)->guided && ( (*min)->exons.size() < 2 || (*i)->exons.size() < 2 ) ) {
//            // join
//            
//            (*min)->join(&(*i).ref());
//            
//            i = keep.erase(i);
//        } else {
//            ++i;
//        }
//    }

    transcripts = keep; 
}
