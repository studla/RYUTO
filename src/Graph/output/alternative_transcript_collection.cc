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

void alternative_transcript_collection::finalize_borders(exon_meta* meta) {
     for (graph_list<lazy<transcript> >::iterator it = transcripts.begin(); it!= transcripts.end(); ++it) {
        (*it)->finalize_borders(meta);
     } 
}

void alternative_transcript_collection::filter_transcripts(int id) {
    
    
     #ifdef ALLOW_DEBUG
      logger::Instance()->debug("Filter Set " + std::to_string(id) + "\n");
     #endif
    
    // this requires FINALIZED transcripts 
    
    std::list<std::pair<rpos, rpos> > regions;
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
    
    graph_list<lazy<transcript> > keep;
    for (std::list<std::pair<rpos, rpos>>::iterator reg_it = regions.begin(); reg_it != regions.end(); ++reg_it) {
    
        
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
    
    // now the nested filter on the rest!
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
