/* 
 * File:   bam_reader.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on January 20, 2016, 1:13 PM
 */

#include <deque>
#include <string>
#include <vector>
#include <math.h>
#include <queue>
#include <boost/unordered/unordered_map.hpp>
#include <tuple>
#include <stdlib.h>

#include <omp.h>

#include "bam_reader.h"
#include "Datatype_Templates/misc_types.h"
#include "Chromosome/exon.h"
#include "../Datatype_Templates/move.h"
#include "Chromosome/connection_iterator.h"
#include "../Options/options.h"
#include "../Logger/logger.h"
#include "Chromosome/read_collection.h"

bam_reader::bam_reader() {
}


bam_reader::~bam_reader() {
}


void bam_reader::finalize(const std::string &chrom_name) {
        
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Finalize " + chrom_name + ".\n");
    #endif
    
    chromosome *chrom_fwd, *chrom_rev; 
    #pragma omp critical(chromosome_map_lock)
    {
        chrom_fwd = &chromosome_map_fwd[chrom_name];
        if (options::Instance()->is_stranded()) {
            chrom_rev = &chromosome_map_rev[chrom_name];
        }        
    }
    
//    split_independent_components(chrom_fwd);
//    //split_independent_components2(chrom_fwd);
//    if (options::Instance()->is_stranded()) {
//        split_independent_components(chrom_rev);
//        //split_independent_components2(chrom_rev);
//    }
    
    #pragma omp critical(iterator_map_lock)
    {
    #pragma omp critical(chromosome_map_lock)
    {
        if (options::Instance()->is_stranded()) {
            
            rcount total_frags =  chromosome_map_fwd[chrom_name].frag_count + chromosome_map_rev[chrom_name].frag_count;
            rcount total_reads =  chromosome_map_fwd[chrom_name].read_count + chromosome_map_rev[chrom_name].read_count;
            
            chromosome_map_fwd[chrom_name].frag_count = total_frags;
            chromosome_map_rev[chrom_name].frag_count = total_frags;
            chromosome_map_fwd[chrom_name].read_count = total_reads;
            chromosome_map_rev[chrom_name].read_count = total_reads;
            
            iterator_map[chrom_name] = connection_iterator(&chromosome_map_fwd[chrom_name], &chromosome_map_rev[chrom_name]);
            
        } else {
            iterator_map[chrom_name] = connection_iterator(&chromosome_map_fwd[chrom_name]);
        }
    }
    }
}


void bam_reader::split_independent_component(connected *conn, greader_list<connected> &all_connected) {
    
    // we want to split conn into a suitable list of connected!
  
    struct reg_save {
        greader_refsorted_list<exon*> exons;
    };
        
    std::list<reg_save> regions;
        
    for(greader_refsorted_list<raw_atom* >::iterator raw_it = conn->atoms->begin(); raw_it != conn->atoms->end(); ++raw_it) {               
        if (regions.empty()) {

            // we add a new on
            regions.push_back(reg_save()); // we make copies, yeay!
            std::copy((*raw_it)->exons->begin(), (*raw_it)->exons->end(), std::inserter(regions.back().exons, regions.back().exons.end()));
        } else {

            greader_list<reg_save* > matched_regions;
            for (std::list<reg_save>::iterator reg_it = regions.begin(); reg_it != regions.end(); ++reg_it) {
                if (!reg_it->exons.is_disjoint((*raw_it)->exons.ref())) {
                    // this is an overlap!
                    matched_regions.push_back(&*reg_it);
                }
            }

            if (matched_regions.size() == 0) {
                regions.push_back(reg_save()); // we make copies, yeay!
                std::copy((*raw_it)->exons->begin(), (*raw_it)->exons->end(), std::inserter(regions.back().exons, regions.back().exons.end()));
            } else if (matched_regions.size() == 1) {
                std::copy((*raw_it)->exons->begin(), (*raw_it)->exons->end(), std::inserter(matched_regions.back()->exons, matched_regions.back()->exons.begin()));
            } else { // matched_regions.size() > 1

                greader_list<reg_save* >::iterator fm = matched_regions.begin();
                greader_list<reg_save* >::iterator match_it = fm;
                ++match_it;

                std::copy((*raw_it)->exons->begin(), (*raw_it)->exons->end(), std::inserter((*fm)->exons, (*fm)->exons.begin()));
                for (; match_it != matched_regions.end(); ++match_it) {
                    std::copy((*match_it)->exons.begin(), (*match_it)->exons.end(), std::inserter((*fm)->exons, (*fm)->exons.begin()));

                    for (std::list<reg_save>::iterator reg_it = regions.begin(); reg_it != regions.end(); ++reg_it) {
                        if (&*reg_it == *match_it) {
                            regions.erase(reg_it);
                            break;
                        }
                    }
                } 
            }    
        }
    }
        
    // now create connected as separate regions!

    if (regions.size() < 2) {
        all_connected.push_back(*conn);
        return;
    }
        
    connected old_one = *conn; // all is lazy, no problem!
        
    for (std::list<reg_save>::iterator reg_it = regions.begin(); reg_it != regions.end(); ++reg_it) {

//            logger::Instance()->info("New Region.\n");

        all_connected.push_back(connected());
        connected* conn_it = &all_connected.back();
        conn_it->intel_count = old_one.intel_count;
        conn_it->avg_split = old_one.avg_split;

        conn_it->start = (*reg_it->exons.begin())->start;
        conn_it->end = (*reg_it->exons.rbegin())->end;
        conn_it->bam_count = old_one.bam_count;
        
        for (greader_refsorted_list<exon*>::iterator fexi = reg_it->exons.begin(); fexi != reg_it->exons.end(); ++fexi) {
            conn_it->fossil_exons->push_back(*fexi);
        }

        // second pass for matches to sort them in
        for(greader_refsorted_list<raw_atom* >::iterator raw_it = old_one.atoms->begin(); raw_it != old_one.atoms->end(); ++raw_it) { 

            if (!reg_it->exons.is_disjoint((*raw_it)->exons.ref())) {
                // this is a match!
                conn_it->atoms->insert(conn_it->atoms->end(), *raw_it);
            }
        }

        // clean out now missing pairs
        for(greader_refsorted_list<raw_atom* >::iterator a_it = conn_it->atoms->begin(); a_it != conn_it->atoms->end(); ++a_it) {
            for (paired_map<raw_atom*, rcount >::iterator p_it = (*a_it)->paired.begin(); p_it!= (*a_it)->paired.end(); ) {
                if (std::find(conn_it->atoms->begin(), conn_it->atoms->end(), p_it->first) == conn_it->atoms->end()) {
                    p_it = (*a_it)->paired.erase(p_it);
                } else {
                    ++p_it;
                }
            }
        }
    }    
}

void bam_reader::discard(const std::string &chrom_name) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->info("Discard finished " + chrom_name + ".\n");
    #endif
    
    #pragma omp critical(iterator_map_lock)
    {
    #pragma omp critical(chromosome_map_lock)
    {
        if (options::Instance()->is_stranded()) {
            chromosome_map_fwd.erase(chrom_name);
            chromosome_map_rev.erase(chrom_name);
            
            iterator_map.erase(chrom_name); 
        } else {
           chromosome_map_fwd.erase(chrom_name);
           iterator_map.erase(chrom_name);
        }
        
    }
    }
}

unsigned int bam_reader::get_num_connected(const std::string &chrom_name) {
    
    unsigned int total;
    #pragma omp critical(iterator_map_lock)
    {
        total = iterator_map[chrom_name].total();
    }
    return total;
}

bool bam_reader::populate_next_group(const std::string &chrom_name, greader_list<connected> &all_connected, exon_meta* meta, unsigned int &bam_count_absolute) {
    
    greader_list<connected>::iterator ob;
    bool exit = false;
    connection_iterator* it;
    #pragma omp critical(iterator_map_lock)
    {
        it = &iterator_map[chrom_name];
        if (!it->next(ob, meta->order_index)) {
            exit = true;
        }
    }
    if (exit) {
        #ifdef ALLOW_DEBUG
        logger::Instance()->info("Exit.\n");
        #endif
        return false;
    }
        
    // general statistics
    meta->chromosome = chrom_name;
    if (it->in_fwd) {
        meta->strand = "+";
    } else {
        meta->strand = "-";
    }
    meta->region_frags = it->fwd->frag_count; // same for both - and + strand in stranded!
    meta->region_reads = it->fwd->read_count; // same for both - and + strand in stranded!
  
    if (options::Instance()->is_stranded()) {
        meta->avrg_read_length = it->fwd->average_read_lenghts / 2 + it->bwd->average_read_lenghts / 2 ;
    } else {
        meta->avrg_read_length = it->fwd->average_read_lenghts / 2;
    }

    bam_count_absolute = it->fwd->bam_count;
    
    // split up components
    split_independent_component(&*ob, all_connected);
    
    logger::Instance()->info("Populate Group " + chrom_name + " " + std::to_string((*ob->fossil_exons.ref().begin())->start  ) + "-" + std::to_string((*ob->fossil_exons.ref().rbegin())->end) +".\n");

    return true;
}


void bam_reader::populate_next_single(const std::string &chrom_name, connected *ob, unsigned int &bam_count_absolute, pre_graph* raw, exon_meta* meta) {
    
    logger::Instance()->info("Populate Single " + chrom_name + " " + std::to_string((*ob->fossil_exons.ref().begin())->start  ) + "-" + std::to_string((*ob->fossil_exons.ref().rbegin())->end) +".\n");

    // set size, and add in mean
    meta->set_size(ob->fossil_exons.ref().size());
    raw->set_size(ob->fossil_exons.ref().size()); 
    
    unsigned int i=0;
    for (greader_list<exon* >::iterator e_it = ob->fossil_exons.ref().begin(); e_it != ob->fossil_exons.ref().end(); ++e_it,++i) {
        (*e_it)->id = i;
        meta->exons[i] = exon_meta::exon_meta_info();
        meta->exons[i].left = (*e_it)->start;
        meta->exons[i].right = (*e_it)->end;
        meta->exons[i].exon_length = meta->exons[i].right - meta->exons[i].left + 1;
        meta->exons[i].id = i; 
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Meta Exon"  + std::to_string( i ) + " " + std::to_string( meta->exons[i].left ) + " "  + std::to_string(meta->exons[i].right) +"\n");
        #endif
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Atoms " + std::to_string(ob->atoms.ref().size())+".\n");
    #endif
    
    raw->bam_count = ob->bam_count;
    raw->bam_count_total = bam_count_absolute;
    
    // now convert all atoms to pregraph types
    // atoms are sorted!
    unsigned int id = 0;
    for (greader_refsorted_list<raw_atom* > ::iterator atom = ob->atoms.ref().begin(); atom != ob->atoms.ref().end(); ++atom) {
        
        if ((*atom)->count == 0 && !(*atom)->reference_atom) {
          
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Raw Omitted " + std::to_string( (long) *atom) + " "  + (*atom)->to_string() +"\n");
            #endif

            continue;
        }
        
        (*atom)->id = id;
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Atom " + std::to_string( (long) *atom) + " "  + (*atom)->to_string() + " " + std::to_string((*atom)->id) +"\n");
        #endif
        
        raw->singled_bin_list.push_back(exon_group(i, (*atom)->exons->size()));
        exon_group* new_group = &raw->singled_bin_list.back();
        
        greader_refsorted_list<exon*>::iterator ae_it = (*atom)->exons.ref().begin();
        new_group->range_start = (*ae_it)->id;
        for (; ae_it != (*atom)->exons.ref().end(); ++ae_it) {
            new_group->set( (*ae_it)->id, true );
        }
        --ae_it;
        new_group->range_end = (*ae_it)->id;

        ++id;
    }
    
    std::set<unsigned int> guide_starts;
    std::set<unsigned int> guide_ends;
    for (greader_refsorted_list<raw_atom* > ::iterator atom = ob->atoms.ref().begin(); atom != ob->atoms.ref().end(); ++atom) {
           
       if ((*atom)->count == 0 && !(*atom)->reference_atom) {

            continue;
        }
       
       #ifdef ALLOW_DEBUG
       logger::Instance()->debug("Next Atom " + std::to_string((*atom)->id) + " " + (*atom)->to_string() + ".\n");
       #endif
       
       exon_group* lr = &raw->singled_bin_list[(*atom)->id];
       lr->read_count = (*atom)->count;
       lr->frag_count = lr->read_count - (*atom)->paired_count;
       lr->length_filterd = (*atom)->length_filtered;
       lr->drain_evidence = (*atom)->drain_evidence;
       lr->source_evidence = (*atom)->source_evidence;
       lr->reference_atom = (*atom)->reference_atom;
        
       if (lr->reference_atom) {
           if (lr->range_start != 0) {
               guide_starts.insert(lr->range_start);
           }
           if (lr->range_end != i-1) {
               guide_ends.insert(lr->range_end);
           }
           lr->reference_name = (*atom)->reference_name;
       }
       
       lr->total_lefts = (*atom)->total_lefts;
       lr->total_rights = (*atom)->total_rights;
       
       lr->lefts = (*atom)->lefts;
       lr->rights = (*atom)->rights;
              
//       if ((*a_it)->exons.ref().size() == 1) {
//           
//           std::map< rpos,rcount >::iterator li = lr->lefts->begin();
//           std::map< rpos,rcount >::iterator ri = lr->rights->begin();
//           std::vector< std::pair<rpos, rcount > > queue ;
//           
//           while (ri != lr->rights->end() ) {
//               if (li != lr->lefts->end() && li->first < = ri->first) {
//                   queue.push_back(*li);
//                   ++li;
//               } else {
//                   
//                    rcount capacity = ri->second;
//                    while (true) {
//                        rcount fc = std::min(queue.back().second, capacity);
//                        lr->base_count_start += (ri->first - li->first + 1) * fc;
//
//                        capacity -= fc;
//                        if (capacity == 0) {
//                            ++ri;
//                            break;
//                        }
//
//                        queue.back().second -= fc;
//                        if (queue.back().second == 0) {
//                            queue.pop_back();
//                        }
//                    }
//               }
//           }
//       } else {
//           for (std::map< rpos,rcount >::iterator li = lr->lefts->begin(); li != lr->lefts->end(); ++li) {
//               lr->base_count_start += (meta->exons[lr->range_start].right - li->first + 1) * li->second;
//           }
//           for (std::map< rpos,rcount >::iterator li = lr->rights->begin(); li != lr->rights->end(); ++li) {
//               lr->base_count_end += ( li->first - meta->exons[lr->range_end].left + 1) * li->second;
//           }           
//       }
       
       std::map< rpos,rcount >::iterator hsi = (*atom)->hole_starts->begin();
       std::map< rpos,rcount >::iterator hei = (*atom)->hole_ends->begin();
       
       unsigned int index = 0;
       for (greader_refsorted_list<exon*>::iterator ae_it = (*atom)->exons.ref().begin(); ae_it != (*atom)->exons.ref().end(); ++ae_it, ++index) {
           while (hsi != (*atom)->hole_starts->end() && hsi->first <= meta->exons[(*ae_it)->id].right && hsi->first >=  meta->exons[(*ae_it)->id].left) {
               lr->hole_starts[index].insert(*hsi);
               lr->hole_start_counts[index] += hsi->second;
               ++hsi;
           }
           while (hei != (*atom)->hole_ends->end() && hei->first <= meta->exons[(*ae_it)->id].right && hei->first >=  meta->exons[(*ae_it)->id].left) {
               lr->hole_ends[index].insert(*hei);
               lr->hole_end_counts[index] += hei->second;
               ++hei;
           } 
       }
       
       paired_map<raw_atom* , rcount>::iterator it = (*atom)->paired.begin();
       while(it != (*atom)->paired.end() && it->first->count == 0 && !it->first->reference_atom) {
            ++it;
        }
       if (it != (*atom)->paired.end() ) {
                      
           graph_list<paired_exon_group>::iterator start;
           
            exon_group* rr = &raw->singled_bin_list[it->first->id]; 

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("RR1 " + std::to_string( (long) it->first) + " " + std::to_string(it->first->id) + " " + it->first->to_string() + ".\n");
            #endif
            
            start = raw->paired_bin_list.insert(raw->paired_bin_list.end(), paired_exon_group(lr, rr, it->second));

            ++it;
            
            for (; it!= (*atom)->paired.end(); ++it) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("RR2 " + std::to_string( (long) it->first) + " " + std::to_string(it->first->id) + " " + it->first->to_string() + ".\n");
                #endif
                
                if (it->first->count > 0 || it->first->reference_atom) {
                    exon_group* rr = &raw->singled_bin_list[it->first->id];
                    start = raw->paired_bin_list.insert(start, paired_exon_group(lr, rr, it->second));
                }
            }
            
            std::sort(start, raw->paired_bin_list.end()); // sorting is done to achieve consistent output
       }
       
       raw->frag_count_region += lr->frag_count;
       raw->read_count_region += lr->read_count;
    }
    
    if (!guide_starts.empty() || !guide_ends.empty()) {
         for(graph_list<exon_group>::iterator sb_it = raw->singled_bin_list.begin(); sb_it != raw->singled_bin_list.end(); ++sb_it) {
             if (guide_starts.find(sb_it->range_start) != guide_starts.end()) {
                 sb_it->source_evidence = true;
             }
             if (guide_ends.find(sb_it->range_end) != guide_ends.end()) {
                 sb_it->drain_evidence = true;
             }
         }
    }
    
    #pragma omp critical(chromosome_map_lock)
    {
        chromosome *chrom_fwd, *chrom_rev;
        
        chrom_fwd = &chromosome_map_fwd[chrom_name];
        raw->frag_count_chrom = chrom_fwd->frag_count;
        raw->read_count_chrom = chrom_fwd->read_count; 
        raw->average_fragment_length = 2*chrom_fwd->average_read_lenghts + ob->avg_split;
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Average Fragsize " + std::to_string(chrom_fwd->average_read_lenghts) + " " + std::to_string(ob->avg_split) + "\n");
        #endif
        
        if (options::Instance()->is_stranded()) {
            chrom_rev = &chromosome_map_rev[chrom_name];
            raw->frag_count_chrom += chrom_rev->frag_count;
            raw->read_count_chrom += chrom_rev->read_count;
            
            rpos flen = 2*chrom_rev->average_read_lenghts + ob->avg_split;
            if (raw->average_fragment_length < flen) {
                raw->average_fragment_length = flen;
            }
        }        
    }

   #ifdef ALLOW_DEBUG
   logger::Instance()->debug("Read Counts " + std::to_string(raw->read_count_region) + " " + std::to_string(raw->read_count_chrom) +"\n");
   logger::Instance()->debug("Frag Counts " + std::to_string(raw->frag_count_region) + " " + std::to_string(raw->frag_count_chrom) +"\n");
   logger::Instance()->debug("Average Fragsize " + std::to_string(raw->average_fragment_length) + "\n");
   #endif
   
   // we can finally kill all contents of the connected
   ob->atoms.ref().clear();
   ob->fossil_exons.ref().clear();
   ob->reads.ref().clear();
    
}

unsigned long bam_reader::return_read_count(const std::string &file_name) {
    
    htsFile* file = hts_open(file_name.c_str(), "r");
    hts_idx_t* idx = sam_index_load(file, file_name.c_str());
    bam_hdr_t *hdr = sam_hdr_read(file);
    
    uint64_t mapped_total = 0;
    
    uint64_t mapped;
    uint64_t unmapped;
        
    for (unsigned int i=0;  i < hdr->n_targets; i++)
    {
        hts_idx_get_stat(idx, i, &mapped, &unmapped);
        mapped_total += mapped;
    }

    hts_close(file);
    return mapped_total;
}

void bam_reader::return_chromosome_names(const std::string &file_name, greader_name_set<std::string> &return_list) {
    
    // open reader
    htsFile* file = hts_open(file_name.c_str(), "r");
    
    // open header
    bam_hdr_t *hdr = sam_hdr_read(file);
    

    for(int i = 0; i < hdr->n_targets; i++) {
        return_list.insert(hdr->target_name[i]);
    }
    
    hts_close(file);
}


void bam_reader::read_chromosome(std::string file_name, std::string chrom_name, const unsigned int input_prefix, bool last_file) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Read chromosome " + file_name + " - "+ chrom_name+" is last: " + std::to_string(last_file) +".\n");
    #endif
    
    std::string id_prefix = std::to_string(input_prefix);
    
    // get or create the chromosome by name
    // we only use fwd if unstranded
    chromosome *chrom_fwd, *chrom_rev; 
    
    #pragma omp critical(chromosome_map_lock)
    {
        chrom_fwd = &chromosome_map_fwd[chrom_name];
        if (options::Instance()->is_stranded()) {
            chrom_rev = &chromosome_map_rev[chrom_name];
        }        
    }
    
    ++chrom_fwd->bam_count;
    if (options::Instance()->is_stranded()) {
         ++chrom_rev->bam_count;
    }
    
    // open reader
    htsFile* file = hts_open(file_name.c_str(), "r");
    // we need an index, because otherwise we need to read the whole file multiple times
    // we can use the same name, it will add the bai automatically
    hts_idx_t* idx = sam_index_load(file, file_name.c_str());
    
    // we also need to the header for indexes
    bam_hdr_t *hdr = sam_hdr_read(file);
    hts_itr_t *iter = NULL;
    
    // iterator and object to read stuff into
    iter = sam_itr_querys(idx, hdr, chrom_name.c_str()); 
    bam1_t *read = bam_init1();
    
    // we need to move through the list of already there exons
    r_border_set<rpos>::iterator ex_start_fwd_it = chrom_fwd->fixed_exon_starts.ref().begin();
    r_border_set<rpos>::iterator ex_end_fwd_it = chrom_fwd->fixed_exon_ends.ref().begin();
    r_border_set<rpos>::iterator ex_start_rev_it, ex_end_rev_it;
    if (options::Instance()->is_stranded()) {
        ex_start_rev_it = chrom_rev->fixed_exon_starts.ref().begin();
        ex_end_rev_it = chrom_rev->fixed_exon_ends.ref().begin();
    }    
    
    // strands can overlap without
    rpos left_border_fwd = 0;
    rpos right_border_fwd = 0;
    rpos left_border_rev = 0;
    rpos right_border_rev = 0;
    
    
    // main loop of this, we always read ALL reads, but compacting is called multiple times
    while ( sam_itr_next(file, iter, read) >= 0) {
	
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Next Sam Iter " + std::string(bam_get_qname(read)) + "\n");
        #endif
        
        if (read->core.flag & BAM_FSECONDARY || read->core.flag & BAM_FUNMAP || read->core.qual < 1) { // this is a secondary alignment
            // only use primary alignments
            continue;
        }
        
        uint32_t* cigar = bam_get_cigar(read);
        bool long_intron = false;
        bool too_small = false;
        for(int i = 0; i < read->core.n_cigar; i++) {
            const int op = bam_cigar_op(cigar[i]);
            const int ol = bam_cigar_oplen(cigar[i]);
            // do we have an N?
            if (op == BAM_CREF_SKIP && ol > options::Instance()->get_maximal_intron_size()) {
                logger::Instance()->info("Skip Read With Overly Long Intron " + std::string(bam_get_qname(read)) + "\n");
                long_intron = true;
                break;
            }
            if (op == BAM_CMATCH && i != 0 && i != read->core.n_cigar - 1 && ol < 3) {
                too_small = true;
                break;
            }
        }
        if(long_intron || too_small) {
            continue;
        }
        
        // process stranded or unstranded
         if (options::Instance()->is_stranded()) {
             // we need to get strand information, use XS tag from mappers
              
            char xs = '.'; 
            uint8_t* ptr = bam_aux_get(read, "XS");
            if (ptr)
            {
                    char src_strand_char = bam_aux2A(ptr);
                    if (src_strand_char == '-') {
                        xs = '-';
                    }
                    else if (src_strand_char == '+') {
                        xs = '+';
                    }         
            } 
            
            if (options::Instance()->get_strand_type() == options::unknown && xs == '.') {
                 // for unknown we test the XS tag
                    logger::Instance()->warning("No strand information on read ");
                    logger::Instance()->warning(bam_get_qname(read));
                    logger::Instance()->warning(". Add XS tag, specify library or turn of stranded option. Read was Skipped. \n");
                    continue;

            } else { 
                char strand = '.'; 
                uint32_t sam_flag = read->core.flag;
                                
                bool antisense_aln = sam_flag & BAM_FREVERSE; //16
                //                1                             64                           1
                if (((sam_flag & BAM_FPAIRED) && (sam_flag & BAM_FREAD1)) || !(sam_flag & BAM_FPAIRED)) // first-in-pair or single-end
                {
                        switch(options::Instance()->get_strand_type())
                        {
                            case options::FF:
                            case options::FR:
                                    (antisense_aln) ? strand = '-' : strand = '+';
                                    break;
                            case options::RF:
                            case options::RR:
                                    (antisense_aln) ? strand = '+' :  strand = '-';
                                    break;
                        }
                }
                else if ((sam_flag & BAM_FPAIRED) && (sam_flag & BAM_FREAD2)) // second-in-pair read
                {
                        switch(options::Instance()->get_strand_type())
                        {
                            case options::FF:
                            case options::RF:
                                    (antisense_aln) ?  strand = '-' : strand = '+';
                                    break;
                            case options::FR:
                            case options::RR:
                                    (antisense_aln) ? strand = '+' : strand = '-';
                                    break;
                        }
                }
                if (strand == '.') {
                    strand = xs;
                }
                
                if (strand == '+') {
                    if (xs != '-')
                        process_read(read, left_border_fwd, right_border_fwd, chrom_fwd, id_prefix, ex_start_fwd_it, ex_end_fwd_it, last_file);
                } else if (strand == '-'){
                    if (xs != '+')
                        process_read(read, left_border_rev, right_border_rev, chrom_rev, id_prefix, ex_start_rev_it, ex_end_rev_it, last_file);
                }
            }
            
         } else {
                process_read(read, left_border_fwd, right_border_fwd, chrom_fwd, id_prefix, ex_start_fwd_it, ex_end_fwd_it, last_file);
         }
    }
    
    if (options::Instance()->is_stranded()) {
        finish_block(chrom_fwd, left_border_fwd, right_border_fwd, ex_start_fwd_it, ex_end_fwd_it, last_file);
        finish_block(chrom_rev, left_border_rev, right_border_rev, ex_start_rev_it, ex_end_rev_it, last_file);
    } else {
        finish_block(chrom_fwd, left_border_fwd, right_border_fwd, ex_start_fwd_it, ex_end_fwd_it, last_file);
    }
    
    if (last_file || options::Instance()->is_parse_heuristic()) {
         if (options::Instance()->is_stranded()) {
            reset_reads(chrom_fwd);
            reset_reads(chrom_rev);
         } else {
            reset_reads(chrom_fwd); 
         }
    }
    
    // destroy everything again
    hts_idx_destroy(idx);
    hts_itr_destroy(iter);
    bam_destroy1(read);
    bam_hdr_destroy(hdr);
    hts_close(file);
    
}




//
// while ( sam_itr_next(file, iter, read) >= 0) {
//	
////        if ( read->core.flag & BAM_FSECONDARY ) { // this is a secondary alignment
////            // only use primary alignments
////            continue;
////        }
//        
//        // process stranded or unstranded
//         if (options::Instance()->is_stranded()) {
//             // we need to get strand information, use XS tag from mappers
//             
//            uint8_t* ptr = bam_aux_get(read, "XS");
//            if (ptr)
//            {
//                    char src_strand_char = bam_aux2A(ptr);
//                    if (src_strand_char == '-') {
//                        process_read(read, left_border_rev, right_border_rev, chrom_rev, id_prefix, ex_start_rev_it, ex_end_rev_it, last_file); 
//                    }
//                    else if (src_strand_char == '+') {
//                        process_read(read, left_border_fwd, right_border_fwd, chrom_fwd, id_prefix, ex_start_fwd_it, ex_end_fwd_it, last_file);
//                    }         
//            } else {
//                if (options::Instance()->get_strand_type() == options::unknown) {
//                 // for unknown we test the XS tag
//                
//                    logger::Instance()->warning("No strand information on read ");
//                    logger::Instance()->warning(bam_get_qname(read));
//                    logger::Instance()->warning(". Add XS tag, specify library or turn of stranded option. Read was Skipped. \n");
//                    continue;
//                } else {    
//                    bool fwd;
//                    uint32_t sam_flag = read->core.flag;
//                    bool antisense_aln = sam_flag & BAM_FREVERSE;
//                    if (((sam_flag & BAM_FPAIRED) && (sam_flag & BAM_FREAD1)) || !(sam_flag & BAM_FPAIRED)) // first-in-pair or single-end
//                    {
//                            switch(options::Instance()->get_strand_type())
//                            {
//                                case options::FF:
//                                case options::FR:
//                                        (antisense_aln) ? fwd = false : fwd = true;
//                                        break;
//                                case options::RF:
//                                case options::RR:
//                                        (antisense_aln) ? fwd = true :  fwd = false;
//                                        break;
//                            }
//                    }
//                    else // second-in-pair read
//                    {
//                            switch(options::Instance()->get_strand_type())
//                            {
//                                case options::FF:
//                                case options::RF:
//                                        (antisense_aln) ?  fwd = false : fwd = true;
//                                        break;
//                                case options::FR:
//                                case options::RR:
//                                        (antisense_aln) ? fwd = true :  fwd = false;
//                                        break;
//                            }
//                    }
//                
//                    if (fwd) {
//                        logger::Instance()->debug("FWD  \n");
//                        process_read(read, left_border_fwd, right_border_fwd, chrom_fwd, id_prefix, ex_start_fwd_it, ex_end_fwd_it, last_file);
//                    } else {
//                        logger::Instance()->debug("REV \n");
//                        process_read(read, left_border_rev, right_border_rev, chrom_rev, id_prefix, ex_start_rev_it, ex_end_rev_it, last_file);
//                    }
//                }
//            }
//        } else {
//                process_read(read, left_border_fwd, right_border_fwd, chrom_fwd, id_prefix, ex_start_fwd_it, ex_end_fwd_it, last_file);
//        }
//    }

void bam_reader::process_read( bam1_t *bread, rpos &left_border, rpos &right_border, chromosome* chrom, const std::string id_prefix,
        r_border_set<rpos>::iterator &ex_start_it, r_border_set<rpos>::iterator &ex_end_it, bool last_file) {
    
    rpos left, right;
    greader_list<interval> junctions;
    greader_list<std::pair<rpos, rpos> > splices;
    rread* prev = NULL;
    if (!chrom->read_queue.empty()) {
        prev = &*chrom->read_queue.rbegin();
    }
    rread* new_read = parse_read(bread, chrom, junctions, splices, left, right, id_prefix, ex_start_it, ex_end_it);

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add read ");
    if (new_read->id_set && ! bread->core.flag & BAM_FSECONDARY) logger::Instance()->debug(new_read->ids.ref().front());
    logger::Instance()->debug( " " + std::to_string(left) + " " + std::to_string(right) + " BORDER " + std::to_string(right_border) + ".\n");
    #endif

    // heuristic early merge down of identical reads!
    if (prev != NULL && new_read->left_limit == prev->left_limit && new_read->right_limit == prev->right_limit) {
        // boundaries are the same, test intervals for compacting!
        bool merge = true;
        greader_list<interval>::reverse_iterator cur_it = junctions.rbegin();
        greader_list<interval>::reverse_iterator prev_it = chrom->interval_queue.rbegin();
        for (; cur_it != junctions.rend(); ++cur_it, ++prev_it) {
            if ( cur_it->left != prev_it->left || cur_it->right != prev_it->right ) {
                merge = false;
                break;
            }
        }
        
        if( merge ) {
            
            for (greader_list<std::pair<rpos, rpos> >::iterator is = splices.begin(); is != splices.end(); ++is) {
                chrom->splice_queue[*is] += 1;
            }
            
            ++prev->count;
            if (new_read->id_set) {
                prev->add_id(new_read->ids.ref().front());
            }
            chrom->read_queue.pop_back();
            return;
        }
        
    }

    rpos ejd = options::Instance()->get_exon_join_distance();    

    // we have a new region, so start new with net read
    if (left_border == right_border && right_border == 0) {
        left_border = left;
        right_border = right;
    }

    if (left <= right_border + 1 + ejd && right >= right_border ) {
        right_border =  right;
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Added read extends border.\n");
        #endif
    } else if (left > right_border + 1 + ejd) {
        
        // not overlapping on chromosome, so we need to finish up so far collected data
        rread last = *new_read;
        chrom->read_queue.pop_back();
        finish_block(chrom, left_border, right_border, ex_start_it, ex_end_it, last_file);
        
        rread* re_add = chrom->addQueuedRead(last);
        for (greader_list<interval>::iterator it = junctions.begin(); it != junctions.end(); ++it) {
            it->parent = re_add;
        }
        
        left_border = left;
        right_border = right;
    } else {
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("No border manipulation.\n");
        #endif
    }
     
    for (greader_list<std::pair<rpos, rpos> >::iterator is = splices.begin(); is != splices.end(); ++is) {
        chrom->splice_queue[*is] += 1;
    }
    std::copy(junctions.begin(), junctions.end(), std::back_inserter(chrom->interval_queue));
}

rread* bam_reader::parse_read( bam1_t *bread,  chromosome* chrom, greader_list<interval> &junctions,
        greader_list<std::pair<rpos, rpos> > &splices,
        rpos &left, rpos &right, const std::string id_prefix,
        r_border_set<rpos>::iterator &ex_start_it, r_border_set<rpos>::iterator &ex_end_it) {
    
    rread* new_read;
    
    ++chrom->read_count;
    // && bread->core.flag & BAM_FPROPER_PAIR Why only proper pair??
    if ( bread->core.flag & BAM_FPAIRED &&  bread->core.tid == bread->core.mtid) { // this is a proper paired read on same ref
        // add this read to the chromosome
        std::string id(bam_get_qname(bread));
        id.append("_").append(id_prefix);

        // new read object and add it to chromosome (we cannot have any name twice)
        new_read = chrom->addQueuedRead(rread(id));
        
        if (bread->core.mpos > bread->core.pos) { // count only left pair!
            ++chrom->frag_count;
        }
    } else {
        new_read = chrom->addQueuedRead(rread());
        ++chrom->frag_count;
    }
    
    bam1_core_t *rcore = &bread->core;
    left = rcore->pos + 1; // original is 0 based, we make this 1 based, much nicer
    rpos start = left;
    rpos offset = 0;
    
    // stepwise averaging
    chrom->average_read_lenghts +=  (rcore->l_qseq - chrom->average_read_lenghts) / chrom->read_count;
    
    // advance the known exon iterators
    // we use those to keep track of new split evidence to avoid double definitions
    ex_start_it = chrom->fixed_exon_starts.ref().lower_bound(ex_start_it, chrom->fixed_exon_starts.ref().end(), left);
    ex_end_it = chrom->fixed_exon_ends.ref().lower_bound(ex_end_it, chrom->fixed_exon_ends.ref().end(), left);
    
    std::deque<rpos> new_ends;
    std::deque<rpos> new_starts;
    greader_list<interval> new_junctions;
    
    // go over the read and it's split info
    // M = 0 I = 1 D = 2 N = 3 S = 4 H = 5 P = 6
    uint32_t* cigar = bam_get_cigar(bread);
    for(int i = 0; i < rcore->n_cigar; i++) {

        const int op = bam_cigar_op(cigar[i]);
        const int ol = bam_cigar_oplen(cigar[i]);

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Operator " + std::to_string(op) + " " + std::to_string(ol) + ".\n");
        #endif
        
        // do we have an N?
        if (op == BAM_CREF_SKIP) {
               
            rpos end = start + offset;
            new_junctions.push_back(interval(new_read) );
            interval* interv = &new_junctions.back();
            interv->left = start;
            interv->right = end-1;
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Add Interval1 " + std::to_string(start) + " "  + std::to_string(end-1)+".\n");
            #endif
                         
            start = end + ol;
            offset = 0;
        } else if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF) {
            offset = offset + ol;
        } 
    }
    
    // last interval
    rpos end = start + offset;
    new_junctions.push_back(interval(new_read) );
    interval* interv = &new_junctions.back();
    interv->left = start;
    interv->right = end-1; // THIS is including obviously
    
    // filter out splits sites that are too small to count
    
    greader_list<interval>::iterator nj = new_junctions.begin();
    greader_list<interval>::iterator njn = nj;
    ++njn;
    while (njn != new_junctions.end()) {
        if (njn->left - nj->right < options::Instance()->get_min_intron_length()) {
            njn->left = nj->left;
        } else {
            junctions.push_back(*nj);
            new_read->add_length(&junctions.back());
            splices.push_back( std::make_pair(nj->right+1, njn->left-1));
            
            new_ends.push_back(nj->right);
            new_starts.push_back(njn->left);
        }
        nj = njn;
        ++njn;
    }
    junctions.push_back(*nj);
    new_read->add_length(&junctions.back());
    
    // this is over filtering
    bool no_start = false;
    if (junctions.front().right - junctions.front().left + 1 < options::Instance()->get_min_junction_anchor() && junctions.size() > 1) {
        no_start = true;
        splices.pop_front();
    }
    bool no_end = false;
    if (junctions.back().right - junctions.back().left + 1 < options::Instance()->get_min_junction_anchor() && junctions.size() > 1) {
        no_end = true;
        if (!splices.empty()) splices.pop_back();
    }
    
    for (std::deque<rpos>::iterator it = new_ends.begin(); it != new_ends.end(); ++it) {
        bool drop = ( it == new_ends.begin() && no_start ) || ( it == new_ends.end()-1 && no_end );
        add_known_end(chrom, *it, ex_end_it, !drop);
    }
    
    for (std::deque<rpos>::iterator it = new_starts.begin(); it != new_starts.end(); ++it) {
        bool drop = ( it == new_starts.begin() && no_start ) || ( it == new_starts.end()-1 && no_end );
        add_known_start(chrom, *it, ex_start_it, !drop);
    }
    
    new_read->set_left_limit(junctions.front().left);
    new_read->set_right_limit(junctions.back().right);
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Interval2 " + std::to_string(start) + " "  + std::to_string(end-1)+".\n");
    #endif
    
    right = end-1;
    
    return new_read;
}


void bam_reader::add_known_start( chromosome* chrom, const rpos pos,
         r_border_set<rpos>::iterator &ex_start_it, bool evidence) {
     
    r_border_set<rpos>::iterator lower = chrom->fixed_exon_starts.ref().lower_bound(ex_start_it, chrom->fixed_exon_starts.ref().end(), pos);
    
    unsigned int max_extend = options::Instance()->get_max_pos_extend();
    if (lower != chrom->fixed_exon_starts.ref().end()) {
        if (pos + max_extend >= *lower && pos <= *lower + max_extend) {
            return; // this already is a known site!
        }
    }
    
    if (lower != chrom->fixed_exon_starts.ref().begin()) { // only test previous element if it exists!
        --lower;
        if (pos + max_extend >= *lower  && pos <= *lower + max_extend) {
            return; // this already is a known site!
        }
    }
    
    // this were all possible hit
    // we found something entirely new!
    chrom->known_starts.push_back(std::make_pair(pos, evidence));
    
}


void bam_reader::add_known_end( chromosome* chrom, const rpos pos,
       r_border_set<rpos>::iterator &ex_end_it, bool evidence) {
    
    r_border_set<rpos>::iterator lower = chrom->fixed_exon_ends.ref().lower_bound(ex_end_it, chrom->fixed_exon_ends.ref().end(), pos);
    
    unsigned int max_extend = options::Instance()->get_max_pos_extend();
    if (lower != chrom->fixed_exon_ends.ref().end()) { 
        if (pos + max_extend >= *lower  && pos <= *lower + max_extend) {
            return; // this already is a known site!
        }
    }
    
    if (lower != chrom->fixed_exon_ends.ref().begin()) { // only test previous element if it exists!
        --lower;
        if (pos + max_extend >= *lower && pos <= *lower + max_extend) {
            return; // this already is a known site!
        }
    }
    
    // this were all possible hit
    // we found something entirely new!
    chrom->known_ends.push_back(std::make_pair(pos, evidence));
    
}


void bam_reader::finish_block(chromosome* chrom,  rpos &left,  rpos &right, r_border_set<rpos>::iterator &ex_start_it, r_border_set<rpos>::iterator &ex_end_it, bool last_file) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Finish block " + std::to_string(left) + " - " + std::to_string(right) + ".\n");
    #endif
     
    // create overlapping regions from new intervals
    greader_list<std::pair<rpos, rpos> > raw;
    create_raw_exons(chrom, raw, right);
    
    // obscure safe condition for very low read count input
    if (raw.empty()) {
        return;
    }
    
//    for(greader_list<std::pair<rpos, rpos> >::iterator ri = raw.begin(); ri != raw.end(); ++ri) {
//        logger::Instance()->debug("RAW EXON "+ std::to_string(ri->first) + " " + std::to_string(ri->second) +"\n");
//    }
//    
//    
//    for (greader_list<connected>::iterator it = chrom->chrom_fragments.begin() ; it != chrom->chrom_fragments.end(); ++it) {
//        logger::Instance()->debug("Pre Connected "+ std::to_string(it->start) + " " + std::to_string(it->end) +"\n");
//    }
    
    //  we update the know connected areas
    connected* conn = &*insert_fragment(chrom, left, right);
    ++conn->bam_count;
    
//    for (greader_list<connected>::iterator it = chrom->chrom_fragments.begin() ; it != chrom->chrom_fragments.end(); ++it) {
//        logger::Instance()->debug("After Connected "+ std::to_string(it->start) + " " + std::to_string(it->end) +"\n");
//        
//    }
//    logger::Instance()->debug("Inside Connected "+ std::to_string(conn->start) + " " + std::to_string(conn->end) +"\n");
//     
//    for(greader_refsorted_list<raw_atom*>::iterator atom_it = conn->atoms.ref().begin(); atom_it != conn->atoms.ref().end(); ++atom_it) {
//        logger::Instance()->debug("Preexisting Atom A " + std::to_string((long) *atom_it) + " " + (*atom_it)->to_string() + ".\n");
//    }
//    for(greader_list<raw_atom>::iterator atom_it = chrom->atoms.begin(); atom_it != chrom->atoms.end(); ++atom_it) {
//        logger::Instance()->debug("Base Atom A " + std::to_string((long) &*atom_it) + " " + atom_it->to_string() + ".\n");
//    }
//    #ifdef ALLOW_DEBUG
//    for (greader_list<exon* >::iterator e_it = conn->fossil_exons.ref().begin(); e_it != conn->fossil_exons.ref().end(); ++e_it) {
//         logger::Instance()->debug("FINAL EXON A: " +  std::to_string((*e_it)->start) + "-" + std::to_string((*e_it)->end)+"\n");
//    }
//    #endif
    
    // if we have new data, cluster for best support!
    // this filters out unsupported clusters by low coverage
    greader_list<rpos> clustered_starts;
    greader_list<std::pair<rpos, bool> >::iterator end_starts = cluster(chrom->known_starts, clustered_starts, right);
    
    
    greader_list<rpos> clustered_ends;
    greader_list<std::pair<rpos, bool> >::iterator end_ends = cluster(chrom->known_ends, clustered_ends, right);
    
    // reset known lists for next round
    chrom->known_starts.erase(chrom->known_starts.begin(), end_starts );
    chrom->known_ends.erase(chrom->known_ends.begin(), end_ends);
     
    // the starts and ends should reflect the the clustered starts and ends!
    solidify_raw_exons_ends(chrom, raw, clustered_starts, clustered_ends);
//    if (!conn->guided && options::Instance()->is_trimming()) trim_exons_1(chrom, raw, clustered_starts, clustered_ends);
    if (!conn->guided && options::Instance()->is_trimming()) trim_exons_2(chrom, raw, clustered_starts, clustered_ends);
    
    // combine them with others
    update_existing_exons(conn, chrom, raw, left, right);
     
//    #ifdef ALLOW_DEBUG
//    for (greader_list<exon* >::iterator e_it = conn->fossil_exons.ref().begin(); e_it != conn->fossil_exons.ref().end(); ++e_it) {
//         logger::Instance()->debug("FINAL EXON B: " +  std::to_string((*e_it)->start) + "-" + std::to_string((*e_it)->end)+"\n");
//    }
//    #endif
//    
//    for(greader_refsorted_list<raw_atom*>::iterator atom_it = conn->atoms.ref().begin(); atom_it != conn->atoms.ref().end(); ++atom_it) {
//        logger::Instance()->debug("Preexisting Atom B " + std::to_string((long) *atom_it) + " " + (*atom_it)->to_string() + ".\n");
//    }
//    for(greader_list<raw_atom>::iterator atom_it = chrom->atoms.begin(); atom_it != chrom->atoms.end(); ++atom_it) {
//        logger::Instance()->debug("Base Atom B " + std::to_string((long) &*atom_it) + " " + atom_it->to_string() + ".\n");
//    }
//    
    // split exons on clustered
    split_exons(conn, chrom, clustered_starts, left, right, 1);
    split_exons(conn, chrom, clustered_ends, left, right, 0);
    
//    #ifdef ALLOW_DEBUG
//    for (greader_list<exon* >::iterator e_it = conn->fossil_exons.ref().begin(); e_it != conn->fossil_exons.ref().end(); ++e_it) {
//         logger::Instance()->debug("FINAL EXON C: " +  std::to_string((*e_it)->start) + "-" + std::to_string((*e_it)->end)+"\n");
//    }
//    #endif
//
//    for(greader_refsorted_list<raw_atom*>::iterator atom_it = conn->atoms.ref().begin(); atom_it != conn->atoms.ref().end(); ++atom_it) {
//        logger::Instance()->debug("Preexisting Atom C " + std::to_string((long) *atom_it) + " " + (*atom_it)->to_string() + ".\n");
//    }
//    for(greader_list<raw_atom>::iterator atom_it = chrom->atoms.begin(); atom_it != chrom->atoms.end(); ++atom_it) {
//        logger::Instance()->debug("Base Atom C " + std::to_string((long) &*atom_it) + " " + atom_it->to_string() + ".\n");
//    }
//    
//    
    // now look at new separators
    add_to_fixed_clusters(chrom->fixed_exon_starts, clustered_starts, ex_start_it);
    add_to_fixed_clusters(chrom->fixed_exon_ends, clustered_ends, ex_end_it);
    
    
//    logger::Instance()->debug("Clustered Starts: ");
//    for (greader_list<rpos>::iterator f = clustered_starts.begin(); f!= clustered_starts.end(); ++f) {
//        logger::Instance()->debug(std::to_string(*f) + ",");
//    }
//    logger::Instance()->debug("\n");
//    logger::Instance()->debug("Clustered Ends: ");
//    for (greader_list<rpos>::iterator f = clustered_ends.begin(); f!= clustered_ends.end(); ++f) {
//        logger::Instance()->debug(std::to_string(*f) + ",");
//    }
//    logger::Instance()->debug("\n");
//    
//    logger::Instance()->debug("Fixed Starts: ");
//    for (r_border_set<rpos>::iterator f = chrom->fixed_exon_starts->begin(); f!= chrom->fixed_exon_starts->end(); ++f) {
//        logger::Instance()->debug(std::to_string(*f) + ",");
//    }
//    logger::Instance()->debug("\n");
//    logger::Instance()->debug("Fixed Ends: ");
//    for (r_border_set<rpos>::iterator f = chrom->fixed_exon_ends->begin(); f!= chrom->fixed_exon_ends->end(); ++f) {
//        logger::Instance()->debug(std::to_string(*f) + ",");
//    }
//    logger::Instance()->debug("\n");
    
    // now assign new reads
    assign_reads(conn, chrom);
    
//    for ( greader_list<rread>::iterator r_it = chrom->read_queue.begin(); r_it != chrom->read_queue.end(); ++r_it) {
//        logger::Instance()->debug("FREAD1: ");
//        for (greader_list<std::string>::iterator id = r_it->ids->begin(); id !=  r_it->ids->end(); ++id) {
//            logger::Instance()->debug(*id + ", ");
//        }
//        logger::Instance()->debug("\nExons1: ");
//        for (greader_refsorted_list<exon*>::iterator e = r_it->atom->exons->begin(); e != r_it->atom->exons->end(); ++e) {
//             logger::Instance()->debug( std::to_string((*e)->start) + "-" + std::to_string((*e)->end) + ", ");
//        }
//        logger::Instance()->debug("\n");
//    }
    
    
    filter_outer_read_junctions(chrom);
    
//    for ( greader_list<rread>::iterator r_it = chrom->read_queue.begin(); r_it != chrom->read_queue.end(); ++r_it) {
//        logger::Instance()->debug("FREAD2: ");
//        for (greader_list<std::string>::iterator id = r_it->ids->begin(); id !=  r_it->ids->end(); ++id) {
//            logger::Instance()->debug(*id + ", ");
//        }
//        logger::Instance()->debug("\nExons2: ");
//        for (greader_refsorted_list<exon*>::iterator e = r_it->atom->exons->begin(); e != r_it->atom->exons->end(); ++e) {
//             logger::Instance()->debug( std::to_string((*e)->start) + "-" + std::to_string((*e)->end) + ", ");
//        }
//        logger::Instance()->debug("\n");
//    }
// 
    // reduce down to minimal amount of atoms
    reduce_atoms(conn, chrom);
//    
//    for(greader_refsorted_list<raw_atom*>::iterator atom_it = conn->atoms.ref().begin(); atom_it != conn->atoms.ref().end(); ++atom_it) {
//        logger::Instance()->debug("Preexisting Atom D " + std::to_string((long) *atom_it) + " " + (*atom_it)->to_string() + ".\n");
//    }
//    for(greader_list<raw_atom>::iterator atom_it = chrom->atoms.begin(); atom_it != chrom->atoms.end(); ++atom_it) {
//        logger::Instance()->debug("Base Atom D " + std::to_string((long) &*atom_it) + " " + atom_it->to_string() + ".\n");
//    }
//    
    chrom->read_queue.clear();
    chrom->interval_queue.clear();
    chrom->splice_queue.clear();
       
    mark_or_reduce_paired_atoms(conn, chrom , conn->atoms.ref().begin(), conn->atoms.ref().end());
        
    if (options::Instance()->is_parse_heuristic() || last_file) {
        // remove all individual reads to avoid memory build up
        if (!conn->guided) filter_bins(conn, chrom);
        reduce_reads(conn);
    }
    
}

// ######### CLUSTERING #########

greader_list<std::pair<rpos, bool> >::iterator bam_reader::cluster( greader_list<std::pair<rpos, bool> > &in,  greader_list<rpos> &out, rpos &right) {
    
    
     #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Cluster Group.\n");
     #endif
            
    // nothing to do here
    if (in.empty()) {
        return in.end();
    }
    
    // we need a sorted one to get found
    in.sort();
    
    // get option
    const unsigned int extend = options::Instance()->get_max_pos_extend();
    
    // init basic values
    unsigned int count = 1;
    greader_list<std::pair<rpos, bool> >::iterator it = in.begin();
    rpos current = it->first;
    bool evidence = it->second;
    ++it;
    
    std::vector<std::pair<rpos,unsigned int> > basic_grouping;
     
    // now search for connected blocks
    for (; it != in.end() && it->first <= right; ++it) {
        if (it->first == current) { // repeated element
            ++count;
            evidence = evidence || it->second;
        } else {   // different element
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Add cluster candidate group: "+ std::to_string(current) + ", " + std::to_string(count)+".\n");
            #endif
            
                basic_grouping.push_back(std::make_pair(current, count));
            
            // do we have combine to real cluster?
            if (current + extend < it->first ) { // *it > current by sorting
                // end of cluster found, process values in temp
                
                if (evidence) {
                    DKMeans(basic_grouping, out, extend);  
                }
                basic_grouping.clear();
                evidence = false;
            }
            
            current = it->first;
            count = 1;
            evidence = evidence || it->second;
            
        }
    }
    
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Add cluster candidate group: "+ std::to_string(current) + ", " + std::to_string(count)+".\n");
    #endif
    
    basic_grouping.push_back(std::make_pair(current, count));
    if (evidence) {
        DKMeans(basic_grouping, out, extend);
    }
  
    return it;
    
}


void bam_reader::DKMeans( std::vector<std::pair<rpos,unsigned int> > &in,  greader_list<rpos> &out, const unsigned int &extend) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Inner Custer Group.\n");
    #endif
    
    // filter low coverage!
    unsigned int count = 0;
    for (std::vector<std::pair<rpos,unsigned int> >::iterator in_it = in.begin(); in_it != in.end(); ++in_it) {
            count += in_it->second;
    }
    if (count < options::Instance()->get_min_junction_coverage()) {
        #ifdef ALLOW_DEBUG
         logger::Instance()->debug("Cluster filtered " + std::to_string(in[0].first)+".\n");
         #endif
        return;
    }
    
    greader_list<rpos>::iterator end_it = out.end();
    
    const unsigned int size = in.size();
    
    if (in.size() == 1) {
         out.push_back(in[0].first);
         #ifdef ALLOW_DEBUG
         logger::Instance()->debug("Clustered Position added " + std::to_string(in[0].first)+".\n");
         #endif
         return;
    }
    
    //################ FORWARD ################
    
    
    // 2D dynamic lookup table arrays
    double *costs = new double[size*size];
    unsigned int *back = new unsigned int[size*size];
    
    // we use unnormalized weighted least square
    // d = sum w_i (x_i - mu)^2
    
    // temporary mean values
    double mean_x1;
    unsigned int n_x1;
    
    costs[0] = 0.0;
    back[0] = 1;
    
    
    // init first row without backvalue k=0
    mean_x1 = in[0].first;
    n_x1 = in[0].second;
    for(unsigned int i = 1; i < size; ++i) {
            back[i] = 0;
               
            if (in[i].first > in[0].first + 2 * extend ) {
                // cannot make this extend, add special cancel value
                costs[i] = -1.0; // all negative will be treated as invalid

            } else {
                double new_mean_x1 = mean_x1 + in[i].second / (double) (n_x1 + in[i].second)
                        * (in[i].first - mean_x1);
                costs[i ] = costs[(i-1)]  + in[i].second * (in[i].first - mean_x1)* (in[i].first - new_mean_x1);
                        
                mean_x1 = new_mean_x1;
                n_x1 += in[i].second;
            }
    }
    
    // now update first diagonal as always 0
//    for(unsigned int k = 1; k < size; ++k) {
//        costs[k*size+k] = 0;
//    }
    

    for(unsigned int k = 1; k < size; ++k) {    
        for(unsigned int i = k+1; i < size; ++i) {
            
            double d = 0.0;
            double mean_xj = 0.0;
            unsigned int n_xj = 0;
            
            int min = -1;
            double min_cost = -1; 
             
            for(unsigned int j = i; j > k; --j) {
            
                if (in[i].first > in[j].first + 2 * extend ) {
                    break;
                }
                
                double new_mean_xj = mean_xj + in[j].second / (double) (n_xj + in[j].second)
                        * (in[j].first - mean_xj);
                d = d + in[j].second * (in[j].first - mean_xj)* (in[j].first - new_mean_xj);

                mean_xj = new_mean_xj;
                n_xj += in[j].second;
               
                if (costs[(k-1) * size + j -1] >= 0) {
                    if (min = -1) {
                        min = j;
                        min_cost = d + costs[(k-1) * size + j -1];
                    } else {
                        if( d + costs[(k-1) * size + j -1] < min_cost) {
                            min = j;
                            min_cost = d + costs[(k-1) * size + j -1];
                        }
                    }
                }
                
            }
            
            if (min < 0) {
                // no minimum found, invalid state!
                costs[ k * size + i] = -1;
                 
            } else {
                
                costs[ k * size + i] = min_cost;
                back[ k * size + i ] = min;
            }
            
        }
    }
    
//    logger::Instance()->debug("Matrix \n");
//    for(unsigned int k = 0; k < size; ++k) {    
//        for(unsigned int i = k; i < size; ++i) {
//            logger::Instance()->debug(" " + std::to_string(costs[ k * size + i]) + ";" + std::to_string(back[ k * size + i]));
//        }
//        logger::Instance()->debug("\n");
//    }
    
    
    //################ Backwards ################
    
    // find smallest number of clusters k that contains all points
    unsigned int end = 0;
    for (; end < size; ++end ) {
        if (costs[ end * size + (size-1)] > 0) {
            break;
        }
    }

    unsigned int range_start;
    unsigned int range_end = size - 1;
    for ( int k = end; k>=0; --k) {
        range_start = back[ k * size + range_end];
        
        double mean = 0.0 ;
        unsigned int n = 0;
        for (unsigned int i = range_start; i <= range_end; ++i) {
            
            mean = mean + in[i].second / (double) (n + in[i].second)
                        * (in[i].first - mean);
            n += in[i].second;  
        }
        
        
        if (n < options::Instance()->get_min_junction_coverage()) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Clustered Position below threshold " + std::to_string(round(mean))+".\n");
            #endif
        } else {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Clustered Position added " + std::to_string(round(mean))+".\n");
            #endif
            end_it = out.insert(end_it, round(mean));
        }
        
        range_end = range_start - 1;
    }
    
    // clean up
    delete [] costs;
    delete [] back;
}

void bam_reader::add_to_fixed_clusters(lazy<r_border_set<rpos> > &fixed, greader_list<rpos> &new_clust, r_border_set<rpos>::iterator &pos_mark) {
    
    // for deque, linear join to new list is fastest, change if you switch type of r_border_set
    
    if (fixed.ref().empty()) { // first round, make this fast
        // unfortunately we need to copy and not move, as new is read individually
        // should not cost too much time though, when in doubt profile
        std::copy(new_clust.begin(), new_clust.end(), std::back_inserter(fixed.ref()) );
        pos_mark = fixed.ref().end();
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Copy full.\n");
        #endif
        
    } else {
        
        lazy<r_border_set<rpos> > new_list;
        
        rpos border;
        if (pos_mark == fixed.ref().end()) {
            border = *(--pos_mark);
        } else {
            border = *pos_mark;
        }
        
        // move through both lists in parallel and copy/ move to new list
        r_border_set<rpos>::iterator i1 = fixed.ref().begin();
        greader_list<rpos>::iterator i2 = new_clust.begin();
        
        while (true) {
            
            if (i1 == fixed.ref().end()) {
                 std::copy(i2, new_clust.end(), std::back_inserter(new_list.ref()) );
                 break;
            }
            
            if (i2 == new_clust.end()) {
                 _MOVE_RANGE(i1, fixed.ref().end(), std::back_inserter(new_list.ref()));
                 break;
            }
                       
            if (*i1 < *i2) {
                new_list.ref().push_back(_MOVE(*i1));
                ++i1;
                
            } else {
                new_list.ref().push_back(_MOVE(*i2));
                ++i2;
            }
           
        }
        
        fixed = new_list;
        pos_mark = fixed.ref().lower_bound(fixed.ref().begin(), fixed.ref().end(), border);
    }
    
}

// ######### EXONS #########

void bam_reader::create_raw_exons( chromosome* chrom, greader_list<std::pair<rpos, rpos> > &out, rpos &right) {
    
    greader_list<std::pair<rpos, rpos> > pre_out;
    
    chrom->interval_queue.sort(); // in bam only reads where sorted
    
    //std::vector<rpos> internal;
    std::priority_queue< std::pair<rpos, rcount>,  std::vector<std::pair<rpos, rcount>>, std::greater<std::pair<rpos, rcount> > > end_queue; //= std::priority_queue<rpos,  std::vector<rpos>, std::greater<rpos> >(internal);
    
    unsigned int count = 0;
    bool in = false;
    rpos start;
    for (greader_list<interval>::iterator it = chrom->interval_queue.begin(); it != chrom->interval_queue.end() && (it)->right <= right; ++it) {
    
//        logger::Instance()->debug("Interval " + std::to_string((it)->left) + " "  + std::to_string((it)->right)+".\n");
        
        count += it->parent->count; // increase count for current interval
        end_queue.push( std::make_pair((it)->right, it->parent->count) );
       
        while ((it)->left > end_queue.top().first + 1) {
            count -= end_queue.top().second;
            // minus count because we have already added counts that are right of end_queue pos
            if (in && count - it->parent->count < options::Instance()->get_min_coverage()) {
                in = false;
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Connected area A " + std::to_string(start) + " "  + std::to_string(end_queue.top().first) +".\n");
//            #endif
                pre_out.push_back(std::make_pair(start, end_queue.top().first));        
            }
            
            end_queue.pop();
        }
        // did we find a start?
        if (!in && count >= options::Instance()->get_min_coverage()) {
            in = true;
            start = (it)->left;
        } 
    }
    
    if (in && count >= options::Instance()->get_min_coverage()) {
            // we need to find the end still!
            rpos last = end_queue.top().first;
            while ( count >= options::Instance()->get_min_coverage()) { // guaranteed to enter it once!
                count -= end_queue.top().second;
                last = end_queue.top().first;
                end_queue.pop();
            }
        
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Connected area B " + std::to_string(start) + " "  + std::to_string(chrom->interval_queue.back().right)+".\n");
//            #endif
            pre_out.push_back(std::make_pair(start, last));        
    }
    
    
      // this removes way to much! There actually are exons smaller than that!
//    if (options::Instance()->get_min_junction_anchor() != 0) {
//        for(greader_list<std::pair<rpos, rpos> >::iterator it = out.begin(); it != out.end(); ) { // outs are disconnected by definition
//            if (it->second - it->first + 1 < options::Instance()->get_min_junction_anchor()) {
//                it = out.erase(it);
//            } else {
//                ++it;
//            }
//        }
//    }
    
    greader_list<std::pair<rpos, rpos> >::iterator it = pre_out.begin();
    if( it == pre_out.end() ) {        
        out = pre_out; // kinda useless
        return;
    }
    greader_list<std::pair<rpos, rpos> >::iterator next = it;        
    ++next;
    while( next != pre_out.end() ) { // connect output
        if (next->first - it->second < options::Instance()->get_exon_join_distance()) {
            next->first = it->first;
        } else {
            out.push_back(*it);
        }
        it = next;
        ++next;   
    }
    out.push_back(*it);
}

void bam_reader::trim_exons_2(chromosome* chrom, greader_list<std::pair<rpos, rpos> > &raw, greader_list<rpos> &starts, greader_list<rpos> &ends) {

    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Trim IN.\n");
    #endif
    
    greader_list<exon> split_raw;
    
    greader_list<rpos>::iterator s_it = starts.begin();
    greader_list<rpos>::iterator e_it = ends.begin();
    
    
    for(greader_list<std::pair<rpos, rpos> >::iterator raw_it = raw.begin(); raw_it != raw.end(); ++raw_it) {
        
        #ifdef ALLOW_DEBUG
           logger::Instance()->debug("RAW it " + std::to_string(raw_it->first) + " "  + std::to_string(raw_it->second) +".\n");
        #endif
        
        rpos pos = raw_it->first;
        bool trimmable_start = true;
        
        while (true) {
            
            if (s_it == starts.end() && e_it == ends.end()) { 
                break;
            }
            
            rpos next;
            bool trimmable_end;
            bool next_trimmable_start;
            if (s_it == starts.end() || (e_it != ends.end() && *e_it + 1 < *s_it ) ) {            
                next = *e_it + 1;
                trimmable_end = false;
                next_trimmable_start = true;
                if (next-1 > raw_it->second) {
                    break;
                }
                ++e_it;
            } else if (e_it == ends.end() || *e_it + 1 > *s_it ) {
                next = *s_it;
                trimmable_end = true;
                next_trimmable_start = false;
                if (next-1 > raw_it->second) {
                    break;
                }
                ++s_it;
            } else {
                next = *s_it;
                trimmable_end = false;
                next_trimmable_start = false;
                if (next-1 > raw_it->second) {
                    break;
                }
                ++s_it;
                ++e_it;
            }
            
            split_raw.push_back(exon(pos, next-1, trimmable_start, trimmable_end));
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Trim A " + std::to_string(pos) + " "  + std::to_string(next-1) + " : " + std::to_string(trimmable_start) + " "  + std::to_string(trimmable_end) +".\n");
            #endif
            
            trimmable_start = next_trimmable_start;
            pos = next;
        }
        
        if (pos != raw_it->second) {
            split_raw.push_back(exon(pos, raw_it->second, trimmable_start, true));
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Trim B " + std::to_string(pos) + " "  + std::to_string(raw_it->second) +  " : " + std::to_string(trimmable_start) + " 1.\n");
            #endif
        }
    }
    
    // here we have the raw exons
    bool modified = false;
    for(greader_list<exon>::iterator sr_it = split_raw.begin(); sr_it != split_raw.end(); ++sr_it) {
        
        #ifdef ALLOW_DEBUG
           logger::Instance()->debug("Test Trim " + std::to_string(sr_it->start) + " "  + std::to_string(sr_it->end) +".\n");
        #endif
        
        std::priority_queue< std::pair<rpos, rcount>,  std::vector<std::pair<rpos, rcount>>, std::greater<std::pair<rpos, rcount> > > end_queue;
        
        unsigned int count = 0;
        unsigned int start_count = 0;
        rcount max = 0;
        rcount min = std::numeric_limits<rcount>::max();
        for (greader_list<interval>::iterator it = chrom->interval_queue.begin(); it != chrom->interval_queue.end() && it->left <= sr_it->end; ++it) {
        
            if (it->right < sr_it->start) {
                continue;
            }
            
            count += it->parent->count; // increase count for current interval
            if (it->left <= sr_it->start) start_count += count;
            end_queue.push( std::make_pair((it)->right, it->parent->count) );
            
            while ((it)->left > end_queue.top().first + 1) {
                count -= end_queue.top().second;
                // minus count because we have already added counts that are right of end_queue pos
                end_queue.pop();
            }
            
            if (count > max) {
                max = count;      
            }
            if (count < min) {
                min = count;
            }
        }
        if (count > max) {
            max = count;      
        }
        if (count < min) {
            min = count;
        }
        
        while (!end_queue.empty()) {
            end_queue.pop();
        }
        
        rcount percentile = max*0.1;
        if (percentile < 15) {
            percentile = 15;
        }
        if (percentile > 60) {
           percentile = 60;
        }
        
        #ifdef ALLOW_DEBUG
           logger::Instance()->debug("Max Min Percentile " + std::to_string(max) + " " + std::to_string(min) + " "  + std::to_string(percentile) +".\n");
        #endif
        
        if (min > percentile || max < 60) {
            continue;
        }
        
        greader_list<std::pair<rpos, rpos> > pre_out;

        unsigned int minimal_report = 40;
        unsigned int cut_length = 100;
        
        count = 0;
        bool in = false;
        rpos start;
        for (greader_list<interval>::iterator it = chrom->interval_queue.begin(); it != chrom->interval_queue.end() && it->left <= sr_it->end; ++it) {
        
            if (it->right < sr_it->start) {
                continue;
            }
            
            count += it->parent->count; // increase count for current interval
            end_queue.push( std::make_pair((it)->right, it->parent->count) );

            while ((it)->left > end_queue.top().first + 1) {
                count -= end_queue.top().second;
                // minus count because we have already added counts that are right of end_queue pos
                if (in && count - it->parent->count < percentile) {
                    in = false;
                    
                    if (!pre_out.empty() && start - pre_out.back().second < minimal_report) {
                        pre_out.back().second = end_queue.top().first;
                        #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Extend Area " + std::to_string(end_queue.top().first) +".\n");
                        #endif
                    } else {
                        if (end_queue.top().first - start + 1 >= minimal_report) {
                            pre_out.push_back(std::make_pair(start, end_queue.top().first));
                            
                             #ifdef ALLOW_DEBUG
                                logger::Instance()->debug("Area " + std::to_string(start) + " " + std::to_string(end_queue.top().first) +".\n");
                             #endif
                        }
                    }        
                }

                end_queue.pop();
            }
            
            if (pre_out.size() > 2) {
                 break;       
            }
            
            // did we find a start?
            if (!in && count >= percentile) {
                in = true;
                start = std::max((it)->left, sr_it->start);
            } 
        }
        if (in && count >= percentile) {
            
            if (!pre_out.empty() && start - pre_out.back().second < minimal_report) {
                pre_out.back().second = sr_it->end;
                #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Extend Area " + std::to_string(sr_it->end) +".\n");
                #endif
            } else {
                if (sr_it->end - start + 1 >= minimal_report) {
                    pre_out.push_back(std::make_pair(start, sr_it->end));
                    
                    #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Area " + std::to_string(start) + " " + std::to_string( sr_it->end) +".\n");
                    #endif
                }
            }        
        }
        if (pre_out.size() > 2) {
                 continue;       
        }
        
        // we missuse the fixed parts for trimming information!
        if (pre_out.size() == 1) {
            if (sr_it->start == pre_out.back().first) {
                
//                greader_list<exon>::iterator sr_it_p = sr_it;
//                ++sr_it_p;
                
                if (pre_out.back().second - pre_out.back().first + 1 >= cut_length && sr_it->end - pre_out.back().second >= cut_length && sr_it->fixed_end
                        && get_max_region(chrom, pre_out.back().first, pre_out.back().second) > 50 ) {
//                        && (sr_it_p == split_raw.end() || sr_it->end+1 != sr_it_p->start || get_max_region(chrom, sr_it_p->start, sr_it_p->end) > count + 20 )) {
                    //    && test_decreasing(chrom, pre_out.back().first, pre_out.back().second)) {
                    // we cut this to the end!
                    sr_it->end = pre_out.back().second;
                    modified = true;
                    #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Trimmed End.\n");
                    #endif
                }
            } else if (sr_it->end == pre_out.back().second) {
                
//                greader_list<exon>::iterator sr_it_p = sr_it;
//                if (sr_it != split_raw.begin()) --sr_it_p;
                
                if (pre_out.back().second - pre_out.back().first + 1 >= cut_length && pre_out.back().first - sr_it->start >= cut_length && sr_it->fixed_start 
                       && get_max_region(chrom, pre_out.back().first, pre_out.back().second) > 50 ){
                    //   && (sr_it == split_raw.begin() || sr_it_p->end+1 != sr_it->start || get_max_region(chrom, sr_it_p->start, sr_it_p->end) > start_count + 20)) {
                    //    && test_increasing(chrom, pre_out.back().first, pre_out.back().second)) {
                    // we cut this to the start!
                    sr_it->start = pre_out.back().first;
                    modified = true;
                    #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Trimmed Start.\n");
                    #endif
                }
            }
        } else {
//            if (sr_it->start == pre_out.front().first && sr_it->end == pre_out.back().second
//                    && pre_out.back().second - pre_out.back().first + 1 >= cut_length 
//                    && pre_out.front().second - pre_out.front().first + 1 >= cut_length
//                    && pre_out.back().first - pre_out.front().second +1 >= cut_length
//                    && get_max_region(chrom, pre_out.front().first, pre_out.front().second) > 50
//                    && get_max_region(chrom, pre_out.back().first, pre_out.back().second) > 50) {
//                   // && test_decreasing(chrom, pre_out.front().first, pre_out.front().second)
//                   // && test_increasing(chrom, pre_out.back().first, pre_out.back().second)) {
//                   // && sr_it->fixed_end && sr_it->fixed_start) {
//             
//                sr_it = split_raw.insert(sr_it, exon(sr_it->start, pre_out.front().second));
//                ++sr_it;
//                sr_it->start = pre_out.back().first;
//                modified = true;
//                #ifdef ALLOW_DEBUG
//                        logger::Instance()->debug("Double Trimmed.\n");
//                #endif
//            }
        }
    }
    
    if (modified) {
        raw.clear();
        for(greader_list<exon>::iterator sr_it = split_raw.begin(); sr_it != split_raw.end(); ++sr_it) {
            if (!raw.empty() && raw.back().second + 1 == sr_it->start) {
                raw.back().second = sr_it->end;
            } else {
                raw.push_back(std::make_pair(sr_it->start, sr_it->end));
            }
        }
    }
    
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Trim OUT.\n");
    #endif
}

void bam_reader::trim_exons_1(chromosome* chrom, greader_list<std::pair<rpos, rpos> > &raw, greader_list<rpos> &starts, greader_list<rpos> &ends) {

    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Trim 1 IN.\n");
    #endif
    
    greader_list<exon> split_raw;
    
    greader_list<rpos>::iterator s_it = starts.begin();
    greader_list<rpos>::iterator e_it = ends.begin();
    
    
    for(greader_list<std::pair<rpos, rpos> >::iterator raw_it = raw.begin(); raw_it != raw.end(); ++raw_it) {
        
        #ifdef ALLOW_DEBUG
           logger::Instance()->debug("RAW it " + std::to_string(raw_it->first) + " "  + std::to_string(raw_it->second) +".\n");
        #endif
        
        rpos pos = raw_it->first;
        bool trimmable_start = true;
        
        while (true) {
            
            if (s_it == starts.end() && e_it == ends.end()) { 
                break;
            }
            
            rpos next;
            bool trimmable_end;
            bool next_trimmable_start;
            if (s_it == starts.end() || (e_it != ends.end() && *e_it + 1 < *s_it ) ) {            
                next = *e_it + 1;
                trimmable_end = false;
                next_trimmable_start = true;
                if (next-1 > raw_it->second) {
                    break;
                }
                ++e_it;
            } else if (e_it == ends.end() || *e_it + 1 > *s_it ) {
                next = *s_it;
                trimmable_end = true;
                next_trimmable_start = false;
                if (next-1 > raw_it->second) {
                    break;
                }
                ++s_it;
            } else {
                next = *s_it;
                trimmable_end = false;
                next_trimmable_start = false;
                if (next-1 > raw_it->second) {
                    break;
                }
                ++s_it;
                ++e_it;
            }
            
            split_raw.push_back(exon(pos, next-1, trimmable_start, trimmable_end));
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Trim A " + std::to_string(pos) + " "  + std::to_string(next-1) + " : " + std::to_string(trimmable_start) + " "  + std::to_string(trimmable_end) +".\n");
            #endif
            
            trimmable_start = next_trimmable_start;
            pos = next;
        }
        
        if (pos != raw_it->second) {
            split_raw.push_back(exon(pos, raw_it->second, trimmable_start, true));
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Trim B " + std::to_string(pos) + " "  + std::to_string(raw_it->second) +  " : " + std::to_string(trimmable_start) + " 1.\n");
            #endif
        }
    }
    
    // here we have the raw exons
    bool modified = false;
    rcount total_max = 0;
    std::deque<std::deque<std::pair<rpos, rpos> > > all_regions;
    std::deque<std::deque<float> > all_aves;
    
    for(greader_list<exon>::iterator sr_it = split_raw.begin(); sr_it != split_raw.end(); ++sr_it) {
        
        #ifdef ALLOW_DEBUG
           logger::Instance()->debug("Test Trim " + std::to_string(sr_it->start) + " "  + std::to_string(sr_it->end) +".\n");
        #endif
        
//        if (! (!sr_it->fixed_start && !sr_it->fixed_end) ) {
//            continue;
//        }  
           
        all_regions.push_back(std::deque<std::pair<rpos, rpos> >());   
        std::deque<std::pair<rpos, rpos> >& regions = all_regions.back();
        all_aves.push_back(std::deque<float>());
        std::deque<float>& aves = all_aves.back();
        
        // find 0 seperated parts
        std::priority_queue< std::pair<rpos, rcount>,  std::vector<std::pair<rpos, rcount>>, std::greater<std::pair<rpos, rcount> > > end_queue;
        
        unsigned int count = 0;
        rcount bases = 0;
        rpos start_pos = sr_it->start;
        for (greader_list<interval>::iterator it = chrom->interval_queue.begin(); it != chrom->interval_queue.end() && it->left <= sr_it->end; ++it) {
        
            if (it->right < sr_it->start) {
                continue;
            }

            end_queue.push( std::make_pair((it)->right, it->parent->count) );

            rpos last_pos = 0;
            while ((it)->left > end_queue.top().first + 1) {
                count -= end_queue.top().second;
                last_pos = end_queue.top().first;
                // minus count because we have already added counts that are right of end_queue pos
                end_queue.pop();
            }
            
            if (count > total_max) {
                total_max = count;      
            }
            
            if (count == 0 && last_pos != 0) {
                bases += (std::min((it)->right, last_pos) - std::max((it)->left, start_pos) + 1 ) * it->parent->count;
                rpos length = last_pos - start_pos + 1;
//                logger::Instance()->debug("B1 " + std::to_string(bases) + " "  + std::to_string(length) + ".\n");
                float average = bases/(float)length;
                
                regions.push_back(std::make_pair(start_pos, last_pos));
                aves.push_back(average);
                
                bases = 0;
                start_pos = (it)->left;
            }

            bases += (std::min((it)->right, sr_it->end) - std::max((it)->left, start_pos) + 1 ) * it->parent->count;
            count += it->parent->count; // increase count for current interval
        }

        rpos length = sr_it->end - start_pos + 1;
//        logger::Instance()->debug("B2 " + std::to_string(bases) + " "  + std::to_string(length) + ".\n");
        float average = bases/(float)length;
           
        regions.push_back(std::make_pair(start_pos, sr_it->end));
        aves.push_back(average);   
        
    }
    
    std::deque<std::deque<std::pair<rpos, rpos> > >::iterator ari = all_regions.begin();
    std::deque<std::deque<float> >::iterator aai = all_aves.begin();
    for(greader_list<exon>::iterator sr_it = split_raw.begin(); sr_it != split_raw.end(); ++sr_it, ++ari, ++aai) {
        
        std::deque<std::pair<rpos, rpos> >& regions = *ari;
        std::deque<float>& aves = *aai;
        
        if (regions.size() == 1) {
            continue;
        }
        
        std::deque<std::pair<rpos, rpos> >::iterator ri_a = regions.begin();
        std::deque<float>::iterator ai_a = aves.begin(); 
        std::deque<std::pair<rpos, rpos> >::iterator ri_b = ri_a;
        ++ri_b;
         std::deque<float>::iterator ai_b = ai_a;
        ++ai_b;
        for (; ri_b != regions.end(); ++ri_b, ++ri_a, ++ai_b, ++ai_a) {
            float max = std::max(*ai_a, *ai_b);
            float min = std::min(*ai_a, *ai_b);

            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Split Test at " + std::to_string(ri_a->second) + " _ " + std::to_string(ri_b->first) + " " + std::to_string(min) + " " + std::to_string(max) + ".\n");
            #endif
            if ( (max * options::Instance()->get_trimming_rate() > min && min < 20) || (min > 25 ) ) {
                // we do want to seperate those
                sr_it = split_raw.insert(sr_it, exon(sr_it->start, ri_a->second));
                ++sr_it;
                sr_it->start = ri_b->first;
                modified = true;
                
                #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Split Gap at " + std::to_string(ri_a->second) + " _ " + std::to_string(ri_b->first) + ".\n");
                #endif
            }
        }
    }
    
    if (modified) {
        raw.clear();
        for(greader_list<exon>::iterator sr_it = split_raw.begin(); sr_it != split_raw.end(); ++sr_it) {
            if (!raw.empty() && raw.back().second + 1 == sr_it->start) {
                raw.back().second = sr_it->end;
            } else {
                raw.push_back(std::make_pair(sr_it->start, sr_it->end));
            }
        }
    }
    
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Trim OUT.\n");
    #endif
}


bool bam_reader::test_increasing(chromosome* chrom, rpos start, rpos end) {
    
    std::priority_queue< std::pair<rpos, rcount>,  std::vector<std::pair<rpos, rcount>>, std::greater<std::pair<rpos, rcount> > > end_queue;
        
    unsigned int count = 0;
    rcount max = 0;
    for (greader_list<interval>::iterator it = chrom->interval_queue.begin(); it != chrom->interval_queue.end() && it->left <= end; ++it) {

        if (it->right < start) {
            continue;
        }

        count += it->parent->count; // increase count for current interval
        end_queue.push( std::make_pair((it)->right, it->parent->count) );

        while ((it)->left > end_queue.top().first + 1) {
            count -= end_queue.top().second;
            // minus count because we have already added counts that are right of end_queue pos
            end_queue.pop();
        }
        if (count > max) {
            max = count;      
        //} else if ( (max - count) * 100 / max > 25 && max - count > 20) {
        } else if ( max - count > 15) {
            return false;
        }
    }
    
    return true;
}

bool bam_reader::test_decreasing(chromosome* chrom, rpos start, rpos end) {
    
    std::priority_queue< std::pair<rpos, rcount>,  std::vector<std::pair<rpos, rcount>>, std::greater<std::pair<rpos, rcount> > > end_queue;
        
    unsigned int count = 0;
    rcount min = std::numeric_limits<rcount>::max();
    for (greader_list<interval>::iterator it = chrom->interval_queue.begin(); it != chrom->interval_queue.end() && it->left <= end; ++it) {

        if (it->right < start) {
            continue;
        }

        count += it->parent->count; // increase count for current interval
        end_queue.push( std::make_pair((it)->right, it->parent->count) );

        while ((it)->left > end_queue.top().first + 1) {
            count -= end_queue.top().second;
            // minus count because we have already added counts that are right of end_queue pos
            end_queue.pop();
        }
        if (count < min) {
            min = count;      
        //} else if ( (count - min) * 100 / count > 25  && count - min > 20 ) {
        } else if ( count - min > 15 ) {
            return false;
        }
    }
    
    return true;
}


rcount bam_reader::get_max_region(chromosome* chrom, rpos start, rpos end) {
    
    std::priority_queue< std::pair<rpos, rcount>,  std::vector<std::pair<rpos, rcount>>, std::greater<std::pair<rpos, rcount> > > end_queue;
        
    unsigned int count = 0;
    rcount max = 0;
    for (greader_list<interval>::iterator it = chrom->interval_queue.begin(); it != chrom->interval_queue.end() && it->left <= end; ++it) {

        if (it->right < start) {
            continue;
        }

        count += it->parent->count; // increase count for current interval
        end_queue.push( std::make_pair((it)->right, it->parent->count) );

        while ((it)->left > end_queue.top().first + 1) {
            count -= end_queue.top().second;
            // minus count because we have already added counts that are right of end_queue pos
            end_queue.pop();
        }
        if (count > max) {
            max = count;      
        }
    }
    
    if (count > max) {
        max = count;      
    }
    
    #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Found MAX " + std::to_string(max) + " " + std::to_string(start) + " " + std::to_string(end)  +".\n");
    #endif
    
    return max;
}


bool bam_reader::get_average_to_first_zero_from_left(chromosome* chrom, rpos start, rpos end, float& average, rpos & length, rpos & new_end, rpos & next_start) {
    
    std::priority_queue< std::pair<rpos, rcount>,  std::vector<std::pair<rpos, rcount>>, std::greater<std::pair<rpos, rcount> > > end_queue;
        
    unsigned int count = 0;
    rcount bases = 0;

    for (greader_list<interval>::iterator it = chrom->interval_queue.begin(); it != chrom->interval_queue.end() && it->left <= end; ++it) {

        if (it->right < start) {
            continue;
        }

        end_queue.push( std::make_pair((it)->right, it->parent->count) );

        rpos last_pos = 0;
        while ((it)->left > end_queue.top().first + 1) {
            count -= end_queue.top().second;
            last_pos = end_queue.top().first;
            // minus count because we have already added counts that are right of end_queue pos
            end_queue.pop();
        }
        if (count == 0 && last_pos != 0) {
            bases += (std::min((it)->right, last_pos) - std::max((it)->left, start) + 1 ) * it->parent->count;
            length = last_pos - start + 1;
            average = bases/(float)length;
            new_end = last_pos;
            next_start = (it)->left;
            return true;
        }
        
        bases += (std::min((it)->right, end) - std::max((it)->left, start) + 1 ) * it->parent->count;
        count += it->parent->count; // increase count for current interval

    }

    length = end - start + 1;
    average = bases/(float)length;
    
    return false; 
}

bool bam_reader::get_average_to_first_zero_from_right(chromosome* chrom, rpos start, rpos end, float& average, rpos & length, rpos & new_start, rpos & next_end) {
    
    std::priority_queue< std::pair<rpos, rcount>,  std::vector<std::pair<rpos, rcount>>, std::greater<std::pair<rpos, rcount> > > end_queue;
        
    unsigned int count = 0;
    rcount bases = 0;

    for (greader_list<interval>::reverse_iterator it = chrom->interval_queue.rbegin(); it != chrom->interval_queue.rend() && it->right >= start; ++it) {

        if (it->left > end) {
            continue;
        }

        end_queue.push( std::make_pair((it)->left, it->parent->count) );

        rpos last_pos = 0;
        while ((it)->right + 1 < end_queue.top().first) {
            count -= end_queue.top().second;
            last_pos = end_queue.top().first;
            // minus count because we have already added counts that are right of end_queue pos
            end_queue.pop();
        }
        if (count == 0 && last_pos != 0) {
            bases += (std::min((it)->right, end) - std::max((it)->left, last_pos) + 1 ) * it->parent->count;
            length = last_pos - start + 1;
            average = bases/(float)length;
            new_start = last_pos;
            next_end = (it)->left;
            return true;
        }
        
        bases += (std::min((it)->right, end) - std::max((it)->left, start) + 1 ) * it->parent->count;
        count += it->parent->count; // increase count for current interval

    }

    length = end - start + 1;
    average = bases/(float)length;
    
    return false; 
}

void bam_reader::solidify_raw_exons_ends(chromosome* chrom, greader_list<std::pair<rpos, rpos> > &raw, greader_list<rpos> &starts, greader_list<rpos> &ends) {
    
    const unsigned int extend = options::Instance()->get_max_pos_extend();
    r_border_set<rpos>::iterator fs_it = chrom->fixed_exon_starts->begin();
    r_border_set<rpos>::iterator fe_it = chrom->fixed_exon_ends->begin();
    
    greader_list<rpos>::iterator cs_it = starts.begin();
    greader_list<rpos>::iterator ce_it = ends.begin();
    
    for(greader_list<std::pair<rpos, rpos> >::iterator raw_it = raw.begin(); raw_it != raw.end(); ) {
        
        rpos start = raw_it->first;
        rpos end = raw_it->second;
        
        fs_it = std::lower_bound(fs_it, chrom->fixed_exon_starts.ref().end(), start);
        fe_it = std::lower_bound(fe_it, chrom->fixed_exon_ends.ref().end(), end);
        cs_it = std::lower_bound(cs_it, starts.end(), start);
        ce_it = std::lower_bound(ce_it, ends.end(), end);
        
        // we look for updated starts
        if (fs_it!= chrom->fixed_exon_starts.ref().end() && start + extend >= *fs_it && start <= *fs_it + extend) {
            // match to fixed
            raw_it->first = *fs_it;
        } else if (cs_it!= starts.end() && start + extend >= *cs_it && start <= *cs_it + extend) {
            // match to cluster
            raw_it->first = *cs_it;
        } else if (fs_it!= chrom->fixed_exon_starts.ref().begin() && start + extend >= *(fs_it-1) && start <= *(fs_it-1) + extend) {
            // match to fixed previous
            raw_it->first = *(fs_it-1);
        } else if (cs_it!= starts.begin() && start + extend >= *(cs_it-1) && start <= *(cs_it-1) + extend) {
            // match to cluster previous
            raw_it->first = *(cs_it-1);
        }
               
        // we look for updated starts
        if (fe_it!= chrom->fixed_exon_ends.ref().end() && end + extend >= *fe_it && end <= *fe_it + extend) {
            // match to fixed
            raw_it->second = *fe_it;
        } else if (ce_it!= ends.end() && end + extend >= *ce_it && end <= *ce_it + extend) {
            // match to cluster
            raw_it->second = *ce_it;
        } else if (fe_it!= chrom->fixed_exon_ends.ref().begin() && end + extend >= *(fe_it-1) && end <= *(fe_it-1) + extend) {
            // match to fixed previous
            raw_it->second = *(fe_it-1);
        } else if (ce_it!= ends.begin() && end + extend >= *(ce_it-1) && end <= *(ce_it-1) + extend) {
            // match to cluster previous
            raw_it->second = *(ce_it-1);
        }
        
        
        if (raw_it->first > raw_it->second || raw_it->second - raw_it->first < options::Instance()->get_min_raw_exon_size()) {
            raw_it = raw.erase(raw_it);
        } else {
            ++raw_it;
        }
    }
}


void bam_reader::update_existing_exons( connected* connected, chromosome* chrom, greader_list<std::pair<rpos, rpos> > &raw, rpos &left,  rpos &right) {
    
    #ifdef ALLOW_DEBUG

            logger::Instance()->debug("Update existing. " + std::to_string(connected->fossil_exons.ref().size()) + "\n");
            for ( greader_list<exon* >::iterator fix_it = connected->fossil_exons.ref().begin(); fix_it != connected->fossil_exons.ref().end(); ++fix_it) {
                logger::Instance()->debug("Fossil " + std::to_string((*fix_it)->start) + " " + std::to_string((*fix_it)->end) + "\n");
            }
           for (  greader_list<std::pair<rpos, rpos> >::iterator raw_it = raw.begin(); raw_it != raw.end(); ++raw_it) {
               logger::Instance()->debug("Raw " + std::to_string(raw_it->first) + " " + std::to_string(raw_it->second) + "\n");
           }
    #endif
    
    greader_list<std::pair<rpos, rpos> >::iterator raw_it = raw.begin();
    greader_list<exon* >::iterator fix_it = connected->fossil_exons.ref().begin();
    
    #ifdef ALLOW_DEBUG
    if (fix_it != connected->fossil_exons.ref().end()) logger::Instance()->debug("Fossil Start " + std::to_string((*fix_it)->start) + " " + std::to_string((*fix_it)->end) + "\n");
    #endif
    
    // move fix to first in range
    while (fix_it != connected->fossil_exons.ref().end() && (*fix_it)->end < left) {
        ++fix_it;
        
        #ifdef ALLOW_DEBUG
        if (fix_it != connected->fossil_exons.ref().end())
            logger::Instance()->debug("Fossil move " + std::to_string((*fix_it)->start) + " " + std::to_string((*fix_it)->end) + "\n");
        #endif
    }
    
    for ( ; raw_it !=  raw.end(); ++raw_it) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Raw it " + std::to_string(raw_it->first) + " " + std::to_string(raw_it->second) + "\n");
        #endif
              
        // find fist possibly intersecting
        while (fix_it != connected->fossil_exons.ref().end() && raw_it->first > (*fix_it)->end) { 
            ++fix_it;
        }
        
        if (fix_it == connected->fossil_exons.ref().end()) {
            chrom->fossil_exons.push_back(exon(raw_it->first, raw_it->second)); 
            connected->fossil_exons.ref().push_back(&chrom->fossil_exons.back());
            fix_it = connected->fossil_exons.ref().end();
            #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Insert without competing.\n");
            #endif
            
            continue;
        }
        
        if (raw_it->second < (*fix_it)->start ) {
            // this means no overlap! so just add it
            
            // add to new exon to the right and stay at new insert
            chrom->fossil_exons.push_back(exon(raw_it->first, raw_it->second));
            fix_it = connected->fossil_exons.ref().insert(fix_it, &chrom->fossil_exons.back());
            ++fix_it;
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Add new Exon in between existing without touch: " + std::to_string(raw_it->first) + " - "+std::to_string(raw_it->second)+".\n");
            #endif
            
            continue;
        }
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Fossil pre " + std::to_string((*fix_it)->start) + " " + std::to_string((*fix_it)->end) + "\n");
        #endif
        if (raw_it->first < (*fix_it)->start && raw_it->second >= (*fix_it)->start ) { // we have an actual overlap top the front
            // this by the iteration is the first such overlap, therefore no previous to consider
            // hence we extend by the given length
            
            if((*fix_it)->fixed_start) { // fix, so add new exon
               
               // add to the left then return to current position
                
               chrom->fossil_exons.push_back(exon(raw_it->first, (*fix_it)->start-1));
               fix_it = connected->fossil_exons.ref().insert(fix_it, &chrom->fossil_exons.back());
               (*fix_it)->fixed_end = true;
                ++fix_it;
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Add new Exon before existing: " + std::to_string(raw_it->first) + " - "+std::to_string((*fix_it)->start-1)+".\n");
                #endif
                
            } else {
                // just extend existing one
                (*fix_it)->start = raw_it->first;
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Extend existing Exon to left: " + std::to_string((*fix_it)->start) + " - "+std::to_string((*fix_it)->end)+".\n");
                #endif
            }
        }
     
        // now loop till end and fix possible holes
        // do we overlap to the end?
        while (fix_it != connected->fossil_exons.ref().end() && raw_it->second > (*fix_it)->end) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Right \n");
            #endif
                        
            // extend to right, see for next element
            greader_list<exon* >::iterator next = fix_it;
            ++next;
            
            if (next == connected->fossil_exons.ref().end() || (*next)->start > raw_it->second) {
                // overlap, but no following exon in the overlap

                if ((*fix_it)->fixed_end) {
                     // add to new exon to the right and stay at new insert
                    chrom->fossil_exons.push_back(exon((*fix_it)->end+1, raw_it->second));
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Add new Exon after existing: " + std::to_string((*fix_it)->end+1) + " - "+std::to_string(raw_it->second)+".\n");
                    #endif
                    fix_it = connected->fossil_exons.ref().insert(next, &chrom->fossil_exons.back());
                    (*fix_it)->fixed_start = true;
                    
                } else {
                    // modify existing
                    (*fix_it)->end = raw_it->second;
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Extend existing Exon to right: " + std::to_string((*fix_it)->start) + " - "+std::to_string((*fix_it)->end)+".\n");
                    #endif
                }
                ++fix_it;
            } else {

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Test Fix " + std::to_string((*fix_it)->fixed_start) + " "  + std::to_string((*fix_it)->fixed_end) + " " + std::to_string((*next)->fixed_start) + " " + std::to_string((*next)->fixed_end) + "\n");
                #endif

                // this means we have to close the gap between two exons
                if ( (*fix_it)->fixed_end && (*next)->fixed_start) {
                    // insert new between two 

                     #ifdef ALLOW_DEBUG
                     logger::Instance()->debug("Add new Exon in between existing: " + std::to_string((*fix_it)->end+1) + " - "+std::to_string((*next)->start-1)+".\n");
                     #endif
                     
                     if ((*fix_it)->end +1 != (*next)->start) {
                        chrom->fossil_exons.push_back(exon((*fix_it)->end+1, (*next)->start-1));
                        fix_it = connected->fossil_exons.ref().insert(next, &chrom->fossil_exons.back());
                        (*fix_it)->fixed_start = true;
                        (*fix_it)->fixed_end = true;
                     }
                     ++fix_it;
                } else if (!(*fix_it)->fixed_end && (*next)->fixed_start) {
                    // extend and fix left
                    (*fix_it)->end = (*next)->start-1;
                    (*fix_it)->fixed_end = true;

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Extend existing Exon to left fixed: " + std::to_string((*fix_it)->start) + " - "+std::to_string((*fix_it)->end)+".\n");
                    #endif

                    ++fix_it;
                } else if ( (*fix_it)->fixed_end && !(*next)->fixed_start) {
                    (*next)->start = (*fix_it)->end+1;
                    (*next)->fixed_start = true;

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Extend existing Exon to right fixed: " + std::to_string((*fix_it)->start) + " - "+std::to_string((*next)->end)+".\n");
                    #endif

                    ++fix_it;
                } else {
                    // we need to merge two separate exons...
                    // get everything in next and erase itr

                    greader_refsorted_list<raw_atom* >::iterator single_first = connected->atoms.ref().end();
                    greader_refsorted_list<raw_atom* >::iterator single_next = connected->atoms.ref().end();
                    greader_refsorted_list<raw_atom* >::iterator both = connected->atoms.ref().end();
                    
                    unsigned int count = 0;
                    for (greader_refsorted_list<raw_atom* >::iterator m = connected->atoms.ref().begin(); m !=  connected->atoms.ref().end() ; ++m) {
                        
                        if( (*m)->exons->size() > 2) {
                            continue;
                        }
                        
                        greader_refsorted_list<exon*>::iterator first_s = (*m)->exons.ref().find(*fix_it);
                        greader_refsorted_list<exon*>::iterator next_s = (*m)->exons.ref().find(*next);
                        
                        if (first_s != (*m)->exons.ref().end() && next_s != (*m)->exons.ref().end()) {
                            both = m;
                            ++count;
                        } else if (first_s != (*m)->exons.ref().end() ) {
                            single_first = m;
                            ++count;
                        } else if (next_s != (*m)->exons.ref().end() ) {
                            single_next = m;
                            ++count;
                        } 
                    }
                    
                    if (single_next != connected->atoms.ref().end() && (single_first != connected->atoms.ref().end() || both != connected->atoms.ref().end())) {
                        // switch out in paired
                        for (greader_refsorted_list<raw_atom* >::iterator m = connected->atoms.ref().begin(); m !=  connected->atoms.ref().end() ; ++m) {
                            
                            // single partner
                            paired_map<raw_atom*, rcount >::iterator ri = (*m)->paired.find(*single_first);
                            if (ri != (*m)->paired.end()) {
                                (*m)->paired[*single_next] += ri->second;
                                (*m)->paired.erase(ri);
                            }
                            // both partner
                            ri = (*m)->paired.find(*both);
                            if (ri != (*m)->paired.end()) {
                                (*m)->paired[*single_next] += ri->second;
                                (*m)->paired.erase(ri);
                            }
                            
                        }
                        if (single_first != connected->atoms.ref().end()) connected->atoms.ref().erase(single_first);
                        if (both != connected->atoms.ref().end()) connected->atoms.ref().erase(both);
                        
                    } else  if (single_first != connected->atoms.ref().end() && both != connected->atoms.ref().end()) {
                        for (greader_refsorted_list<raw_atom* >::iterator m = connected->atoms.ref().begin(); m !=  connected->atoms.ref().end() ; ++m) {
                            
                            // both partner
                            paired_map<raw_atom*, rcount >::iterator ri = (*m)->paired.find(*both);
                            if (ri != (*m)->paired.end()) {
                                (*m)->paired[*single_first] += ri->second;
                                (*m)->paired.erase(ri);
                            }
                            
                        }
                        connected->atoms.ref().erase(both);
                    }
                        
                    (*next)->start = (*fix_it)->start;    
                        
                    // sadly we need to resort...
                    lazy<greader_refsorted_list<raw_atom* > > new_atom_order;
                    
                    for (greader_refsorted_list<raw_atom* >::iterator m = connected->atoms.ref().begin(); m !=  connected->atoms.ref().end() ; ++m) {
                        
                        greader_refsorted_list<exon*>::iterator old_left = (*m)->exons.ref().find(*fix_it);
                        if (old_left != (*m)->exons.ref().end()) {
                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Switch " + std::to_string( (*fix_it)->start ) + "-" + std::to_string((*fix_it)->end) + " ; " + std::to_string( (*old_left)->start ) + "-" + std::to_string((*old_left)->end) + " ; " + std::to_string( (*next)->start ) + "-" + std::to_string((*next)->end) + ".\n");
                            #endif
                            (*m)->exons.ref().erase(*old_left); // erase if found
                            (*m)->exons.ref().insert(*next); // insert other instead
                        }
                        
                        new_atom_order->insert(*m);
                    }
                    connected->atoms = new_atom_order;
                     
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Joined two existing Exons: " + std::to_string((*fix_it)->start) + " - "+std::to_string((*fix_it)->end)+".\n");
                    #endif
                    // now erase fix_it from real list
                    greader_refsafe_list<exon>::iterator rem = std::find(chrom->fossil_exons.begin(), chrom->fossil_exons.end(), **fix_it); 
                    fix_it = connected->fossil_exons.ref().erase(fix_it);                
                    chrom->fossil_exons.erase(rem);
                      
                }

            }
            
        } 
    }
    
}

void bam_reader::split_exons( connected* connected, chromosome* chrom, greader_list<rpos> &splits, rpos &left,  rpos &right, int correction) {
    
    // move linear through both and split exons
    
    greader_list<exon* >::iterator e_it = connected->fossil_exons.ref().begin();
    greader_list<rpos>::iterator s_it = splits.begin();
    
    
    while (e_it != connected->fossil_exons.ref().end() && (*e_it)->end < left ) {
        ++e_it;
    }
    
    unsigned int max_extend = options::Instance()->get_max_pos_extend();
    
    for (; e_it != connected->fossil_exons.ref().end() && (*e_it)->start <= right && s_it != splits.end(); ++e_it ) {
        while (*s_it <= (*e_it)->end && s_it != splits.end()) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test " + std::to_string(*s_it) + " " + std::to_string((*e_it)->start) + ":" + std::to_string((*e_it)->fixed_start) + "-" + std::to_string((*e_it)->end) + ":" + std::to_string((*e_it)->fixed_end) + ".\n");
            #endif
            
            if ( !(*e_it)->fixed_start && correction == 1 &&
                    ((*s_it >= (*e_it)->start  && *s_it - (*e_it)->start <= max_extend) || (*s_it < (*e_it)->start  &&  (*e_it)->start - *s_it <= max_extend))) {
                
                 if (*s_it > (*e_it)->start) {
                    for (greader_refsorted_list<raw_atom* >::iterator a = connected->atoms.ref().begin(); a !=  connected->atoms.ref().end() ; ++a) {
                        for (std::map< rpos,rcount >::iterator li = (*a)->lefts->begin(); li != (*a)->lefts->end(); ) {
                            if ( li->first < *s_it && li->first >= (*e_it)->start) {
                                (*a)->lefts.ref()[*s_it] += li->second;
                                li = (*a)->lefts->erase(li);
                            } else {
                                ++li;
                            }
                        }

                        for (std::map< rpos,rcount >::iterator ri = (*a)->rights->begin(); ri != (*a)->rights->end(); ) {
                            if ( ri->first < *s_it && ri->first >= (*e_it)->start) {
                                (*a)->rights.ref()[*s_it] += ri->second;
                                ri = (*a)->rights->erase(ri);
                            } else {
                                ++ri;
                            }
                        }

                        for (std::map< rpos,rcount >::iterator li = (*a)->hole_starts->begin(); li != (*a)->hole_starts->end(); ) {
                            if ( li->first < *s_it && li->first >= (*e_it)->start) {
                                (*a)->hole_starts.ref()[*s_it] += li->second;
                                li = (*a)->hole_starts->erase(li);
                            } else {
                                ++li;
                            }
                        }
                        
                        for (std::map< rpos,rcount >::iterator ri = (*a)->hole_ends->begin(); ri != (*a)->hole_ends->end(); ) {
                            if ( ri->first < *s_it && ri->first >= (*e_it)->start) {
                                (*a)->hole_ends.ref()[*s_it] += ri->second;
                                ri = (*a)->hole_ends->erase(ri);
                            } else {
                                ++ri;
                            }
                        }
                    }
                }
                
                (*e_it)->fixed_start = true;
                (*e_it)->start = *s_it;
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Set Fixed Start.\n");
                #endif
            }
            if ( !(*e_it)->fixed_end && correction == 0 &&
                    ((*s_it >= (*e_it)->end  && *s_it - (*e_it)->end <= max_extend) || (*s_it < (*e_it)->end  &&  (*e_it)->end - *s_it <= max_extend))) {
                
                if (*s_it < (*e_it)->end) {
                    for (greader_refsorted_list<raw_atom* >::iterator a = connected->atoms.ref().begin(); a !=  connected->atoms.ref().end() ; ++a) {
                        for (std::map< rpos,rcount >::iterator li = (*a)->lefts->begin(); li != (*a)->lefts->end(); ) {
                            if ( li->first > *s_it && li->first <= (*e_it)->end) {
                                (*a)->lefts.ref()[*s_it] += li->second;
                                li = (*a)->lefts->erase(li);
                            } else {
                                ++li;
                            }
                        }

                        for (std::map< rpos,rcount >::iterator ri = (*a)->rights->begin(); ri != (*a)->rights->end(); ) {
                            if ( ri->first > *s_it && ri->first <= (*e_it)->end) {
                                (*a)->rights.ref()[*s_it] += ri->second;
                                ri = (*a)->rights->erase(ri);
                            } else {
                                ++ri;
                            }
                        }

                        for (std::map< rpos,rcount >::iterator li = (*a)->hole_starts->begin(); li != (*a)->hole_starts->end(); ) {
                            if ( li->first > *s_it && li->first <= (*e_it)->end) {
                                (*a)->hole_starts.ref()[*s_it] += li->second;
                                li = (*a)->hole_starts->erase(li);
                            } else {
                                ++li;
                            }
                        }
                        
                        for (std::map< rpos,rcount >::iterator ri = (*a)->hole_ends->begin(); ri != (*a)->hole_ends->end(); ) {
                            if ( ri->first > *s_it && ri->first <= (*e_it)->end) {
                                (*a)->hole_ends.ref()[*s_it] += ri->second;
                                ri = (*a)->hole_ends->erase(ri);
                            } else {
                                ++ri;
                            }
                        }
                    }
                }
                 
                (*e_it)->fixed_end = true;
                (*e_it)->end = *s_it;
                 
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Set Fixed End.\n");
                #endif
            }
                      
            if (*s_it >= (*e_it)->start && *s_it+1-correction - (*e_it)->start > max_extend && (*e_it)->end - *s_it > max_extend) {
                // we have an overlap AND need to split 
                // NOTE: if points are outside of exon, they are rejected for splitting as they be construction lie within max_extend
                
                // insert new exon to the left
                chrom->split_exon(*s_it-correction, e_it, connected);
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Split Exon at " + std::to_string(*s_it) + ".\n");
                #endif
                
            } else if (*s_it < (*e_it)->start && (*e_it)->start - *s_it > max_extend 
                    && ( e_it == connected->fossil_exons.ref().begin() || *s_it - (*(e_it-1))->end > max_extend)) {
                // in between, remove this one
                s_it = splits.erase(s_it);
                continue;
            }
            ++s_it;
        }
    }
    --e_it;
    for ( ; s_it != splits.end(); ) {
        
        // TODO: simplify
        if ( !(*e_it)->fixed_start && 
                    ((*s_it >= (*e_it)->start  && *s_it - (*e_it)->start <= max_extend) || (*s_it < (*e_it)->start  &&  (*e_it)->start - *s_it <= max_extend))) {
                (*e_it)->fixed_start = true;
        }
        if ( !(*e_it)->fixed_end && 
                    ((*s_it >= (*e_it)->end  && *s_it - (*e_it)->end <= max_extend) || (*s_it < (*e_it)->end  &&  (*e_it)->end - *s_it <= max_extend))) {
                (*e_it)->fixed_end = true;
         }
        
        if (*s_it >= (*e_it)->end && *s_it - (*(e_it))->end > max_extend) {
           s_it = splits.erase(s_it);
        } else {
            ++s_it;
        }
    }
}

//############ Fragments ##############

greader_list<connected>::iterator bam_reader::insert_fragment(chromosome* chrom,  rpos &left,  rpos &right) {
    
    greader_list<connected>::iterator merge_start, merge_end;
    bool found = false;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Test Frag " + std::to_string(left) + " " + std::to_string(right) + "\n");
    #endif

    for(greader_list<connected>::iterator it = chrom->chrom_fragments.begin();it != chrom->chrom_fragments.end(); ++it){
         #ifdef ALLOW_DEBUG
         logger::Instance()->debug("InFRAG " + std::to_string(it->start) + " " + std::to_string(it->end) + "\n");
         #endif
    }
    
    greader_list<connected>::iterator it = chrom->chrom_fragments.begin();
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Fraga " + std::to_string(it->start) + " " + std::to_string(it->end) + "\n");
    #endif
    for ( ; it != chrom->chrom_fragments.end() && right >= it->start; ++it) {
        // mark all overlap, connected areas 
        
        #ifdef ALLOW_DEBUG
         logger::Instance()->debug("Fragb " + std::to_string(it->start) + " " + std::to_string(it->end) + "\n");
        #endif

        if ( (it->end >= left && it->end <= right) || (it->start >= left && it->start <= right) || (it->start < left && it->end > right )) {
           
            if (!found) {
                found = true;
                merge_start = it;
            }
            merge_end = it;
        }
        
    }
    
    greader_list<connected>::iterator ret;
    if (!found) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Insert new Fragment " + std::to_string(left) + " - "+std::to_string(right)+".\n");
        #endif
        
        greader_list<connected>::iterator el = chrom->chrom_fragments.insert(it, connected());
         
        el->start = left;
        el->end = right;
        
        lazy< std::deque<read_collection> > new_inner = chrom->reads.add_inner();
        el->reads.ref().push_back(new_inner);
        
        ret = el;
              
    } else if (found && merge_start==merge_end) {
        // just modify this one, exons are changed later
        
        ret = merge_start;
        
        if (left < merge_start->start) {
            merge_start->start = left;
        }
        if (right > merge_start->end) {
            merge_start->end = right;
        }
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Extend Region with " + std::to_string(left) + " - "+std::to_string(right) + " to " + std::to_string(merge_start->start) + " - " + std::to_string(merge_start->end) +".\n");
        #endif
        
    } else {
        
        // modify first one to merge
        greader_list<connected>::iterator m = merge_start;
        
//         logger::Instance()->debug("======== " + std::to_string(m->atoms.ref().size()) + "\n");
//        logger::Instance()->debug("Connected " + std::to_string(m->start) + " " + std::to_string(m->end) + "\n");
//        for ( greader_list<exon* >::iterator fix_it = m->fossil_exons.ref().begin(); fix_it != m->fossil_exons.ref().end(); ++fix_it) {
//            logger::Instance()->debug("Fossil " + std::to_string((*fix_it)->start) + " " + std::to_string((*fix_it)->end) + "\n");
//        }
        
        ++m;
        greader_list<connected>::iterator end = merge_end;
        ++end;
        
        for (; m != end; ++m) {
//             logger::Instance()->debug("========" + std::to_string(m->atoms.ref().size()) + "\n");
//            logger::Instance()->debug("Connected " + std::to_string(m->start) + " " + std::to_string(m->end) + "\n");
//            for ( greader_list<exon* >::iterator fix_it = m->fossil_exons.ref().begin(); fix_it != m->fossil_exons.ref().end(); ++fix_it) {
//                logger::Instance()->debug("Fossil " + std::to_string((*fix_it)->start) + " " + std::to_string((*fix_it)->end) + "\n");
//            }
            
            _MOVE_RANGE(m->fossil_exons.ref().begin(), m->fossil_exons.ref().end(), std::inserter(merge_start->fossil_exons.ref(), merge_start->fossil_exons.ref().end()));
            _MOVE_RANGE(m->atoms.ref().begin(), m->atoms.ref().end(), std::inserter(merge_start->atoms.ref(), merge_start->atoms.ref().end()));
            
            _MOVE_RANGE(m->reads.ref().begin(), m->reads.ref().end(), std::back_inserter(merge_start->reads.ref()));
            
            merge_start->avg_split = ((merge_start->avg_split * merge_start->intel_count) + (m->avg_split * m->intel_count)) / (merge_start->intel_count + m->intel_count);
            merge_start->intel_count += m->intel_count;
            merge_start->bam_count = std::max(merge_start->bam_count, m->bam_count);
        }     
     
        if (left < merge_start->start) {
            merge_start->start = left;
        }
        if (right > merge_end->end) {
            merge_start->end = right;
        } else {
            merge_start->end = merge_end->end;
        }
        
       // mark_or_reduce_paired_atoms(&*merge_start, chrom, merge_start->atoms.ref().begin(), merge_start->atoms.ref().end()); done later anyway right now
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Join two Regions with " + std::to_string(left) + " - "+std::to_string(right) + " to " + std::to_string(merge_start->start) + " - " + std::to_string(merge_start->end) + " " + std::to_string(merge_start->atoms.ref().size()) +".\n");
        #endif

        ++merge_start;
        ret = chrom->chrom_fragments.erase(merge_start, end);
        --ret;
         
    }
    
    for(greader_list<connected>::iterator it = chrom->chrom_fragments.begin();it != chrom->chrom_fragments.end(); ++it){
         #ifdef ALLOW_DEBUG
         logger::Instance()->debug("OUTFRAG " + std::to_string(it->start) + " " + std::to_string(it->end) + "\n");
         #endif
    }
    
    return ret;
    
}


// ######### ATOMS #########

void bam_reader::assign_reads( connected* conn, chromosome* chrom) {
    
    // loop over exons (sorted)
    greader_list<interval >::iterator i_it_start = chrom->interval_queue.begin();
    
    for (greader_list<exon* >::iterator e_it = conn->fossil_exons.ref().begin(); e_it != conn->fossil_exons.ref().end(); ++e_it) {
        
        // loop over (still sorted) intervals
        for (greader_list<interval >::iterator i_it = i_it_start; i_it != chrom->interval_queue.end(); ++i_it) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("START " + std::to_string((i_it)->left)+ "-" + std::to_string((i_it)->right)+ ";" + std::to_string((*e_it)->start) + " - " + std::to_string((*e_it)->end) + "\n");
            
                    if ((i_it)->parent->id_set)  {
                        logger::Instance()->debug( "ID " + *(i_it)->parent->ids->begin() + "\n");
                    }
            #endif
            
            // test for starts
            if ((i_it)->right < (*e_it)->start && i_it == i_it_start) { // cannot overlap anymore, therefore increase start counter to one ahead
                ++i_it_start;
                #ifdef ALLOW_DEBUG
                 logger::Instance()->debug("Skip Start. \n");
                #endif
                continue;
            }
            
            if ((i_it)->left > (*e_it)->end) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Break. \n");
                 #endif
                break;
            }
            
//            if (    ( (i_it)->right >= (*e_it)->end  && (
//                                    (i_it)->left <= (*e_it)->start  // encased
//                                    ||  ( (i_it)->left <= (*e_it)->end && (i_it)->left >= (*e_it)->start 
//                                        && (*e_it)->end - (i_it)->left > options::Instance()->get_max_pos_extend()) ) 
//                    )  // overlap bigger merge
//                    || ( (i_it)->left < (*e_it)->start && (i_it)->right >= (*e_it)->start &&  (i_it)->right - (*e_it)->start> options::Instance()->get_max_pos_extend() )
//                    || ( (i_it)->right <= (*e_it)->end && (i_it)->left >= (*e_it)->start)
//                ) {
                
            if ((i_it)->right >= (*e_it)->start && (i_it)->left <= (*e_it)->end) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Overlap \n");
                #endif

                raw_atom* atom;
                if ( (i_it)->parent->atom == NULL) {

                    atom = (i_it)->parent->create_atom();
                            
//                    atom->reads.ref().push_back(read_collection((i_it)->parent->get_left_limit(), (i_it)->parent->right_limit, (i_it)->parent->length));
//                    if ((i_it)->parent->id_set) {
//                       atom->reads.ref().back().add_id((i_it)->parent->id);
//                    }
                } else {
                    atom = (i_it)->parent->atom;
                }
                
                // test for uncomplete overlaps due to unsupported moves
                if ( (i_it)->parent->get_left_limit() != (i_it)->left && (i_it)->left > (*e_it)->start + options::Instance()->get_min_raw_exon_size()
                        || (i_it)->parent->get_right_limit() != (i_it)->right && (i_it)->right + options::Instance()->get_min_raw_exon_size() < (*e_it)->end) {
                    (i_it)->parent->block = true;
                }
                atom->exons.ref().insert(*e_it);
            }
        }
    }
}


void bam_reader::filter_outer_read_junctions(chromosome* chrom) {

    lazy<r_border_set<rpos> > fs = chrom->fixed_exon_starts;
    lazy<r_border_set<rpos> > fe = chrom->fixed_exon_ends;
    
    for ( greader_list<rread>::iterator r_it = chrom->read_queue.begin(); r_it != chrom->read_queue.end(); ++r_it) {

        if (r_it->atom == NULL || r_it->atom->exons->size() < 2) {
            // this read was capped by a filter !
            continue;
        }
         
//      for (std::map< std::pair<rpos, rpos>, unsigned int >::iterator i =  chrom->splice_queue.begin(); i != chrom->splice_queue.end(); ++i) {
//            logger::Instance()->debug("Known Splice " + std::to_string(i->first.first ) + ":" + std::to_string(i->first.second) + " " + std::to_string(i->second) + "\n" );
//      }
        
        greader_refsorted_list<exon*>::iterator ei = r_it->atom->exons->begin();
        greader_refsorted_list<exon*>::iterator ein = ei;
        ++ein; 
        for ( ; ein != r_it->atom->exons->end() && !r_it->block; ++ei, ++ein ) {
            
//            logger::Instance()->debug("Test SPlice " + std::to_string((*ei)->end+1) + ":" + std::to_string((*ein)->start-1) + "\n" );

            if ((*ei)->end+1 == (*ein)->start) {
                continue;
            }
            
            std::map< std::pair<rpos, rpos>, unsigned int >::iterator splice_count = chrom->splice_queue.find(std::make_pair( (*ei)->end+1, (*ein)->start-1) );
            if ( splice_count == chrom->splice_queue.end() && splice_count->second < options::Instance()->get_min_junction_coverage()) {
//                logger::Instance()->debug("block\n" );
                r_it->block = true;
            }
        }
        if (r_it->block) {
            continue;
        }

        // test first and last atom wheter it is long enough atom
        rpos begin_start = (*r_it->atom->exons->begin())->start;
        rpos begin_end = (*r_it->atom->exons->begin())->end;

        rpos end_start = (*r_it->atom->exons->rbegin())->start;
        rpos end_end = (*r_it->atom->exons->rbegin())->end;

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("TEST START " + std::to_string(begin_start) + ":" + std::to_string(begin_end) + " to " + std::to_string(r_it->left_limit) + "\n" );
        logger::Instance()->debug("TEST END " + std::to_string(end_start) + ":" + std::to_string(end_end) + " to " + std::to_string(r_it->right_limit) + "\n" );
        #endif
        
        if (r_it->left_limit < begin_start) { // in case other one was missed
            r_it->left_limit = begin_start;
        }
        bool cut_first;
        if (begin_end - begin_start + 1 > options::Instance()->get_min_raw_exon_size()) { // exon bigger than cutoff
            cut_first = begin_end - r_it->left_limit + 1 < options::Instance()->get_min_raw_exon_size();
            if (!cut_first) {
                if(!fe->sorted_find(begin_end) && !fs->sorted_find(begin_end+1)) { // unsupported
                    cut_first = true;
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Not in Ends\n");
                    #endif
                }
            } 
        } else {
            cut_first = begin_start != r_it->left_limit;
        } 
        

        if (r_it->right_limit > end_end) { // in case other one was missed
            r_it->right_limit = end_end;
        }
        bool cut_last;
        if (end_end - end_start + 1 > options::Instance()->get_min_raw_exon_size()) {
            cut_last = r_it->right_limit - end_start + 1 < options::Instance()->get_min_raw_exon_size();
            if (!cut_last) {
                if(!fs->sorted_find(end_start) && !fe->sorted_find(end_start-1)) { // unsupported
                    cut_last = true;
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Not in Starts\n");
                    #endif
                }
            }
        } else {
            cut_last = end_end != r_it->right_limit;
        } 
        
        if (cut_first && cut_last && r_it->atom->exons->size() == 2) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("ERASE READ\n");
            #endif
            r_it->block = true;
            continue;
        }

        if (cut_first) {
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Cut First\n");
            #endif
            r_it->atom->exons->erase(r_it->atom->exons->begin());   // erase and set to new begin
            r_it->left_limit = (*r_it->atom->exons->begin())->start; 
        }

        if (cut_last) {
            #ifdef ALLOW_DEBUG
             logger::Instance()->debug("Cut Last\n");
            #endif
            r_it->atom->exons->erase(std::prev(r_it->atom->exons->end()));   // erase and set to new end
            r_it->right_limit = (*r_it->atom->exons->rbegin())->end; 
        }

    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("DONE\n");
    #endif
}


void bam_reader::reduce_atoms(connected* conn, chromosome* chrom) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Reduce connected .\n");
    #endif
     
    // these need to be inserted
    for ( greader_list<rread>::iterator r_it = chrom->read_queue.begin(); r_it != chrom->read_queue.end(); ++r_it) {
        
        if (r_it->atom == NULL || r_it->block) {
            // this read was capped by a filter !
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Filtered Read " + std::to_string(r_it->left_limit) + " " + std::to_string(r_it->right_limit) + " - " + std::to_string(r_it->count) + " Blocked " + std::to_string(r_it->block) + ".\n");
            if (r_it->atom != NULL) {
                    logger::Instance()->debug("Blocked Atom " + r_it->atom->to_string() + ".\n");
                 
            }
            #endif
            continue;
        }
        
        #ifdef ALLOW_DEBUG
         logger::Instance()->debug("Search Existing Atom " + r_it->atom->to_string() + ".\n");
        #endif
        
        // into the existing atoms marked by the 
        raw_atom* atom;
        greader_refsorted_list<raw_atom*>::iterator atom_it = conn->atoms.ref().find( r_it->atom );
        if (atom_it == conn->atoms.ref().end()) {
            
            // atom does not exist yet, so just add
            chrom->atoms.push_back(*r_it->atom);
            atom = &chrom->atoms.back();
            conn->atoms.ref().insert(atom);
            // we can take the atom as is, but read has still not been added
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("New Atom " + (atom)->to_string() + ".\n");
            #endif
            
        } else {
           
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Existing Atom " + (*atom_it)->to_string() + ".\n");
            #endif
            
          // we found the correct one, so merge into it
          atom =  *atom_it; 
        }
          
        // try and find read collection // min max used for filtered exon boundaries
        read_collection* rc_joined = new read_collection(std::max(r_it->left_limit, (*atom->exons->begin())->start), std::min(r_it->right_limit, (*atom->exons->rbegin())->end), r_it->length, atom);
        greader_refsorted_list<read_collection*>::iterator rc_it = atom->reads.ref().find( rc_joined );
        read_collection* rc;
        if (rc_it == atom->reads.ref().end()) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("New Collection.\n");
            #endif
            
            conn->reads.ref().push_to_end(*rc_joined);
            rc = &conn->reads.ref().get_end();
            atom->reads->insert(rc);
        } else {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Existing Collection .\n");
            #endif
            
            rc = *rc_it;
        }
      
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Collection: " + std::to_string(rc->left_limit) + " " + std::to_string(rc->right_limit) + ".\n");
        #endif

        // now add info of current read!
        if (r_it->id_set) {
           _MOVE_RANGE( r_it->ids.ref().begin(), r_it->ids.ref().end(), std::back_inserter(rc->open_ids.ref()));
        }
        rc->count += r_it->count;
        rc->length_filtered = rc->length_filtered && rc_joined->length_filtered;
        
        delete rc_joined;
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Final Graph ++++++++++++++++++++++++++++++++\n");
    logger::Instance()->debug("Reduced to " + std::to_string(conn->atoms.ref().size())+" atoms.\n"); 
    #endif
}

void bam_reader::mark_or_reduce_paired_atoms( connected* conn, chromosome* chrom , const greader_refsorted_list<raw_atom*>::iterator &atom_start, const greader_refsorted_list<raw_atom*>::iterator &atom_end) {
    
    boost::unordered_map<std::string, std::tuple<raw_atom*, read_collection*, std::string* > > id_map;
    
    // set current indices, makes joining faster!
    unsigned int i=0;
    for (greader_list<exon* >::iterator e_it = conn->fossil_exons.ref().begin(); e_it != conn->fossil_exons.ref().end(); ++e_it,++i) {
        (*e_it)->id = i;
    }
    
    for (greader_refsorted_list<raw_atom*>::iterator a_it = atom_start; a_it != atom_end; ++a_it) {
        
        #ifdef ALLOW_DEBUG
          logger::Instance()->debug("RAW: " + std::to_string((*(*a_it)->exons->begin())->id) + " " + std::to_string((*(*a_it)->exons->rbegin())->id) + "\n");
        #endif
        
        for (greader_refsorted_list<read_collection*>::iterator c_it = (*a_it)->reads.ref().begin(); c_it != (*a_it)->reads.ref().end(); ++c_it) {
             for (greader_list<std::string>::iterator i_it = (*c_it)->open_ids.ref().begin(); i_it != (*c_it)->open_ids.ref().end(); ++i_it) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("ID: " + *i_it + "\n");
                #endif
                
                boost::unordered_map<std::string, std::tuple<raw_atom*, read_collection*, std::string* > > ::iterator find = id_map.find(*i_it);
            
                 if ( find == id_map.end()) {
                    id_map[*i_it] = std::make_tuple(&*(*a_it), &*(*c_it), &*i_it);
                } else {
                    // we found a pair
                    raw_atom* first = std::get<0>(find->second);
                    read_collection* fcol = std::get<1>(find->second);
                    
                    // test if we want to join these two
                    exon* lastexon = *first->exons.ref().rbegin();
                    exon* firstexon = *(*a_it)->exons.ref().begin();
                    
                    exon* ll = *first->exons.ref().begin();
                    exon* rr = *(*a_it)->exons.ref().rbegin();
                    
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Test exon join " + std::to_string(ll->id) + "-" + std::to_string(lastexon->id) + ":" + std::to_string(firstexon->id) + "-" + std::to_string(rr->id) + " " + std::to_string(first->exons.ref().size()) + "-" + std::to_string((*a_it)->exons.ref().size()) + "\n" );
                    logger::Instance()->debug("L " + first->to_string() + "\n" );
                    logger::Instance()->debug("R " + (*a_it)->to_string() + "\n" );
                    #endif
                  
                    if ( ( (first->exons.ref().size() ==1 || (*a_it)->exons.ref().size() == 1) && lastexon->id == firstexon->id )
                            || ll->id == firstexon->id || rr->id == lastexon->id ) {
                        // this happens when left partner is a subset, i.e. only one exon before a split
                        // do nothing but remove after if

                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Subset \n" );
                        #endif
                        
                        
                    } else if (lastexon->id == firstexon->id || lastexon->id +1 == firstexon->id ) {
                        // join :)

                        if (lastexon->id +1 == firstexon->id) {
                            // there is a gap!, we can only do this if we have a proven split here! one of my major pains
                            
                            bool junction_found = false;;
                            for (greader_refsorted_list<raw_atom*>::iterator ra = atom_start; ra != atom_end; ++ra) {
                                
                                if ( (*(*ra)->exons.ref().begin())->id > lastexon->id ) {
                                    break; // we are sorted and this one is bigger!
                                }
                                
                                // we need to find one with firstexon and lastexon in same raw!
                                if ( (*ra)->exons.ref().find(firstexon) != (*ra)->exons.ref().end() && (*ra)->exons.ref().find(lastexon) != (*ra)->exons.ref().end() ) {
                                    // we found it
                                    junction_found = true;
                                    break;
                                }
                            }
                            if (!junction_found) {
                                continue;
                            }
                        }
                        
                        
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("1Join \n" );
                        #endif
                        
                        raw_atom* joined_atom = new raw_atom();

                        std::copy(first->exons.ref().begin(), first->exons.ref().end(),std::inserter( joined_atom->exons.ref(), joined_atom->exons.ref().end()) );
                        std::copy((*a_it)->exons.ref().begin(), (*a_it)->exons.ref().end(),std::inserter( joined_atom->exons.ref(), joined_atom->exons.ref().end()) );

                        // try to find atom in connected
                        raw_atom* atom;
                        greader_refsorted_list<raw_atom*>::iterator atom_it = conn->atoms.ref().find( joined_atom );
                        if (atom_it == conn->atoms.ref().end()) {

                            // atom does not exist yet, so just add
                            chrom->atoms.push_back(*joined_atom);
                            atom = &chrom->atoms.back();
                            conn->atoms.ref().insert(atom);
                            // we can take the atom as is
 
                        } else {
                          // we found the correct one, so merge into it
                          atom =  *atom_it; 
                        }
                        
                        // try and find read collection
                        read_collection* rc_joined = new read_collection(fcol->left_limit, (*c_it)->right_limit, false, atom);
                        greader_refsorted_list<read_collection*>::iterator rc_it = atom->reads.ref().find( rc_joined );
                        read_collection* rc;
                        if (rc_it == atom->reads.ref().end()) {
                            conn->reads.ref().push_to_end(*rc_joined);
                            rc = &conn->reads.ref().get_end();
                            atom->reads->insert(rc);
                        } else {
                            rc = *rc_it;
                        }
                        
                        rc->count+=2;
                        ++rc->paired_count;
                        
                        --fcol->count;
                        --(*c_it)->count;
                           
                        rc->holes->push_back(std::make_pair(fcol->right_limit, (*c_it)->left_limit));
                        
                        delete joined_atom;
                        delete rc_joined;
                        
//                        ++fcol->paired[&*(*c_it)];
                        
                    } else if (lastexon->id > firstexon->id && rr->id > lastexon->id) {
                        
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Overlap \n" );
                        #endif
                        
                        bool matching = true;
                        
                        greader_refsorted_list<exon*>::iterator new_right_it = (*a_it)->exons.ref().begin();
                        greader_refsorted_list<exon*>::iterator old_left_it =  first->exons->find(*new_right_it);
                        
                        if (old_left_it == first->exons->end()) { // if first one cannot be found, we already lost!
                            matching = false;
                        }
                        // now loop over all remaining in tandem!
                        for (; old_left_it != first->exons->end(); ++new_right_it, ++old_left_it) {
                            
                            #ifdef ALLOW_DEBUG
                            logger::Instance()->debug("Test " + std::to_string((*new_right_it)->id) +"\n" );
                            #endif
                            
                            if ( (*new_right_it)->id !=  (*old_left_it)->id ) {
                                #ifdef ALLOW_DEBUG
                                logger::Instance()->debug("Unfound\n" );
                                #endif
                                matching = false;
                                break;
                            }
                        }
                                               
                        if (matching) {
                            
                            raw_atom* joined_atom = new raw_atom();
                            raw_atom* overlap_atom = new raw_atom();

                            std::copy(first->exons.ref().begin(), first->exons.ref().end(),std::inserter( joined_atom->exons.ref(), joined_atom->exons.ref().end()) );

                            greader_refsorted_list<exon*>::iterator r_it = (*a_it)->exons.ref().begin();
                            for (; r_it !=(*a_it)->exons.ref().end() && (*r_it)->id <= lastexon->id; ++r_it )  {
                                overlap_atom->exons->insert(*r_it);
                            }
                            for (; r_it !=(*a_it)->exons.ref().end(); ++r_it )  {
                                joined_atom->exons->insert(*r_it);
                            }

                            ////// JOINED

                            // try to find JOINED atom in connected
                            raw_atom* atom_joined;
                            greader_refsorted_list<raw_atom*>::iterator atom_it = conn->atoms.ref().find( joined_atom );
                            if (atom_it == conn->atoms.ref().end()) {

                                // atom does not exist yet, so just add
                                chrom->atoms.push_back(*joined_atom);
                                atom_joined = &chrom->atoms.back();
                                conn->atoms.ref().insert(atom_joined);
                                // we can take the atom as is

                            } else {
                                // we found the correct one, so merge into it
                                atom_joined =  *atom_it; 
                            }

                            // try and find JOINED read collection
                            read_collection* rc_joined = new read_collection(fcol->left_limit, (*c_it)->right_limit, false, atom_joined);
                            greader_refsorted_list<read_collection*>::iterator rc_it = atom_joined->reads.ref().find( rc_joined );
                            read_collection* rcj;
                            if (rc_it == atom_joined->reads.ref().end()) {
                                conn->reads.ref().push_to_end(*rc_joined);
                                rcj = &conn->reads.ref().get_end();
                                atom_joined->reads->insert(rcj);
                            } else {
                                rcj = *rc_it;
                            }

                            ////// OVERLAP
                            if (options::Instance()->is_create_overlap_merge_read()) {
                                // try to find OVERLAP atom in connected
                                raw_atom* atom_overlap;
                                atom_it = conn->atoms.ref().find( overlap_atom );
                                if (atom_it == conn->atoms.ref().end()) {

                                    // atom does not exist yet, so just add
                                    chrom->atoms.push_back(*overlap_atom);
                                    atom_overlap = &chrom->atoms.back();
                                    conn->atoms.ref().insert(atom_overlap);
                                    // we can take the atom as is

                                } else {
                                  // we found the correct one, so merge into it
                                  atom_overlap =  *atom_it; 
                                }

                                // try and find OVERLAP read collection
                                read_collection* rc_overlap = new read_collection((*c_it)->left_limit, fcol->right_limit, false, atom_overlap);
                                rc_it = atom_overlap->reads.ref().find( rc_overlap );
                                read_collection* rco;
                                if (rc_it == atom_overlap->reads.ref().end()) {
                                    conn->reads.ref().push_to_end(*rc_overlap);
                                    rco = &conn->reads.ref().get_end();
                                    atom_overlap->reads->insert(rco);
                                } else {
                                    rco = *rc_it;
                                }

                                ++rco->count;


                                delete rc_overlap;
                            }

                            rcj->count+=2;
                            ++rcj->paired_count;

                            --fcol->count;
                            --(*c_it)->count;

                            delete overlap_atom;
                            delete joined_atom;
                            delete rc_joined;
                        }
                        
//                        ++fcol->paired[&*(*c_it)];
                    } else {
                        
                        ++fcol->paired[&*(*c_it)];
                        
                    }
                    
                    (*c_it)->flag_id(*i_it);
                     fcol->flag_id(*std::get<2>(find->second));
                    
                    id_map.erase(find);
                    
                    // do insert size statistic // heuristic estimating as if all exons in between are part of the atom
                    if (fcol->right_limit <= (*c_it)->left_limit) {
                        ++conn->intel_count;
                        
                        rpos len = 0;
                        if (lastexon->id == firstexon->id) {
                            len = (*c_it)->left_limit - fcol->right_limit;
                        } else { 
                            len += lastexon->end - fcol->right_limit + 1;
                            len += (*c_it)->left_limit - firstexon->start + 1;
                            for (unsigned int li = lastexon->id + 1; li < firstexon->id; ++li) {
                                len += conn->fossil_exons->at(li)->end - conn->fossil_exons->at(li)->start + 1; ;
                            }
                        }

                        conn->avg_split += ( len - conn->avg_split) / conn->intel_count;                        
                    }
                } 
            }
        }
    }
    for (greader_refsorted_list<raw_atom*>::iterator a_it = atom_start; a_it != atom_end; ++a_it) {  
        for (greader_refsorted_list<read_collection*>::iterator c_it = (*a_it)->reads.ref().begin(); c_it != (*a_it)->reads.ref().end(); ++c_it) {
             (*c_it)->clean_flagged_ids();
        }
    }
}


void bam_reader::filter_bins(connected* conn, chromosome* chrom) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("------------ Filter Bins\n");
    #endif

    for (greader_refsorted_list<raw_atom*>::iterator a_it = conn->atoms->begin(); a_it != conn->atoms->end(); ++a_it) {
                
        // test if this atom reads overlapping far enough to left and right
        
        if ((*a_it)->exons->size() < 2) {
            continue;
        }
        
        if ((*a_it)->reads.ref().empty()) {
            continue;
        }
        
//        if ((*a_it)->reads->size() == 1 && (*a_it)->exons->size() > 2 ) {
//            (*a_it)->reads.ref().clear();
//            (*a_it)->count = 0;
//            (*a_it)->paired_count = 0;
//            continue;
//        }
        
        rpos begin_end = (*(*a_it)->exons->begin())->end;

        rpos end_start = (*(*a_it)->exons->rbegin())->start;
                
        bool cut_start = true;
        bool cut_end = true;
        
        int i = 1;
        
        for (greader_refsorted_list<read_collection*>::iterator c_it = (*a_it)->reads.ref().begin(); c_it != (*a_it)->reads.ref().end() && ( cut_start || cut_end) ; ++c_it, ++i) {
            
            rpos left = (*c_it)->left_limit;
            rpos right = (*c_it)->right_limit;
            
//            logger::Instance()->info("RC " + std::to_string(left) + "-" + std::to_string(right) + " " + std::to_string(begin_end) + "-" + std::to_string(end_start) + "\n");
//            logger::Instance()->info(" " + std::to_string(i) + "/" + std::to_string((*a_it)->reads->size()) + "\n");
//            
            if (begin_end - left + 1 >= options::Instance()->get_min_junction_anchor()) {
                cut_start = false;
            } 
            if (right - end_start + 1  >= options::Instance()->get_min_junction_anchor()) {
                cut_end = false;
            }          
        }
        
//        logger::Instance()->info("----- Atom " + (*a_it)->to_string() + " " + std::to_string(cut_start) + std::to_string(cut_end) + "\n");

        raw_atom* cut_atom = new raw_atom();
        if (cut_start && cut_end) {
            
            if ((*a_it)->exons->size() == 2) {
                
                (*a_it)->reads.ref().clear();
                (*a_it)->count = 0;
                (*a_it)->paired_count = 0;
                        
                delete cut_atom;
                continue;
            }
            
            std::copy(std::next((*a_it)->exons->begin()), std::prev((*a_it)->exons->end()),std::inserter( cut_atom->exons.ref(), cut_atom->exons.ref().end()) );
        } else if (cut_start) {
            std::copy(std::next((*a_it)->exons->begin()), (*a_it)->exons->end(),std::inserter( cut_atom->exons.ref(), cut_atom->exons.ref().end()) );
        } else if (cut_end) {
            std::copy((*a_it)->exons->begin(), std::prev((*a_it)->exons->end()),std::inserter( cut_atom->exons.ref(), cut_atom->exons.ref().end()) );
        } else {
            delete cut_atom;
            continue;
        }
        
//        logger::Instance()->info("New " + cut_atom->to_string() +  "\n");
                
        raw_atom* atom;
        greader_refsorted_list<raw_atom*>::iterator atom_it = conn->atoms.ref().find( cut_atom );
        if (atom_it == conn->atoms.ref().end()) {

            // atom does not exist yet, so just add
            chrom->atoms.push_back(*cut_atom);
            atom = &chrom->atoms.back();
            conn->atoms.ref().insert(atom);
            // we can take the atom as is

        } else {
          // we found the correct one, so merge into it
          atom =  *atom_it; 
        }
        
        for (greader_refsorted_list<read_collection*>::iterator c_it = (*a_it)->reads.ref().begin(); c_it != (*a_it)->reads.ref().end(); ++c_it) {
            // move all read collections over to shorter form
            
            read_collection* rc_joined = new read_collection(std::max((*c_it)->left_limit, (*atom->exons->begin())->start), std::min((*c_it)->right_limit, (*atom->exons->rbegin())->end), false, atom);
            greader_refsorted_list<read_collection*>::iterator rc_it = atom->reads.ref().find( rc_joined );
            read_collection* rc;
            if (rc_it == atom->reads.ref().end()) {
                conn->reads.ref().push_to_end(*rc_joined);
                rc = &conn->reads.ref().get_end();
                atom->reads->insert(rc);
            } else {
                rc = *rc_it;
            }

            // now add info of current read!
            rc->count += (*c_it)->count;
            rc->paired_count += (*c_it)->paired_count;
            rc->length_filtered = rc->length_filtered && rc_joined->length_filtered;

            for (paired_map<read_collection*, rcount >::iterator pi = (*c_it)->paired.begin(); pi != (*c_it)->paired.end(); ++pi) {
                rc->paired[pi->first] += pi->second;
            }
            for (std::deque<std::pair<rpos, rpos> >::iterator hi = (*c_it)->holes->begin(); hi != (*c_it)->holes->end(); ++hi) {
                rc->holes->push_back(std::make_pair( std::max(hi->first, (*atom->exons->begin())->start), std::min(hi->second, (*atom->exons->rbegin())->end) ));
            }
            
        }
        
        (*a_it)->reads.ref().clear();
        (*a_it)->count = 0;
        (*a_it)->paired_count = 0;
        
        delete cut_atom;
        
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("------------ Filter Bins Out\n");
    #endif
    
}


void bam_reader::reduce_reads(connected* conn) {
    
    for (greader_refsorted_list<raw_atom*>::iterator a_it = conn->atoms->begin(); a_it != conn->atoms->end(); ++a_it) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("------------new ATOM\n");
        #endif
          
        for (greader_refsorted_list<read_collection*>::iterator c_it = (*a_it)->reads.ref().begin(); c_it != (*a_it)->reads.ref().end(); ++c_it) {
          //  logger::Instance()->debug("Add collection \n");
            // basic info to break down
            (*a_it)->paired_count += (*c_it)->paired_count;
            (*a_it)->count += (*c_it)->count;
            (*a_it)->length_filtered = (*a_it)->length_filtered && (*c_it)->length_filtered;
                    
            rcount fragcount = (*c_it)->count - (*c_it)->paired_count;
            
            std::map< rpos,rcount >::iterator fl = (*a_it)->lefts->find((*c_it)->left_limit);
            if (fl == (*a_it)->lefts->end()) {
                (*a_it)->lefts->insert(std::make_pair((*c_it)->left_limit, fragcount));
            } else {
                fl->second += fragcount;
            }
            std::map< rpos,rcount >::iterator fr = (*a_it)->rights->find((*c_it)->right_limit);
            if (fr == (*a_it)->rights->end()) {
                (*a_it)->rights->insert(std::make_pair((*c_it)->right_limit, fragcount));
            } else {
                fr->second += fragcount;
            }

            // transfer down coverage info
            if ((*a_it)->exons.ref().size() == 1) {
                // if we just have one exon, start to end are bases, we also just use the start
//                (*a_it)->base_count_start += ((*c_it)->right_limit - (*c_it)->left_limit + 1) * fragcount;  
            } else {
                // more than two exons, so there is a fixed junction guaranteed in the middle
                // left
//                rpos first_right = (*(*a_it)->exons->begin())->end;
//                rpos last_left = (*(*a_it)->exons->rbegin())->start;
//                
//                (*a_it)->base_count_start += (first_right - (*c_it)->left_limit + 1) * fragcount;
//                (*a_it)->base_count_end += ((*c_it)->right_limit - last_left + 1 ) * fragcount ;
                                
                (*a_it)->total_rights += fragcount;
                (*a_it)->total_lefts += fragcount;                    
            }
            
            for (std::deque<std::pair<rpos, rpos> >::iterator hole_it = (*c_it)->holes->begin(); hole_it != (*c_it)->holes->end();++hole_it) {
                // unfortunately we need to find the actual exons for this
                    
                    std::map< rpos,rcount >::iterator fl = (*a_it)->hole_starts->find(hole_it->first);
                    if (fl == (*a_it)->hole_starts->end()) {
                        (*a_it)->hole_starts->insert(std::make_pair(hole_it->first, 1));
                    } else {
                        fl->second += 1;
                    }
                    std::map< rpos,rcount >::iterator fr = (*a_it)->hole_ends->find(hole_it->second);
                    if (fr == (*a_it)->hole_ends->end()) {
                        (*a_it)->hole_ends->insert(std::make_pair(hole_it->second, 1));
                    } else {
                        fr->second += 1;
                    }  
            }
                    
            // now transfer actual paired
            for(paired_map<read_collection*, rcount >::iterator pair_it = (*c_it)->paired.begin(); pair_it != (*c_it)->paired.end(); ++pair_it) {
              //  logger::Instance()->debug("Add paired \n");
                (*a_it)->paired[pair_it->first->parent] +=  pair_it->second;
            }
        }
        (*a_it)->reads.ref().clear();
    }
    
    for (double_deque_ref<read_collection>::iterator rc_it = conn->reads->begin(); rc_it != conn->reads->end(); ++rc_it) {
        rc_it->ref().clear();
    }
    conn->reads->begin()->ref().clear();

}

void bam_reader::reset_reads(chromosome* chrom) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Reset reads \n");
    #endif
    
    // when this is called ALL read collections should be removed down
    chrom->reads.clear();
    for( greader_list<connected>::iterator con =  chrom->chrom_fragments.begin(); con!=  chrom->chrom_fragments.end(); ++con) {
         lazy< std::deque<read_collection> > new_inner = chrom->reads.add_inner();
         con->reads.ref().push_back(new_inner);
    }
}
