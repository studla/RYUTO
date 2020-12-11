/* 
 * File:   transcript.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 11, 2016, 1:30 PM
 */

#include "transcript.h"

transcript::transcript() : cycle_id_in(0), cycle_id_out(0), guided(false), guide_grouped(false), ignore_guide_gene(false), post_filter_regional_group(0) {
}


transcript::~transcript() {
}

void transcript::join(transcript* o) {
    
    found_edge.join_edge(o->found_edge);

    for (gmap<int, series_struct>::iterator iss = series.begin(); iss != series.end(); ++iss) {
        int id = iss->first;
        iss->second.flow = (iss->second.flow + o->series[id].flow) / 2;
        iss->second.mean = (iss->second.mean + o->series[id].mean) / 2;
        iss->second.score = (iss->second.score + o->series[id].score) / 2;
    }
    length += o->length;

    for(std::deque<transcript_unsecurity>::iterator ui = o->unsecurity_id.begin(); ui != o->unsecurity_id.end(); ++ui) {
        unsecurity_id.push_back(*ui);
    }
    
    std::deque<std::pair<rpos, rpos> > exons;
        
    std::deque<std::pair<rpos, rpos> >::iterator ei = o->exons.begin();
    exons.back().second = ei->second;
    ++ei;
    for (; ei !=  o->exons.end(); ++ei) {
        exons.push_back(*ei);
    }
    
    // TODO
//    unsigned int cycle_id_in;
//    unsigned int cycle_id_out;
    
    // finalized values
;
//    float fpkm;
//    std::string chromosome;
//    std::string strand;
}

void transcript::print(std::ostream &os) {
    
//    os << std::to_string(cycle_id_in) << "\t" << std::to_string(cycle_id_out) << "\t";
//    os << found_edge.to_string() << "\t" << std::to_string(flow) << "\t";
//    
//    std::deque<transcript_unsecurity>::iterator it = unsecurity_id.begin();
//    
//    if (it == unsecurity_id.end()) {
//    
//        os << "-";
//       
//    } else {
//        
//        os << std::to_string(it->position);
//        ++it;
//        for(; it!= unsecurity_id.end(); ++it) {
//            os << ", " << std::to_string(it->position) << "(" << it->evidenced << ")" ;
//        }
//    }
//    
//     os << "\n";
    
}

void transcript::finalize_borders(exon_meta* meta) {
    
    boost::dynamic_bitset<>::size_type index = found_edge.id.find_first();
    unsigned int pos = index;
    length = 0;
    
    index = found_edge.id.find_next(index);
    rpos region_start =  meta->exons[pos].left;
    rpos region_end =  meta->exons[pos].right;
    
    while(index <= found_edge.id.size()) {
        
        unsigned int pos = index;
        if (region_end + 1 == meta->exons[pos].left) { // do we join these exons, because thy are directly adjourning
            
            region_end = meta->exons[pos].right;
            
        } else {
            exons.push_back(std::make_pair(region_start, region_end));
            length += region_end - region_start + 1;
            
            region_start =  meta->exons[pos].left;
            region_end = meta->exons[pos].right; 
        }
       
        index = found_edge.id.find_next(index);
    } 
    exons.push_back(std::make_pair(region_start, region_end));
    length += region_end - region_start + 1;
    
    for (gmap<int, series_struct>::iterator iss = series.begin(); iss != series.end(); ++iss) {
        iss->second.fpkm = iss->second.mean / (float)(2 * meta->avrg_read_length) / (float) (meta->absolute_reads[iss->first]) * 1000000000;
    }
    chromosome = meta->chromosome;
    strand = meta->strand;

    avrg_read_length = meta->avrg_read_length;

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Finalize: "+found_edge.to_string() + "\n");
    logger::Instance()->debug("Meta: "+std::to_string(avrg_read_length) + "\n");
    #endif
}

void transcript::print_count_matrix_entry(std::ostream &os, std::string &gene_id, unsigned int trans_id, int main_input_id, std::set<int> &ids) {

    if (guided) {
        if (ignore_guide_gene) {
            os << "ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group) << "\t"; 
        } else {
            os << guide_gene << "\t";
        }
        os << guide_reference << "\t";

    } else {
        if (guide_grouped) {
              os << guide_gene << "\t";
        } else {
             os << "ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group)  << "\t";
        }
        os  << "ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group) << "." << std::to_string(trans_id) << "\t" ;
    }

    os << std::to_string(length);
        
    for (std::set<int>::iterator iid = ids.begin(); iid != ids.end(); ++iid) {
        
        if (ids.size() > 1 && *iid == main_input_id) {
            continue;
        }
        
        gmap<int, series_struct>::iterator iss = series.find(*iid);
        if (iss == series.end()) {
            os << "\t" << std::to_string(0);
        } else {
            //avrg_read_length = 75;
            long count = (long) (series[iss->first].mean * ((series[iss->first].effective_length)  / (float)(2 * avrg_read_length)));
            if (count == 0 && series[iss->first].mean > 0) {
                 count = 1;
            }
            os << "\t" << std::to_string( count );// << "(" << std::to_string(series[iss->first].effective_length) << "," << std::to_string(series[iss->first].mean) << "," << std::to_string(avrg_read_length) << ")" ;
        }
        
    }
    os << "\n";
}


void transcript::print_gtf_entry(std::ostream &os, std::string &gene_id, unsigned int trans_id, int main_input_id) {
        
    // print transcript line first
    os << chromosome << "\t" << "ryuto" << "\t" << "transcript" << "\t";
    os << std::to_string(exons.begin()->first) << "\t" << std::to_string(exons.rbegin()->second) << "\t";
    os << "0" << "\t"; // TODO: create a SCORE?
    os << strand << "\t" << "." << "\t";
    os << "gene_id \"ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group) << "\"; ";
    os << "transcript_id \"ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group) << "." << std::to_string(trans_id) << "\"; ";
    os << "FPKM \"" << std::to_string(series[main_input_id].fpkm) << "\"; ";
    os << "cov \"" << std::to_string(series[main_input_id].mean) << "\"; ";
    os << "uniform_cov \"" << std::to_string(series[main_input_id].flow) << "\";";
    
    for (gmap<int, series_struct>::iterator iss = series.begin(); iss != series.end(); ++iss) {
        int id = iss->first;
        
        if (id == main_input_id) {
            continue;
        }
        
        os << "FPKM_" << std::to_string(id) << " \"" << std::to_string(series[id].fpkm) << "\"; ";
        os << "cov_" << std::to_string(id) << " \"" << std::to_string(series[id].mean) << "\"; ";
        os << "uniform_cov_" << std::to_string(id) << " \"" << std::to_string(series[id].flow) << "\";";
    }

    os << "unsecurities \"";
    std::deque<transcript_unsecurity>::iterator u_it = unsecurity_id.begin();
    if (u_it == unsecurity_id.end()) {
        os << "-";
    } else {
        os << *u_it;
        ++u_it;
        for(; u_it!= unsecurity_id.end(); ++u_it) {
            os << ", " << *u_it ;
        }
    }
    os << "\";\n";
    
    // then second print individual exons
    for(  std::deque<std::pair<rpos, rpos> >::iterator it = exons.begin(); it != exons.end(); ++it) {
        
        os << chromosome << "\t" << "ryuto" << "\t" << "exon" << "\t";
        os << std::to_string(it->first) << "\t" << std::to_string(it->second) << "\t";
        os << "0" << "\t"; // TODO: create a SCORE?
        os << strand << "\t" << "." << "\t";
        
        os << "gene_id \"ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group) << "\"; ";
        os << "transcript_id \"ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group) << "." << std::to_string(trans_id) << "\"; ";
//        os << "FPKM \"" << std::to_string(fpkm) << "\"; ";

        if (cycle_id_in != 0) {
           os << "cycle_in \"ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group) << "." << std::to_string(trans_id) << "/" << std::to_string(cycle_id_in) << "\"; ";
        }
        if (cycle_id_out != 0) {
           os << "cycle_out \"ryuto." << gene_id << "_" << std::to_string(post_filter_regional_group) << "." << std::to_string(trans_id) << "/" << std::to_string(cycle_id_out) << "\"; ";
        }
//        os << "cov \"" << std::to_string(mean) << "\"; ";
//        os << "uniform_cov \"" << std::to_string(flow) << "\"; ";
        
        if (guided) {
            os << "reference_id \"" << guide_reference << "\"; ";
        }
        
        os << "\n";
    
    }
    
}
