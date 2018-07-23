/* 
 * File:   transcript.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 11, 2016, 1:30 PM
 */

#include "transcript.h"

transcript::transcript() : flow(0), mean(0), cycle_id_in(0), cycle_id_out(0), guided(false) {
}


transcript::~transcript() {
}

void transcript::print(std::ostream &os) {
    
    os << std::to_string(cycle_id_in) << "\t" << std::to_string(cycle_id_out) << "\t";
    os << found_edge.to_string() << "\t" << std::to_string(flow) << "\t";
    
    std::deque<transcript_unsecurity>::iterator it = unsecurity_id.begin();
    
    if (it == unsecurity_id.end()) {
    
        os << "-";
       
    } else {
        
        os << std::to_string(it->position);
        ++it;
        for(; it!= unsecurity_id.end(); ++it) {
            os << ", " << std::to_string(it->position) << "(" << it->evidenced << ")" ;
        }
    }
    
     os << "\n";
    
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
    
    fpkm = mean / (float)(2 * meta->avrg_read_length) / (float) (meta->absolute_reads) * 1000000000;
    
    chromosome = meta->chromosome;
    strand = meta->strand;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Finalize: "+found_edge.to_string()+ " l " + std::to_string(length) + " f " + std::to_string(flow) + " m " + std::to_string(mean) + " s " + std::to_string(score) +"\n");
    #endif
}

void transcript::print_gtf_entry(std::ostream &os, std::string &gene_id, unsigned int trans_id) {
        
    // print transcript line first
    os << chromosome << "\t" << "ryuto" << "\t" << "transcript" << "\t";
    os << std::to_string(exons.begin()->first) << "\t" << std::to_string(exons.rbegin()->second) << "\t";
    os << "0" << "\t"; // TODO: create a SCORE?
    os << strand << "\t" << "." << "\t";
    os << "gene_id \"fres." << gene_id << "\"; ";
    os << "transcript_id \"fres." << gene_id << "." << std::to_string(trans_id) << "\"; ";
    os << "FPKM \"" << std::to_string(fpkm) << "\"; ";
    os << "cov \"" << std::to_string(mean) << "\"; ";
    os << "uniform_cov \"" << std::to_string(flow) << "\";";
    
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
        
        os << "gene_id \"fres." << gene_id << "\"; ";
        os << "transcript_id \"fres." << gene_id << "." << std::to_string(trans_id) << "\"; ";
        os << "FPKM \"" << std::to_string(fpkm) << "\"; ";

        if (cycle_id_in != 0) {
           os << "cycle_in \"fres." << gene_id << "." << std::to_string(trans_id) << "/" << std::to_string(cycle_id_in) << "\"; ";
        }
        if (cycle_id_out != 0) {
           os << "cycle_out \"fres." << gene_id << "." << std::to_string(trans_id) << "/" << std::to_string(cycle_id_out) << "\"; ";
        }
        os << "cov \"" << std::to_string(mean) << "\"; ";
        os << "uniform_cov \"" << std::to_string(flow) << "\"; ";
        
        if (guided) {
            os << "reference_id \"" << guide_reference << "\"; ";
        }
        
        os << "\n";
    
    }
    
}
