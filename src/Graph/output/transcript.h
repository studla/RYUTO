/* 
 * File:   transcript.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 11, 2016, 1:30 PM
 */

#ifndef TRANSCRIPT_H
#define	TRANSCRIPT_H

#include <deque> 
#include "transcript_unsecurity.h"
#include "../flow_graph/exon_edge.h"


class transcript {
public:
    transcript();
    virtual ~transcript();
    
    void print(std::ostream &os);
    void print_gtf_entry(std::ostream &os, std::string &gene_id, unsigned int trans_id);
    
    void finalize_borders(exon_meta* meta);
    
    void join(transcript *o);
    
    exon_edge found_edge;

    capacity_type flow;
    float mean;
    float score;
    std::deque<transcript_unsecurity> unsecurity_id;

    unsigned int cycle_id_in;
    unsigned int cycle_id_out;
    
    bool guided;
    std::string guide_reference;
    
    // finalized values
    std::deque<std::pair<rpos, rpos> > exons;
    rpos length;
    float fpkm;
    std::string chromosome;
    std::string strand;
    
private:
    
    

};

#endif	/* TRANSCRIPT_H */

