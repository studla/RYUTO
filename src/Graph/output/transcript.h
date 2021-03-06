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
#include "../../Datatype_Templates/maps.h"
#include "../flow_graph/coverage/flow_series.h"

class transcript {
public:
    transcript();
    virtual ~transcript();
    
    void print(std::ostream &os);
    void print_gtf_entry(std::ostream &os, std::string &gene_id, unsigned int trans_id, int main_input_id);
    void print_count_matrix_entry(std::ostream &os, std::string &gene_id, unsigned int trans_id, int main_input_id, std::set<int> &ids);
    
    void finalize_borders(exon_meta* meta);
    
    void join(transcript *o);
    
    exon_edge found_edge;

    struct series_struct {
        capacity_type flow;
        rpos effective_length;
        float mean;
        float score;
        float fpkm;
    };
    gmap<int, series_struct> series;
        
    std::deque<transcript_unsecurity> unsecurity_id;

    unsigned int cycle_id_in;
    unsigned int cycle_id_out;
    
    bool guided;
    bool guide_grouped;
    bool ignore_guide_gene;
    std::string guide_reference;
    std::string guide_gene;

    unsigned int post_filter_regional_group;
    
    // finalized values
    std::deque<std::pair<rpos, rpos> > exons;
    rpos length;
    std::string chromosome;
    std::string strand;
    unsigned int avrg_read_length;
    
private:
    
    

};

#endif	/* TRANSCRIPT_H */

