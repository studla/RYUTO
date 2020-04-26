/* 
 * File:   alternative_transcript_collection.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 13, 2016, 9:51 AM
 */

#ifndef ALTERNATIVE_TRANSCRIPT_COLLECTION_H
#define	ALTERNATIVE_TRANSCRIPT_COLLECTION_H

#include "transcript.h"
#include "../../Datatype_Templates/lazy.h"
#include "../../Datatype_Templates/graph_list.h"
#include "transcript_unsecurity.h"


class alternative_transcript_collection {
public:
    alternative_transcript_collection();

    virtual ~alternative_transcript_collection();
    
    graph_list<lazy<transcript> > transcripts;
    int input_main_id;
    
    std::map<int, capacity_type> total_flow_error;
    capacity_type total_filtered_flow;
    
    void print(std::ostream &os);
    void print_gtf(std::ostream &os, std::string &gene_id);
    
    void print_count_matrix(std::ostream &os, std::string &gene_id, std::set<int> &ids);
    
    void finalize_borders(exon_meta* meta); 
    void filter_transcripts(int id);
    void filter_transcripts(std::set<int> &ids);
    
    bool empty();
    
    int size();
    
    void join(alternative_transcript_collection &transcripts);
    
private:

    void compute_region(std::list<std::pair<rpos, rpos> > &regions);
    void vote(int id, graph_list<lazy<transcript> > &keep, std::list<std::pair<rpos, rpos> > &regions);
    void multi_vote(std::set<int> &ids, graph_list<lazy<transcript> > &keep, std::list<std::pair<rpos, rpos> > &regions);
    void filter_nested(int id, graph_list<lazy<transcript> > &keep);
};

#endif	/* ALTERNATIVE_TRANSCRIPT_COLLECTION_H */

