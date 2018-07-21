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
    
    void print(std::ostream &os);
    void print_gtf(std::ostream &os, std::string &gene_id);
    
    void finalize_borders(exon_meta* meta); 
    void filter_transcripts();
    
    bool empty();
    
    int size();
    
    void join(alternative_transcript_collection &transcripts);
    
private:

};

#endif	/* ALTERNATIVE_TRANSCRIPT_COLLECTION_H */

