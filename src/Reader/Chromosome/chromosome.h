/* 
 * File:   chromosome.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 2, 2015, 4:56 PM
 */

#ifndef CHROMOSOME_H
#define	CHROMOSOME_H

#include <set>

#include "../../Datatype_Templates/reader_list.h"
#include "../../Datatype_Templates/reader_refsafe_list.h"
#include "../../Datatype_Templates/reader_border_set.h"
#include "../../Datatype_Templates/double_deque.h"
#include "../../Datatype_Templates/reader_sorted_list.h"

#include "connected.h"
#include "raw_atom.h"
#include "interval.h"
#include "read.h"
#include "exon.h"

class chromosome {
public:
    chromosome();
    virtual ~chromosome();

    rcount frag_count;
    rcount read_count;
    float average_read_lenghts;
     
    // fixed structures from all runs, all share same data from last run
    greader_refsafe_list<exon> fossil_exons; // not ordered
    
    lazy<r_border_set<rpos> > fixed_exon_starts;
    lazy<r_border_set<rpos> > fixed_exon_ends;
    
    greader_list<connected> chrom_fragments; // connected regions
    greader_list<raw_atom> atoms; // atoms used by all connected regions, careful to avoid doubling!
    double_deque<read_collection> reads;  // collection of all reads used by all connected, connected store references to sublists 
    
    // block queues
    greader_list<rread > read_queue;
    greader_list<interval> interval_queue;
    std::map< std::pair<rpos, rpos>, unsigned int > splice_queue;
    
    greader_list<std::pair<rpos, bool> > known_starts;
    greader_list<std::pair<rpos, bool> > known_ends; 
    
    // has guide
    bool has_guide;
    
    rread* addQueuedRead(const rread& r);
    
    // split the given exon at position of the iterator
    void split_exon(rpos pos, greader_list<exon* >::iterator &exon,  connected* connected);
    
    void split_atom_start(raw_atom* atom, raw_atom* new_atom, const rpos &pos,  exon* left);
    void split_atom_end(raw_atom* atom, raw_atom* new_atom, const rpos &pos,  exon* left,  exon* right);
    void split_atom_singleton(raw_atom* atom, raw_atom* new_atom_left, raw_atom* new_atom_split, const rpos &pos,  exon* left);
private:


    
};

#endif	/* CHROMOSOME_H */

