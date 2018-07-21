/* 
 * File:   connected.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 9, 2016, 11:05 AM
 */

#ifndef CONNECTED_H
#define	CONNECTED_H

#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/reader_list.h"
#include "../../Datatype_Templates/lazy.h"

#include "raw_atom.h"
#include "read_collection.h"
#include "../../Datatype_Templates/double_deque.h"

class exon;

class connected {
public:
    connected();
    virtual ~connected();
    
    
    rpos start;
    rpos end;
    
    unsigned int bam_count;
    
    rcount intel_count;
    float avg_split;
    
    bool guided;
    
    lazy<greader_list<exon* > > fossil_exons;
   //  lazy<greader_list<read* > > reads; // we probably don't need these
    
    lazy<greader_refsorted_list<raw_atom* > > atoms; // but we need to join atoms, for obvious reasons!
    
    lazy<double_deque_ref<read_collection> > reads; // store them here because they can change between atoms on splits
    
private:

};

#endif	/* CONNECTED_H */

