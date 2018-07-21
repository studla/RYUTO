/* 
 * File:   raw_atom.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 3, 2016, 4:20 PM
 */

#ifndef RAW_ATOM_H
#define	RAW_ATOM_H

#include "../../Datatype_Templates/reader_sorted_list.h"
#include "../../Datatype_Templates/reader_list.h"
#include "../../Datatype_Templates/lazy.h"

//class exon;
#include "exon.h"
//#include "read_fwd.h"
#include "read.h"
#include "../../Datatype_Templates/maps.h"
#include "read_collection.h"

class raw_atom {
public:
    raw_atom();
    virtual ~raw_atom();

    void remove_read(const read_collection* r);
    
    lazy<greader_refsorted_list<exon*> > exons; // needs < operator on real exon objects
    
    lazy<greader_refsorted_list<read_collection*> > reads; // save reads for later splits/ further processing
    
    // pairing info
    paired_map<raw_atom*, rcount > paired; 
    
    // solidified data
    
    rcount count; // additional to counts in reads!
    rcount paired_count; // additional to info in reads / paired included
    bool length_filtered; // additional to info in reads
    bool source_evidence; // additional to info in reads
    bool drain_evidence; // additional to info in reads
    
    bool reference_atom; // was this created as a reference atom
    bool killed; // was this created as a reference atom
    std::string reference_name;
    
    // all middle exons are fully covered by design, except for merging holes

    // extension to have full controll over coverage!
    lazy<std::map< rpos,rcount > > lefts;
    rcount total_lefts;
    lazy<std::map< rpos,rcount > > rights;
    rcount total_rights;
    
    lazy<std::map< rpos,rcount > > hole_starts;
    lazy<std::map< rpos,rcount > > hole_ends;
    
    unsigned int id; // id by sorting
    
    std::string to_string();
};


// we need to sort atoms unfortunately for compacting them
// since we have not yet any fixation of features, go by exons positions alone
// possibly slow?
bool operator< ( const raw_atom& a1, const  raw_atom& a2);

bool operator== (const raw_atom& a1, const raw_atom& a2);

#endif	/* RAW_ATOM_H */

