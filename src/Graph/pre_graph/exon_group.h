/* 
 * File:   exon_group.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on November 16, 2015, 12:16 PM
 */

#ifndef EXON_GROUP_H
#define	EXON_GROUP_H

#include "paired_exon_group.h"
#include "../flow_graph/exon_edge.h"
#include "../../Datatype_Templates/lazy.h"
#include <boost/cstdint.hpp>

using namespace boost; 

class exon_group {

public:
        
    exon_group(const unsigned int size, const unsigned int holes);
    
    virtual ~exon_group();
    
    unsigned int size;
    unsigned int set_exons;
    
    unsigned int range_start, range_end;
    
    rcount read_count;
    rcount frag_count;
    
    // all middle exons are fully covered by design
    rcount base_count_start;
    rcount base_count_end;
    
    lazy<std::map< rpos,rcount > > lefts;
    rcount total_lefts;
    lazy<std::map< rpos,rcount > > rights;
    rcount total_rights;
    
    std::vector<std::map< rpos,rcount > > hole_starts;
    std::vector<rcount > hole_start_counts;
    std::vector<std::map< rpos,rcount > > hole_ends;
    std::vector<rcount > hole_end_counts;
    
    bool length_filterd;
    bool filtered;
    bool extended;
    bool source_evidence;
    bool drain_evidence;
     
    bool reference_atom; // was this created as a reference atom
    std::string reference_name;
    
    bool backlink_fragment;
    
    exon_edge bin_mask;
    
public:
    
    bool operator[] (const unsigned int pos);
    void update_mask(const exon_group& eg);
    void set(const unsigned int pos, const bool val);
    void reset_maps();
    
};


// objects should be finalized before using these!
bool operator== ( exon_group& g1,  exon_group& g2);

// test for containment,  This is not a sorting!
// g1 contained in g2
bool operator< ( exon_group& g1,  exon_group& g2);

// g2 contained in g1
bool operator> ( exon_group& g1,  exon_group& g2);

// test for overlap! This is not a sorting!
// start of g2 overlaps with end of g1
bool operator>= ( exon_group& g1,  exon_group& g2);


#endif	/* SPLIT_BIN_H */

