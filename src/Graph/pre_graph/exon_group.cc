/* 
 * File:   exon_group.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on November 16, 2015, 12:16 PM
 */        

#include <map>

#include "exon_group.h"
#include "Datatype_Templates/move.h"

exon_group::exon_group(const unsigned int size, const unsigned int se) 
    : size(size), set_exons(se),
      length_filterd(false), has_coverage(false), extended(false), source_evidence(false), drain_evidence(false), reference_atom(false), backlink_fragment(false), bin_mask(size) {
     
}

exon_group::~exon_group() {
}

bool exon_group::operator[] (const unsigned int pos) {
    return bin_mask[pos];
}

// this function only works when the update object is offset to the right
// calling function need to make sure this is the case!
void exon_group::update_mask(const exon_group& eg) {
       
    bin_mask.join_edge(eg.bin_mask);
    range_end = eg.range_end;
    set_exons = bin_mask.id.count();
    reset_maps();
    
}

void exon_group::set(const unsigned int pos, const bool val) {
    bin_mask.set(pos, val);
}

void exon_group::reset_maps() {
    
    count_series.clear();
    
}

// objects should be finalized before using these!
bool operator== ( exon_group& g1,  exon_group& g2){
    
    // are the ranges actually a speedup? TODO: test
    if (g1.range_start != g2.range_start || g1.range_end != g2.range_end) {
        return false;
    }
    
    return g1.bin_mask.id == g2.bin_mask.id;
}

// test for containment,  This is not a sorting!
// g1 contained in g2
bool operator< ( exon_group& g1,  exon_group& g2){
    
    if ( g1.range_start < g2.range_start || g1.range_end > g2.range_end || (g1.range_start == g2.range_start && g1.range_end == g2.range_end) ) {
        return false;
    }
    
    return g1.bin_mask.is_contained_in(g2.bin_mask, g1.range_start, g1.range_end);
}

// g2 contained in g1
bool operator> ( exon_group& g1,  exon_group& g2){
    return g2 < g1;
}

// test for overlap! This is not a sorting!
// start of g2 overlaps with end of g1
bool operator>= ( exon_group& g1,  exon_group& g2){
    if (g1.range_end < g2.range_start || g1.range_end > g1.range_end 
            || g1.range_start > g2.range_start ) {
        return false;
    }
    // test overlap, creating block filter would likely cost more
    for(unsigned int i = g2.range_start; i <= g1.range_end; ++i) {
        if (g1.bin_mask[i] != g2[i]) {
            return false;
        }
    }
    return true;
}

bool operator<(paired_exon_group& p1, paired_exon_group &p2)
{
     // this is only used because of right side unordered map inconsistencies in bam_reader, hence incomplete sorter
    return p1.right_read->bin_mask.true_smaller(p2.right_read->bin_mask);
}

bool pointer_comp_paired_exon_group::operator ()(paired_exon_group* p1, paired_exon_group* p2)
{
    if (p1->left_read->bin_mask == p2->left_read->bin_mask) {
        return p1->right_read->bin_mask.true_smaller(p2->right_read->bin_mask);
    }

    return p1->left_read->bin_mask.true_smaller(p2->left_read->bin_mask);
}