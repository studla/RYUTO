/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   exon_group_count.cc
 * Author: thomas
 * 
 * Created on October 22, 2018, 10:04 AM
 */

#include "exon_group_count.h"

exon_group_count::exon_group_count() : read_count(0), frag_count(0), total_lefts(0), total_rights(0) {
}

exon_group_count::~exon_group_count() {
}

void exon_group_count::init(unsigned int exon_size) {
    
    hole_starts.assign(exon_size, std::map< rpos,rcount >());
    hole_ends.assign(exon_size, std::map< rpos,rcount >());
    hole_start_counts.assign(exon_size,0);
    hole_end_counts.assign(exon_size,0);
}