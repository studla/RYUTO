/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   exon_group_count.h
 * Author: thomas
 *
 * Created on October 22, 2018, 10:04 AM
 */

#ifndef EXON_GROUP_COUNT_H
#define EXON_GROUP_COUNT_H

#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/lazy.h"
#include <map>
#include <vector>

class exon_group_count {
    
public:
    exon_group_count();
    virtual ~exon_group_count();
    
    void init(unsigned int size);
    
    rcount read_count;
    rcount frag_count;
    
    lazy<std::map< rpos,rcount > > lefts;
    rcount total_lefts;
    lazy<std::map< rpos,rcount > > rights;
    rcount total_rights;
    
    std::vector<std::map< rpos,rcount > > hole_starts;
    std::vector<rcount > hole_start_counts;
    std::vector<std::map< rpos,rcount > > hole_ends;
    std::vector<rcount > hole_end_counts;
    
private:

};

#endif /* EXON_GROUP_COUNT_H */

