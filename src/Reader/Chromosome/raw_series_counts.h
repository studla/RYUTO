/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   raw_series_counts.h
 * Author: thomas
 *
 * Created on October 18, 2018, 3:24 PM
 */

#ifndef RAW_SERIES_COUNTS_H
#define RAW_SERIES_COUNTS_H

#include "../../Datatype_Templates/lazy.h"
#include "../../Datatype_Templates/misc_types.h"
#include <map>


class raw_series_counts {
public:
    raw_series_counts();
    virtual ~raw_series_counts();
    
    void add_other_max_min(raw_series_counts & rsc, rpos min, rpos max);
    
    lazy<std::map< rpos,rcount > > lefts;
    rcount total_lefts;
    lazy<std::map< rpos,rcount > > rights;
    rcount total_rights;
    
    lazy<std::map< rpos,rcount > > hole_starts;
    lazy<std::map< rpos,rcount > > hole_ends;
    
    rcount count;
    rcount paired_count;
    
};

#endif /* RAW_SERIES_COUNTS_H */

