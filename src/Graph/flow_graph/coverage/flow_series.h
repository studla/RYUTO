/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   flow_series.h
 * Author: thomas
 *
 * Created on October 22, 2018, 1:42 PM
 */

#ifndef FLOW_SERIES_H
#define FLOW_SERIES_H

#include "capacity_mean.h"
#include "../../../Datatype_Templates/maps.h"
#include "../../../Datatype_Templates/misc_types.h"


class flow_series {
public:
    flow_series();
    virtual ~flow_series();
    
    struct series_struct {
        
        capacity_type capacity = 0;
        capacity_type flow = 0;
        
        float region_average;
        float deviation;
        rcount region_max;
        rcount region_min;
        rcount left;
        rcount right;
        rpos length;
        
        rpos length_to_first_zero_left;
        rpos length_to_first_zero_right;
        
        float average_to_first_zero_from_left;
        float average_to_first_zero_from_right;
        
        capacity_mean mean;
    };
    
    gfmap<int, series_struct> series;
    
    capacity_mean combined_mean;
    capacity_type combined_capacity;
    capacity_type combined_flow;
    
    capacity_mean& get_mean(int id);
    capacity_type& get_flow(int id);
    
private:

};

std::ostream& operator<<(std::ostream& os, const flow_series& fs);

#endif /* FLOW_SERIES_H */

