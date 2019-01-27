/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   flow_series.cpp
 * Author: thomas
 * 
 * Created on October 22, 2018, 1:42 PM
 */

#include "flow_series.h"

flow_series::flow_series() : combined_capacity(0), combined_flow(0) {
    
}

flow_series::~flow_series() {
}

capacity_mean& flow_series::get_mean(int id) {
    if (id == -1) {
        return combined_mean;
    }
    return series[id].mean;
}

capacity_type& flow_series::get_flow(int id) {
    if (id == -1) {
        return combined_flow;
    }
    return series[id].flow;
}


std::ostream& operator<<(std::ostream& os, const flow_series& fs) {
    
    for (gfmap<int, flow_series::series_struct>::const_iterator ssi = fs.series.begin(); ssi != fs.series.end(); ++ssi) {
        capacity_mean cc = ssi->second.mean;
        os << std::to_string(ssi->first) << ":" << std::to_string(ssi->second.capacity) << "/" << std::to_string(ssi->second.flow) << "-" << std::to_string(ssi->second.mean.mean) << "(" << std::to_string(ssi->second.mean.weight) << ")-" << std::to_string(cc.compute_score()) << ";"; 
    }
    capacity_mean cc = fs.combined_mean;
    os << "c:" << std::to_string(fs.combined_capacity) << "/" << std::to_string(fs.combined_flow) << "-" << std::to_string(fs.combined_mean.mean) << "(" << std::to_string(fs.combined_mean.weight) << ")-" << std::to_string(cc.compute_score()) << ";";
    
}