/* 
 * File:   region.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 17, 2017, 11:59 AM
 */

#ifndef REGION_H
#define	REGION_H

#include "../../../Datatype_Templates/misc_types.h"

class region {
public:
    region();
    virtual ~region();
    
    // each subregion is an increasing or decreasing region
    // maximally 50(?) bp long for accuracy
    struct subregion {

        subregion(rpos s, rpos e, rcount sc, rcount ec, rcount bc) : start(s), end(e), start_count(sc), end_count(ec), basecount(bc) {
        };

        rpos start;
        rpos end;
        rcount start_count;
        rcount end_count;
        rcount basecount;
    };
    
    struct average_region {
	average_region(rpos s, rpos e, rcount a) : start(s), end(e), average(a) {};
	rpos start;
	rpos end;
	rcount average;
    };

    rpos total_length;
    std::deque<subregion> subregions; 

    rcount get_max();
    rcount get_min();
    float get_average();
    float get_deviation();
    rcount get_pos_max(rpos length);
    rcount get_pos_max_reverse(rpos length);
    rcount get_right();
    rcount get_left();
    void get_split_average(std::deque<average_region> &ret);
    
    bool is_increasing_region(rpos length);
    bool is_decreasing_region(rpos length);

    void set_new_left(rpos left);
    void set_new_right(rpos right);
    
    rcount get_base_count();
    rpos get_length();
    
    rpos get_length_to_first_zero_from_left();
    rpos get_length_to_first_zero_from_right();
    float get_average_to_first_zero_from_left();
    float get_average_to_first_zero_from_right();
private:

};

#endif	/* REGION_H */

