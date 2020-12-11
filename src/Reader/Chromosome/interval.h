/* 
 * File:   junction_raw.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 2, 2015, 4:10 PM
 */

#ifndef JUNCTION_RAW_H
#define	JUNCTION_RAW_H



#include "../../Datatype_Templates/misc_types.h"
//#include "read_fwd.h"
#include "read.h"

struct interval {
    
    interval( rread* p) {
        parent = p;
    }
    
    rread* parent; // links to original read
    rpos left, right; 
};


bool operator< (const interval& e1, const interval& e2);


const struct junction_sorter
{
    junction_sorter() {};
    
    inline bool operator() (const interval& e1, const interval& e2) {
        
        if (e1.left < e2.left) {
            return true;
        }
        
        if (e1.left > e2.left) {
            return false;
        }
        
        return e1.right < e2.right;
    }
} junction_sorter_singleton;


const struct junction_sorter_left
{
    junction_sorter_left() {};
    inline bool operator() (const interval& e1, const interval& e2) {
        
        return e1.left < e2.left;
    }
} junction_sorter_left_singleton;

const struct junction_sorter_right
{
    junction_sorter_right() {};
    inline bool operator() (const interval& e1, const interval& e2) {
        
        return e1.right < e2.right;
    }
} junction_sorter_right_singleton;

#endif	/* JUNCTION_RAW_H */

