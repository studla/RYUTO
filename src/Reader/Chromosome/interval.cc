/* 
 * File:   interval.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 3, 2016, 4:13 PM
 */

#include "interval.h"


bool operator< (const interval& e1, const interval& e2) {
        
        if (e1.left < e2.left) {
            return true;
        }
        
        if (e1.left > e2.left) {
            return false;
        }
        
        return e1.right < e2.right;
}

