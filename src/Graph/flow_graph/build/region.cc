/* 
 * File:   region.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 17, 2017, 11:59 AM
 */

#include <deque>

#include "region.h"
#include <math.h>
#include "../../../Logger/logger.h"

region::region() : total_length(0) {
}

region::~region() {
}

rcount region::get_max() {
    
    if (subregions.empty()) {
        return 0;
    }
    
    std::deque<subregion>::iterator it = subregions.begin();
    
    rcount max = std::max(it->start_count, it->end_count);
    ++it;
    for (; it != subregions.end(); ++it) {
        
        rcount mval = std::max(it->start_count, it->end_count);
        if (mval > max) {
            max = mval;
        } 
    }
    return max;
}

rcount region::get_min() {
    
    if (subregions.empty()) {
        return 0;
    }
    
    std::deque<subregion>::iterator it = subregions.begin();
    
    rcount min = std::min(it->start_count, it->end_count);
    ++it;
    for (; it != subregions.end(); ++it) {
        rcount mval = std::min(it->start_count, it->end_count);
        if (mval < min) {
            min = mval;
        } 
    }
    return min;
}

float region::get_average() {
    
    rpos total_len = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        total_len += it->end - it->start +1;
    }
    float coverage = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        coverage += it->basecount / float(total_len);
    }
    return coverage;
}

float region::get_deviation() {
    
    rpos total_len = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        total_len += it->end - it->start +1;
    }
    float avr = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        avr += it->basecount / float(total_len);
    }
   
    float dev = 0;
    unsigned int count = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        
        dev += ((float) it->start_count - avr) * ((float) it->start_count - avr)  + ((float) it->end_count - avr) * ((float) it->end_count - avr);
        count += 2;
    }
   
    dev = sqrt(dev / count);
    return dev;
}

rcount region::get_base_count() {
    
    rcount base_count = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        base_count += it->basecount;
    }
    return base_count;
}

rpos region::get_length() {
    
    rpos total_len = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        total_len += it->end - it->start +1;
    }
    return total_len;
}

rpos region::get_length_to_first_zero_from_left() {
    
    rpos total_len = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        total_len += it->end - it->start +1;
        if (it->end_count == 0) {
            break;
        }
    }
    return total_len;
}

rpos region::get_length_to_first_zero_from_right() {
    
    rpos total_len = 0;
    for (std::deque<subregion>::reverse_iterator it = subregions.rbegin(); it != subregions.rend(); ++it) {
        total_len += it->end - it->start +1;
        if (it->start_count == 0) {
            break;
        }
    }
    return total_len;
}

float region::get_average_to_first_zero_from_left() {
    
    rpos total_len = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        total_len += it->end - it->start +1;
        if (it->end_count == 0) {
            break;
        }
    }
    float coverage = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {
        coverage += it->basecount / float(total_len);
        if (it->end_count == 0) {
            break;
        }
    }
    return coverage;
}

float region::get_average_to_first_zero_from_right() {
    
    rpos total_len = 0;
    for (std::deque<subregion>::reverse_iterator it = subregions.rbegin(); it != subregions.rend(); ++it) {
        total_len += it->end - it->start +1;
        if (it->start_count == 0) {
            break;
        }
    }
    float coverage = 0;
    for (std::deque<subregion>::reverse_iterator it = subregions.rbegin(); it != subregions.rend(); ++it) {
        coverage += it->basecount / float(total_len);
        if (it->start_count == 0) {
            break;
        }
    }
    return coverage;
}

void region::get_split_average(std::deque<average_region> &ret) {
    
    rcount basecount = 0;
    rpos total_len = 0;
    rpos start = 0;
    rpos end = 0;
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end(); ++it) {

	if (it->basecount == 0) {
	 	if (total_len > 0) {
			ret.push_back(average_region(start, end, basecount / float(total_len)));
		}
		basecount = 0;
		total_len = 0;
		continue;
	}

	if (total_len == 0) {
		start = it->start;
	}

	basecount += it->basecount;
	total_len += it->end - it->start +1;
	end = it->end; 
    }
    if (total_len > 0) {
	ret.push_back(average_region(start, end, basecount / float(total_len)));
    }
}

rcount region::get_pos_max(rpos length) {
    
    if (subregions.empty()) {
        return 0;
    }
    
    std::deque<subregion>::iterator it = subregions.begin();
    
    rcount max = std::max(it->start_count, it->end_count);
    rcount collected_length = it->end - it->start + 1;
    ++it;
    for (; it != subregions.end() && collected_length < length; ++it) {
                
        rcount mval = std::max(it->start_count, it->end_count);
        if (mval > max) {
            max = mval;
        }
        
        collected_length += it->end - it->start + 1;
    }
    return max;
    
}

rcount region::get_pos_max_reverse(rpos length) {
    
    if (subregions.empty()) {
        return 0;
    }
    
    std::deque<subregion>::reverse_iterator it = subregions.rbegin();
    
    rcount max = std::max(it->start_count, it->end_count);
    rcount collected_length = it->end - it->start + 1;
    ++it;
    for (; it != subregions.rend() && collected_length < length; ++it) {
                
        rcount mval = std::max(it->start_count, it->end_count);
        if (mval > max) {
            max = mval;
        }
        
        collected_length += it->end - it->start + 1;
    }
    return max;
    
}

rcount region::get_left() {
    if (subregions.empty()) {
        return 0;
    }
    return subregions.begin()->start_count;
}

rcount region::get_right() {
    if (subregions.empty()) {
        return 0;
    }
    return subregions.rbegin()->end_count;
}

bool region::is_increasing_region(rpos length) {
    
    if (total_length < 10 || length < 10) {
        return true;
    }
    
    rcount ups = 1; // set to one cause of division
    rcount down = 1;
    rpos upr = 1;
    rpos downr = 1;
    
    rcount collected_length = 0;
            
    // now we have the monotone regions, so far test without length correction // TODO
    for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end() && collected_length < length; ++it)  {
        
        if (it->start_count > it->end_count) {
            down += it->start_count - it->end_count;
            downr += it->end - it->start;
        } else  {
            ups += it->end_count - it->start_count;
            upr += it->end - it->start;
        }
        
        collected_length += it->end - it->start + 1;
    }
    
    // if we have massive one-sided increase
    return  upr >= downr * 1.5 || ups >= down * 1.5;
}

bool region::is_decreasing_region(rpos length) {
    
    if (total_length < 10 || length < 10) {
        return true;
    }
    
    rcount ups = 1; // set to one cause of division
    rcount down = 1;
    rpos upr = 1;
    rpos downr = 1;
    
    rcount collected_length = 0;
    
    // now we have the monotone regions, so far test without length correction // TODO
    for (std::deque<subregion>::reverse_iterator it = subregions.rbegin(); it != subregions.rend() && collected_length < length; ++it)  {
        
        if (it->start_count > it->end_count) {
            down += it->start_count - it->end_count;
            downr += it->end - it->start;
        } else  {
            ups += it->end_count - it->start_count;
            upr += it->end - it->start;
        }
        
        collected_length += it->end - it->start + 1;
    }

    // if we have massive one-sided increase
    return  downr >= upr * 1.5 || down >= ups * 1.5;
}

void region::set_new_left(rpos left) {
	for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end();)  {
		if (it->start >= left) {
			++it;
		} else {
			it = subregions.erase(it);
		}
	}
}
void region::set_new_right(rpos right) {
	for (std::deque<subregion>::iterator it = subregions.begin(); it != subregions.end();)  {
		if (it->end <= right) {
			++it;
		} else {
			it = subregions.erase(it);
		}
	}
}

