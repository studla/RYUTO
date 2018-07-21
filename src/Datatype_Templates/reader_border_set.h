/* 
 * File:   graph_list.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on November 16, 2015, 6:05 PM
 */

#ifndef READER_BORDER_SET_H
#define	READER_BORDER_SET_H

#include <deque>
#include <iterator>
#include "../Options/options.h"

// list type for keeping exon borders in a set that is easy to search through
template <class T>
class r_border_set : public std::deque<T> {
    
    public:
    
    typedef typename std::deque<T>::iterator iterator;    
        
    void setify() {
        // sort
        std::sort(this->begin(), this->end());
        
        // remove doubles
        typename std::deque<T>::iterator it = this->begin();
        T* last = &*it;
        ++it;
        for (; it!=this->end(); ) {
            if (*last == *it) {
                it = this->erase(it);
            } else {
                last = &*it;
                ++it;
            }
        }
    }
    
    iterator lower_bound(iterator start, iterator end, const T& value) {
        return std::lower_bound(start, end, value);
    }
    
    bool sorted_find(T element) {
            unsigned int max_extend = options::Instance()->get_max_pos_extend();
            for (typename std::deque<T>::iterator it = this->begin(); it!= this->end() && element +  max_extend >= *it; ++it) {
                if ((*it >= element  && *it - element <= max_extend) || (*it < element  &&  element - *it <= max_extend)) {
                    return true;
                }
            }
            return false;
        }
    
};




#endif	/* READER_BORDER_SET_H */

