/* 
 * File:   graph_list.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on November 16, 2015, 6:05 PM
 */

#ifndef READER_LIST_H
#define	READER_LIST_H

#include <deque>
#include <algorithm>
#include "../Options/options.h"
#include "misc_types.h"

template <class T>
struct PtrCompg
{
  bool operator()(const T lhs, const T rhs) const  {
      return *lhs < *rhs;  
  }
};

// list type for reading in exons
template <class T>
class greader_list : public std::deque<T> {
    
    public:
        
        void sort() {
            std::sort(this->begin(), this->end());
        }
        
        void sort_ref() {
            std::sort(this->begin(), this->end(), PtrCompg<T>() );
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





#endif	/* GRAPH_LIST_H */

