/* 
 * File:   graph_list.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on November 16, 2015, 6:05 PM
 */

#ifndef READER_SORTED_LIST_H
#define	READER_SORTED_LIST_H

#include <set>

// keeping the exon lists of each raw atom
template <class T>
class greader_sorted_list : public std::set<T> {
    
    
};


template <class T>
struct PtrComp
{
  bool operator()(const T &lhs, const T &rhs) const  {
      return *lhs < *rhs;  
  }
};

template <class T>
class greader_refsorted_list : public std::set<T, PtrComp<T> > {
    
    public:
        bool is_disjoint(std::set<T, PtrComp<T> > &set2)
        {
            if(this->empty() || set2.empty()) return true;

            typename std::set<T, PtrComp<T> >::iterator 
                it1 = this->begin(), 
                it1End = this->end();
            typename std::set<T, PtrComp<T> >::iterator 
                it2 = set2.begin(), 
                it2End = set2.end();


            // double * for refsorted!
            if(**set2.rbegin() < **it1  || **this->rbegin() < **it2 ) return true;

            while(it1 != it1End && it2 != it2End)
            {
                if(**it1 == **it2) return false;
                if(**it1 < **it2) { it1++; }
                else { it2++; }
            }

            return true;
        }
    
};



#endif	/* READER_SORTED_LIST_H */

