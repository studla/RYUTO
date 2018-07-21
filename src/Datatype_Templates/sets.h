/* 
 * File:   maps.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 8, 2015, 5:45 PM
 */

#ifndef SETS_H
#define	SETS_H

#include <boost/unordered_set.hpp>
#include <unordered_set> 
#include <set>

template <class T>
class pair_set : public std::unordered_set<T> {
    
};

template <class T>
class sorted_set : public std::set<T> {
    
};

template <class T>
class path_evidence_set : public boost::unordered_set<T> {
    
};


template <class T>
bool empty_intersection(const std::set<T>& x, const std::set<T>& y)
{
    typename std::set<T>::const_iterator i = x.begin();
    typename std::set<T>::const_iterator j = y.begin();
    while (i != x.end() && j != y.end())
    {
      if (*i == *j)
        return false;
      else if (*i < *j)
        ++i;
      else
        ++j;
    }
    return true;
}

#endif	/* SETS_H */

