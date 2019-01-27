/* 
 * File:   maps.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 8, 2015, 5:45 PM
 */

#ifndef MAPS_H
#define	MAPS_H

#include <boost/unordered_map.hpp> 
#include <boost/functional/hash/extensions.hpp>
#include <functional>
#include <map>
#include <unordered_map>

template <class Key, class Mapped>
class gmap : public boost::unordered_map<Key, Mapped, std::hash<Key> > {
    
};

template <class Key, class Mapped>
class gfmap : public std::map<Key, Mapped > {
    
};

template <class Key, class Mapped>
class reader_chromosome_map : public std::map<Key, Mapped > {
    
};

template <class Key, class Mapped>
class paired_map : public boost::unordered_map<Key, Mapped, std::hash<Key> > {
    
};


template <class Key, class Mapped>
class path_evidence_map : public boost::unordered_map<Key, Mapped, std::hash<Key> >{
    
};

namespace std {
    template <class A, class B> struct hash<std::pair<A, B> >
    {
        size_t operator()(std::pair<A, B> const& v) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, v.first);
            boost::hash_combine(seed, v.second);
            return seed;
        }
    };
}

#endif	/* MAPS_H */

