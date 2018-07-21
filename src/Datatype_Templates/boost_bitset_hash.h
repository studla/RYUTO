/* 
 * File:   boost_bitset_hash.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 8, 2015, 1:45 PM
 */

#ifndef BOOST_BITSET_HASH_H
#define	BOOST_BITSET_HASH_H

# define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <functional>

namespace std {
    template <> struct hash<boost::dynamic_bitset<> >
    {
        size_t operator()(const boost::dynamic_bitset<>& bitset) const {
            return boost::hash_value(bitset.m_bits);
        }
    };
}

namespace boost {
    template <typename B, typename A>
    std::size_t hash_value(const boost::dynamic_bitset<B, A>& bs) {            
        return boost::hash_value(bs.m_bits);
    }
}

#endif	/* BOOST_BITSET_HASH_H */

