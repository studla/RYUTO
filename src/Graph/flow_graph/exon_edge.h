/* 
 * File:   exon_edge.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 7, 2015, 5:01 PM
 */

#ifndef EXON_EDGE_H
#define	EXON_EDGE_H

#include <vector>
#include <string>

#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/boost_bitset_hash.h"
#include "../pre_graph/exon_meta.h"
#include "../../Logger/logger.h"

// This object is mapped to the graph structure to define exons included
// in an edge. Exon lists DO include the nodes the edge connects, as they are 
// essential in uniquely defining it! Therefore edges on the same node WILL
// overlap

class exon_edge {
public:
    exon_edge();
    exon_edge(const exon_edge &e);
    exon_edge(unsigned int maxlength);
    virtual ~exon_edge();
    
     boost::dynamic_bitset<> id;
     int node_index;
     bool left_consecutive;
     bool right_consecutive;
     
     void add_exon(const unsigned int exon);
     bool is_subset_of(const exon_edge &edge);
     bool is_contained_in(const exon_edge &edge, unsigned int range_start, unsigned int range_end);
     void join_edge(const exon_edge& e2);
     void reserve(unsigned int maxlength);
     
     void left_split(const unsigned int exon, exon_edge& left);
     void right_split(const unsigned int exon, exon_edge& right);
  
     bool operator[] (const unsigned int pos);
     
     void set(const unsigned int pos, const bool val);
     
     std::string to_string();
     
     bool true_smaller(const exon_edge& e2);
};

// define new operator for faster comparison! No need to actually loop through
// the set other than to get specific data
bool operator== (const exon_edge& e1, const exon_edge& e2);
bool operator!= (const exon_edge& e1, const exon_edge& e2);
bool operator< (const exon_edge& e1, const exon_edge& e2);

exon_edge& operator |=(exon_edge& e1, exon_edge& e2);

std::ostream& operator<<(std::ostream& os, const exon_edge& dt);

namespace std {
    template <> struct hash<exon_edge>
    {
        size_t operator()(const exon_edge& edge) const {
            return std::hash<boost::dynamic_bitset<> >()(edge.id);
        }
    };
}

namespace boost {
    //template <>
    std::size_t hash_value(const exon_edge& edge);
}

#endif	/* EXON_EDGE_H */

