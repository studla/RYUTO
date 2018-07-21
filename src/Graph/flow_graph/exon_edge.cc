/* 
 * File:   exon_edge.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 7, 2015, 5:01 PM
 * 
 * 
 * This class has multiple purposes. It is used to represent edges in a graph.
 * It is also used to store split information of paired end data. 
 */

#include <algorithm>

#include "exon_edge.h"
#include "Datatype_Templates/move.h"


exon_edge::exon_edge() : node_index(-1), left_consecutive(false), right_consecutive(false) {
}

exon_edge::exon_edge(const exon_edge& e) : id(e.id), node_index(e.node_index), left_consecutive(e.left_consecutive), right_consecutive(e.right_consecutive)
{ 
   //  logger::Instance()->debug("Exon created "+ to_string()  +"\n");
}

// initial constructor
exon_edge::exon_edge(unsigned int maxlength) : id(maxlength), node_index(-1), left_consecutive(false), right_consecutive(false) {
}

// join e2 edge to this
// exons flow from left to right and cannot on last/first element !
void exon_edge::join_edge(const exon_edge& e2)
{ 
    id |= e2.id; // OR them for best results
}

void exon_edge::reserve(unsigned int maxlength) {
    id.resize(maxlength);
}

exon_edge::~exon_edge() {
 //   logger::Instance()->debug("Exon destroyed "+ to_string()  +"\n");
}

void exon_edge::add_exon(const unsigned int exon) {
    id.set(exon, true);
}

bool exon_edge::is_subset_of(const exon_edge& edge) {
    // this uses a blockwise comparison with  1 AND NOT 2
    // should be fast cause cpu near!
    return id.is_subset_of(edge.id);
}

bool exon_edge::is_contained_in(const exon_edge &a, unsigned int range_start, unsigned int range_end) {
    
    boost::dynamic_bitset<> xored = id ^ a.id; 
    // rangestart is not tested itself with find, but cannot always decrease it to before!
    if (xored[range_start]) {
        return false;
    }
        
    unsigned int i = xored.find_next(range_start);
    if (i <= range_end) {
        return false;
    }
    
    return true;
}

 void exon_edge::left_split(const unsigned int exon, exon_edge& left) {
    
    left = exon_edge(*this); // make a clone
    
    // remove unwanted bits
    // we grab into the datatype to make it faster updating on blocks
    unsigned int block = exon / left.id.bits_per_block;
    unsigned int bit_index = left.id.bits_per_block -  exon % left.id.bits_per_block;
    
    //logger::Instance()->debug("Block L " + std::to_string(block) + " " + std::to_string(bit_index) + "\n");
    
    // all 1 block
    boost::dynamic_bitset<>::block_type mask = ~boost::dynamic_bitset<>::block_type(0);
    mask >>= bit_index-1;
    
    // mask for current block
    left.id.m_bits[block] &= mask;
    // all other blocks directly to 0
    for (boost::dynamic_bitset<>::size_type i = block + 1; i < left.id.num_blocks(); ++i) {
      //  logger::Instance()->debug("i " + std::to_string(i) + "\n");
        left.id.m_bits[i] = 0;
    }
    
}

void exon_edge::right_split(const unsigned int exon, exon_edge& right) {
   
    right = exon_edge(*this); // make a clone
    
    // remove unwanted bits
    // we grab into the datatype to make it faster updating on blocks
    unsigned int block = exon / right.id.bits_per_block;
    unsigned int bit_index = exon % right.id.bits_per_block;
    
  //  logger::Instance()->debug("Block R " + std::to_string(block) + " " + std::to_string(bit_index) + "\n");
    
    // all 1 block
    boost::dynamic_bitset<>::block_type mask = ~boost::dynamic_bitset<>::block_type(0);
    mask <<= bit_index;
    
    // mask for current block
    right.id.m_bits[block] &= mask;
    // all other blocks directly to 0
    if (block == 0) {
        return;
    }
    for (boost::dynamic_bitset<>::size_type i = block - 1; i >= 0 ; --i) {
       // logger::Instance()->debug("i " + std::to_string(i) + "\n");
        right.id.m_bits[i] = 0;
        if (i == 0 ){
            break;
        }
    }
}

bool exon_edge::operator[] (const unsigned int pos) {
    return id[pos];
}

void exon_edge::set(const unsigned int pos, const bool val) {
    id.set(pos, val);
}

std::string exon_edge::to_string() {
    std::ostringstream buffer;
    buffer << id;
    std::string str = buffer.str();
    std::reverse(str.begin(), str.end());
    return str;
}

bool exon_edge::true_smaller(const exon_edge& e2) {
    
    for (unsigned i = 0; i < id.size(); ++i) {
        if (id[i] > e2.id[i]) {
            return true;
        } else if (id[i] < e2.id[i]) {
            return false;
        }
    }
    return false;
}

// define new operator for faster comparison! No need to actually loop through
// the set other than to get specific data
bool operator== (const exon_edge& e1, const exon_edge& e2){
    return e1.id == e2.id;
}

bool operator!= (const exon_edge& e1, const exon_edge& e2){
    return ! (e1 == e2);
}

bool operator< (const exon_edge& e1, const exon_edge& e2){
    return e1.id < e2.id;
}

exon_edge& operator |=(exon_edge& e1, exon_edge& e2) {
    e1.id != e2.id;
    return e1;
}


namespace boost {
    //template <>
    std::size_t hash_value(const exon_edge& edge) {            
        return boost::hash_value(edge.id);
    }
}


std::ostream& operator<<(std::ostream& os, const exon_edge& dt)
{
    std::ostringstream buffer;
    buffer << dt.id;
    std::string str = buffer.str();
    std::reverse(str.begin(), str.end());
    os << str;
    os << " " << std::to_string(dt.node_index) << " " << std::to_string(dt.left_consecutive) << " " << std::to_string(dt.right_consecutive);
    return os;
}