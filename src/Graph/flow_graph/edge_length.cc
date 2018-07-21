/* 
 * File:   edge_length.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 18, 2015, 1:14 PM
 */

#include "edge_length.h"

edge_length::edge_length() : first_exon(0), middle(0), last_exon(0) {
}

edge_length::~edge_length() {
}

rpos edge_length::without_first() {
    return middle + last_exon;
}

rpos edge_length::without_last() {
    return first_exon + middle;
}

std::ostream& operator<<(std::ostream& os, const edge_length& dt)
{
    os << std::to_string(dt.first_exon) << "," << std::to_string(dt.middle)<< "," << std::to_string(dt.last_exon);
    return os;
}