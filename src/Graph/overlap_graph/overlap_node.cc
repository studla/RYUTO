/* 
 * File:   overlap_node.cpp
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 1, 2015, 6:14 PM
 */

#include "overlap_node.h"
#include "contained_node.h"


overlap_node::overlap_node( exon_group* ex) : exons(ex), marker(VACANT), activated(true) {
}


overlap_node::~overlap_node() { 
}

void overlap_node::push_link(overlap_node* link) {
    links.push_back(link);
    edge_marker.push_back(false);
}

void overlap_node::push_backlink(overlap_node* link) {
    back_links.push_back(link);
}

void overlap_node::erase_backlink(overlap_node* link) {
    
    back_links.erase( std::find(back_links.begin(), back_links.end(), link));
}

void overlap_node::erase_link(overlap_node* link) {
    
    links.erase( std::find(links.begin(), links.end(), link));
}


void overlap_node::replace_link( overlap_node* old, overlap_node* nn) {
    
    *std::find(links.begin(), links.end(), old) = nn;
}

void overlap_node::replace_backlink( overlap_node* old, overlap_node* nn) {
    
    *std::find(back_links.begin(), back_links.end(), old) = nn;
}

bool operator== (const overlap_node& g1, const overlap_node& g2) {
    return g1.exons == g2.exons;
}