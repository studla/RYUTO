/* 
 * File:   arc_range.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 16, 2015, 5:50 PM
 */

#include "arc_range.h"

#include "range_helper.h"

arc_range::arc_range(unsigned int start, unsigned int end, ListDigraph::Arc arc) :
start(start), end(end), arc(arc), left(NULL), right(NULL){
}

arc_range::~arc_range() {
}

arc_range* arc_range::begin(range_helper* queue) {
    
//    logger::Instance()->warning("In Begin1 " + std::to_string((long) this) + " " + std::to_string(start) + " " + std::to_string(end) + "\n");
            
    if (left == NULL) {
        return this;
    }
    queue->add_parent(this);
//    logger::Instance()->warning("In Begin2 " + std::to_string(queue->parent_queue.size())+ " " + std::to_string((long) this) + " " + std::to_string((long) left) + " " + std::to_string(left->start) + " " + std::to_string(left->end) +"\n");
    return left->begin(queue);
}

arc_range* arc_range::last(range_helper* queue) {
    
//    logger::Instance()->warning("In Last1 " + std::to_string((long) this) + " " + std::to_string(start) + " " + std::to_string(end) + "\n");
    
    if (right == NULL) {
        return this;
    }
    queue->add_parent(this);
    
//    logger::Instance()->warning("In Last2 " + std::to_string(queue->parent_queue.size())+ " " + std::to_string((long) this) + " " + std::to_string((long) left) + " " + std::to_string(left->start) + " " + std::to_string(left->end) +"\n");
    
    return right->last(queue);
}


arc_range* arc_range::next(range_helper* queue) {
    
//    logger::Instance()->warning("In Next " + std::to_string((long) this) + " " + std::to_string(start) + " " + std::to_string(end) + "\n");
    
    arc_range* parent = queue->get_current_parent();
    if( parent == NULL) {
        return NULL;
    } 
    
//    logger::Instance()->warning("In Next 1 Parent  " + std::to_string(parent->start) + " " + std::to_string(parent->end) + "\n");
    
    if (parent->right == this) {
        queue->remove_parent();
        return parent->next(queue);
    }
    
//    logger::Instance()->warning("In Next 2\n");
    
    return parent->right->begin(queue);
}

arc_range* arc_range::previous(range_helper* queue) {
    
//    logger::Instance()->warning("In Prev " + std::to_string((long) this) + "\n");
    
    arc_range* parent = queue->get_current_parent();
    if(parent == NULL) {
        return NULL;
    } 
    
//    logger::Instance()->warning("In Prev 1\n");
    
    if (parent->left == this) {
        queue->remove_parent();
        return parent->previous(queue);
    }
    
//    logger::Instance()->warning("In Prev 1\n");
    
    return parent->left->last(queue);
}