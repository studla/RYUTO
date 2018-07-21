/* 
 * File:   range_helper.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 15, 2015, 1:51 PM
 */

#include <deque>
#include <queue>

#include "range_helper.h"
#include "Datatype_Templates/graph_list.h"

range_helper::range_helper() : iterator_stat(NULL) {
}

range_helper::~range_helper() {
}

void range_helper::addRange(arc_range* a) {
    ranges.push_back(a);
}

arc_range* range_helper::get_current_parent() {
    
    if (parent_queue.empty()) {
        return NULL;
    }
    
    return parent_queue.back();
}

void range_helper::add_parent(arc_range* p) {
    
//    logger::Instance()->warning("Range " + std::to_string((long) this) +" add Queue "+ std::to_string( (long) p ) +" \n");
    parent_queue.push_back(p);
}

void range_helper::remove_parent() {
    
//    logger::Instance()->warning("Range " + std::to_string((long) this) +" remove Queue \n");
    
    parent_queue.pop_back();
}


arc_range* range_helper::begin() {
    
//    logger::Instance()->warning("Range " + std::to_string((long) this) +" Begin \n");
    
    if (ranges.empty()) {
        return NULL;
    }
    
    parent_queue.clear();
    
    list_stat = ranges.begin();
    iterator_stat = (*list_stat)->begin(this);
    
//     logger::Instance()->warning("Range " + std::to_string((long) this) +" Ret " + std::to_string(iterator_stat->start) +" " + std::to_string(iterator_stat->end) +" \n");
    return iterator_stat;
}

arc_range* range_helper::next() {
    
//    logger::Instance()->warning("Range " + std::to_string((long) this) +" Next \n");
    
    if (iterator_stat == NULL) { // guard only!
//        logger::Instance()->warning("Range " + std::to_string((long) this) +" Ret NULL \n");
        return NULL;
    }
    
    iterator_stat = iterator_stat->next(this);
    
    if (iterator_stat == NULL) {
        ++list_stat;
        if (list_stat == ranges.end()) {
//            logger::Instance()->warning("Range " + std::to_string((long) this) +" Ret NULL \n");
            return NULL; // global end
        }
        iterator_stat = (*list_stat)->begin(this);
    }
    
//    logger::Instance()->warning("Range " + std::to_string((long) this) +" Ret " + std::to_string(iterator_stat->start) +" " + std::to_string(iterator_stat->end) +" \n");
    return iterator_stat;
}

arc_range* range_helper::previous() {
    
//    logger::Instance()->warning("Range " + std::to_string((long) this) +" Previous \n");
    
    iterator_stat = iterator_stat->previous(this);
    
    if (iterator_stat == NULL) {
        if (list_stat == ranges.begin()) {
            return NULL; // global end
        }
        --list_stat;
        iterator_stat = (*list_stat)->last(this);
    }
    
    return iterator_stat;
}

// split up the CURRENT split on iterator
arc_range* range_helper::split(unsigned int i, arc_range* left, arc_range* right) {
    
 //    logger::Instance()->warning("Split AT (" + std::to_string(iterator_stat->start)+ "," + std::to_string(iterator_stat->end) + ")" +" " + std::to_string((long) this) + " l " + std::to_string((long) left) + " " + std::to_string(left->start) + "," + std::to_string(left->end)  + " r " + std::to_string((long) right)  + " " + std::to_string(right->start) + "," + std::to_string(right->end) + ".\n");
    
//    left->parent = iterator_stat;
//    right->parent = iterator_stat;
    iterator_stat->left = left;
    iterator_stat->right = right;
    
    add_parent(iterator_stat);
    iterator_stat = right;
    return iterator_stat;
}