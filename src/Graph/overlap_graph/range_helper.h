/* 
 * File:   range_helper.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 15, 2015, 1:51 PM
 */

#ifndef RANGE_HELPER_H
#define	RANGE_HELPER_H

#include <lemon/list_graph.h>
#include "arc_range.h"
#include "../../Datatype_Templates/graph_list.h"
#include "../../Logger/logger.h"

using namespace lemon;

class range_helper {
public:
    
    range_helper();
    virtual ~range_helper();
    
    void addRange(arc_range* a);
    
    graph_list<arc_range*> ranges;
    
    std::deque<arc_range*> parent_queue;
    
    arc_range* get_current_parent();
    void add_parent(arc_range* p);
    void remove_parent();
    
    arc_range* begin();
    arc_range* next();
    arc_range* previous();
    arc_range* split(unsigned int i, arc_range* left, arc_range* right);
    
private:
    
    typename graph_list<arc_range*>::iterator list_stat;
    arc_range* iterator_stat;
     

};

#endif	/* RANGE_HELPER_H */

