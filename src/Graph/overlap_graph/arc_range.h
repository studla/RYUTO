/* 
 * File:   arc_range.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 16, 2015, 5:50 PM
 */

#ifndef ARC_RANGE_H
#define	ARC_RANGE_H

#include <lemon/list_graph.h>

using namespace lemon;

class range_helper;

class arc_range {
public:
    arc_range(unsigned int start, unsigned int end, ListDigraph::Arc arc);
    virtual ~arc_range();
    
    unsigned int start, end;
    ListDigraph::Arc arc; // arc is deleted from graph if children exist!
        
    arc_range* left;
    arc_range* right;
    //arc_range* parent; // parents are missleading, because they are reset if we have multiple starts
    
    arc_range* begin(range_helper* queue);
    arc_range* last(range_helper* queue);
    
    arc_range* next(range_helper* queue);
    arc_range* previous(range_helper* queue);

    
private:

};

#endif	/* ARC_RANGE_H */

