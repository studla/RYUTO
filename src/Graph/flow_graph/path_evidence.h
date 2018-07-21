/* 
 * File:   path_evidence.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on May 24, 2016, 10:41 AM
 */

#ifndef PATH_EVIDENCE_H
#define	PATH_EVIDENCE_H

#include "../../Datatype_Templates/sets.h"
#include "../../Datatype_Templates/maps.h"
#include "../../Datatype_Templates/misc_types.h"
#include <lemon/list_graph.h>

class path_evidence {
public:
    
    path_evidence();
    
    virtual ~path_evidence();
    
    void add_evidence(int id, rcount count);
    void add_blocked(int id);
    
    bool is_blocked(int id);
    bool is_evidence(int id);
    
    void copy_evidence(path_evidence &other);
    
    bool no_evidence();
    
    path_evidence_map<int, rcount>::iterator begin();
    path_evidence_map<int, rcount>::iterator end();
    
private:
    
   path_evidence_set<int> blocked;
   path_evidence_map<int, rcount> evidence;

};

#endif	/* PATH_EVIDENCE_H */

