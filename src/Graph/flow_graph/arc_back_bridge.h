/* 
 * File:   arc_back_bridge.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on June 15, 2016, 2:59 PM
 */

#ifndef ARC_BACK_BRIDGE_H
#define	ARC_BACK_BRIDGE_H

#include "../../Datatype_Templates/sets.h"
#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/lazy.h"

class arc_back_bridge {
public:
    
    typedef path_evidence_set<int>::iterator iterator;
    
    arc_back_bridge();
    virtual ~arc_back_bridge();
    
    lazy< path_evidence_set<int> > bridges;
    
    bool is_evidenced_path();
    bool has_path_evidence(int id);
    void replace_evidence( int old_id, int new_id);
    void add_evidence_if(int old_id, int new_id);
    void remove_id(int id);
    
    iterator begin();
    iterator end();
    
    int size();
    void clear();
    
private:

};

#endif	/* ARC_BACK_BRIDGE_H */

