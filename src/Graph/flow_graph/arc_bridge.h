/* 
 * File:   arc_bridge.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on June 7, 2016, 8:56 AM
 */

#ifndef ARC_BRIDGE_H
#define	ARC_BRIDGE_H

#include "../../Datatype_Templates/maps.h"
#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/lazy.h"

class arc_bridge {
public:
    
    typedef path_evidence_map<int, rcount>::iterator iterator;
    
    arc_bridge();
    virtual ~arc_bridge();
    
    lazy<path_evidence_map<int, rcount> > bridges;
    
    bool is_evidenced_path();
    bool has_path_evidence(int id);
    rcount get_path_evidence(int id);
    void replace_evidence(int old_id, int new_id);
    void add_evidence_if(int old_id, int new_id);
    void remove_id(int id);
        
    int size();
    
    iterator begin();
    iterator end();
    
    void clear();
    
    rcount& operator[] (int id);
   
private:

};

#endif	/* ARC_BRIDGE_H */

