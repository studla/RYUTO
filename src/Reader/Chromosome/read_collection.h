/* 
 * File:   read_collection.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on April 14, 2016, 12:32 PM
 */

#ifndef READ_COLLECTION_H
#define	READ_COLLECTION_H

#include <string>
#include "../../Options/options.h"
#include "../../Logger/logger.h"
#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/reader_list.h"
#include "../../Datatype_Templates/lazy.h"
#include "../../Datatype_Templates/maps.h"

class raw_atom;
class read_collection {
public:
    
    read_collection(rpos left_limit, rpos right_limit, bool length_filtered, raw_atom* parent );
    read_collection(rpos left_limit, rpos right_limit, rpos length, raw_atom* parent );

    virtual ~read_collection();
    
    void flag_id(std::string &id);
    void clean_flagged_ids();
    
    raw_atom* parent;
    
    lazy<gmap<int, greader_list<std::string> > > open_ids; // these are included in count
    
    struct raw_count {
        
        raw_count() : count(0), paired_count(0) {};
        
        rcount count;
        rcount paired_count;
        lazy<std::deque<std::pair<rpos, rpos> > > holes; // these appear within paired joines that have a gap between paired reads
        paired_map<read_collection*, rcount > paired; // not directly registering as paired count as target in other atom!
        
    };
    
    lazy<gmap<int, raw_count > > counts;
    
    rpos left_limit, right_limit;
    bool length_filtered;
    
private:

};


bool operator< ( const read_collection& r1, const  read_collection& r2);

bool operator== (const read_collection& r1, const read_collection& r2);

#endif	/* READ_COLLECTION_H */

