/* 
 * File:   read_collection.cpp
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on April 14, 2016, 12:32 PM
 */

#include "read_collection.h"

read_collection::read_collection(rpos ll, rpos rl, bool length_filtered, raw_atom* parent ) : parent(parent), left_limit(ll), right_limit(rl), length_filtered(length_filtered) { 
    
}


read_collection::read_collection(rpos ll, rpos rl, rpos length, raw_atom* parent ) : parent(parent), left_limit(ll), right_limit(rl) {
    
    length_filtered = options::Instance()->is_min_readlength_set();
    
    if (options::Instance()->is_min_readlength_set() && length > options::Instance()->get_min_readlength()) {
        length_filtered = false;
    }
     
}

read_collection::~read_collection() {
}

// string MUST be from within list!
void read_collection::flag_id(std::string &id) {
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug( "RM ID " + id  + "\n" );
    #endif
  //  logger::Instance()->debug( "RM ID FOUND " + std::to_string(std::find(open_ids.ref().begin(), open_ids.ref().end(), id) != open_ids.ref().end())  + "\n" );
   // open_ids.ref().erase(std::find(open_ids.ref().begin(), open_ids.ref().end(), id));
    id = "";
}

void read_collection::clean_flagged_ids() {
    for (gmap<int, greader_list<std::string> >::iterator oi_it = open_ids.ref().begin(); oi_it != open_ids.ref().end(); oi_it++) {
        for (greader_list<std::string>::iterator i_it = oi_it->second.begin(); i_it != oi_it->second.end(); ) {
            if (*i_it == "") {
                i_it = oi_it->second.erase(i_it);
            } else {
                ++i_it;
            }
        }
    }
}

bool operator< ( const read_collection& r1, const  read_collection& r2) {
    
    if (r1.left_limit != r2.left_limit ) {
        return r1.left_limit < r2.left_limit;
    } 
    
    return r1.right_limit < r2.right_limit;
}

bool operator== (const read_collection& r1, const read_collection& r2) {
    return r1.left_limit == r2.left_limit && r1.right_limit == r2.right_limit;
}