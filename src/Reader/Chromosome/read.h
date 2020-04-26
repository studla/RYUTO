/* 
 * File:   read.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 2, 2015, 4:12 PM
 */

#ifndef READ_H
#define	READ_H


#include <string>
#include "../../Datatype_Templates/reader_list.h"
#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/lazy.h"
#include "../../Datatype_Templates/maps.h"

class interval;
class raw_atom;

class rread {
public:
    rread();
    rread(const std::string &id, int index);
    rread(const rread &r);
    virtual ~rread();
    
    // for finding pairs again
    lazy<gmap<int, greader_list<std::string> > >  ids;
    unsigned int global_count;
    gmap<int, unsigned int> count;
    
    bool id_set;
     
    raw_atom* atom; // the atom this is transformed to
    
    rpos left_limit, right_limit, length;
    
    bool block;
    bool primary;
    
    raw_atom* create_atom();
    
    rpos get_left_limit();
    rpos get_right_limit();
    
    void set_left_limit(rpos ll);
    void set_right_limit(rpos rl);
    
    void add_length(interval* i) ;
    rpos get_length();
    
    void add_count(int index);
    void add_id(std::string &s, int index);
    
private:

};

#endif	/* READ_H */

