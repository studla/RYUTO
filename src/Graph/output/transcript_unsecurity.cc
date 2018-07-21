/* 
 * File:   transcript_unsecurity.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 11, 2016, 1:55 PM
 */

#include "transcript_unsecurity.h"

transcript_unsecurity::transcript_unsecurity(unsigned int id, evidence e) : position(-1), id(id), evidenced(e){
}

transcript_unsecurity::~transcript_unsecurity() {
}

std::ostream& operator<<(std::ostream& os, const transcript_unsecurity& tu)
{
    os << std::to_string(tu.left) << "-" << std::to_string(tu.right) << "ID" << std::to_string(tu.id) << "S" << std::to_string(tu.evidenced);
    return os;
}

bool operator< (const transcript_unsecurity& t1, const transcript_unsecurity& t2){
        
        return t1.id < t2.id;
}

bool operator== (const transcript_unsecurity& t1, const transcript_unsecurity& t2) {
    return t1.id == t2.id && t1.position == t2.position && t1.evidenced == t2.evidenced;
}
