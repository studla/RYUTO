/* 
 * File:   unsecurity_id.cpp
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 7, 2018, 11:26 AM
 */

#include "unsecurity_id.h"

unsecurity_id::unsecurity_id() : id(0), resolvable(true) {
}

unsecurity_id::unsecurity_id(unsigned int id, bool r) : id(id), resolvable(r) {

}


unsecurity_id::~unsecurity_id() {
}

