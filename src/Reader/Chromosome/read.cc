/* 
 * File:   read.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 2, 2015, 4:12 PM
 */

#include <deque>

#include "read.h"
#include "interval.h"
#include "raw_atom.h"

rread::rread() : count(1), id_set(false), atom(NULL), left_limit(0), right_limit(0), length(0), block(false) {
}

rread::rread(const std::string &id) :count(1), id_set(true), atom(NULL), left_limit(0), right_limit(0), length(0), block(false) {
    ids.ref().push_back(id);
}

//rread::rread(rread* r) : ids(r->ids), count(r->count), id_set(r->id_set), atom(r->atom), left_limit(r->left_limit), right_limit(r->right_limit), length(r->length) {
//    
//}

rread::rread(const rread& r) : ids(r.ids), count(r.count), id_set(r.id_set), atom(r.atom), left_limit(r.left_limit), right_limit(r.right_limit), length(r.length), block(r.block) {
    
}

rread::~rread() {
    if (atom != NULL) {
        delete atom;
    }
    
}

void rread::add_id(std::string &s) {
    ids.ref().push_back(s);
    id_set = true;
}

raw_atom* rread::create_atom() {
    atom =  new raw_atom();
    return atom;
}

rpos rread::get_left_limit() {
    return left_limit;
}

rpos rread::get_right_limit() {
    return right_limit;
}

void rread::set_left_limit(rpos ll) {
    left_limit = ll;
}

void rread::set_right_limit(rpos rl){
    right_limit = rl;
}

void rread::add_length(interval* i) {
        length += i->right - i->left + 1;
}

rpos rread::get_length() {
    return length;
}