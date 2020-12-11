/* 
 * File:   paired_exon_group.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 1, 2015, 8:44 PM
 */

#ifndef PAIRED_EXON_GROUP_H
#define	PAIRED_EXON_GROUP_H

#include "../../Datatype_Templates/misc_types.h"
#include "../../Datatype_Templates/maps.h"

class exon_group;

struct paired_exon_group {
public:
    paired_exon_group( exon_group* lr,  exon_group* rr, gmap<int, rcount> &c) : left_read(lr), right_read(rr), count(c) {};
    virtual ~paired_exon_group() {};
    
    exon_group* left_read;
    exon_group* right_read;
    gmap<int, rcount> count; // how often was this encountered
    
private:

};

struct pointer_comp_paired_exon_group {
    bool operator()(paired_exon_group* p1, paired_exon_group* p2);
};


bool operator<(const paired_exon_group& p1, const paired_exon_group &p2);

#endif	/* PAIRED_EXON_GROUP_H */

