/* 
 * File:   unsecurity_id.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 7, 2018, 11:26 AM
 */

#ifndef UNSECURITY_ID_H
#define	UNSECURITY_ID_H

class unsecurity_id {
public:
    unsecurity_id();
    unsecurity_id(unsigned int id, bool r);
    virtual ~unsecurity_id();
    
    unsigned int id;
    bool resolvable;

};

#endif	/* UNSECURITY_ID_H */

