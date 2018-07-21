/* 
 * File:   bam_controller.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 16, 2016, 2:51 PM
 */

#ifndef BAM_CONTROLLER_H
#define	BAM_CONTROLLER_H

#include "base_controller.h"
#include <string>

class bam_controller : public base_controller {
public:
    bam_controller();
    virtual ~bam_controller();
    
    void execute( std::vector<std::string>  &files);
    

    
private:

};

#endif	/* BAM_CONTROLLER_H */

