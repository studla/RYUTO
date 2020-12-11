#include <iostream>
#include <string>

#include <sam.h>
#include <hts.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Options/options.h"
#include "Reader/bam_reader.h"
#include "Logger/logger.h"
#include "Controller/bam_controller.h"

int main(int argc, char **argv)
{
    
    #ifdef _OPENMP
        omp_set_nested(1); // allow nesting
        omp_set_dynamic(0); // no automatic adjustments to thread counts
    #endif
    
    options::Instance()->init(argc, argv);        

    bam_controller* controller = new bam_controller(); 
    controller->execute(options::Instance()->get_bam_files());
   
    delete controller;
    return 0;
}
