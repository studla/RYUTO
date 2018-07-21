/* 
 * File:   bam_controller.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 16, 2016, 2:51 PM
 */

#include "bam_controller.h"
#include "../Reader/bam_reader.h"
#include "../Reader/bam_guided_reader.h"

bam_controller::bam_controller() {
}

bam_controller::~bam_controller() {
}

void bam_controller::execute( std::vector<std::string>& files) {
    
    reader_base* reader;
    if (options::Instance()->get_gtf_file().empty()) {
        reader = new bam_reader();
    } else {
        reader = new bam_guided_reader(options::Instance()->get_gtf_file());
    }
    
    base_controller::execute(*reader, files);
    
    delete reader;
}