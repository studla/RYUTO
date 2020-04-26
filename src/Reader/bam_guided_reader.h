/* 
 * File:   bam_guided_reader.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 16, 2017, 2:25 PM
 */

#ifndef BAM_GUIDED_READER_H
#define	BAM_GUIDED_READER_H

#include "gffReader.h"
#include "bam_reader.h"

class bam_guided_reader : public bam_reader {
public:
    bam_guided_reader(const std::string &gff_path);

    virtual ~bam_guided_reader();
    
    virtual void read_chromosome(std::vector<std::string> file, std::string chromosome);
    
private:
    
    gffReader* gff_reader;

};

#endif	/* BAM_GUIDED_READER_H */

