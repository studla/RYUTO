/* 
 * File:   bam_guided_reader.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 16, 2017, 2:25 PM
 */

#include "bam_guided_reader.h"

bam_guided_reader::bam_guided_reader(const std::string &gff_path)  {
    gff_reader = new gffReader(gff_path);
}


bam_guided_reader::~bam_guided_reader() {
    delete gff_reader; 
}


void bam_guided_reader::read_chromosome(std::string file_name, std::string chrom_name, const unsigned int input_prefix, bool last_file) {
    

    // at this point we needs to fill in our info from gff_reader before we call the original function without guiding
    chromosome *chrom_fwd, *chrom_rev; 
    #pragma omp critical(chromosome_map_lock)
    {
        chrom_fwd = &chromosome_map_fwd[chrom_name];
        if (options::Instance()->is_stranded()) {
            chrom_rev = &chromosome_map_rev[chrom_name];
        }        
    }
        
    if (!chrom_fwd->has_guide) {
        gff_reader->initialize_chromosome(chrom_name, chrom_fwd, chrom_rev);
        chrom_fwd->has_guide = true;
        if (options::Instance()->is_stranded()) {
            chrom_rev->has_guide = true;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Region Initialized to \n");
        logger::Instance()->debug("FORWARD \n");
        for(greader_list<connected>::iterator c_it = chrom_fwd->chrom_fragments.begin(); c_it != chrom_fwd->chrom_fragments.end(); ++c_it) {
            for(greader_list<exon* >::iterator f_it = c_it->fossil_exons->begin(); f_it!=c_it->fossil_exons->end(); ++f_it) {
                logger::Instance()->debug("Exon " + std::to_string((*f_it)->start) + " - "+ std::to_string((*f_it)->end) +" \n");
             }
//            for(greader_refsorted_list<raw_atom* >::iterator r_it = c_it->atoms->begin(); r_it != c_it->atoms->end(); ++r_it) {
//                logger::Instance()->debug("Atom \n");
//                for(greader_refsorted_list<exon*>::iterator f_it = (*r_it)->exons->begin(); f_it != (*r_it)->exons->end(); ++f_it) {
//                    logger::Instance()->debug("Exon " + std::to_string((*f_it)->start) + " - "+ std::to_string((*f_it)->end) +" \n");
//                }
//            }
        }
        for(greader_list<raw_atom>::iterator a_it = chrom_fwd->atoms.begin(); a_it != chrom_fwd->atoms.end(); ++a_it) {
            logger::Instance()->debug("Atom \n");
            for(greader_refsorted_list<exon*>::iterator f_it = (a_it)->exons->begin(); f_it != (a_it)->exons->end(); ++f_it) {
                logger::Instance()->debug("Exon " + std::to_string((*f_it)->start) + " - "+ std::to_string((*f_it)->end) +" \n");
            }
        }
        if (options::Instance()->is_stranded()) {
            logger::Instance()->debug("BACKWARD \n");
            for(greader_list<connected>::iterator c_it = chrom_rev->chrom_fragments.begin(); c_it != chrom_rev->chrom_fragments.end(); ++c_it) {
                for(greader_list<exon* >::iterator f_it = c_it->fossil_exons->begin(); f_it!=c_it->fossil_exons->end(); ++f_it) {
                    logger::Instance()->debug("FExon " + std::to_string((*f_it)->start) + " - "+ std::to_string((*f_it)->end) +" \n");
                 }
                for(greader_refsorted_list<raw_atom* >::iterator r_it = c_it->atoms->begin(); r_it != c_it->atoms->end(); ++r_it) {
                    logger::Instance()->debug("Atom \n");
                    for(greader_refsorted_list<exon*>::iterator f_it = (*r_it)->exons->begin(); f_it != (*r_it)->exons->end(); ++f_it) {
                        logger::Instance()->debug("Exon " + std::to_string((*f_it)->start) + " - "+ std::to_string((*f_it)->end)  +" \n");
                    }
                }
            }
        }
        #endif
    }
    // call original function
    bam_reader::read_chromosome(file_name, chrom_name, input_prefix, last_file);
}