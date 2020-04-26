/* 
 * File:   gffReader.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 15, 2017, 4:06 PM
 */

#ifndef GFFREADER_H
#define	GFFREADER_H

#include <fstream>
#include <iostream>
#include <ios>
#include <string>
#include <set>
#include <deque>
#include <unordered_map>
#include "Chromosome/chromosome.h"

class gffReader {
public:
    gffReader(const std::string &gff_path);

    virtual ~gffReader();
    
    void initialize_chromosome(std::string chrom_name, chromosome *chrom_fwd, chromosome *chrom_rev );
    
private:
   
    struct transcript_info {
    
        public:
            
            transcript_info()  {
            };
            
            std::string name;
            std::string gene;
            rpos start, end;
            // id in the flow graph 
            // start at 0 and keep it tight!
            std::set<std::pair<rpos, rpos> > exons;
            
            bool operator< ( const transcript_info& t2);
            
    };
    
    void parseAttribute(std::string &all_attributes, std::vector<std::string> &transcript_ids, const std::string& attribute);
    void insertExon(std::deque<transcript_info> &trans, std::unordered_map<std::string, transcript_info *> &id_to_trans, std::vector<std::string> &transcript_ids, std::vector<std::string> &gene_ids, rpos &exon_start, rpos &exon_end);
    void add_all_into_chromosome(std::deque<transcript_info> &trans, chromosome *chrom);
    void add_one_into_chromosome(std::deque<transcript_info *> &region, chromosome *chrom, rpos start, rpos end);
    
    bool gff3;
    std::string gff_path;
};



#endif	/* GFFREADER_H */

