/* 
 * File:   reader_base.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 2, 2015, 5:51 PM
 */

#ifndef READER_BASE_H
#define	READER_BASE_H

#include "Chromosome/chromosome.h"
#include "Chromosome/connection_iterator.h"
#include "../Datatype_Templates/reader_list.h"
#include "../Datatype_Templates/maps.h"
#include "../Datatype_Templates/reader_name_set.h"
#include "../Graph/pre_graph/pre_graph.h"
#include "../Graph/pre_graph/exon_meta.h"
#include <string>

typedef reader_chromosome_map<std::string, chromosome> rmap;
typedef reader_chromosome_map<std::string, connection_iterator> c_it_map;

class reader_base {
public:
    reader_base();
    virtual ~reader_base();
    
    virtual unsigned long return_read_count(const std::string &file) = 0;
    
    virtual void return_chromosome_names(const std::string &file, greader_name_set<std::string> &return_list) = 0;
    
    virtual void read_chromosome(std::vector<std::string > file_names, std::string chromosome) = 0;
    
    
    virtual void finalize(const std::string &chromosome) = 0;
    
    virtual unsigned int get_num_connected(const std::string &chromosome) = 0;
    
    virtual void populate_next_single(const std::string &chrom_name, connected *ob, pre_graph* raw, exon_meta* meta) = 0;
    virtual bool populate_next_group(const std::string &chrom_name, greader_list<connected> &all_connected, exon_meta* meta) = 0;
    
    virtual void discard(const std::string &chromosome) = 0;
     
private:

};

#endif	/* READER_BASE_H */

