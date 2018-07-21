/* 
 * File:   bam_reader.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on January 20, 2016, 1:13 PM
 */

#ifndef BAM_READER_H
#define	BAM_READER_H

#include "reader_base.h"

#include <sam.h>
#include <hts.h>
#include "Chromosome/read_collection.h"
#include "../Datatype_Templates/double_deque.h"

class bam_reader : public reader_base {
public:
    bam_reader();
    virtual ~bam_reader();
    
    unsigned long return_read_count(const std::string &file);
    void return_chromosome_names(const std::string &file, greader_name_set<std::string> &return_list);
    virtual void read_chromosome(std::string file, std::string chromosome, const unsigned int input_prefix, bool last_file);
   
    void finalize(const std::string &chromosome);
    unsigned int get_num_connected(const std::string &chromosome);
    void populate_next_single(const std::string &chrom_name, connected *ob, unsigned int &bam_count_absolute, pre_graph* raw, exon_meta* meta);
    bool populate_next_group(const std::string &chrom_name, greader_list<connected> &all_connected, exon_meta* meta, unsigned int &bam_count_absolute);
            
    void discard(const std::string &chrom_name);
    
protected:

    void finish_block( chromosome* chrom,  rpos &left,  rpos &right, r_border_set<rpos>::iterator &ex_start_it,  r_border_set<rpos>::iterator &ex_end_it, bool last_file);
    
    void process_read( bam1_t *bread, rpos &left_border, rpos &right_border,  chromosome* chrom, const std::string id_prefix,
        r_border_set<rpos>::iterator &ex_start_it, r_border_set<rpos>::iterator &ex_end_it, bool last_file);
    rread* parse_read( bam1_t *bread, chromosome* chrom, greader_list<interval> &junctions, greader_list<std::pair<rpos, rpos> > &splices, rpos &left,  rpos &right, const std::string id_prefix,
         r_border_set<rpos>::iterator &ex_start_it,  r_border_set<rpos>::iterator &ex_end_it);
    void add_known_start( chromosome* chrom, const rpos pos,
         r_border_set<rpos>::iterator &ex_start_it, bool evidence);
    void add_known_end( chromosome* chrom, const rpos pos,
         r_border_set<rpos>::iterator &ex_end_it, bool evidence);
    
    greader_list<std::pair<rpos, bool> >::iterator cluster( greader_list<std::pair<rpos, bool> > &in,  greader_list<rpos> &out,  rpos &right);
    void DKMeans( std::vector<std::pair<rpos,unsigned int> > &in,  greader_list<rpos> &out, const unsigned int &extend);
    void add_to_fixed_clusters(lazy<r_border_set<rpos> > &fixed, greader_list<rpos> &new_clust, r_border_set<rpos>::iterator &pos_mark);
    
    void create_raw_exons( chromosome* chrom, greader_list<std::pair<rpos, rpos> > &out, rpos &right);
    void trim_exons_1(chromosome* chrom, greader_list<std::pair<rpos, rpos> > &raw, greader_list<rpos> &starts, greader_list<rpos> &ends);
    void trim_exons_2(chromosome* chrom, greader_list<std::pair<rpos, rpos> > &raw, greader_list<rpos> &starts, greader_list<rpos> &ends);
    bool test_increasing(chromosome* chrom, rpos start, rpos end);
    bool test_decreasing(chromosome* chrom, rpos start, rpos end);
    rcount get_max_region(chromosome* chrom, rpos start, rpos end);
    bool get_average_to_first_zero_from_left(chromosome* chrom, rpos start, rpos end, float& average, rpos & length, rpos & n1, rpos & n2);
    bool get_average_to_first_zero_from_right(chromosome* chrom, rpos start, rpos end, float& average, rpos & length, rpos & n1, rpos & n2);
    
    void update_existing_exons( connected* connected, chromosome* chrom, greader_list<std::pair<rpos, rpos> > &raw, rpos &left,  rpos &right);
    void split_exons( connected* connected, chromosome* chrom, greader_list<rpos> &splits, rpos &left,  rpos &right, int correction);
    
    void solidify_raw_exons_ends(chromosome* chrom, greader_list<std::pair<rpos, rpos> > &raw, greader_list<rpos> &starts, greader_list<rpos> &ends);
    
    greader_list<connected>::iterator  insert_fragment( chromosome* chrom,  rpos &left,  rpos &right);
    void assign_reads( connected* conn, chromosome* chrom);
    void filter_outer_read_junctions(chromosome* chrom);
    void reduce_atoms( connected* conn, chromosome* chrom);
    void mark_or_reduce_paired_atoms( connected* conn, chromosome* chrom , const greader_refsorted_list<raw_atom*>::iterator &atom_start, const greader_refsorted_list<raw_atom*>::iterator &atom_end);
    
    void filter_bins(connected* conn, chromosome* chrom);
    void reduce_reads(connected* conn);
    void reset_reads(chromosome* chrom);
    
    void split_independent_component(connected *conn, greader_list<connected> &connected) ;
    
    rmap chromosome_map_fwd, chromosome_map_rev;
    c_it_map iterator_map;
    
};

#endif	/* BAM_READER_H */

