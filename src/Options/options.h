/* 
 * File:   options.h
 * Author: thomas
 *
 * Created on November 5, 2015, 5:35 PM
 */

#ifndef OPTIONS_H
#define	OPTIONS_H

#include <string>
#include <vector>
#include <boost/cstdint.hpp>
#include <sys/stat.h>
#include <map>

class options{
    
private:   
    // hide block cause singleton
    options();
    virtual ~options();
    options(options const&){};
    options& operator=(options const&){};  // assignment operator is private

public:
    
    enum strandedness : uint8_t { 
        unknown, FF, FR, RF, RR
    };
    
    enum delete_on : uint8_t { 
        single, group_vote_low, group_vote_high, all
    };
    
    
private:
    
    // the singleton :)
    static options* instance;
    
    // actual option values
    
    // general options
    std::vector<std::string> bam_files;             
    std::string gtf_file = "";                      
    
    int pool = 1;
    
    std::map<unsigned int, unsigned int> input_to_id;
    //std::map< unsigned int, std::vector<unsigned int> > group_to_inputs;
    bool grouped = false;
    
    bool parse_heuristic = true; 
    bool compute_all_singles = false;
    bool print_graphs = false;           
    bool debug = false; 
    int parallel_chromosomes;
    int parallel_graphs;
    std::string outdir;
    
    // exon building options
    unsigned int max_pos_extend = 1;        
    unsigned int min_coverage = 1;      
    bool stranded = false;
    // if no XS is given, which library do we look at?
    strandedness strand_type = unknown;
    unsigned int min_intron_length = 5; // base pairs           
    unsigned int maximal_intron_size = 1000000;                
    
    // Graph Building and flow options
    float guide_trust_level = 0.2;
    bool graph_denoising = true;                                
    bool keep_introns = false;                                  
    unsigned int intron_retention_threshold = 30; // 30% of overlap to both sides!
    unsigned int broken_intron_retention_threshold = 15; 
    
    unsigned int exon_join_distance = 5; //0; // base pairs    
    bool trimming = true;
    float trimming_rate = 0.2;
    
    // Graph Filter options
    unsigned int minimal_edge_capacity = 1; // filter out edges below this before flow computations
   // unsigned int arc_filter = 8;  // adjourning edges lower than % are killed  // was  (6)
    unsigned int min_junction_coverage = 1; // filter out junctions with less coverage than this (excluding those < min_junction_anchor)! // ADD
    unsigned int min_junction_anchor = 10; // filter out junctions with less length evidence than this!  (10) // ADD
    float low_edge_mark = 4;
    
//    delete_on delete_vote = delete_on::group_vote;
    unsigned int vote_percentage_low = 30;
    unsigned int vote_percentage_high = 60;

    // Output Filter options
    unsigned int percent_filter = 0; // percent minimal value to filter final transcripts
    unsigned int capacity_filter = 0; // absolute minimal value to filter final transcripts    //internal
    float mean_filter = 4;
    unsigned int min_single_coverage = 20; // absolute minimal value to filter final transcripts that have no introns
    unsigned int min_single_length = 300; // absolute minimal value to filter final transcripts that have no introns
    unsigned int min_transcript_length_base = 200; // absolute minimal length on assembled transcripts // ADD
    unsigned int min_transcript_length_extension = 50; // absolute minimal length on assembled transcripts // ADD
    float minimal_score = 3;
    float group_mean_min = 5.5;
    
   // technical Graph and Flow Options
    bool min_readlength_set = false;
    unsigned int min_readlength = 30;
    unsigned int min_paired_distance = 0;
    unsigned int max_paired_distance = 100000;
    unsigned int max_arc_split = 500;
    unsigned int min_raw_exon_size = 3; // should be +1 to max_pos_extend!
    unsigned int max_path_enumeration = 10000;
    
    // INTERNAL
    // global constants for now
    bool create_overlap_merge_read = false;
    unsigned int coverage_change_limit = 20;
    
    bool disjoint_input = true; // extra filter not useful so far!
    
    bool secret_exon_mode = false;

    // UNUSED
    unsigned int max_flow_iteration = 5;
    unsigned int coverage_bias_limit = 20; // for potential adding more sources (!later)

    
    
    
    inline bool file_exists (const std::string& name) {
        struct stat buffer;   
        return (stat (name.c_str(), &buffer) == 0); 
    }
    
public:
    // retrieve instance and init function for main
    static options* Instance();
    void init(int argc, char **argv);
    
    // data access
    std::vector<std::string>& get_bam_files() {
        return bam_files;
    }
    
    std::string& get_gtf_file() {
        return gtf_file;
    }

    unsigned int get_pooling() {
        return pool;
    }
    
    bool is_grouped() {
        return grouped;
    }
    
    std::map<unsigned int, unsigned int>&  get_input_to_id() {
        return input_to_id;
    }
    //std::map< unsigned int, std::vector<unsigned int> >& get_group_to_inputs(){
    //    return group_to_inputs;
    //}
    
    unsigned int get_max_pos_extend() {
        return max_pos_extend;
    }
    
    unsigned int get_min_coverage() {
        return min_coverage;
    }
    
    unsigned int get_min_readlength() {
        return min_readlength;
    }
    
    bool is_min_readlength_set() {
        return min_readlength_set;
    }
    
    bool is_stranded() {
        return stranded;
    }
    
    bool is_print_graphs() {
        return print_graphs;
    }
    
    bool is_debug() {
        return debug;
    }
    
    int get_parallel_chromosomes() {
        return parallel_chromosomes;
    }
    
    int get_parallel_graphs() {
        return parallel_graphs;
    }
    
    std::string& get_output_dir() {
        return outdir;
    }

    bool is_parse_heuristic() {
        return parse_heuristic;
    }
    
    bool is_compute_all_singles() {
        return compute_all_singles;
    }
    
    int get_max_enumerated_paths() {
        return max_path_enumeration;
    }
    
    unsigned int get_min_paired_distance() {
        return min_paired_distance;
    }
    
    unsigned int get_max_paired_distance() {
        return max_paired_distance;
    }
    
    unsigned int get_percentage_filter() {
        return percent_filter;
    }
    
    unsigned int get_capacity_filter() {
        return capacity_filter;
    }
    
    float get_mean_filter() {
        return mean_filter;
    }
    
    unsigned int get_minimal_edge_capacity() {
        return minimal_edge_capacity;
    }
    
    unsigned int get_min_raw_exon_size() {
        return min_raw_exon_size;
    }
    
    bool is_create_overlap_merge_read() {
        return create_overlap_merge_read;
    }
    
    bool is_graph_denoising() {
        return graph_denoising;
    }
    
    unsigned int get_max_flow_iteration() {
        return max_flow_iteration;
    }
    
    unsigned int get_max_arc_split() {
        return max_arc_split;
    }
    
    strandedness get_strand_type() {
        return strand_type;
    }
    
    unsigned int get_min_single_coverage() {
        return min_single_coverage;
    }
    
    unsigned int get_min_single_length() {
        return min_single_length;
    }
    
    unsigned int get_min_junction_anchor() {
        return min_junction_anchor;
    }
    
    unsigned int get_min_junction_coverage() {
        return min_junction_coverage;
    }
    
    unsigned int get_min_intron_length() {
        return min_intron_length;
    }
    
    bool is_keep_introns(){
        return keep_introns;
    }
    unsigned int get_intron_retention_threshold() {
        return intron_retention_threshold;
    }
    unsigned int get_broken_intron_retention_threshold() {
        return broken_intron_retention_threshold;
    }
 
    unsigned int get_exon_join_distance() {
        return exon_join_distance;
    }
    
    unsigned int get_maximal_intron_size() {
        return maximal_intron_size;
    }
    
    float get_guide_trust_level() {
        return guide_trust_level;
    }
    
    bool is_disjoint_input() {
        return disjoint_input;
    }
    
    bool is_trimming(){
        return trimming;
    }
    
    float get_low_edge_mark() {
        return low_edge_mark;
    }
    
    float get_minimal_score() {
        return minimal_score;
    }
    
    float get_group_mean_min() {
        return group_mean_min;
    }
    
    float get_trimming_rate() {
        return trimming_rate;
    }
    
    unsigned int get_min_transcript_length_base() {
        return min_transcript_length_base;
    }
    
    unsigned int get_min_transcript_length_extension() {
        return min_transcript_length_extension;
    }
    
    unsigned int get_vote_percentage_low() {
        return vote_percentage_low;
    }
    
    unsigned int get_vote_percentage_high() {
        return vote_percentage_high;
    }
       
    // internal
    unsigned int get_coverage_change_limit() {
        return coverage_change_limit;
    }

    bool is_secret_exon_counting() {
        return secret_exon_mode;
    }
    
    // unused
    unsigned int get_coverage_bias_limit() {
        return coverage_bias_limit;
    }
    
    bool vote(unsigned int vote_count, unsigned int total_count, delete_on delete_vote) {
        bool result;
        switch (delete_vote){
            case all:
                result = vote_count == total_count;
                break;
            case group_vote_high:
                result = vote_count * 100 / total_count >= vote_percentage_high;
                break;
            case group_vote_low:
                result = vote_count * 100 / total_count >= vote_percentage_low;
                break;
            case single:
                result = vote_count > 0;
                break;    
        } 
        return result;
    }    
    
//    delete_on get_delete_vote_option() {
//        return delete_vote;
//    }
    
//    unsigned int get_delete_vote_percentage() {
//        return delete_vote_percentage;
//    }
    

};

#endif	/* OPTIONS_H */

