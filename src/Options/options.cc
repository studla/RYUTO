/* 
 * File:   options.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on November 5, 2015, 5:35 PM
 */

#include "options.h"
#include "Logger/logger.h"

#include <iostream>
#include <boost/program_options.hpp>

#include <sstream>

namespace po = boost::program_options;

options* options::instance = NULL;

options::options() {}

options::~options() {
    if (instance) {
        delete instance;
    }
}

options* options::Instance()
{
   if (!instance)   // Only allow one instance of class to be generated.
      instance = new options;

   return instance;
}

void options::init(int argc, char **argv) {
    po::variables_map vm;
    
    po::options_description basic_o("Basic options");
    basic_o.add_options()
    ("gtf-guide,g", po::value<std::string>(), "Path to a gtf/gff file containing reference genes. (default: none)")
    ("guide-trustlevel", po::value<unsigned int>(), "Level of trust in the given guides. 0 no trust. 100 maximal trust. Abstract non-linear variable! (default: 20)")
    ("pool,w", po::value<unsigned int>(), "Number of consecutive inputs to merge internally in groups.")
    //("diff-exp-groups,e", po::value<std::string>(), "Input grouping for differential expression by index of inputs, e.g. '1-2;3-4'; or '1;2;' .")
    ("nc", po::value<unsigned int>(), "Number of chromosomes processed in parallel. (# of threads = nc * ng; default: 1)")
    ("ng", po::value<unsigned int>(), "Number of graph per chromosome processed in parallel. (# of threads = nc * ng; default: 1)")
    ("out-dir,o", po::value<std::string>(), "The output directory for all data. (default: ./)")
    ("library-type,l", po::value<std::string>(), "fr-unstranded, fr-firststrand, fr-secondstrand, ff-unstranded, ff-firststrand, or ff-secondstrand according to TopHat and Cufflinks definitions. (default: fr-unstranded)")
    ("unstranded,s", "Add this option to ignore strand information and call everything on + strand. (default: none)")
    ("base-coverage-filter,c", po::value<float>(), "Minimum average coverage for features to be filtered. Sets low-edge-mark, mean-filter, and score-filter (default: 4.0)")
    #ifdef ALLOW_DEBUG
    ("debug,d", "Show debug information.")
    #endif
    ("help,h", "Show this help message.")
    ;
        
    po::options_description general_o("General options");
    general_o.add_options()
    ("no-parse-heuristic", "If several input files are used, use this to disable the heuristic to merge fragments. Increases needed RAM. (default: use heuristic)")
    ("compute-singles", "If several input files are used, use this to compute transcripts for all individual transcripts. (default: only joined output)")
    //("print-graphs", "Produce files with intermediate graphs of all genes for each chromosome. (default: no)")
    ;
    
    po::options_description exon_o("Exon building options");
    exon_o.add_options()
    ("pos-merge", po::value<unsigned int>(), "Positions within this range of bases will be clustered together. (default: 1)")
    //("min-cov", po::value<unsigned int>(), "Minimum coverage to support an exon per dataset. (default: 1)")
    ("min-intron-length", po::value<unsigned int>(), "Minimal length of predicted introns. (default: 5)")
    ("max-intron-length", po::value<unsigned int>(), "Maximal length of predicted introns. (default: 1,000,000)")
    ;
    
    po::options_description graph_o("Graph building and flow options");
    graph_o.add_options()
    ("no-denoising", "Add this option to disable denoising of flow graph capacities.")
    ("keep-introns", "Add this option to keep all introns with read coverage intact.")
    ("intron-threshold", po::value<unsigned int>(), "Percentage of coverage to existing splices needed to keep introns intact. (default: 30%)")
    ("broken-intron-threshold", po::value<unsigned int>(), "Percentage of coverage to existing splices needed to keep fragmented intron pieces intact. (default: 15%)")
    ("exon-join-distance", po::value<unsigned int>(), "Join exons no more distant than this amount of bp into one. (default: 5)")
    ("no-trimming", "Disable trimming of predicted transcripts based on coverage.")
    ("trimming-rate", po::value<unsigned int>(), "Don't join exons if one is smaller than this percentage of the bigger (default: 20%)")
    ;

    po::options_description graph_filter_o("Graph filter options");
    graph_filter_o.add_options()
    ("low-edge-mark", po::value<float>(), "Capacity of edges in the splitgraph considered unlikely. (default: 4.0)")
    //("minimal-capacity", po::value<unsigned int>(), "Minimum capacity of edges in the splitgraph. (default: 1)")
    //("arc-filter-percentage", po::value<unsigned int>(), "Arcs accounting for less than this percentage of capacity on a node are removed. (default: 8%)" )
    ("min-junction-coverage", po::value<unsigned int>(), "Minimal coverage to consider intron junctions in the graph. (default: 1)" )
    ("min-junction-anchor", po::value<unsigned int>(), "Minimal length of bp evidence needed left and right of a junction to consider it. (default: 10)" )
    ;
    
    po::options_description filter_o("Transcript output filter options");
    filter_o.add_options()
    ("percent-filter,p", po::value<unsigned int>(), "Minimum capacity to report retrieved transcripts compared to the largest transcript found. (default: 0%)")
    ("mean-filter", po::value<float>(), "Minimum average capacity to report retrieved transcripts. (default: 4.0)")
    ("score-filter", po::value<float>(), "Minimum score to report retrieved transcripts. (default: 3.0)")
    ("min-single-coverage", po::value<unsigned int>(), "Minimal average coverage off single exon transcripts to report. (default: 20)" )
    ("min-single-length", po::value<unsigned int>(), "Minimal length off single exon transcripts to report. (default: 300)" )
    ("min-transcript-length-base", po::value<unsigned int>(), "Minimal length of a transcripts to report in bp. (default: 200)" )
    ("min-transcript-length-extension", po::value<unsigned int>(), "Additional length in bp required per exon. (default: 50)" )
    ("region-group-filter", po::value<float>(), "Regions where all transcripts are below this coverage are discarded. (default: 5.5)")
    ;
    
    po::options_description tech_o("Advanced options");
    tech_o.add_options()
    ("min-readlength,r", po::value<unsigned int>(), "Reads below this length will not be used for atom evidence. (default: 30)")
    ("min-paired-distance,x", po::value<unsigned int>(), "Minimal distance between two paired reads in order to be considered. (default: 0)")
    ("max-paired-distance,y", po::value<unsigned int>(), "Maximal distance between two paired reads in order to be considered. (default: 100000)")
    ("arc-granularity", po::value<unsigned int>(), "Maximum number of arcs per original arc on pseudo-polynomial convex cost flows. (default: 500)")
    ("min-raw-exon-length", po::value<unsigned int>(), "Minimal number of basepairs needed to constitute a raw exon. (default: 3)")
    ("max-path-enumeration", po::value<unsigned int>(), "Maximal number of distinct paths enumerated on each side of an ambiguous node. (default: 10,000)") 
    ;
    
    po::options_description hidden("");
    hidden.add_options()
        ("exon-counting", "Secret mode to output only exon counts.")
        ("sam", po::value<std::vector<std::string> >(), "Input alignment Bam/Sam")
        ;
    
    po::positional_options_description pd;
    pd.add("sam", -1);
    
    po::options_description opts;
    opts.add(basic_o).add(general_o).add(exon_o).add(graph_o).add(graph_filter_o).add(filter_o).add(tech_o).add(hidden);
    
    po::store(po::command_line_parser(argc, argv).options(opts).positional(pd).run(),vm);
    
    if(vm.count("help")) {
        std::cout << "Use as: ./ryuto [options]* <input1.bam> <input2.bam> ..." << std::endl << std::endl;
        std::cout << basic_o << std::endl << std::endl;
        std::cout << " ----- Advanced Options ----- " << std::endl << std::endl;
        std::cout << general_o << std::endl;
        std::cout << exon_o << std::endl;
        std::cout << graph_o << std::endl;
        std::cout << graph_filter_o << std::endl;
        std::cout << filter_o << std::endl;
        std::cout << tech_o << std::endl;
        std::cout << "(Version 1.5m)" << std::endl;
        std::exit(0);
    }
    
//    if(vm.count("print-graphs")) {
//        print_graphs = true;
//    } 
    
    if(vm.count("debug")) {
        debug = true;
    } else {
        debug = false;
    }
    
    if(vm.count("exon-counting")) {
        secret_exon_mode = true;
    }



    if(!vm.count("sam")) {
        logger::Instance()->error("No Input BAM/SAM file.");
        std::exit(1);
    }    
    bam_files = vm["sam"].as<std::vector<std::string>>();
    for (std::vector<std::string>::iterator bf_it = bam_files.begin(); bf_it != bam_files.end(); ++bf_it) {
        if (!file_exists(*bf_it)) {
            logger::Instance()->error("One or more specified BAM/SAM files do not exist.");
            std::exit(1);
        }
    }
    
    if (vm.count("guide-trustlevel")) {
         guide_trust_level = vm["guide-trustlevel"].as<unsigned int>() / 100.0;
    }
    
    if (vm.count("gtf-guide")) {
        gtf_file = vm["gtf-guide"].as<std::string>();
        if (!file_exists(gtf_file)) {
            logger::Instance()->error("GTF guide specified but does not exist.");
            std::exit(1);
        }
        if (guide_trust_level < 0.6) {
            percent_filter = guide_trust_level * 10 + 4;// EXTRA SECURITY
        } else {
            percent_filter = 10;
        }
    } else {
        gtf_file = "";
    }
    
    if (vm.count("pool")) {
         pool = vm["pool"].as<unsigned int>();
    }
    
    /*if (vm.count("diff-exp-groups")) {
        
         // assume grouping is correct
        unsigned int group_counter = 1;
        unsigned int id_counter = 0;
        
        std::string dfs = vm["diff-exp-groups"].as<std::string>();
        std::stringstream ss(dfs);
        std::string item;
        while (std::getline(ss, item, ';'))
        {
            if (item.empty()) {
                continue;
            }
            
            std::vector<unsigned int> new_group;
            
            std::size_t found = item.find("-");
            if (found == std::string::npos) {  // no range, assume single index
                int i = std::stoi(item);
                new_group.push_back(i);
                input_to_id[i-1] = id_counter;
                id_counter++;
            } else {
                int is = std::stoi(item.substr(0, found-1));
                int ie = std::stoi(item.substr(found+1));
                for ( unsigned int i = is; i <= ie ; ++i) {
                    new_group.push_back(i);
                    input_to_id[i-1] = id_counter + (i - is) / pool;
                }
                id_counter += (ie - is) / pool;
            }
            group_to_inputs[group_counter] = new_group;
            group_counter++;
        }
        
        for(unsigned int i = 0; i< bam_files.size(); ++i) {
            if ( input_to_id.find(i) == input_to_id.end()){
                // id was not found in grouping, kill this!
                logger::Instance()->error("Inputs and Groups do not match.");
                std::exit(1);
            }
        }
        
        grouped = true;
    } else { */
        grouped = false;
        
        for(unsigned int i = 0; i< bam_files.size(); ++i) {
            input_to_id[i] = i / pool;
        }
        
   // }
    
    if (vm.count("nc")){
        parallel_chromosomes = vm["nc"].as<unsigned int>();
        
        if (parallel_chromosomes < 1) {
            logger::Instance()->error("Number of threads nc must be > 0.");
            std::exit(1);
        }
        
    } else {
        parallel_chromosomes = 1;
    }
    if (vm.count("ng")){
        parallel_graphs = vm["ng"].as<unsigned int>();
        
        if (parallel_graphs < 1) {
            logger::Instance()->error("Number of threads ng must be > 0.");
            std::exit(1);
        }
        
    } else {
        parallel_graphs = 1;
    }
    
    if (vm.count("out-dir")) {
        outdir = vm["out-dir"].as<std::string>();
    } else {
        outdir = "./";
    }
    
    if(vm.count("unstranded")) {
        stranded = false;
    }
    
    if (vm.count("library-type")) {
        std::string type = vm["library-type"].as<std::string>();
        
        if (type == "fr-unstranded") {
            stranded = true;
            strand_type = unknown;
        } else if (type == "fr-firststrand") {
            stranded = true;
            strand_type = RF;
        } else if (type == "fr-secondstrand") {
            stranded = true;
            strand_type = FR;
        } else if (type == "ff-unstranded") {   
            stranded = true;
            strand_type = FF;
        } else if (type == "ff-firststrand") {
            stranded = true;
            strand_type = FF;
        } else if (type == "ff-secondstrand") {
            stranded = true;
            strand_type = RR;
        }
    }
    
    if(vm.count("base-coverage-filter")) {
        low_edge_mark = vm["base-coverage-filter"].as<float>();
        mean_filter = low_edge_mark;
        minimal_score = low_edge_mark - 1;
    }
    
    if (vm.count("low-edge-mark")) {
        low_edge_mark = vm["low-edge-mark"].as<float>();
    }

    if (vm.count("mean-filter")) {
        mean_filter = vm["mean-filter"].as<float>();
    }
    
    if (vm.count("score-filter")) {
         minimal_score = vm["score-filter"].as<float>();
    }
    
    if (vm.count("no-parse-heuristic")) {
         parse_heuristic = false;
    }
    
    if (vm.count("compute-singles")) {
         compute_all_singles = true;
    }
    
    if(vm.count("pos-merge")){
        max_pos_extend = vm["pos-merge"].as<unsigned int>();
    }    
    
//    if(vm.count("min-cov")){
//        min_coverage = vm["min-cov"].as<unsigned int>();
//    }    
    
    if (vm.count("min-intron-length")) {
         min_intron_length = vm["min-intron-length"].as<unsigned int>();
    }
    if (vm.count("max-intron-length")) {
         maximal_intron_size = vm["max-intron-length"].as<unsigned int>();
    }    
    
    if(vm.count("no-denoising")) {
        graph_denoising = false;
    }
    if (vm.count("keep-introns")) {
         keep_introns = true;
    }
    
    if (vm.count("intron-threshold")) {
         intron_retention_threshold = vm["intron-threshold"].as<unsigned int>();
    }
    if (vm.count("broken-intron-threshold")) {
         broken_intron_retention_threshold = vm["broken-intron-threshold"].as<unsigned int>();
    }    
    
    if (vm.count("exon-join-distance")) {
         exon_join_distance = vm["exon-join-distance"].as<unsigned int>();
    }
    
    if (vm.count("no-trimming")) {
         trimming = false;
    }

    if (vm.count("trimming-rate")) {
         trimming_rate = vm["trimming-rate"].as<unsigned int>() / 100.0;
    }
    
    if (vm.count("min-junction-coverage")) {
         min_junction_coverage = vm["min-junction-coverage"].as<unsigned int>();
    }
    if (vm.count("min-junction-anchor")) {
         min_junction_anchor = vm["min-junction-anchor"].as<unsigned int>();
    }
    
    if (vm.count("percent-filter")) {
         percent_filter = vm["percent-filter"].as<unsigned int>();
    }    
    
    if (vm.count("min-single-coverage")) {
        min_single_coverage = vm["min-single-coverage"].as<unsigned int>();
    }
    
    if (vm.count("min-single-length")) {
        min_single_length = vm["min-single-length"].as<unsigned int>();
    }
    
    if (vm.count("min-transcript-length-base")) {
         min_transcript_length_base = vm["min-transcript-length-base"].as<unsigned int>();
    }
    if (vm.count("min-transcript-length-extension")) {
         min_transcript_length_extension= vm["min-transcript-extension"].as<unsigned int>();
    }
 
    if (vm.count("region-group-filter")) {
         group_mean_min = vm["region-group-filter"].as<float>();
    }    
 
    if (vm.count("min-readlengt")) {
        min_readlength = vm["min-readlengt"].as<unsigned int>();
        min_readlength_set = true;
    } 
    
    if (vm.count("min-paired-distance")) {
         min_paired_distance = vm["min-paired-distance"].as<unsigned int>();
    }
    
    if (vm.count("max-paired-distance")) {
         max_paired_distance = vm["max-paired-distance"].as<unsigned int>();
    }  
    
    if (vm.count("arc-granularity")) {
        max_arc_split = vm["arc-granularity"].as<unsigned int>();
    }
    
    if (vm.count("min-raw-exon-length")) {
        min_raw_exon_size = vm["min-raw-exon-length"].as<unsigned int>();
    }
    
    if (vm.count("max-path-enumeration")) {
         max_path_enumeration = vm["max-path-enumeration"].as<unsigned int>();
    }
    
//    if(vm.count("disjoint-input")) {
//        disjoint_input = true;
//    }

//    if (vm.count("capacity-filter")) {
//         capacity_filter = vm["capacity-filter"].as<unsigned int>();
//    }
    
//    if (vm.count("minimal-capacity")) {
//         minimal_edge_capacity = vm["minimal-capacity"].as<unsigned int>();
//    }

    if (vm.count("max-flow-iteration")) {
         max_flow_iteration = vm["max-flow-iteration"].as<unsigned int>();
    }
//    if (vm.count("arc-filter-percentage")) {
//        arc_filter = vm["arc-filter-percentage"].as<unsigned int>();
//    }

//    if (!disjoint_input) {
//        minimal_edge_capacity = minimal_edge_capacity * bam_files.size();
//        capacity_filter = capacity_filter * bam_files.size();
//        min_single_coverage = min_single_coverage * bam_files.size();
//    }
    
}

