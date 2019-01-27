/* 
 * File:   base_controller.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 2, 2015, 9:39 AM
 */

#include "base_controller.h"
#include <omp.h>
#include <set>
#include <fstream>
#include "../Logger/logger.h"
#include "../Options/options.h"
#include <cstdio>
#include "../Datatype_Templates/move.h"
#include <sys/stat.h>


#include <sys/types.h>
#include <sys/sysinfo.h>

#include <sam.h>
#include <hts.h>
#include <vector>

void print_memory() {
    
    struct sysinfo memInfo;

    sysinfo (&memInfo);
    long long totalVirtualMem = memInfo.totalram;
    //Add other values in next statement to avoid int overflow on right hand side...
    totalVirtualMem += memInfo.totalswap;
    totalVirtualMem *= memInfo.mem_unit;
    
    long long physMemUsed = memInfo.totalram - memInfo.freeram;
    //Multiply in next statement to avoid int overflow on right hand side...
    physMemUsed *= memInfo.mem_unit;
    
    logger::Instance()->info("VirtualTotal " + std::to_string(totalVirtualMem) + " PhysUsed " + std::to_string(physMemUsed) + "\n");
}


base_controller::base_controller() {
}

base_controller::~base_controller() {
}

void base_controller::execute( reader_base& reader,  std::vector<std::string> &files) {
    
    // test we can create the output
    if (!create_path(0777, options::Instance()->get_output_dir()) ) {
        return;
    }
    
    // first we need to get a list of ALL chromosomes
    std::set<int> ids;
    unsigned int i = 0;
    
    greader_name_set<std::string> names;
    std::vector<unsigned long> read_counts;
    for (std::vector<std::string>::iterator str = files.begin(); str != files.end(); ++str, ++i) {
        
        // check that all files are present!
        
        htsFile* file = hts_open(str->c_str(), "r");
        if (file == NULL) {
            logger::Instance()->error("Input unreadable: " + *str +"\n");
            return;
        }
        
        hts_idx_t* idx = sam_index_load(file, str->c_str());
        if (idx == NULL) {
            logger::Instance()->error("No index file not found for " + *str +"\n");
            return;
        }
        
        reader.return_chromosome_names(*str, names);
        read_counts.push_back(reader.return_read_count(*str));
        
        ids.insert(i);
    }
    std::set<int> base_ids = ids;
    int main_id = 0;
    if (ids.size() > 1) {
        ids.insert(-1);
        main_id = -1;
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Thread Set " + std::to_string(options::Instance()->get_parallel_chromosomes()) + " " + std::to_string(options::Instance()->get_parallel_graphs()) + "\n");
    #endif

    // now that we know what to look for, the real fun begins
    // run create chromosomes in parallel
    
    std::vector<const std::string *> name_access;
    name_access.reserve(names.size());
    for (greader_name_set<std::string>::iterator chrom = names.begin(); chrom != names.end(); ++chrom) {
        name_access.push_back(&*chrom);
    }
    
    #pragma omp parallel for num_threads(options::Instance()->get_parallel_chromosomes())
    for (unsigned int c = 0; c < names.size() ; ++c) {
        
        const std::string *chrom = name_access[c];
        
        logger::Instance()->info("Started Reading Chromosome "+ *chrom + " ========================= ");
             #ifdef _OPENMP
                logger::Instance()->info(std::to_string(omp_get_thread_num()+1) + "/" + std::to_string(omp_get_num_threads()));
             #endif
        logger::Instance()->info("\n");
        
        
        // this code is called for each thread individually, with no fixed order
        unsigned int i = 0;
        for (std::vector<std::string>::iterator file = files.begin(); file != files.end(); ++file, ++i) {
            
            logger::Instance()->info("Reading "+ *chrom + " from file "+ *file + "\n");
            
            std::vector<std::string>::iterator next = file;
            ++next;
            
            reader.read_chromosome(*file, *chrom, i, next==files.end() );           
        }
        reader.finalize(*chrom);

        
        logger::Instance()->info("Finalized "+ *chrom + " ========================= \n");
        
        
        // second level of parallelization
        unsigned int total =  reader.get_num_connected(*chrom);
        #pragma omp parallel for num_threads(options::Instance()->get_parallel_graphs())
        for (unsigned int j = 0; j < total; ++j) {
            
            
            logger::Instance()->info("Next Region Reading Chromosome "+ *chrom + " ========================= ");
                #ifdef _OPENMP
                    logger::Instance()->info(std::to_string(omp_get_thread_num()+1) + "/" + std::to_string(omp_get_num_threads()));
                #endif
            logger::Instance()->info("\n");
            
            exon_meta* meta = new exon_meta();
            greader_list<connected> sub_connected;
            reader.populate_next_group(*chrom, sub_connected, meta);
            meta->absolute_reads = read_counts;
            
            unsigned int index = meta->order_index;
            
            std::map<int, alternative_transcript_collection > group_transcripts;
            
            for (greader_list<connected>::iterator ci = sub_connected.begin(); ci != sub_connected.end(); ++ci) {
                pre_graph* graph = new pre_graph();
                
                reader.populate_next_single(*chrom, &*ci, graph, meta);

                #ifdef ALLOW_DEBUG
                print_connected_region(graph, meta);
                #endif

                base_manager * manager = create_flow_handler(graph, meta, *chrom, base_ids);
                bool has_graph;
				if (graph->singled_bin_list.size() < 50000) { // less than 50k atoms
					has_graph = manager->build_extended_splitgraph(); 
				} else { 
					has_graph = manager->build_basic_splitgraph(); 
				}
                if (has_graph) {
                    manager->extract_transcripts_from_flow(); 
                }
                if (manager->has_single_exons()) {
                    manager->add_single_exons();
                }
                if (has_graph || manager->has_single_exons()) {
                    for (std::set<int>::iterator iid = ids.begin(); iid != ids.end(); ++iid) {
                        manager->transcripts[*iid].finalize_borders(meta);
                        group_transcripts[*iid].join(manager->transcripts[*iid]);
                    }
                }
                
                delete graph;
                delete manager;
            }
            
            delete meta;
            
            for (std::set<int>::iterator iid = ids.begin(); iid != ids.end(); ++iid) {
                
                group_transcripts[*iid].input_main_id = main_id;
                group_transcripts[*iid].filter_transcripts(*iid);

                std::ofstream tmp_b(options::Instance()->get_output_dir()+"/tmp_" + *chrom + "_input" + std::to_string(*iid) + "_" + std::to_string(index)+ ".gtf");
                std::string gene_id = std::to_string(c) + "_" + std::to_string(index);
                group_transcripts[*iid].print_gtf(tmp_b, gene_id);
                tmp_b.close();      
            }
        }
        
        reader.discard(*chrom);

        logger::Instance()->info("Join " + *chrom + "\n");
        
        // concat all result files
        
//        if (options::Instance()->is_print_graphs()) {
//            std::ofstream res_a(options::Instance()->get_output_dir()+"/"+*chrom + ".graph") ;
//            for (unsigned int j = 0; j < total; ++j) {
//
//                std::string s_a = options::Instance()->get_output_dir()+"/tmp_" + *chrom + "_" + std::to_string(j)+".graph";
//
//                std::ifstream tmp_a(s_a);
//                if (!tmp_a.fail() && tmp_a.peek() != std::ifstream::traits_type::eof()) {
//                    res_a << tmp_a.rdbuf();
//                }
//                tmp_a.close();
//
//                std::remove(s_a.c_str());
//            }
//            res_a.close();   
//        }
        
        for (std::set<int>::iterator iid = ids.begin(); iid != ids.end(); ++iid) {
            std::ofstream res_b(options::Instance()->get_output_dir()+"/"+*chrom + "_input" + std::to_string(*iid) + ".gtf") ; 
            for (unsigned int j = 0; j < total; ++j) {

                std::string s_b = options::Instance()->get_output_dir()+"/tmp_" + *chrom + "_input" + std::to_string(*iid) + "_" + std::to_string(j)+".gtf";

                std::ifstream tmp_b(s_b) ;

                if (!tmp_b.fail() && tmp_b.peek() != std::ifstream::traits_type::eof()) {
                    res_b << tmp_b.rdbuf();
                }
                tmp_b.close();

                std::remove(s_b.c_str());
            }
            res_b.close();
        }
    }
    
    for (std::set<int>::iterator iid = ids.begin(); iid != ids.end(); ++iid) {
        
        std::string full_name = options::Instance()->get_output_dir()+"/";
        if (*iid == -1 || ids.size() == 1) {
            full_name += "transcripts.gtf";
        } else {
            full_name += "transcripts_" + std::to_string(*iid) + ".gtf";
        }
        
        std::ofstream full(full_name) ; 
        
        for (greader_name_set<std::string>::iterator name_it = names.begin(); name_it != names.end(); ++name_it) {

            std::string chrp = options::Instance()->get_output_dir()+"/"+*name_it + "_input" + std::to_string(*iid) + ".gtf";

            std::ifstream chrg(chrp); 

            if (chrg.fail() || chrg.peek() == std::ifstream::traits_type::eof()) {
                std::remove(chrp.c_str());
                continue;
            }

            full << chrg.rdbuf();
            chrg.close();

           std::remove(chrp.c_str());

        }
        full.close();
    }
}

base_manager * base_controller::create_flow_handler(pre_graph* raw, exon_meta* meta, const std::string &chrom, std::set<int> &ids) {
    // TODO which one / handling
    //return new basic_coverage_flow(raw, meta, chrom);
    //return new mincost_flow_base(raw, meta, chrom);
    //return new mincost_flow_linear_limit(raw, meta, chrom);
    return new mincost_flow_square(raw, meta, chrom, ids);
    //return new mincost_flow_square_hcost(raw, meta, chrom);
    //return new mincost_flow_square_limit(raw, meta, chrom);
}

void base_controller::print_connected_region( pre_graph* graph, exon_meta* meta) {
    
    logger::Instance()->debug("META ----------- \n");
            
    for (int i = 0; i < meta->size; ++i) {
        logger::Instance()->debug("Exon " + std::to_string(meta->exons[i].left) + " - " + std::to_string(meta->exons[i].right) + " Len " + std::to_string(meta->exons[i].exon_length) + ".\n");
    }
    
    logger::Instance()->debug("ATOMS ----------- \n");
    
    for ( graph_list<exon_group>::iterator it = graph->singled_bin_list.begin(); it != graph->singled_bin_list.end(); ++it) {
        logger::Instance()->debug("Atom " + it->bin_mask.to_string() + " ");
//        for( gmap<int, exon_group_count>::iterator egci = it->count_series.begin(); egci != it->count_series.end(); ++egci) {
//            logger::Instance()->debug( std::to_string(egci->first) + " : " + std::to_string(egci->second.read_count) + " " + std::to_string(egci->second.frag_count) + " ; ");
//            logger::Instance()->debug( " ( " + std::to_string(egci->second.total_lefts) + "-" + std::to_string(egci->second.total_rights) + "  ");
//            for ( std::vector<rcount >::iterator hi = egci->second.hole_start_counts.begin(); hi != egci->second.hole_start_counts.end(); ++hi) {
//                logger::Instance()->debug( " S " + std::to_string(*hi) + " ");
//            }
//            for ( std::vector<rcount >::iterator hi = egci->second.hole_end_counts.begin(); hi != egci->second.hole_end_counts.end(); ++hi) {
//                logger::Instance()->debug( " E " + std::to_string(*hi) + " ");
//            }
//        }
        logger::Instance()->debug("\n");
    }
    
    logger::Instance()->debug("PAIRED ATOMS ----------- \n");
    
    for ( graph_list<paired_exon_group>::iterator it = graph->paired_bin_list.begin(); it != graph->paired_bin_list.end(); ++it) {         
        logger::Instance()->debug("Atom Left  " + it->left_read->bin_mask.to_string()+"\n");
        logger::Instance()->debug("Atom Right " + it->right_read->bin_mask.to_string()+"\n");
    }
    
    logger::Instance()->debug("----------- \n");
}


bool base_controller::create_path( mode_t mode, std::string& path )
{
    struct stat st;

    for( std::string::iterator iter = path.begin() ; iter != path.end(); )
    {
         std::string::iterator newIter = std::find( iter, path.end(), '/' );
         std::string newPath = std::string( path.begin(), newIter);

         if( stat( newPath.c_str(), &st) != 0)
         {           
             if( mkdir( newPath.c_str(), mode) != 0 && errno != EEXIST )
             {
                logger::Instance()->error("Cannot create output folder [" + newPath + "] : " + strerror(errno)+"\n");
                return false;
             }
         }
         else
            if( !S_ISDIR(st.st_mode) )
             {
                 errno = ENOTDIR;
                logger::Instance()->error("Output path [" + newPath + "] is not a dir.\n");
                 return false;
             }


         iter = newIter;
         if( newIter != path.end() )
             ++ iter;
    }
    return true;
}
