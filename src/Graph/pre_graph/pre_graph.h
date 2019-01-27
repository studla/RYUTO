/* 
 * File:   pre_graph.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on November 17, 2015, 9:48 AM
 */

#ifndef PRE_GRAPH_H
#define	PRE_GRAPH_H

#include "../../Datatype_Templates/graph_list.h"
#include "../overlap_graph/overlap_node.h"
#include "../overlap_graph/contained_node.h"
#include "exon_group.h"


class pre_graph {
public:
    // set with no paired data!
    pre_graph();

    virtual ~pre_graph();
    
  //  void reduce_to_bins();
    
   void set_size(const unsigned int size);
   
   void initialize_exon_gaps_single_raw();
   void initialize_exon_gaps_paired_raw();
   
private:
  
//    void update_parent(const exon_group* base, const exon_group* merged);
    
public:
    
    unsigned int size;
    
    // basic trans splicing can be inserted into these as well
    // as at this point this difference is mood for computations
    // we need DAG for some operations however, so backlinks 
    // (from circular or repeating transcripts) have to be stored extra
    // this can be done by pairing
    
    // all data broken apart to single read parts!
    graph_list<exon_group> singled_bin_list;
    // paired info on single reads
    graph_list<paired_exon_group> paired_bin_list;
        
    // stored circular backlinks
    graph_list<paired_exon_group> chim_circle_bin_list;
    
    //  {     read part 2     } {    read part 1      }
    //     |----| |----| |--  ...  --| |----| |----|
    // read 1/2 have to included IN single bins
    // read 1 to be stored first  in pair (left_pair)
    // pairs define direct backlinks 
    
    
    // we extract paired data and single data that can be used to disambiguate nodes
    // after contracting, nodes always represent an exon unsecurity 
    graph_list<paired_exon_group *>*  pairs_for_exon;
    graph_list<std::pair<exon_group *, gmap<int, rcount> > >*  single_for_exon;
    
    //meta-data // TODO, WHAT DO WE NEED?
    bool paired;
    
    rpos average_fragment_length;
      
    //we keep the info transcripts from annotation
    graph_list<exon_group *> guide_transcripts;
    
    
};

#endif	/* PRE_GRAPH_H */

