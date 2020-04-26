/* 
 * File:   base_controller.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 2, 2015, 9:39 AM
 */

#ifndef BASE_CONTROLLER_H
#define	BASE_CONTROLLER_H

#include <vector>
#include <string>
#include "../Reader/reader_base.h"
#include "../Graph/flow_manager/base_manager.h"
#include "../Graph/flow_manager/mincost_flow_base.h"
#include "../Graph/flow_manager/mincost_flow_square.h"
#include "../Graph/flow_manager/mincost_flow_square_hcost.h"

class base_controller {
public:
    base_controller();
    virtual ~base_controller();
    
    virtual void execute( std::vector<std::string> &files) = 0;

protected:
   
    void execute( reader_base& reader,  std::vector<std::string> &files);
    
private:

    base_manager *  create_flow_handler(pre_graph* raw, exon_meta* meta, const std::string &chrom, std::set<int> &ids);
    void print_connected_region( pre_graph* graph, exon_meta* meta);
    
    bool create_path( mode_t mode, std::string& path);
    std::string getFileName(std::string filePath, bool withExtension = true, char seperator = '/');
};

#endif	/* BASE_CONTROLLER_H */

