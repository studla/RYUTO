/* 
 * File:   edge_types.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 8, 2015, 4:20 PM
 */

#ifndef EDGE_TYPES_H
#define	EDGE_TYPES_H

#include <boost/cstdint.hpp>

namespace edge_types
{
    
    enum edge_type : uint8_t { 
        BASE_TYPE, EXON, NODE, BACKLINK,  HELPER
    };
    
}

#endif	/* EDGE_TYPES_H */

