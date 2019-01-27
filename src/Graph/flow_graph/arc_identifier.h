/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   arc_identifier.h
 * Author: thomas
 *
 * Created on November 26, 2018, 2:03 PM
 */

#ifndef ARC_IDENTIFIER_H
#define ARC_IDENTIFIER_H

#include "exon_edge.h"
#include "edge_types.h"
#include "edge_length.h"

class arc_identifier {
public:
    arc_identifier();
    virtual ~arc_identifier();
    
    exon_edge edge_specifier;
    edge_types::edge_type edge_type;
    edge_length edge_lengths;
    unsigned int cycle_id_in;
    unsigned int cycle_id_out;
        
private:

};

std::ostream& operator<<(std::ostream& os, const arc_identifier& ai);


#endif /* ARC_IDENTIFIER_H */

