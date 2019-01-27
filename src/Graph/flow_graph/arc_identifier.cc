/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   arc_identifier.cc
 * Author: thomas
 * 
 * Created on November 26, 2018, 2:03 PM
 */

#include "arc_identifier.h"

arc_identifier::arc_identifier() : edge_type(edge_types::BASE_TYPE), cycle_id_in(0), cycle_id_out(0) {
}

arc_identifier::~arc_identifier() {
}

std::ostream& operator<<(std::ostream& os, const arc_identifier& ai)
{
    os << ai.edge_specifier << ", " << std::to_string(ai.edge_type) << ", " << ai.edge_lengths;
    return os;
}