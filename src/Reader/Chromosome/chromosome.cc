/* 
 * File:   chromosome.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 2, 2015, 4:56 PM
 */

#include <deque>
#include "../../Datatype_Templates/move.h"
#include "chromosome.h"

chromosome::chromosome() : bam_count(0), frag_count(0), read_count(0), average_read_lenghts(0), has_guide(false)  {
}


chromosome::~chromosome() {
}

//TODO:: TODO split with new type

rread* chromosome::addQueuedRead(const rread& r){
    read_queue.push_back(r);
    return &read_queue.back();
}

//interval* chromosome::addInterval(const interval& i) {
//    intervals.push_back(i);
//    return &intervals.back();
//}

void chromosome::split_exon(rpos pos, greader_list<exon* >::iterator &it,  connected* connected) {
    
    // create new exon and insert it, set position data
    
    exon* right = *it;
    fossil_exons.push_back(exon(right->start, pos));
    it = connected->fossil_exons.ref().insert(it, &fossil_exons.back());
    exon* left = *it;
    ++it; // return it to previous state on right element
    
    left->fixed_end = true;
    left->fixed_start = right->fixed_start;
    right->start = pos+1;
    right->fixed_start = true;
    
    
    for ( greader_refsorted_list<raw_atom* >::iterator at  = connected->atoms.ref().begin(); at != connected->atoms.ref().end(); ) {
       
        raw_atom* atom = *at;
        
        exon* first = *atom->exons.ref().begin();
        exon* last = *atom->exons.ref().rbegin();
        
         if (first == right && first == last ) { // single atoms without split are possibly split up!
            
            atoms.push_back(raw_atom()); // create the left atom
            raw_atom* new_atom_left = &atoms.back();
            
            atoms.push_back(raw_atom()); // create the right atom
            raw_atom* new_atom_right =  &atoms.back();
            
            split_atom_singleton(atom, new_atom_left, new_atom_right , pos, left);
            
            at = connected->atoms.ref().erase(at);
            connected->atoms.ref().insert(new_atom_left);
            connected->atoms.ref().insert(new_atom_right); 
            connected->atoms.ref().insert(atom);
             
        } else if (first == right) {
            
            atoms.push_back(raw_atom()); // create the shorter atom
            raw_atom* new_atom =  &atoms.back();
            
            split_atom_start(atom, new_atom, pos, left);
           
            at = connected->atoms.ref().erase(at);
            connected->atoms.ref().insert(new_atom);
            connected->atoms.ref().insert(atom);
            
        } else if (last == right) {
            
            atoms.push_back(raw_atom()); // create the shorter atom
            raw_atom* new_atom = &atoms.back();
            
            split_atom_end(atom, new_atom, pos, left, right);
            
            at = connected->atoms.ref().erase(at);
            connected->atoms.ref().insert(new_atom);
            connected->atoms.ref().insert(atom);

        } else if ( atom->exons.ref().find(right) != atom->exons.ref().end() ) {
        
            // the atom stays unchanged
           atom->exons.ref().insert(left);     // just add in new left exon
           ++at;
           
        } else {
           ++at; 
        }
       
    }
    
}


void chromosome::split_atom_start(raw_atom* atom, raw_atom* new_atom, const rpos &pos, exon* left) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Split Start  " + std::to_string(pos) + "\n");
    #endif
    
    // first process collected reads, their info on pairing stays intact!
    
    lazy<greader_refsorted_list<read_collection*> > left_reads; // extended atoms!
    lazy<greader_refsorted_list<read_collection*> > right_reads;
    
    for (greader_refsorted_list<read_collection*>::iterator r_it = atom->reads->begin(); r_it != atom->reads->end(); ++r_it) {    
        if ( (*r_it)->left_limit <= pos ) {
            // we are left of cutting point, so this reaches into new left exon
            left_reads->insert(left_reads->end(), *r_it);
        } else {
            right_reads->insert(right_reads->end(), *r_it);
        }
    }
    
    rcount lefts_moved = 0;
    for (std::map< rpos,rcount >::iterator l_it = atom->lefts->begin(); l_it != atom->lefts->end();) {
       if ( l_it->first > pos ) {
            // we are left of cutting point, so this reaches into new left exon       
            lefts_moved += l_it->second;
            // we "missuse" the holes for moved starts and ends
            std::map< rpos,rcount >::iterator fl = atom->hole_ends->find(l_it->first);
            if (fl == atom->hole_ends->end()) {
                atom->hole_ends->insert(*l_it);
            } else {
                fl->second += l_it->second;
            }
            l_it = atom->lefts->erase(l_it);
        } else {
            ++l_it;
        }
    }
    atom->total_lefts -= lefts_moved;

    atom->reads = _MOVE(left_reads); // we use the old atom, with reads extending into the new exon removed
    // empties are possible, but since we have a cut, with new reads this will be filled up again
    
    new_atom->reads = _MOVE(right_reads); //assorted reads
    
    // if the parse heuristic is used, we have already eradicated reads unfortunately
    
     // TODO: for now use overestimation
    new_atom->count = atom->count;
    new_atom->paired_count = atom->paired_count;
    new_atom->length_filtered = atom->length_filtered;
    
    // set exons to split
    
    new_atom->exons = lazy<greader_refsorted_list<exon*> >();
    new_atom->exons.ref() = atom->exons.ref(); // linear COPY of RB tree, therefore safe
    // right only
    
    atom->exons.ref().insert(left); // all
   
}


void chromosome::split_atom_end(raw_atom* atom, raw_atom* new_atom, const rpos &pos, exon* left, exon* right) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Split End  " + std::to_string(pos) + "\n");
    #endif

     // first process collected reads, their info on pairing stays intact!
    
    lazy<greader_refsorted_list<read_collection*> > left_reads; // extended atoms!
    lazy<greader_refsorted_list<read_collection*> > right_reads;
    
    for (greader_refsorted_list<read_collection*>::iterator r_it = atom->reads->begin(); r_it != atom->reads->end(); ++r_it) {
        if ( (*r_it)->right_limit <= pos ) {
            // we are left of cutting point, so this reaches into new left exon
            left_reads->insert(left_reads->end(), *r_it);
        } else {
            right_reads->insert(right_reads->end(), *r_it);
        }
    }
    

    rcount rights_moved = 0;
    for (std::map< rpos,rcount >::iterator r_it = atom->rights->begin(); r_it != atom->rights->end();) {
       if ( r_it->first <= pos ) {
                
            rights_moved += r_it->second;
            // we "missuse" the holes for moved starts and ends
            std::map< rpos,rcount >::iterator fr = atom->hole_starts->find(r_it->first);
            if (fr == atom->hole_starts->end()) {
                atom->hole_starts->insert(*r_it);
            } else {
                fr->second += r_it->second;
            }
            r_it = atom->rights->erase(r_it);
        } else {
            ++r_it;
        }
    }
    atom->total_rights -= rights_moved;
    
    
    atom->reads = _MOVE(right_reads); 
    new_atom->reads = _MOVE(left_reads); //assorted reads
    
    // if the parse heuristic is used, we have already eradicated reads unfortunately
    
    // TODO: for now use overestimation
    new_atom->count = atom->count;
    new_atom->paired_count = atom->paired_count;
    new_atom->length_filtered = atom->length_filtered;
    
    // set exons to split
    
    atom->exons.ref().insert(left); // left and right
    
    new_atom->exons = lazy<greader_refsorted_list<exon*> >();
    new_atom->exons.ref() = atom->exons.ref(); // linear COPY of RB tree, therefore safe
    new_atom->exons.ref().erase(right); // just left, cause insert before
    
    // update exons with new atom

}



void chromosome::split_atom_singleton(raw_atom* atom, raw_atom* new_atom_left, raw_atom* new_atom_right , const rpos &pos, exon* left) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Split Singleton  " + std::to_string(pos) + "\n");
    #endif
    
    // first process collected reads, their info on pairing stays intact!
    
    lazy<greader_refsorted_list<read_collection*> > left_reads; // extended atoms!
    lazy<greader_refsorted_list<read_collection*> > right_reads;
    lazy<greader_refsorted_list<read_collection*> > split_reads;
    
    for (greader_refsorted_list<read_collection*>::iterator r_it = atom->reads->begin(); r_it != atom->reads->end(); ++r_it) {
        
        if ((*r_it)->right_limit <= pos ) {
            left_reads->insert(left_reads->end(), *r_it);
        } else if ((*r_it)->left_limit > pos ) {
            right_reads->insert(right_reads->end(), *r_it);
        } else {
            split_reads->insert(split_reads->end(), *r_it);
        }
        
    }
        
    atom->total_rights = 0;
    for (std::map< rpos,rcount >::iterator r_it = atom->rights->begin(); r_it != atom->rights->end();) {

       if ( r_it->first <= pos ) {
                
            // we "missuse" the holes for moved starts and ends
            std::map< rpos,rcount >::iterator fr = atom->hole_starts->find(r_it->first);
            if (fr == atom->hole_starts->end()) {
                atom->hole_starts->insert(*r_it);
            } else {
                fr->second += r_it->second;
            }
            r_it = atom->rights->erase(r_it);
        } else {
           atom->total_rights += r_it->second;
            ++r_it;
        }
    }
    atom->total_lefts = 0;
    for (std::map< rpos,rcount >::iterator l_it = atom->lefts->begin(); l_it != atom->lefts->end();) {

       if ( l_it->first > pos ) {

            // we "missuse" the holes for moved starts and ends
            std::map< rpos,rcount >::iterator fl = atom->hole_ends->find(l_it->first);
            if (fl == atom->hole_ends->end()) {
                atom->hole_ends->insert(*l_it);
            } else {
                fl->second += l_it->second;
            }
            l_it = atom->lefts->erase(l_it);
        } else {
           atom->total_lefts += l_it->second;
            ++l_it;
        }
    }

    // if the parse heuristic is used, we have already eradicated reads unfortunately
    // TODO: for now use overestimation
    new_atom_left->count = atom->count;
    new_atom_left->paired_count = atom->paired_count;
    new_atom_left->length_filtered = atom->length_filtered;
    new_atom_right->count = atom->count;
    new_atom_right->paired_count = atom->paired_count;
    new_atom_right->length_filtered = atom->length_filtered;
    
    // exon manipulation       
    atom->reads = _MOVE(split_reads); 
    new_atom_left->reads = _MOVE(left_reads); //assorted reads
    new_atom_right->reads = _MOVE(right_reads);
    
    new_atom_left->exons.ref().insert(left);  // just left
    
    new_atom_right->exons = lazy<greader_refsorted_list<exon*> >();
    new_atom_right->exons.ref() = atom->exons.ref(); // linear COPY of RB tree, therefore safe
    // just right
    
    atom->exons.ref().insert(left); // both
}