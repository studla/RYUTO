/* 
 * File:   gffReader.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 15, 2017, 4:06 PM
 */

#include "gffReader.h"
#include <boost/algorithm/string.hpp>
#include <algorithm>    // std::sort


gffReader::gffReader(const std::string &gff_path) : gff_path(gff_path) {
    
    std::ifstream gff(gff_path);
    gff3 = false;
    std::string line;
    if ( std::getline(gff, line)) {
        // we do actually have a first line
        if ( line.length() >= 15  && line.compare(0,15, "##gff-version 3")==0 ) {
            gff3 = true;
        }
    }
    gff.close();
}


gffReader::~gffReader() {
    
}

void gffReader::initialize_chromosome(std::string chrom_name, chromosome *chrom_fwd, chromosome *chrom_rev ) {
 
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Read GFF Guide " + gff_path + " - "+ chrom_name+" .\n");
    #endif
    
    // get file handle from start
    std::ifstream gff(gff_path);
    
    // each entry is one transcript, we find transcript overlaps later    
    std::deque<transcript_info> fwd_trans;
    std::deque<transcript_info> rev_trans;

    std::unordered_map<std::string, transcript_info *> id_to_trans_fwd;
    std::unordered_map<std::string, transcript_info *> id_to_trans_rev;
    
    // we need to loop through whole file
    std::string line;
    while (std::getline(gff, line))
    {
        if (line.length() > 0 && line[0] == '#') {
            continue;
        }
   
        // we do some splitting on tab
        std::vector<std::string> tokens;
        std::istringstream iss(line);
        std::string token;
        while(std::getline(iss, token, '\t')) {
            tokens.push_back(token);
        }
        
        if (tokens.size() < 9 ) {
            logger::Instance()->warning("Format error in guiding input");
            continue;
        }
        
        if (tokens[0] != chrom_name || tokens[2]!= "exon") {
            continue; // we are only interested in actual exons with their parent/transc ids on right chromosome
        }
        
        bool minus_strand = false;
        std::vector<std::string> transcript_ids;
        std::vector<std::string> gene_ids;
        rpos exon_start, exon_end;
        
        // now it depends on format
        if (gff3) {   // gff 3 format
            parseAttribute(tokens[8], transcript_ids, "Parent"); // TODO add Gene ID support
        } else { // gtf or gff 2
            parseAttribute(tokens[8], transcript_ids, "transcript_id"); 
            parseAttribute(tokens[8], gene_ids, "gene_id"); 
        }
        
        if ( tokens[6] == "-") {
            minus_strand = true;
        }
        
        exon_start = std::stoi( tokens[3] );
        exon_end = std::stoi( tokens[4] );
        
        if (exon_start > exon_end) {
            // order can be switched on minus sometimes
            rpos tmp = exon_start;
            exon_start = exon_end;
            exon_end = tmp;
        }
        
        if (options::Instance()->is_stranded() && minus_strand) {
            insertExon(rev_trans, id_to_trans_rev, transcript_ids, gene_ids, exon_start, exon_end);
        } else {
            insertExon(fwd_trans, id_to_trans_fwd, transcript_ids, gene_ids, exon_start, exon_end);
        }
        
    }  
    gff.close();
    
    // the file is completely read in, now create the real base structure
    if (options::Instance()->is_stranded()) {
        add_all_into_chromosome(fwd_trans, chrom_fwd);
        add_all_into_chromosome(rev_trans, chrom_rev);
    } else {
        add_all_into_chromosome(fwd_trans, chrom_fwd);
    }
    
}


void gffReader::parseAttribute(std::string &all_attributes, std::vector<std::string> &transcript_ids, const std::string& attribute) {
    
    std::string token;
    std::istringstream all_attributes_str;
    all_attributes_str.str(all_attributes);
    while(std::getline(all_attributes_str, token, ';')) { // split on ; first
        boost::algorithm::trim( token );

        if (token.length() > attribute.length()+2 && token.compare(0, attribute.length(), attribute)==0 ) {
            // we found the actual attribute
            token.erase(0, attribute.length()+1); // cut of attribute and separator
            if (token[0]=='"') {
                token.erase(0, 1);
            }
            if (token[token.length()-1]=='"') {
                token.pop_back();
            }
                        
            // here we could still have comma seperated list 
            std::string id;
            std::istringstream token_str;
            token_str.str(token);
            while(std::getline(token_str, id, ',')) {
                transcript_ids.push_back(id);
            }
        }
    }
    
}

void gffReader::insertExon(std::deque<transcript_info> &trans, std::unordered_map<std::string, transcript_info *> &id_to_trans, std::vector<std::string> &transcript_ids, std::vector<std::string> &gene_ids, rpos &exon_start, rpos &exon_end) {
    
    std::vector<std::string>::iterator g_id_it = gene_ids.begin();
    for (std::vector<std::string>::iterator t_id_it = transcript_ids.begin(); t_id_it != transcript_ids.end(); ++t_id_it, ++g_id_it) { 
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Insert Exon " + *t_id_it + " " + std::to_string(exon_start) + " - "+ std::to_string(exon_end) +" .\n");
        #endif
        
        transcript_info *info;
        
        std::unordered_map<std::string, transcript_info *>::iterator t = id_to_trans.find(*t_id_it);
        if (t == id_to_trans.end()) {
            
            // we need to insert it new
            trans.push_back(transcript_info());
            info = &trans.back();
            id_to_trans.insert(std::make_pair(*t_id_it, info));
            info->start = exon_start;
            info->end = exon_end;
            info->name = *t_id_it;
            info->gene = *g_id_it;
        } else {
            // we have found it
            info = t->second;
        }
        
        // we have inserted/found trans, now insert exon
        info->exons.insert(std::make_pair(exon_start, exon_end));
        if (exon_end > info->end) {
            info->end = exon_end;
        }
        if (exon_start < info->start) {
            info->start = exon_start;
        }
    }
}


void gffReader::add_all_into_chromosome(std::deque<transcript_info> &trans, chromosome *chrom) {
    
    if (trans.empty()) {
        return;
    }
    
    if (!std::is_sorted(trans.begin(), trans.end())) {
        std::sort(trans.begin(), trans.end());
    }
    
    // we assume transcripts are properly ordered by exon position
    
    std::deque<transcript_info>::iterator t_it = trans.begin();
    std::deque<transcript_info *> region;
    region.push_back(&*t_it);
    rpos start = t_it->start;
    rpos end = t_it->end;
    ++t_it;
    for (; t_it != trans.end(); ++t_it) {

        // test if this is an overlapping region to current
        if (t_it->end >= start && t_it->start <= end) {
            // this is an overlap, join
            start = std::min(start, t_it->start);
            end= std::max(end, t_it->end);
            region.push_back(&*t_it);
            
        } else {

            // insert is solo
            add_one_into_chromosome(region, chrom, start, end);
            region.clear();
            region.push_back(&*t_it);
            start = t_it->start;
            end = t_it->end;
        }
    }
    add_one_into_chromosome(region, chrom, start, end);
    
}

void gffReader::add_one_into_chromosome(std::deque<transcript_info *> &region, chromosome *chrom, rpos start, rpos end) {
    
//    logger::Instance()->debug("ADD INTO \n");
//    for(std::deque<transcript_info *>::iterator r_it = region.begin(); r_it != region.end(); ++r_it) {
//        logger::Instance()->debug("GTF: ");
//        for (std::set<std::pair<rpos, rpos> >::iterator e_it = (*r_it)->exons.begin(); e_it != (*r_it)->exons.end(); ++e_it) {
//             logger::Instance()->debug( "("  + std::to_string(e_it->first) + "," + std::to_string(e_it->second) + ") - ");
//        }
//        logger::Instance()->debug("\n");
//    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Connected " + std::to_string((*region.begin())->start) + " " + std::to_string((*region.rbegin())->end) + ".\n");
    logger::Instance()->debug("Add first Atom " + std::to_string((*region.begin())->start) + " " + std::to_string((*region.begin())->end) + ".\n");
    #endif
    
    std::set<rpos> starts;
    std::set<rpos> ends;
    
    // we create a new connected region
    chrom->chrom_fragments.push_back(connected());
    connected * con = &chrom->chrom_fragments.back();
    
    con->guided = true;
    
    chrom->atoms.push_back(raw_atom());
    raw_atom * at = &chrom->atoms.back();
    con->atoms->insert(at);
    at->reference_atom = true;
    at->reference_name = (*region.begin())->name;
    at->reference_gene = (*region.begin())->gene;

    // loop over exons and create the exon structure 
    std::deque<transcript_info *>::iterator r_it = region.begin();
    for (std::set<std::pair<rpos, rpos> >::iterator e_it = (*r_it)->exons.begin(); e_it != (*r_it)->exons.end(); ++e_it ) {
        
        chrom->fossil_exons.push_back(exon(e_it->first, e_it->second));
        exon * ex = &chrom->fossil_exons.back();
        ex->fixed_start = true;
        ex->fixed_end = true;
        con->fossil_exons->push_back( ex);
        
        at->exons->insert(ex);
        
        starts.insert(e_it->first);
        ends.insert(e_it->second);
    }
    
    con->start = start;
    con->end = end;
    lazy< std::deque<read_collection> > new_inner = chrom->reads.add_inner();
    con->reads.ref().push_back(new_inner);
    
    // the rest on already initialised
    ++r_it;
    // TODO: I left the code in for non-fixed edges for now, do we need it?
    for (; r_it != region.end(); ++r_it) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Add next Atom " + std::to_string((*r_it)->start) + " " + std::to_string((*r_it)->end) + ".\n");
        #endif
        
        chrom->atoms.push_back(raw_atom());
        raw_atom * at = &chrom->atoms.back();
        at->reference_atom = true;
        at->reference_name = (*r_it)->name;
        at->reference_gene = (*r_it)->gene;

        // we need to do this on tandem;
        greader_list<exon* >::iterator fe_it = con->fossil_exons->begin();
        for (std::set<std::pair<rpos, rpos> >::iterator e_it = (*r_it)->exons.begin(); e_it != (*r_it)->exons.end(); ++e_it ) {
            
            starts.insert(e_it->first);
            ends.insert(e_it->second);
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Exons " +  std::to_string(e_it->first) + " " +  std::to_string(e_it->second)  + "\n");
            #endif
            
            // while we are definitely smaller, advance reference
            while ( fe_it != con->fossil_exons->end() &&  (*fe_it)->end < e_it->first) {
                ++fe_it;
            }
            #ifdef ALLOW_DEBUG
            if (fe_it != con->fossil_exons->end()) logger::Instance()->debug("FE1 " +  std::to_string((*fe_it)->start) + " " +  std::to_string((*fe_it)->end)  + "\n");
            #endif
            if (fe_it == con->fossil_exons->end()) {
            
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Add to End\n");
                #endif
                
                chrom->fossil_exons.push_back(exon(e_it->first, e_it->second)); 
                con->fossil_exons->push_back(&chrom->fossil_exons.back());
                fe_it = --con->fossil_exons->end();
                
                (*fe_it)->fixed_start = true;
                (*fe_it)->fixed_end = true;
                
                at->exons->insert(&chrom->fossil_exons.back());
                
                continue;
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("FE2 " +  std::to_string((*fe_it)->start) + " " +  std::to_string((*fe_it)->end)  + "\n");
            #endif
            if (e_it->second < (*fe_it)->start ) {
                // exon in before next fixed, no overlap 
            
                // add to new exon to the right and stay at new insert
                chrom->fossil_exons.push_back(exon(e_it->first, e_it->second));
                fe_it = con->fossil_exons.ref().insert(fe_it, &chrom->fossil_exons.back());
                
                (*fe_it)->fixed_start = true;
                (*fe_it)->fixed_end = true;
                
                at->exons->insert(&chrom->fossil_exons.back());
                ++fe_it;
                continue;
            }
            
            // new the rest is all overlap
            at->exons->insert(*fe_it); // we have guaranteed overlap already
            if (e_it->first < (*fe_it)->start && e_it->second >= (*fe_it)->start ) { // we have an actual overlap to the front
                // this by the iteration is the first such overlap, therefore no previous to consider
                // hence we extend by the given length

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("FE3 " +  std::to_string((*fe_it)->start) + " " +  std::to_string((*fe_it)->end)  + "\n");
                #endif
                
                if((*fe_it)->fixed_start) { // fix, so add new exon

                   // add to the left then return to current position
                   chrom->fossil_exons.push_back(exon(e_it->first, (*fe_it)->start-1));
                   fe_it = con->fossil_exons.ref().insert(fe_it, &chrom->fossil_exons.back());
                   at->exons->insert(&chrom->fossil_exons.back());
                   
                   (*fe_it)->fixed_start = true;
                   (*fe_it)->fixed_end = true;
                    ++fe_it;

                } else {
                    // just extend existing one
                    (*fe_it)->start = e_it->first;
                }
            }
            
            // we need to split up the original fossil for some cases
            if (e_it->first > (*fe_it)->start && e_it->first < (*fe_it)->end) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("FE4 " +  std::to_string((*fe_it)->start) + " " +  std::to_string((*fe_it)->end)  + "\n");
                #endif
                
                // we add in a new exon
                chrom->fossil_exons.push_back(exon((*fe_it)->start, e_it->first-1));                
                fe_it = con->fossil_exons.ref().insert(fe_it, &chrom->fossil_exons.back());                
                greader_list<exon* >::iterator new_fe = fe_it;
                ++fe_it;
                                
                (*new_fe)->fixed_start = true;
                (*new_fe)->fixed_end = true;
                (*fe_it)->fixed_start = true;
                (*fe_it)->start = e_it->first;
                
                for (greader_refsorted_list<raw_atom* >::iterator m = con->atoms.ref().begin(); m !=  con->atoms.ref().end() ; ++m) {
                    greader_refsorted_list<exon*>::iterator old_left = (*m)->exons.ref().find(*fe_it);
                    if (old_left != (*m)->exons.ref().end()) {
                        (*m)->exons->insert(*new_fe); // insert other instead
                    }
                }
            }
            
            // now loop till end and fix possible holes
            // do we overlap to the end?
            while (fe_it != con->fossil_exons.ref().end() && e_it->second > (*fe_it)->end ) {

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("FE5 " +  std::to_string((*fe_it)->start) + " " +  std::to_string((*fe_it)->end)  + "\n");
                #endif
                
                // extend to right, see for next element
                greader_list<exon* >::iterator next = fe_it;
                ++next;

                if (next == con->fossil_exons.ref().end() || (*next)->start > e_it->second) {
                    // overlap, but no following exon in the overlap
                 
                    if ((*fe_it)->fixed_end) {
                         // add to new exon to the right and stay at new insert
                        chrom->fossil_exons.push_back(exon((*fe_it)->end+1, e_it->second));
                        fe_it = con->fossil_exons.ref().insert(next, &chrom->fossil_exons.back());
                        
                        at->exons->insert(&chrom->fossil_exons.back());
                        (*fe_it)->fixed_start = true;
                        (*fe_it)->fixed_end = true;
                    } else {
                        // modify existing
                        (*fe_it)->end = e_it->second;
                    }
                    
                    at->exons->insert(*fe_it);
                    
                    ++fe_it;                 
                } else {
                                   
                    // this means we have to close the gap between two exons
                    if ( (*fe_it)->fixed_end && (*next)->fixed_start) {
                        // insert new between two 

                        at->exons->insert(*fe_it);
                        at->exons->insert(*next);
                        
                        if ((*fe_it)->end +1 != (*next)->start) {
                            chrom->fossil_exons.push_back(exon((*fe_it)->end+1, (*next)->start-1));
                            fe_it = con->fossil_exons.ref().insert(next, &chrom->fossil_exons.back());
                            at->exons->insert(&chrom->fossil_exons.back());

                            (*fe_it)->fixed_start = true;
                            (*fe_it)->fixed_end = true;
                        }
                        
                        ++fe_it;

                    } else if (!(*fe_it)->fixed_end && (*next)->fixed_start) {
                        // extend and fix left
                        
                        at->exons->insert(*fe_it);
                        at->exons->insert(*next);
                        
                        (*fe_it)->end = (*next)->start-1;
                        (*fe_it)->fixed_end = true;
                        ++fe_it;
                    } else if ( (*fe_it)->fixed_end && !(*next)->fixed_start) {
                        
                        at->exons->insert(*fe_it);
                        at->exons->insert(*next);
                        
                        (*next)->start = (*fe_it)->end+1;
                        (*next)->fixed_start = true;
                        ++fe_it;
                    } else {
                        // we need to merge two separate exons...
                        // get everything in next and erase itr
                        (*next)->start = (*fe_it)->start;
                        at->exons->insert(*next);
                        for (greader_refsorted_list<raw_atom* >::iterator m = con->atoms.ref().begin(); m !=  con->atoms.ref().end() ; ++m) {

                            greader_refsorted_list<exon*>::iterator old_left = (*m)->exons.ref().find(*fe_it);
                            if (old_left != (*m)->exons.ref().end()) {
                                (*m)->exons.ref().erase(*old_left); // erase if found
                                (*m)->exons.ref().insert(*next); // insert other instead
                            }

                        }
                        // now erase fix_it from real list
                        greader_refsafe_list<exon>::iterator rem = std::find(chrom->fossil_exons.begin(), chrom->fossil_exons.end(), **fe_it); 
                        fe_it = con->fossil_exons.ref().erase(fe_it);
                        chrom->fossil_exons.erase(rem);
                    }
                }
            }
            
            
            // we need to split up the original fossil for some cases
            if (fe_it != con->fossil_exons.ref().end() && e_it->second < (*fe_it)->end && e_it->second > (*fe_it)->start) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("FE6 " +  std::to_string((*fe_it)->start) + " " +  std::to_string((*fe_it)->end)  + "\n");
                #endif
                
                // we add in a new exon
                chrom->fossil_exons.push_back(exon((*fe_it)->start, e_it->second));                
                fe_it = con->fossil_exons.ref().insert(fe_it, &chrom->fossil_exons.back());                
                greader_list<exon* >::iterator new_fe = fe_it;
                ++fe_it;
                
                at->exons->insert(*new_fe);
                at->exons->erase(*fe_it);
                
                (*new_fe)->fixed_start = true;
                (*new_fe)->fixed_end = true;
                (*fe_it)->fixed_start = true;
                (*fe_it)->start = e_it->second + 1;
                
                for (greader_refsorted_list<raw_atom* >::iterator m = con->atoms.ref().begin(); m !=  con->atoms.ref().end() ; ++m) {
                    greader_refsorted_list<exon*>::iterator old_left = (*m)->exons.ref().find(*fe_it);
                    if (old_left != (*m)->exons.ref().end()) {
                        (*m)->exons->insert(*new_fe); // insert other instead
                    }
                }
            }
            
        }
        
        con->atoms->insert(at);     
    }

    // add in neighbouring fragments!
    for (greader_refsorted_list<raw_atom* >::iterator ra = con->atoms->begin(); ra !=  con->atoms->end(); ++ra) {

//        logger::Instance()->debug("Test Next Raw Atom Region \n");

        greader_refsorted_list<exon*>::iterator prev = (*ra)->exons->begin();
        greader_refsorted_list<exon*>::iterator next = prev;
        ++next;

        greader_refsorted_list<exon*> collected;

        while (next != (*ra)->exons->end()) {

//            logger::Instance()->debug("RAExon " + std::to_string((*prev)->start) + "-" + std::to_string((*prev)->end) + " : " + std::to_string((*next)->start) + "-" + std::to_string((*next)->end)  +" \n");            

            if ( (*prev)->end + 1 == (*next)->start ) {
                if (collected.empty()) {
                    collected.insert(*prev);
                }
                collected.insert(*next);
//                logger::Instance()->debug("Connected \n");

            } else {
                if (!collected.empty()) {
//                    logger::Instance()->debug("Write out \n");

                    // we add a new minimal evidence atom!
                    chrom->atoms.push_back(raw_atom());
                    raw_atom * at = &chrom->atoms.back();
                    for(greader_refsorted_list<exon*>::iterator ci = collected.begin(); ci != collected.end(); ++ci) {
                        at->exons->insert(*ci);
                    }
                    
                    raw_atom* atom;
                    greader_refsorted_list<raw_atom*>::iterator atom_it = con->atoms.ref().find( at );
                    if (atom_it == con->atoms.ref().end()) {
                        con->atoms->insert(at);
                        atom = at;
                    } else {
                        atom = *atom_it;
                    }
                    atom->has_coverage = true;
                    
                    collected.clear();
                }
            }
            prev = next;
            ++next;
        }
        if (!collected.empty()) {
//            logger::Instance()->debug("Write out \n");

            // we add a new minimal evidence atom!
            chrom->atoms.push_back(raw_atom());
            raw_atom * at = &chrom->atoms.back();
            for(greader_refsorted_list<exon*>::iterator ci = collected.begin(); ci != collected.end(); ++ci) {               
                at->exons->insert(*ci);
            }
            
            raw_atom* atom;
            greader_refsorted_list<raw_atom*>::iterator atom_it = con->atoms.ref().find( at );
            if (atom_it == con->atoms.ref().end()) {
                con->atoms->insert(at);
                atom = at;
            } else {
                atom = *atom_it;
            }
            atom->has_coverage = true;

            collected.clear();
        }
    }

    
    std::copy(starts.begin(), starts.end(), std::back_inserter(chrom->fixed_exon_starts.ref()));
    std::copy(ends.begin(), ends.end(), std::back_inserter(chrom->fixed_exon_ends.ref()));
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("GFF ADDED STARTS: ");
    for (std::set<rpos>::iterator f = starts.begin(); f!= starts.end(); ++f) {
        logger::Instance()->debug(std::to_string(*f) + ",");
    }
    logger::Instance()->debug("\n");
    logger::Instance()->debug("GFF ADDED ENDS: ");
    for (std::set<rpos>::iterator f = ends.begin(); f!= ends.end(); ++f) {
        logger::Instance()->debug(std::to_string(*f) + ",");
    }
    logger::Instance()->debug("\n");
    #endif
    
}


bool gffReader::transcript_info::operator< ( const transcript_info& t2) const {
    
    if (start != t2.start ) {
        return start < t2.start;
    } 
    
    return end < t2.end;
}
