/* 
 * File:   capacity_mean.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on February 23, 2018, 2:39 PM
 */

#include "capacity_mean.h"
#include "../../../Logger/logger.h"
#include <math.h>
#include <cmath>
#include <deque>

capacity_mean::capacity_mean() : mean(0), hidden_score(0), weight(0) {
}

capacity_mean::capacity_mean(float m, rpos w) : mean(m), hidden_score(0), weight(w) {
    scores.push_back(m);
}


void capacity_mean::update(capacity_mean &o) {
    
   // logger::Instance()->debug("Update " + std::to_string(o.mean) + "/" + std::to_string(mean) + " " + std::to_string(o.weight) + "/" + std::to_string(weight) + "\n");
    
    mean = mean + o.weight/(float) (weight+o.weight) * (o.mean -  mean);
    weight += o.weight;

    std::copy(std::begin(o.scores), std::end(o.scores), std::back_inserter(scores));
}

void capacity_mean::sum_up(capacity_mean& rhs, unsigned int i) {
    mean = mean + 1/(float) (i) * (rhs.mean -  mean);
}


void capacity_mean::reduce(float percentage) {
    //logger::Instance()->debug("Reduce by " + std::to_string(percentage) + "\n");
    
    if (std::isnan(percentage) || percentage < 0.0001) {
        percentage = 0.001;
    }
    
    mean = mean * percentage;
    if (mean < 0.01 ) {
        mean = 0.01;
    }
    for (std::deque<float>::iterator si = scores.begin(); si != scores.end(); ++si) {
        *si = *si * percentage;
        if (*si < 0.01 ) {
            *si = 0.01;
        }
    }
}

void capacity_mean::assign_mean(float m) {

    if (std::isnan(m) || m < 0.01 ) {
        m = 0.01;
    }
    mean = m;
}

float capacity_mean::compute_score() {
    
//    logger::Instance()->debug("Compute Score \n");
    
    if (scores.empty()) {
        return hidden_score;
    }
 
    std::deque<float> sc = scores;
    
    while(sc.size() > 1) {
        
//        logger::Instance()->debug("SC ");
//        for(std::deque<float>::iterator it = sc.begin(); it != sc.end(); ++it) {
//            logger::Instance()->debug(std::to_string(*it) + ", ");
//        }
//        logger::Instance()->debug("\n");
        
        float ratio = 0;
        std::deque<float>::iterator p1, p2;
        
        std::deque<float>::iterator si1 = sc.begin();
        std::deque<float>::iterator si2 = si1;
        ++si2;
        for (; si2 != sc.end(); ++si1, ++si2) {
            float w1 = *si1;
            float w2 = *si2;
            
            float r = 0;    
            if(w1 > w2) {
                r = w1/w2;
            } else {
                r = w2/w1;
            }
            
            if (r < ratio || si1 == sc.begin()) {
                ratio = r;
                p1 = si1;
                p2 = si2;
            }
        }
        
//        logger::Instance()->debug("Combine " + std::to_string(*p1) + " " + std::to_string(*p2) + "\n");
        
        *p1 = sqrt(*p1 * *p2);
        sc.erase(p2);
    }
    return sc[0];
}

void capacity_mean::reduce_score() {
    float sc = compute_score();
    scores.clear();
    scores.push_back(sc);
}

std::string capacity_mean::to_string() {
    std::string ret = "Scores ";
    for(std::deque<float>::iterator si = scores.begin(); si != scores.end(); ++si) {
        ret+= std::to_string(*si) + " ";
    }
    return ret;
}


std::ostream& operator<<(std::ostream& os, const capacity_mean& cm)
{
    capacity_mean cmc = cm; // copy over scores cause they are destroyed
    os << std::to_string(cm.mean) << "," << std::to_string(cm.weight) << " : " << std::to_string(cmc.compute_score()) << "(" << cm.scores.size() << ")";
    return os;
}
