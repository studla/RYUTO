/* 
 * File:   capacity_mean.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 23, 2018, 2:39 PM
 */

#ifndef CAPACITY_MEAN_H
#define	CAPACITY_MEAN_H

#include <iostream>
#include <deque>
#include "../../Datatype_Templates/misc_types.h"

class capacity_mean {
public:
        capacity_mean();
        capacity_mean(float m, rpos w);
        
        void update(capacity_mean o);
        void reduce(float percentage);
        
        float compute_score();
        void reduce_score();
        std::string to_string();
        
        float mean;
        std::deque<float> scores;
        float hidden_score;
        rpos weight;
                
private:

};

std::ostream& operator<<(std::ostream& os, const capacity_mean& dt);


#endif	/* CAPACITY_MEAN_H */

