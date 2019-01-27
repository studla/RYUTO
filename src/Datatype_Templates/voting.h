/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   voting.h
 * Author: thomas
 *
 * Created on January 3, 2019, 11:12 AM
 */

#ifndef VOTING_H
#define VOTING_H

template<class InputIt>
int most_common_mapped(InputIt begin, InputIt end)
{
    std::unordered_map<int, int> counts;
    for (InputIt it = begin; it != end; ++it) {
            ++counts[it->second];
    }
    return std::max_element(counts.begin(), counts.end(),
            [] (const std::pair<int, int>& pair1, const std::pair<int, int>& pair2) {
            return pair1.second < pair2.second;})->first;
}

template<class InputIt, class Set>
void most_common_multi(InputIt begin, InputIt end, Set &ret)
{
    std::unordered_map<int, int> counts;
    for (InputIt it = begin; it != end; ++it) {
            ++counts[*it];
    }
    int max = std::max_element(counts.begin(), counts.end(),
            [] (const std::pair<int, int>& pair1, const std::pair<int, int>& pair2) {
            return pair1.second < pair2.second;})->second;
    for (std::unordered_map<int, int>::iterator ci = counts.begin(); ci != counts.end(); ++ci) {
        if (ci->second == max) {
            ret.insert(ci->first);
        }
    }        
}

#endif /* VOTING_H */

