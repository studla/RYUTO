/* 
 * File:   logger.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on January 29, 2016, 11:51 AM
 */

#include "logger.h"
#include "../Options/options.h"

#include <iostream>

logger* logger::instance = NULL;

logger::logger() {
}

logger::~logger() {
    if (instance)
        delete instance;
}

logger* logger::Instance()
{
   if (!instance)   // Only allow one instance of class to be generated.
      instance = new logger;

   return instance;
}

void logger::debug(const std::string& text) {
    if (options::Instance()->is_debug()) {
        std::cout << text;
    }
}

void logger::error(const std::string& text) {
    std::cerr << text;
}

void logger::warning(const std::string& text) {
    std::cout << text;
}

void logger::info(const std::string& text) {
    std::cout << text;
}