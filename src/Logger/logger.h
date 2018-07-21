/* 
 * File:   logger.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on January 29, 2016, 11:51 AM
 */

#ifndef LOGGER_H
#define	LOGGER_H

#include <string>

class logger {
public:
    
    static logger* Instance();
    void debug(const std::string &text);
    void error(const std::string &text);
    void warning(const std::string &text);
    void info(const std::string& text);
    
private:
    logger();
    logger(const logger& orig) {};
    virtual ~logger() ;
    logger& operator=(logger const&){};
    
    static logger* instance;

};

#endif	/* LOGGER_H */

