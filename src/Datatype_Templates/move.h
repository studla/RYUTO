/* 
 * File:   move.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on December 2, 2015, 6:29 PM
 */

#ifndef MOVE_H
#define	MOVE_H

#if __cplusplus >= 201103L
    #include <utility>
    #include <iterator>

    #define _MOVE(__val) std::move(__val)
    #define _MOVE_RANGE(__it1, __it2, __res ) std::move(__it1, __it2, __res )
    #define _MOVE_IT(__val) std::make_move_iterator(__val)
#else
    #define _MOVE(__val) (__val)
    #define _MOVE_RANGE(__it1, __it2, __res ) std::copy(__it1, __it2, __res )
    #define _MOVE_IT(__val) (__val)
#endif

#endif	/* MOVE_H */

