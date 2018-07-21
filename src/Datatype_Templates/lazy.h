/* 
 * File:   lazy.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on February 4, 2016, 3:14 PM
 */

#ifndef LAZY_H
#define	LAZY_H


#include <boost/shared_ptr.hpp>

// this class is used to lazy invoke objects
// this is especially useful for copying around lists as references instead of deep copies

template<class T>
    class lazy 
    {
      public:
      
      boost::shared_ptr<T> l;
      
      protected:

        void lazy_int()
        {
          if (!l.use_count()) l.reset(new T());
        }
        
        void copy_in(const lazy &r)
        {
          l = r.l;
        }
      public:
          
        typedef T Type;

        lazy()
        {
        }

        lazy(const lazy &r)
        {
          copy_in(r);
        }

        ~lazy()
        {
          
        }

        lazy &operator=(const lazy &r)
        {
          copy_in(r);
          return *this;
        }

        T *operator->() { lazy_int(); return l.get(); }
        
        T &ref() { lazy_int(); return *l.get(); }

        bool has_ref() {return l.use_count();};

        const T &const_ref() const { assert(l.use_count());  return *l.get(); } 
    };

#endif	/* LAZY_H */

