/* 
 * File:   double_deque.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on April 15, 2016, 1:52 PM
 */

#ifndef DOUBLE_DEQUE_H
#define	DOUBLE_DEQUE_H

#include <deque>
#include <list>
#include "lazy.h"

template <class T>
class double_deque : public std::list< lazy< std::deque<T> > > {
    public:
        
        void push_to_end(const T& el) {
            this->back().ref().push_back(el);
        }
        
        T & get_end() {
            this->back().ref().back();
        }
        
        lazy< std::deque<T> > add_inner() {
            this->push_back( lazy< std::deque<T> >() ); 
            return this->back();
        }
        
        void remove_list(lazy< std::deque<T> > el) {
            this->erase(std::find(this->begin(), this->end(), el));
        }
};


template <class T>
class double_deque_ref : public std::deque< lazy< std::deque<T> > >{
public:
    
    // TODO: this class will fail for empty sublists!
    class full_iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
        
        friend class double_deque_ref;
        
    protected:
    
        double_deque_ref<T>* container;
        typename std::deque< lazy< std::deque<T> > >::iterator outer;
        typename std::deque<T>::iterator inner;
        
    public:
        
        full_iterator(typename std::deque< lazy< std::deque<T> > >::iterator &outer, typename std::deque<T>::iterator &inner, const double_deque_ref<T> &c) : outer(outer), inner(inner), container(&c) {
            
        }
        
        bool operator==(full_iterator const& rhs) const 
        {
            return (inner==rhs.inner);
        }

        bool operator!=(full_iterator const& rhs) const 
        {
            return !(*this==rhs);
        }

        full_iterator& operator++() 
        {
            
            ++inner;
            if (inner == (*outer)->end()) {    
                ++outer;
                if (outer != container->end())
                    inner = (*outer)->begin();
            }
            return *this;
        }   

        full_iterator operator++(T) 
        {
            full_iterator tmp (outer, inner, *container);
            ++(*this);
            return tmp;
        }

        full_iterator& operator--() 
        {
            if (inner == (*outer)->begin()) {
                --outer;
                inner = (*outer)->end();
            }
            --inner;
            
            return *this;
        }

        full_iterator operator--(T) 
        {
            full_iterator tmp (outer, inner, *container);
            --(*this);
            return tmp;
        }

        T& operator* () const 
        {
            return *inner;
        }
    };
    
    
    full_iterator dbegin() {
        return full_iterator(this->begin(), this->front()->begin(), this);
    }
    full_iterator dend() {
        return full_iterator(this->end(), this->back()->end(), this);
    }
    
    void push_to_end(const T& el) {
        this->back().ref().push_back(el);
    }
    T & get_end() {
        return this->back().ref().back();
    }
    
private:
    
    
    
};

#endif	/* DOUBLE_DEQUE_H */

