/* 
 * File:   topsort.h
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 *
 * Created on November 24, 2016, 1:02 PM
 */

#ifndef TOPSORT_H
#define	TOPSORT_H

#include <lemon/list_graph.h>
#include <lemon/dfs.h>

namespace pff {

    template <typename Digraph, typename List>
    class TopologicalListVisitor : public lemon::DfsVisitor<Digraph> {
    public:

      typedef typename Digraph::Node Node;
      typedef typename Digraph::Arc edge;

      TopologicalListVisitor(List& order)
        : _order(order) {}

      void leave(const Node& node) {
         _order.push_front(node);
      }

    private:
      List& _order;

    };

    template <typename List, typename NodeMap>
    void order_to_nodemap(List& list, NodeMap& nodemap) {

        unsigned int i = 0;
        for (typename List::iterator it = list.begin(); it != list.end(); ++it, ++i) {
            nodemap.set(*it, i);
        }

    }
    
    
    template <typename Digraph, typename List>
    class TopologicalArcListVisitor : public lemon::DfsVisitor<Digraph> {
    public:

      typedef typename Digraph::Node Node;
      typedef typename Digraph::Arc edge;

      TopologicalArcListVisitor(List& order)
        : _order(order) {}

      void backtrack(const edge& arc) {
         _order.push_front(arc);
      }

    private:
      List& _order;

    };
    
    
}

#endif	/* TOPSORT_H */

