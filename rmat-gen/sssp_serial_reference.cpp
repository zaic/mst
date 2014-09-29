#include <iostream>
#include <vector>
#include <algorithm>
#include "defs.h"

using namespace std;

/* Reference SSSP implementation --- Dijkstra algorithm 
 * it is needed to carefully count the number of traversed edges
 * dist must be initialized as -1
*/
extern "C" void sssp(vertex_id_t root, graph_t *G, weight_t *dist, uint64_t *traversed_edges)
{
    /* (distance,node) */
    vector< pair<weight_t, vertex_id_t> > pq(G->n);
    uint64_t nedges = 0;
    pq.at(0) = make_pair(0, root); 
    int len = 1; 
    dist[root] = 0;
    while ( len ) {
        pair<weight_t, vertex_id_t> s = pq[0];
        pop_heap(pq.begin(), pq.begin()+len, greater<pair<weight_t, vertex_id_t> >()); 
        len--;
        if ( (dist[s.second] >= s.first) || (dist[s.second] == -1) ) {
            for ( unsigned int i = G->rowsIndices[s.second]; i < G->rowsIndices[s.second+1]; ++i ) {
                if ( (dist[G->endV[i]] > dist[s.second] + G->weights[i]) || (dist[G->endV[i]] == -1) ) {
                    dist[G->endV[i]] = dist[s.second] + G->weights[i];
                    pq.at(len) = make_pair(dist[G->endV[i]], G->endV[i]); 
                    len++;
                    if ((unsigned)len+1 > pq.size()) { pq.resize(2*len); }
                    push_heap(pq.begin(), pq.begin()+len, greater<pair<weight_t, vertex_id_t> >());
                }
	            /* we count every traversed edge */
	            ++nedges;
            }
        }
    }
    *traversed_edges = nedges;
}

