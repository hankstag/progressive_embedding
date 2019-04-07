#include "mst.h"
#include <igl/adjacency_list.h>
#include <igl/adjacency_matrix.h>
#include <stdio.h>

// geneate a minimum spanning tree for the mesh [V,F], taking s as root
// use prim algorithm [https://www.geeksforgeeks.org/prims-minimum-spanning-tree-mst-greedy-algo-5/]

int minKey(int key[], bool mstSet[], int nV) 
{ 
    // Initialize min value 
    int min = INT_MAX, min_index; 
    
    for (int v = 0; v < nV; v++) 
        if (mstSet[v] == false && key[v] < min) 
            min = key[v], min_index = v; 
    
    return min_index; 
}

void mst(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    int s,
    std::vector<int>& T,
    std::set<int>& mst_set,
    Eigen::SparseMatrix<int>& graph
){
    
    int nV = V.rows();
    int parent[nV];
    int key[nV];
    bool mstSet[nV];  
  
    // Initialize all keys as INFINITE 
    for (int i = 0; i < nV; i++) 
        key[i] = INT_MAX, mstSet[i] = false; 
  
    // Always include first 1st vertex in MST. 
    // Make key 0 so that this vertex is picked as first vertex. 
    key[0] = 0;      
    parent[0] = -1; // First node is always root of MST  
  
    // The MST will have nV vertices 
    for (int count = 0; count < nV-1; count++) 
    { 
        // Pick the minimum key vertex from the  
        // set of vertices not yet included in MST 
        int u = minKey(key, mstSet, nV); 
  
        // Add the picked vertex to the MST Set 
        mstSet[u] = true; 
  
        // Update key value and parent index of  
        // the adjacent vertices of the picked vertex.  
        // Consider only those vertices which are not  
        // yet included in MST 
        for (int v = 0; v < nV; v++) 
  
        // graph(u,v) is non zero only for adjacent vertices of m 
        // mstSet[v] is false for vertices not yet included in MST 
        // Update the key only if graph[u][v] is smaller than key[v] 
        if (graph.coeff(u,v) && mstSet[v] == false && graph.coeff(u,v) < key[v]) 
            parent[v] = u, key[v] = graph.coeff(u,v); 
    } 
    T.resize(nV);
    for(int i=0;i<nV;i++)
        T[i] = parent[i];
    //printMST(parent, nV, graph); 
}