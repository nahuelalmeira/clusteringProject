#ifndef GUARD_graphFunctions_hpp
#define GUARD_graphFunctions_hpp

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm> // set_intersection
#include <assert.h>
#include <fstream>
#include <cmath>

#define Edge std::pair<int,int>
#define source first
#define target second


extern long seed;

class Graph {

    public:

        int N;                // Number of nodes
        int M;                // Number of edges
        int n_triangles;      // Number of triangles
        int n_triads;         // Number of triads in the graph
        double C;             // Newman clustering coefficient
        double Cws;           // Watts-Strogatz clustering coefficient

        std::vector<std::set<int> > neighbors;       // Neighbor set per node
        std::vector<int> degree;                    // Degree sequence
        std::vector<std::vector<bool> > adjacency;  // Adjacency matrix
        std::vector<int> triangles;                 // Number of triangles per node
        std::vector<int> triads;                    // Number of triads per node
        std::vector<double> clustering;             // Local clustering coefficient
        std::vector<int> stubs;                     // Each node is repeated a number of
                                                    // times equal to its degree

        /* Construct empty graph */
        Graph() {}

        /* Construct from edge list */
        Graph(std::vector<Edge > edgeList);

        /* Compute number of triads (K2 subgraphs) */
        void computeTriads(); 

        /* Compute number of triangles in the graph */
        void computeTriangles();

        /* Compute Newmann clustering coefficient */
        void computeC();

        /* Compute Watts-Strogatz clustering coefficient */
        void computeCws(); 

        /* Print basic information of graph */
        void summary(std::ostream& out) const; 

        /* Print edge list */
        void printEdgeList(std::ostream& out) const;

        /* Print list of stubs (half edges) */
        void printStubs(std::ostream& out) const;

        /* Print local clustering coefficient */
        void printCi(std::ostream& out) const; 

        /* Writes an edgelist file */
        void writeEdgeList(std::string file_name) const; 

        /* Check if graph is connected */
        bool isConnected() const;

    private:

        /* Auxiliar function for BFS search */
        void DFSUtil(int v, bool visited[]) const;

};






double getCNorm(double C, double Cmax, double Cran);
void printEdge(const Edge e, std::ostream& out);
Edge makeEdge(int a, int b);
bool isEdge(const Graph& G, int a, int b);
void testConsistency(const Graph& G);
Graph createGraph(const std::string file_name, bool verbose=false);
void breakEdge(Graph &G, Edge e, double *delta_C, double *delta_Cws);
void createEdge(Graph &G, Edge e, double *delta_C, double *delta_Cws);
bool areValidEdges(const Graph& G, Edge e1, Edge e2);
void selectPairOfLinks(const Graph& G, Edge *p_e1, Edge *p_e2, long *seed);
int doubleEdgeSwap(Graph &G, Edge e1, Edge e2, double *delta_C, 
                   double *delta_Cws, std::string mode, bool connected, long *seed, double mu);

#endif