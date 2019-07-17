#include "../headers/graphFunctions.hpp"
#include "../headers/ran2.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

using std::vector;
using std::set;
using std::string;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::make_pair;
using std::ios;
using std::accumulate;

/* Construct from edge list */
Graph::Graph(vector<Edge > edgeList) {
    int source, target;
    int max_node = 0;
    for(Edge e : edgeList) {
        source = e.source;
        target = e.target;
        if(source > max_node) max_node = source;
        if(target > max_node) max_node = target;            
    }
    N = max_node + 1;
    M = edgeList.size();

    adjacency.resize(N);
    for(int i=0; i<N; i++) adjacency[i].resize(N, false);
    degree.resize(N, 0);
    neighbors.resize(N);
    for(int i=0; i<N; i++) neighbors[i].clear();

    for(Edge e : edgeList) {
        source = e.source;
        target = e.target;
        neighbors[source].insert(target);
        neighbors[target].insert(source);
        adjacency[source][target] = true;
        adjacency[target][source] = true;
        degree[source]++; degree[target]++;
        stubs.push_back(source);
        stubs.push_back(target);
    }
    computeTriads();
    computeTriangles();
    computeC();
    computeCws();
}

/* Compute number of triads (K2 subgraphs) */
void Graph::computeTriads() {
    n_triads = 0;
    int t_i = 0;
    triads.resize(N, 0);
    for(int i=0; i<N; i++) {
        t_i = degree[i] * (degree[i] - 1) / 2;
        triads[i] = t_i;
        n_triads += t_i;
    }
}

/* Compute number of triangles in the graph */
void Graph::computeTriangles() {
    n_triangles = 0;
    triangles.resize(N, 0);
    for(int i=0; i<N; i++) {
        for(int j : neighbors[i])
            for(int k: neighbors[i])
                if(adjacency[j][k]) {
                    triangles[i]++;
                    n_triangles++;
                }
        triangles[i] /= 2;
    }
    n_triangles /= 6;
}

/* Compute Newmann clustering coefficient */
void Graph::computeC() {
    if(n_triads) C = 3. * n_triangles / n_triads;
    else         C = 0;
}

/* Compute Watts-Strogatz clustering coefficient */
void Graph::computeCws() {
    double c_i;
    Cws = 0;
    clustering.resize(N, 0.);
    for(int i=0; i<N; i++) {
        if(degree[i] < 2) c_i = 0;
        else              c_i = float(triangles[i]) / triads[i];
        clustering[i] = c_i;
        Cws += c_i;
    }
    Cws = Cws / N;
}

/* Print basic information of graph */
void Graph::summary(ostream& out) const {
    out << "Graph info"                    << endl;
    out << "N           = " << N           << endl;
    out << "M           = " << M           << endl;
    out << "n_triangles = " << n_triangles << endl;
    out << "n_triads    = " << n_triads    << endl;
    out << "C           = " << C           << endl;
    out << "Cws         = " << Cws         << endl;
}

/* Print edge list */
void Graph::printEdgeList(ostream& out) const {
    for(int a=0; a<N; a++)
        for(int b : neighbors[a])
            out << "(" << a << "," << b << ")" << endl;
}

/* Print list of stubs (half edges) */
void Graph::printStubs(ostream& out) const {
    out << "Stub list: ";
    for(int i : stubs) out << i << " ";
    out << endl;
}

/* Print local clustering coefficient */
void Graph::printCi(ostream& out) const {
    for(int i=0; i<N; i++)
        out << "C[" << i << "] = " << clustering[i] << endl;
}

/* Writes an edgelist file */
void Graph::writeEdgeList(string file_name) const {
    ofstream file;
    file.open(file_name);
    
    for(int source=0; source<N; source++) {
        for(int target : neighbors[source])
            if(source<target) 
                file << source << " " << target << endl;	
    }
}

/* Check if graph is connected */
bool Graph::isConnected() const {
    bool visited[N];
    for(int i=0; i<N; i++)
        visited[i] = false;
    DFSUtil(0, visited);
    for(int i=0; i<N; i++) {
        if(!visited[i]) return false;
    }
    return true;
}

/* Auxiliar function for BFS search */
void Graph::DFSUtil(int v, bool visited[]) const {
    visited[v] = true;
    for(int w : neighbors[v]) {
        if(!visited[w])
            DFSUtil(w, visited);
    }
}


double getCNorm(double C, double Cmax, double Cran) {
    return (C - Cran) / (Cmax - Cran);
}

/* Print an edge as a tuple */
void printEdge(const Edge e, ostream& out) {
    out << "(" << e.source << "," << e.target << ")" << endl;
}

/* Create an edge a---b, where a < b */
Edge makeEdge(int a, int b) {
	if (a > b) return make_pair(b, a);
	return make_pair(a, b);
}

/* Check if there is an edge between nodes a and b */
bool isEdge(const Graph& G, int a, int b) {
    return G.adjacency[a][b];
    }

/* Perform different consistency tests over graph*/
void testConsistency(const Graph& G) {

    /* Size of edgelist must equal property M */
    //assert(G.edges.size() == G.M);

    /* Stubs list size must be twice the number of edges */
    assert(G.stubs.size() == 2*G.M);

    /* Degree sequence must be constant */
    for(int i=0; i<G.N; i++)
        assert(G.degree[i] == G.neighbors[i].size());
    

    /* Graph cannnot have isolated nodes */
    for(vector<set<int > >::const_iterator it = G.neighbors.begin();
        it != G.neighbors.end(); ++it)
        assert(it->size());

    /* Total degree must be equal to 2*M */
    int total_degree = 0;
    for(vector<set<int > >::const_iterator it = G.neighbors.begin();
        it != G.neighbors.end(); ++it)
        total_degree += it->size();
    assert(total_degree == 2*G.M);

}

/* Create a Graph from an adjacency list file.
   Nodes must be numbered from 0 to N-1 */
Graph createGraph(const string file_name, bool verbose) {

    ifstream file;
    
    Graph G;
    int a, b;
    int N = 0; // Number of nodes
    int M = 0; // Number of edges

    
    /* First loop to get N and M */
    file.open(file_name);
    if(!file.is_open()) {
        cerr << "Could not open file '" << file_name << "'" << endl;
        assert(0);
    }

    if(verbose) cout << "Reading file '"<< file_name << "'" << endl;

    while(true) {
        file >> a; file >> b;
        if(file.eof()) break;
        if(a==b) continue;

        if(a>N) N = a;
        if(b>N) N = b;
        M++;
    }
    file.close();
    N++;
    
    if(verbose) cout << "Initializing data structures" << endl;
    /* Initialize data structures */
    G.N = N;
    G.M = M;
    G.adjacency.resize(N);
    for (int i=0; i<G.adjacency.size(); i++)
		G.adjacency[i].resize(N,false);
    G.degree.resize(N, 0);
    G.neighbors.resize(N);
	for (int i=0; i < N; i++)
		G.neighbors[i].clear();

    if(verbose) cout << "Adding data to graph" << endl;
    /* Second loop to add data to graph */
    file.open(file_name, ios::in);
    while(true) {
        file >> a >> b;
        if(file.eof()) break;
        if(a==b) continue;

        G.adjacency[a][b] = true; G.adjacency[b][a] = true;
        G.degree[a]++; G.degree[b]++;
        G.neighbors[a].insert(b); G.neighbors[b].insert(a);
        G.stubs.push_back(a);
        G.stubs.push_back(b);
    }
    file.close();

    if(verbose) cout << "Computing triads" << endl;
    G.computeTriads();

    if(verbose) cout << "Computing triangles" << endl;
    G.computeTriangles();

    if(verbose) cout << "Computing clustering coefficients" << endl;
    G.computeCws(); 
    G.computeC(); 

    testConsistency(G);

    return G;
}

/* Break an edge of the graph */
void breakEdge(Graph &G, Edge e, double *delta_C, double *delta_Cws) {
	int a, b, broken_triangles;
    double delta_c_i;
	a = e.source; b = e.target;

	/* Remove edge from all data structures */
	G.neighbors[a].erase(b); G.neighbors[b].erase(a);
	G.adjacency[a][b] = false; G.adjacency[b][a] = false;
	
	/* Count number of broken triangles */
	set<int> common_neighbors;
	set_intersection(G.neighbors[a].begin(), G.neighbors[a].end(),
	                 G.neighbors[b].begin(), G.neighbors[b].end(),
                     inserter(common_neighbors, common_neighbors.begin()));
	broken_triangles = common_neighbors.size();
    
    /* Update quantities */
    G.n_triangles -= broken_triangles;

    G.triangles[a] -= broken_triangles;
    G.triangles[b] -= broken_triangles;

    *delta_Cws = 0.;

    if (G.degree[a] > 1) {
        delta_c_i = 1. * broken_triangles / G.triads[a];
        G.clustering[a] -= delta_c_i;
        *delta_Cws -= delta_c_i / G.N;
    }
    if (G.degree[b] > 1) {
        delta_c_i = 1. * broken_triangles / G.triads[b];
        G.clustering[b] -= delta_c_i;
        *delta_Cws -= delta_c_i / G.N;
    }

    /* Triangles centered in common neighbors */
    for(int c : common_neighbors) {
        G.triangles[c]--;
        delta_c_i = 1. / G.triads[c];
        G.clustering[c] -= delta_c_i;
        *delta_Cws -= delta_c_i / G.N;
    }

    *delta_C = -3. * broken_triangles / G.n_triads;
    //G.C += *delta_C;
    //G.Cws += *delta_Cws;
}

/* Creates an edge in the graph */
void createEdge(Graph &G, Edge e, double *delta_C, double *delta_Cws) {
	int a, b, new_triangles;
    double delta_c_i;
	a = e.source; b = e.target;

	/* Add edge to all data structures */
	G.neighbors[a].insert(b); G.neighbors[b].insert(a);
	G.adjacency[a][b] = true; G.adjacency[b][a] = true;
	
	/* Count number of broken triangles */
	set<int> common_neighbors;
	set_intersection(G.neighbors[a].begin(), G.neighbors[a].end(),
	                 G.neighbors[b].begin(), G.neighbors[b].end(),
                     inserter(common_neighbors, common_neighbors.begin()));
	new_triangles = common_neighbors.size();
    
    /* Update quantities */
    G.n_triangles += new_triangles;

    G.triangles[a] += new_triangles;
    G.triangles[b] += new_triangles;

    *delta_Cws = 0.;
    if (G.degree[a] > 1) {
        delta_c_i = 1. * new_triangles / G.triads[a];
        G.clustering[a] += delta_c_i;
        *delta_Cws += delta_c_i / G.N;
    }
    if (G.degree[b] > 1) {
        delta_c_i = 1. * new_triangles / G.triads[b];
        G.clustering[b] += delta_c_i;
        *delta_Cws += delta_c_i / G.N;
    }

    /* Triangles centered in common neighbors */
    for(int c : common_neighbors) {
        G.triangles[c]++;
        delta_c_i = 1. / G.triads[c];
        G.clustering[c] += delta_c_i;
        *delta_Cws += delta_c_i / G.N;
    }

    *delta_C = 3. * new_triangles / G.n_triads;
    //G.C += *delta_C;
    //G.Cws += *delta_Cws;
}

/* Check that edges e1 and e2 are not the same or adjacent */
bool areValidEdges(const Graph& G, Edge e1, Edge e2) {
	int a = e1.source, b = e1.target;
	int c = e2.source, d = e2.target;

	if(e1 == e2) return false;
	
	if(a == b) return false;
	if(a == c) return false;
	if(a == d) return false;
	if(b == c) return false;
	if(b == d) return false;
	if(c == d) return false;
	
	//if(isLink(G, a, c) || isLink(G, b, d)) return false;
	
	return true;
}

/* Selects a pair of non-adjacent links with uniform probability 
   Each link is selected in the following way: first, a node is 
   randomly selected with probability proportional to its degree. 
   Then, a neighbor of the node is randomly selected with uniform 
   probability among the neighbors.*/
void selectPairOfLinks(const Graph& G, Edge *p_e1, Edge *p_e2, long *seed) {

	Edge e1, e2;
	int source, target, idx;
	set<int >::iterator it; 
	set<int > neighbors;
    int M = G.M;
	do {
		idx = floor(ran2(seed)*2*M);
		source = G.stubs[idx];
		neighbors = G.neighbors[source];
		it = neighbors.begin();
		idx = floor(ran2(seed)*neighbors.size());
		advance(it, idx);
		target = *it;
		e1 = makeEdge(source, target);
	    //assert(G.edges.find(e1) != G.edges.end());

		idx = floor(ran2(seed)*2*M);
		source = G.stubs[idx];
		neighbors = G.neighbors[source];
		it = neighbors.begin();
		idx = floor(ran2(seed)*neighbors.size());
		advance(it, idx);
		target = *it;
		e2 = makeEdge(source, target);
		//assert(G.edges.find(e2) != G.edges.end());
			
	} while(!areValidEdges(G, e1, e2));
	
	//assert(e1.source!=e1.target);
	//assert(e2.source!=e2.target);
		
	*p_e1 = e1;
	*p_e2 = e2;
}

int doubleEdgeSwap(Graph &G, Edge e1, Edge e2, 
                   double *delta_C, double *delta_Cws, 
                   string mode, bool connected, long *seed, double mu) {

    Edge f1, f2;
    int valid_swap = 0;
    double sub_delta_C = 0;
    double sub_delta_Cws = 0;
    double delta, prob;
    int reverse_swap = false;

    *delta_C = 0;
    *delta_Cws = 0;

    /* This part is necesary for uniformity in the rewiring preocess.
       With equal probability it tries to connect (source, source), 
       (target, target) and (source, target), (target, source)  */
    if(ran2(seed)<0.5) {
        if( !(isEdge(G, e1.source, e2.source) || isEdge(G, e1.target, e2.target)) ) {
            f1 = makeEdge(e1.source, e2.source);
            f2 = makeEdge(e1.target, e2.target);
            valid_swap = 1;
        }
    }
    else {
        if( !(isEdge(G, e1.source, e2.target) || isEdge(G, e1.target, e2.source)) ) {
            f1 = makeEdge(e1.source, e2.target);
            f2 = makeEdge(e1.target, e2.source);
            valid_swap = 1;
        }
    }
    if(!valid_swap) return 0;

    /* Perform double-edge swap */
    if(valid_swap) {
        breakEdge(G, e1, &sub_delta_C, &sub_delta_Cws);
        *delta_C += sub_delta_C;
        *delta_Cws += sub_delta_Cws;

        breakEdge(G, e2, &sub_delta_C, &sub_delta_Cws);
        *delta_C += sub_delta_C;
        *delta_Cws += sub_delta_Cws;

        createEdge(G, f1, &sub_delta_C, &sub_delta_Cws);
        *delta_C += sub_delta_C;
        *delta_Cws += sub_delta_Cws;

        createEdge(G, f2, &sub_delta_C, &sub_delta_Cws);
        *delta_C += sub_delta_C;
        *delta_Cws += sub_delta_Cws;
    } 

    /* Check wether to reverse swap or not */
    if(mode != "random") {
        if(mode == "maxC")        delta = G.n_triads * (*delta_C) / 3.;
        else if(mode == "maxCws") delta = G.n_triads * (*delta_Cws) / 3.;

        prob = ran2(seed);
        if(exp(mu*delta) < (1 - prob)) {
            reverse_swap = true;

        }           
    }

    if(!reverse_swap && connected) {
        if(!G.isConnected()) reverse_swap = true;
    }

    /* Perform reverse double-edge swap if necesary */
    if(reverse_swap) {
        breakEdge(G, f1, &sub_delta_C, &sub_delta_Cws);
        *delta_C += sub_delta_C;
        *delta_Cws += sub_delta_Cws;

        breakEdge(G, f2, &sub_delta_C, &sub_delta_Cws);
        *delta_C += sub_delta_C;
        *delta_Cws += sub_delta_Cws;

        createEdge(G, e1, &sub_delta_C, &sub_delta_Cws);
        *delta_C += sub_delta_C;
        *delta_Cws += sub_delta_Cws;

        createEdge(G, e2, &sub_delta_C, &sub_delta_Cws);
        *delta_C += sub_delta_C;
        *delta_Cws += sub_delta_Cws;

        /* After reverse swap, C and Cws must remain constant */
        assert(fabs(*delta_C) < 10E-8);
        assert(fabs(*delta_Cws) < 10E-8);
    } 

    return valid_swap;
}