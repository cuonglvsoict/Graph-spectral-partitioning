#ifndef _GRAPH_HEADER_
#define _GRAPH_HEADER_

#include<iostream>
#include<vector>

using namespace std;

typedef vector<vector<double>> matrix;

class Graph
{
public:
	int n;									// # vertices
	int m;									// # edges
	matrix adjacency;
	matrix laplacian;
	matrix degree;							// diagonal matrix, degree[v][v] = degree(v)
	
	vector<double> laplacian_spectrum;
	matrix laplacian_eigenvectors;

	vector<int> A;							// 1st set of vertices
	vector<int> B;							// 2nd set of vertices
	vector<int> S;							// vertex separator
	vector<pair<int, int>> edge_separator;

	// auxiliary sets
	vector<int> A_prime;					// initial 1st set
	vector<int> B_prime;					// initial 2nd set
	

	int load_adjacency_matrix(string file_path);
	void calculate_laplacian_matrix();
	void calculate_laplacian_spectrum();
	void perform_initial_partition();
	void calculate_vertex_separator();
	void calculate_partitions();
};


#endif // !_GRAPH_HEADER_