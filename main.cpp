#include"Graph.h"

int main() {

	Graph graph;
	graph.load_adjacency_matrix("data/small_test.mtx");
	// graph.load_adjacency_matrix("data/bcspwr09.mtx");

	graph.calculate_laplacian_matrix();
	graph.calculate_laplacian_spectrum();

	// display the eigenvalues of the laplacian matrix
	/*
	for (int i = 0; i < graph.n; i++) {
		cout << graph.laplacian_spectrum[i] << "\t";
	}
	cout << endl;
	
	// dislay the second eigenvector
	for (int i = 0; i < graph.n; i++) {
		cout << graph.laplacian_eigenvectors[1][i] << "\t";
	}
	cout << endl;
	*/

	graph.perform_initial_partition();
	
	// display edge separator
	/*
	cout << "Intial partitions' sizes: " << graph.A_prime.size() << "\t" << graph.B_prime.size() << endl;
	cout << "Size of edge separator: " << graph.edge_separator.size() << endl;

	for (pair<int, int> edge : graph.edge_separator) {
		cout << edge.first << "\t" << edge.second << endl;
	}

	for (int i : graph.A_prime) {
		cout << i << ", ";
	}
	cout << endl;

	for (int i : graph.B_prime) {
		cout << i << ", ";
	}
	cout << endl;
	*/

	graph.calculate_vertex_separator();
	graph.calculate_partitions();

	cout << "-----------------------------------------" << endl;
	cout << endl << "Results:" << endl;
	cout << "Size of 1st partition: " << graph.A.size() << endl;
	cout << "Size of 2nd partition: " << graph.B.size() << endl;
	cout << "Size of vertex separator: " << graph.S.size() << endl;
	cout << "Size of edge separator: " << graph.edge_separator.size() << endl;

	/*
	cout << "Vertex separator: ";
	for (int u : graph.S) {
		cout << u << " ";
	}
	cout << endl;

	cout << "Edge separator: " << endl;
	for (pair<int, int> edge : graph.edge_separator) {
		cout << "\t(" << edge.first << " - " << edge.second << ")" << endl;
	}
	*/

	return 0;
}