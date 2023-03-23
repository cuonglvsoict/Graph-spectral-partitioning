#include<fstream>
#include<string>
#include<set>
#include<algorithm>

#include"Graph.h"

// function for the spectral decomposition
void tred2(int n, matrix& V, vector<double>& d, vector<double>& e);
void tql2(int n, vector<double> &d, vector<double> &e, matrix &V);

int Graph::load_adjacency_matrix(string file_path) {
	fstream fin(file_path);

	if (!fin.is_open()) {
		cout << "Failed to open the coordinate file..." << endl;
		return 0;
	}

	string file_comment;
	getline(fin, file_comment);

	fin >> n >> n >> m;

	adjacency.clear();
	for (int i = 0; i < n; i++) {
		vector<double> row_i(n, 0);
		adjacency.push_back(row_i);
	}

	int u, v;
	int count_edge = 0;
	for (int i = 0; i < m; i++) {
		fin >> u >> v;

		if (u == v) continue;

		count_edge++;
		adjacency[u-1][v-1] = 1;
		adjacency[v-1][u-1] = 1;
	}
	m = count_edge;

	fin.close();

	cout << "Input graph loaded successfully..." << endl;
	cout << "\tNumber of nodes: " << n << endl;
	cout << "\tNumber of edges: " << count_edge << endl;

	return 1;
}

void Graph::calculate_laplacian_matrix() {
	// calculate degree
	degree.clear();
	for (int i = 0; i < n; i++) {
		vector<double> d_i(n, 0);

		for (int j = 0; j < n; j++) {
			d_i[i] += adjacency[i][j];
		}

		degree.push_back(d_i);
	}

	laplacian.clear();
	for (int i = 0; i < n; i++) {
		vector<double> lap_i(n, 0);

		for (int j=0; j<n; j++) {
			lap_i[j] = degree[i][j] - adjacency[i][j];
		}

		laplacian.push_back(lap_i);
	}

	cout << "Laplacian matrix calculated..." << endl;
}

void Graph::calculate_laplacian_spectrum() {

	// Tridiagonalize
	laplacian_eigenvectors.clear();
	for (int i = 0; i < n; i++) {
		vector<double> eigenvector_i(n, 0);
		laplacian_eigenvectors.push_back(eigenvector_i);

		for (int j = 0; j < n; j++) {
			laplacian_eigenvectors[i][j] = laplacian[i][j];
		}
	}

	vector<double> offdiag(n, 0);
	laplacian_spectrum = vector<double>(n, 0);

	tred2(n, laplacian_eigenvectors, laplacian_spectrum, offdiag);

	// Diagonalize
	tql2(n, laplacian_spectrum, offdiag, laplacian_eigenvectors);

	cout << "Laplacian matrix's spectrum calculated..." << endl;
}

void Graph::perform_initial_partition() {
	// vector<double> sorted_2nd_eig(laplacian_eigenvectors[1].begin(), laplacian_eigenvectors[1].end());
	vector<double> sorted_2nd_eig;
	for (int i = 0; i < n; i++) {
		sorted_2nd_eig.push_back(laplacian_eigenvectors[i][1]);
	}
	sort(sorted_2nd_eig.begin(), sorted_2nd_eig.end());
	double median = sorted_2nd_eig[n / 2];

	A_prime.clear();
	B_prime.clear();

	vector<int> tie;
	for (int i = 0; i < n; i++) {
		if (laplacian_eigenvectors[i][1] < median) {
			A_prime.push_back(i);
		}
		else if (laplacian_eigenvectors[i][1] > median) {
			B_prime.push_back(i);
		}
		else {
			tie.push_back(i);
		}
	}

	if (tie.size() == 1) {
		if (A_prime.size() <= B_prime.size()) {
			A_prime.push_back(tie[0]);
		}
		else {
			B_prime.push_back(tie[0]);
		}
	}
	else if (tie.size() > 1) {
		int k = tie.size() / 2;
		A_prime.insert(A_prime.end(), tie.begin(), tie.begin() + k - 1);
		B_prime.insert(B_prime.end(), tie.begin() + k, tie.end());
	}

	// extract edge separator
	edge_separator.clear();
	for (int u : A_prime) {
		for (int v : B_prime) {
			if (adjacency[u][v]) {
				edge_separator.push_back(make_pair(u, v));
			}
		}
	}

	cout << "Intial partition completed..." << endl;
}

int dfs(matrix& H, int left_node, vector<int>& right_nodes, vector<int>& left_matching, vector<int>& mark) {
	for (int v : right_nodes) {
		if (H[left_node][v] && !mark[v]) {
			mark[v] = 1;

			if (left_matching[v] < 0 || dfs(H, left_matching[v], right_nodes, left_matching, mark)) {
				left_matching[v] = left_node;
				return 1;
			}
		}
	}

	return 0;
}

void trace_alternating_paths(vector<int> &match_with_A1, vector<int> &match_with_B1, matrix& H, int u, set<int>& Z, vector<int> &mark) {
	if (Z.count(u) > 0) {
		return;
	}

	Z.insert(u);
	mark[u] = 1;

	for (int v = 0; v < H.size(); v++) {
		if (H[u][v] && match_with_B1[u] != v && !mark[v]) {
			Z.insert(v);
			mark[v] = 1;

			trace_alternating_paths(match_with_A1, match_with_B1, H, match_with_A1[v], Z, mark);
		}
	}
}

void Graph::calculate_vertex_separator() {
	// calcualte the bipartite graph H = (A1, B1, edge_separator)
	set<int> tmpA, tmpB;
	for (pair<int, int> edge : edge_separator) {
		tmpA.insert(edge.first);
		tmpB.insert(edge.second);
	}
	vector<int> A1(tmpA.begin(), tmpA.end());
	vector<int> B1(tmpB.begin(), tmpB.end());

	matrix H;
	for (int i = 0; i < n; i++) {
		vector<double> row_i(n, 0);
		H.push_back(row_i);
	}

	for (pair<int, int> edge : edge_separator) {
		H[edge.first][edge.second] = 1;
	}

	// find the maximum matching of H
	vector<int> match_with_B1(n, -1); // i in B1, matching[i] in A1
	for (int u : A1) {
		vector<int> mark(H.size(), 0);
		dfs(H, u, B1, match_with_B1, mark);
	}

	vector<int> match_with_A1(n, -1);
	for (int i = 0; i < n; i++) {
		if (match_with_B1[i] >= 0) {
			match_with_A1[match_with_B1[i]] = i;
		}
	}

	// calculate the minimum vertex cover from the maximum matching
	set<int> Z;
	vector<int> mark(n, 0);
	for (pair<int, int> edge : edge_separator) {
		H[edge.second][edge.first] = 1;
	}

	for (int u : B1) {
		if (match_with_B1[u] < 0) {
			trace_alternating_paths(match_with_A1, match_with_B1, H, u, Z, mark);
		}
	}

	// S = (B1 \ Z) U (A1 intersect Z)
	S.clear();
	for (int u : B1) {
		if (Z.count(u) < 1) {
			S.push_back(u);
		}
	}

	for (int u : Z) {
		if (tmpA.count(u) > 0) {
			S.push_back(u);
		}
	}

	cout << "Vertex separator calculated... " << endl;
}

void Graph::calculate_partitions() {
	vector<int> mark(n, 0);
	for (int u : S) {
		mark[u] = 1;
	}

	A.clear();
	B.clear();
	for (int u : A_prime) {
		if (!mark[u]) {
			A.push_back(u);
		}
	}

	for (int u : B_prime) {
		if (!mark[u]) {
			B.push_back(u);
		}
	}
}

/**
	 * Symmetric Householder reduction to tridiagonal form.
	 * This implementation is based on a java implementation in JAMA package.
	 * 
	 * This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch,
	 * and Wilkinson, Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the
	 * corresponding Fortran subroutine in EISPACK.
*/
void tred2(int n, matrix &V, vector<double> &d, vector<double> &e) {
	for (int j = 0; j < n; j++) {
		d[j] = V[n - 1][j];
	}

	// Householder reduction to tridiagonal form.
	for (int i = n - 1; i > 0; i--) {

		// Scale to avoid under/overflow.
		double scale = 0.0;
		double h = 0.0;
		for (int k = 0; k < i; k++) {
			scale = scale + abs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i - 1];
			for (int j = 0; j < i; j++) {
				d[j] = V[i - 1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		}
		else {
			// Generate Householder vector.
			for (int k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i - 1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i - 1] = f - g;
			for (int j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.
			for (int j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j + 1; k <= i - 1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for (int j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for (int j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for (int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i - 1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i - 1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.
	for (int i = 0; i < n - 1; i++) {
		V[n - 1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i + 1];
		if (h != 0.0) {
			for (int k = 0; k <= i; k++) {
				d[k] = V[k][i + 1] / h;
			}
			for (int j = 0; j <= i; j++) {
				double g = 0.0;
				for (int k = 0; k <= i; k++) {
					g += V[k][i + 1] * V[k][j];
				}
				for (int k = 0; k <= i; k++) {
					V[k][j] -= g * d[k];
				}
			}
		}
		for (int k = 0; k <= i; k++) {
			V[k][i + 1] = 0.0;
		}
	}
	for (int j = 0; j < n; j++) {
		d[j] = V[n - 1][j];
		V[n - 1][j] = 0.0;
	}
	V[n - 1][n - 1] = 1.0;
	e[0] = 0.0;
}

/**
	 * Symmetric tridiagonal QL algorithm.
	 * This implementation is based on a java implementation in JAMA package.
	 * 
	 * This is derived from the Algol procedures tql2, by Bowdler, Martin, Reinsch,
	 * and Wilkinson, Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the
	 * corresponding Fortran subroutine in EISPACK.
 */
void tql2(int n, vector<double>& d, vector<double>& e, matrix& V) {
	for (int i = 1; i < n; i++) {
		e[i - 1] = e[i];
	}
	e[n - 1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	double eps = pow(2.0, -52.0);
	for (int l = 0; l < n; l++) {
		// Find small subdiagonal element
		tst1 = max(tst1, abs(d[l]) + abs(e[l]));
		int m = l;
		while (m < n) {
			if (abs(e[m]) <= eps * tst1) {
				break;
			}
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.
		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1; // (Could check iteration count here.)

				// Compute implicit shift
				double g = d[l];
				double p = (d[l + 1] - g) / (2.0 * e[l]);
				double r = hypot(p, 1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l + 1] = e[l] * (p + r);
				double dl1 = d[l + 1];
				double h = g - d[l];
				for (int i = l + 2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.
				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l + 1];
				double s = 0.0;
				double s2 = 0.0;
				for (int i = m - 1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot(p, e[i]);
					e[i + 1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i + 1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.
					for (int k = 0; k < n; k++) {
						h = V[k][i + 1];
						V[k][i + 1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.
			} while (abs(e[l]) > eps * tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.
	for (int i = 0; i < n - 1; i++) {
		int k = i;
		double p = d[i];
		for (int j = i + 1; j < n; j++) {
			if (d[j] < p) { // NH find smallest k>i
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i]; // swap k and i
			d[i] = p;
			for (int j = 0; j < n; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
}

/*
void check_spetral_decomposition(matrix& V, vector<double>& eigenvalues, matrix& eigenvectors) {
	// Vx = lambda x
	for (int i = 0; i < eigenvalues.size(); i++) {
		double lambda = eigenvalues[i];

		vector<double> a(eigenvectors[0].size(), 0);
		for (int j = 0; j < a.size(); j++) {
			for (int k = 0; k < V[j].size(); k++) {
				a[j] += V[j][k] * eigenvectors[i][k];
			}
		}

		for (int j = 0; j < a.size(); j++) {
			cout << a[j] << "\t";
		}
		cout << endl;

		for (int j = 0; j < eigenvectors[i].size(); j++) {
			cout << lambda * eigenvectors[i][j] << "\t";
		}
		cout << endl;
		cout << "--------------------" << endl;
	}
}
*/