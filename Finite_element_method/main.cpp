#include "inmost.h"
#include <stdio.h>

#include <cmath>

using namespace INMOST;
using namespace std;

enum BoundCondType
{
	BC_DIR = 1,
	BC_NEUM = 2
};

const double dx = 1.0;
const double dy = 1.0;
const double dxy = 0.0;
const double a = 1;
const double pi = 3.1415926535898;

double conc_an(double x, double y)
{
    //return (sin(pi * x) * sinh(pi * (1 - y)) + sin(pi * y) * sinh(pi * (1 - x))) / sinh(pi);
    return (1 / sinh(pi)) * (sin(pi * x) * sinh(pi * (1 - y))) + (1 / sinh(2 * pi)) * (sin(2 * pi * y) * sinh(2 * pi * x)) + (1 / sinh(3 * pi)) * (sin(3 * pi * x) * sinh(3 * pi * y)) + (1 / sinh(4 * pi)) * (sin(4 * pi * y) * sinh(4 * pi * (1 - x))) - (sin(5 * pi * x) * sin(6 * pi * y)) / (61 * (pi * pi));
}

double source(double x, double y)
{
    //return 0;
    return sin(5 * pi * x) * sin(6 * pi * y);
}

// Class including everything needed
class Problem
{
private:
	/// Mesh
	Mesh &m;
	// =========== Tags =============
	/// Solution tag: 1 real value per node
	Tag tagConc;
	/// Diffusion tensor tag: 3 real values (Dx, Dy, Dxy) per cell
	Tag tagD;
	/// Boundary condition type tag: 1 integer value per node
	Tag tagBCtype;
	/// Boundary condition value tag: 1 real value per node, sparse on nodes
	Tag tagBCval;
	/// Right-hand side tag: 1 real value per node, sparse on nodes
	Tag tagSource;
	/// Analytical solution tag: 1 real value per node
	Tag tagConcAn;
	/// Global index tag: 1 integer value per node
	Tag tagGlobInd;

	// =========== Tag names ===========
	const string tagNameConc = "Concentration";
	const string tagNameD = "Diffusion_tensor";
	const string tagNameBCtype = "BC_type";
	const string tagNameBCval = "BC_value";
	const string tagNameSource = "Source";
	const string tagNameConcAn = "Concentration_analytical";
	const string tagNameGlobInd = "Global_Index";

	// =========== Markers
	/// Marker for Dirichlet nodes
	MarkerType mrkDirNode;
	/// Number of Dirichlet nodes
	unsigned numDirNodes;

public:
	Problem(Mesh &m_);
	~Problem();
	void initProblem();
	void assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs);
	void assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc);
	void run();
	double get_L2_norm();
	double linear_approx_tri(double x, double y, const Cell &);
};

Problem::Problem(Mesh &m_) : m(m_)
{
}

Problem::~Problem()
{
}

void Problem::initProblem()
{
	// Init tags
	tagConc = m.CreateTag(tagNameConc, DATA_REAL, NODE, NONE, 1);
	tagD = m.CreateTag(tagNameD, DATA_REAL, CELL, NONE, 3);
	tagBCtype = m.CreateTag(tagNameBCtype, DATA_INTEGER, NODE, NODE, 1);
	tagBCval = m.CreateTag(tagNameBCval, DATA_REAL, NODE, NODE, 1);
	tagSource = m.CreateTag(tagNameSource, DATA_REAL, NODE, NONE, 1);
	tagConcAn = m.CreateTag(tagNameConcAn, DATA_REAL, NODE, NONE, 1);
	tagGlobInd = m.CreateTag(tagNameGlobInd, DATA_INTEGER, NODE, NONE, 1);

	// Cell loop
	// 1. Check that cell is a triangle
	// 2. Set diffusion tensor values
	for (Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++)
	{
		Cell c = icell->getAsCell();
		ElementArray<Node> nodes = c.getNodes();
		if (nodes.size() != 3)
		{
			cout << "Cell is not a triangle!!!";
			exit(1);
		}

		c.RealArray(tagD)[0] = dx;
		c.RealArray(tagD)[1] = dy;
		c.RealArray(tagD)[2] = dxy;
	}

	// Node loop
	// 1. Write analytical solution and source
	// 2. Check if node is boundary, set BC type and value, mark it
	// 3. Assign global indices
	mrkDirNode = m.CreateMarker();
	numDirNodes = 0;
	int glob_ind = 0;
	for (Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++)
	{
		Node n = inode->getAsNode();
		double xn[3];
		n.Centroid(xn);
		n.Real(tagConcAn) = conc_an(xn[0], xn[1]);
		n.Real(tagSource) = source(xn[0], xn[1]);

		if (n.Boundary())
		{
			n.SetMarker(mrkDirNode);
			n.Integer(tagBCtype) = BC_DIR;
			numDirNodes++;
			n.Real(tagBCval) = n.Real(tagConcAn);
		}
		else
		{
			n.Integer(tagGlobInd) = glob_ind;
			glob_ind++;
		}
	}
	printf("Number of Dirichlet nodes: %d\n", numDirNodes);
}

double basis_func(double x_, double y_, const Cell &c, const Node &n)
{
	ElementArray<Node> nodes = c.getNodes();
	unsigned n_ind = 0;
	double x[3];
	double y[3];
	for (unsigned i = 0; i < 3; i++)
	{
		if (n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}

	if (n_ind == 0)
	{
		return ((x_ - x[2]) * (y[1] - y[2]) - (x[1] - x[2]) * (y_ - y[2])) /
			   ((x[0] - x[2]) * (y[1] - y[2]) - (x[1] - x[2]) * (y[0] - y[2]));
	}
	else if (n_ind == 1)
	{
		return ((x_ - x[2]) * (y[0] - y[2]) - (x[0] - x[2]) * (y_ - y[2])) /
			   ((x[1] - x[2]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[1] - y[2]));
	}
	else if (n_ind == 2)
	{
		return ((x_ - x[0]) * (y[1] - y[0]) - (x[1] - x[0]) * (y_ - y[0])) /
			   ((x[2] - x[0]) * (y[1] - y[0]) - (x[1] - x[0]) * (y[2] - y[0]));
	}
	else
	{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
}

rMatrix basis_func_grad(const Cell &c, const Node &n)
{
	ElementArray<Node> nodes = c.getNodes();
	double x[3];
	double y[3];
	// gradient of the basis function
	rMatrix grad(2, 1);
	unsigned n_ind = 0;
	for (unsigned i = 0; i < 3; i++)
	{
		if (n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}

	if (n_ind == 0)
	{
		grad(0, 0) = (y[1] - y[2]);
		grad(1, 0) = -(x[1] - x[2]);
		grad /= ((x[0] - x[2]) * (y[1] - y[2]) - (x[1] - x[2]) * (y[0] - y[2]));
	}
	else if (n_ind == 1)
	{
		grad(0, 0) = (y[0] - y[2]);
		grad(1, 0) = -(x[0] - x[2]);
		grad /= ((x[1] - x[2]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[1] - y[2]));
	}
	else if (n_ind == 2)
	{
		grad(0, 0) = (y[1] - y[0]);
		grad(1, 0) = -(x[1] - x[0]);
		grad /= ((x[2] - x[0]) * (y[1] - y[0]) - (x[1] - x[0]) * (y[2] - y[0]));
	}
	else
	{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
	return grad;
}

void coords_from_barycentric(double *node_x, double *node_y, double *eta, double *x, double *y)
{
	*x = node_x[0] * eta[0] + node_x[1] * eta[1] + node_x[2] * eta[2];
	*y = node_y[0] * eta[0] + node_y[1] * eta[1] + node_y[2] * eta[2];
}

double Problem::linear_approx_tri(double x, double y, const Cell &c)
{
	ElementArray<Node> nodes = c.getNodes();
	double res = 0.0;
	for (unsigned i = 0; i < 3; i++)
	{
		res += nodes[i].Real(tagConc) * basis_func(x, y, c, nodes[i]);
	}
	return conc_an(x, y) - res;
}

double integrate_over_triangle(const Cell &c, const Node n, double (*f)(double, double, const Cell &, const Node &))
{
	double res = 0.0;
	double w3 = 0.205950504760887;
	double w6 = 0.063691414286223;
	double eta3[3] = {0.124949503233232, 0.437525248383384, 0.437525248383384};
	double eta6[3] = {0.797112651860071, 0.165409927389841, 0.037477420750088};

	ElementArray<Node> nodes = c.getNodes();
	if (nodes.size() != 3)
	{
		printf("Cell is not a triangle, has %lu nodes!\n", nodes.size());
		exit(1);
	}
	// Coordinates of triangle nodes
	double node_x[3], node_y[3];
	// Set them
	for (unsigned i = 0; i < 3; i++)
	{
		double c[3];
		nodes[i].Centroid(c);
		node_x[i] = c[0];
		node_y[i] = c[1];
	}

	// Add contribution from all combinations in eta3
	double x, y, val;
	double eta[3];
	eta[0] = eta3[0];
	eta[1] = eta3[1];
	eta[2] = eta3[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	// printf("x = %e, y = %e, val = %e\n", x, y, val);
	res += w3 * val;
	eta[0] = eta3[1];
	eta[1] = eta3[2];
	eta[2] = eta3[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w3 * val;
	eta[0] = eta3[2];
	eta[1] = eta3[0];
	eta[2] = eta3[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w3 * val;

	// Add contribution from all combinations in eta6
	eta[0] = eta6[0];
	eta[1] = eta6[1];
	eta[2] = eta6[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[0];
	eta[1] = eta6[2];
	eta[2] = eta6[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[1];
	eta[1] = eta6[0];
	eta[2] = eta6[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[1];
	eta[1] = eta6[2];
	eta[2] = eta6[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[2];
	eta[1] = eta6[0];
	eta[2] = eta6[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;
	eta[0] = eta6[2];
	eta[1] = eta6[1];
	eta[2] = eta6[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = f(x, y, c, n);
	res += w6 * val;

	res *= c.Volume();
	return res;
}

double integrate_over_triangle1(const Cell &c, Problem &p)
{
	double res = 0.0;
	double w3 = 0.205950504760887;
	double w6 = 0.063691414286223;
	double eta3[3] = {0.124949503233232, 0.437525248383384, 0.437525248383384};
	double eta6[3] = {0.797112651860071, 0.165409927389841, 0.037477420750088};

	ElementArray<Node> nodes = c.getNodes();
	if (nodes.size() != 3)
	{
		printf("Cell is not a triangle, has %lu nodes!\n", nodes.size());
		exit(1);
	}
	// Coordinates of triangle nodes
	double node_x[3], node_y[3];
	// Set them
	for (unsigned i = 0; i < 3; i++)
	{
		double c[3];
		nodes[i].Centroid(c);
		node_x[i] = c[0];
		node_y[i] = c[1];
	}

	// Add contribution from all combinations in eta3
	double x, y, val;
	double eta[3];
	eta[0] = eta3[0];
	eta[1] = eta3[1];
	eta[2] = eta3[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	// printf("x = %e, y = %e, val = %e\n", x, y, val);
	res += w3 * val;
	eta[0] = eta3[1];
	eta[1] = eta3[2];
	eta[2] = eta3[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	res += w3 * val;
	eta[0] = eta3[2];
	eta[1] = eta3[0];
	eta[2] = eta3[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	res += w3 * val;

	// Add contribution from all combinations in eta6
	eta[0] = eta6[0];
	eta[1] = eta6[1];
	eta[2] = eta6[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	res += w6 * val;
	eta[0] = eta6[0];
	eta[1] = eta6[2];
	eta[2] = eta6[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	res += w6 * val;
	eta[0] = eta6[1];
	eta[1] = eta6[0];
	eta[2] = eta6[2];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	res += w6 * val;
	eta[0] = eta6[1];
	eta[1] = eta6[2];
	eta[2] = eta6[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	res += w6 * val;
	eta[0] = eta6[2];
	eta[1] = eta6[0];
	eta[2] = eta6[1];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	res += w6 * val;
	eta[0] = eta6[2];
	eta[1] = eta6[1];
	eta[2] = eta6[0];
	coords_from_barycentric(node_x, node_y, eta, &x, &y);
	val = p.linear_approx_tri(x, y, c);
	res += w6 * val;

	return res * res * c.Volume();
}

double fphi(double x, double y, const Cell &c, const Node &n)
{
	return source(x, y) * basis_func(x, y, c, n);
}

double dot_D(const Cell &c, const Node &n1, const Node &n2, const rMatrix &D)
{
	rMatrix Dgradu(2, 1), gradu = basis_func_grad(c, n1), gradv = basis_func_grad(c, n2);
	Dgradu(0, 0) = D(0, 0) * gradu(0, 0) + D(0, 1) * gradu(1, 0);
	Dgradu(1, 0) = D(1, 0) * gradu(0, 0) + D(1, 1) * gradu(1, 0);
	return Dgradu(0, 0) * gradv(0, 0) + Dgradu(1, 0) * gradv(1, 0);
}

double basis_func_dot_source(double x, double y, const Cell &c, const Node &n)
{
	return basis_func(x, y, c, n) * source(x, y);
}

// [ [phi1,phi1], [phi1,phi2]...   ]
// [   ]
// [   ]
// energetic scalar product: [u,v] = (D*grad(u), grad(v)) = (grad(v))^T * D * grad(u)
void Problem::assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc)
{
	rMatrix D(2, 2);
	D(0, 0) = c.RealArray(tagD)[0];
	D(1, 1) = c.RealArray(tagD)[1];
	D(0, 1) = D(1, 0) = c.RealArray(tagD)[2];

	ElementArray<Node> nodes = c.getNodes();

	A_loc = rMatrix(3, 3);
	rhs_loc = rMatrix(3, 1);

	for (unsigned i = 0; i < 3; i++)
	{
		for (unsigned j = 0; j < 3; ++j)
		{
			A_loc(i, j) = dot_D(c, nodes[i], nodes[j], D) * c.Volume();
		}
		rhs_loc(i, 0) = integrate_over_triangle(c, nodes[i], &basis_func_dot_source);
	}
}

void Problem::assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs)
{
	// Cell loop
	// For each cell assemble local system
	// and incorporate it into global
	for (Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++)
	{
		Cell c = icell->getAsCell();
		rMatrix A_loc, rhs_loc;
		assembleLocalSystem(c, A_loc, rhs_loc);
		// Now A_loc is 3x3, rhs_loc is 3x1
		//

		ElementArray<Node> nodes = c.getNodes();
		unsigned glob_ind[3];
		for (unsigned loc_ind = 0; loc_ind < 3; loc_ind++)
		{
			// How to determine global index?
			glob_ind[loc_ind] = nodes[loc_ind].Integer(tagGlobInd);
		}

		for (unsigned loc_ind = 0; loc_ind < 3; loc_ind++)
		{
			// Consider node with local index 'loc_ind'

			// Check if this is a Dirichlet node
			if (nodes[loc_ind].GetMarker(mrkDirNode))
				continue;

			for (unsigned j = 0; j < 3; j++)
			{
				if (nodes[j].GetMarker(mrkDirNode))
				{
					rhs[glob_ind[loc_ind]] -= A_loc(loc_ind, j) * nodes[j].Real(tagBCval);
				}
				else
					A[glob_ind[loc_ind]][glob_ind[j]] += A_loc(loc_ind, j);
			}
			rhs[glob_ind[loc_ind]] += rhs_loc(loc_ind, 0);
		}
	}
}

double get_c_norm(Mesh &m)
{
	double normC = 0.0;
	for (Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++)
	{
		Node n = inode->getAsNode();
		double x[3];
		n.Centroid(x);
		double diff_n = abs(n.Real(m.GetTag("Concentration_analytical")) - n.Real(m.GetTag("Concentration")));
		if (diff_n > normC)
		{
			normC = diff_n;
		}
	}
	return normC;
}

double Problem::get_L2_norm()
{
	double normL2 = 0.0;
	for (Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++)
	{
		Cell c = icell->getAsCell();
		normL2 += integrate_over_triangle1(c, *this);
	}
	normL2 = sqrt(normL2);
	return normL2;
}

void Problem::run()
{
	// Matrix size
	unsigned N = static_cast<unsigned>(m.NumberOfNodes()) - numDirNodes;
	// Global matrix called 'stiffness matrix'
	Sparse::Matrix A;
	// Solution vector
	Sparse::Vector sol;
	// Right-hand side vector
	Sparse::Vector rhs;

	A.SetInterval(0, N);
	sol.SetInterval(0, N);
	rhs.SetInterval(0, N);

	assembleGlobalSystem(A, rhs);

	A.Save("A.mtx");
	rhs.Save("rhs.mtx");

	string solver_name = "inner_mptiluc";
	Solver S(solver_name);

	S.SetParameter("absolute_tolerance", "1e-14");
	S.SetParameter("relative_tolerance", "1e-14");

	S.SetMatrix(A);
	bool solved = S.Solve(rhs, sol);
	if (!solved)
	{
		printf("Linear solver failed: %s\n", S.GetReason().c_str());
		printf("Number of iterations: %d\n", S.Iterations());
		printf("Residual:             %e\n", S.Residual());
		exit(1);
	}

	for (Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++)
	{
		Node n = inode->getAsNode();
		if (n.GetMarker(mrkDirNode))
		{
			n.Real(tagConc) = n.Real(tagBCval);
			continue;
		}
		unsigned ind = static_cast<unsigned>(n.Integer(tagGlobInd));
		n.Real(tagConc) = sol[ind];
	}
	m.Save("res.vtk");
}

// Compute distance between two nodes
double nodeDist(const Node &n1, const Node &n2)
{
	double x1[3], x2[3];
	n1.Centroid(x1);
	n2.Centroid(x2);
	return sqrt(
		(x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[1] - x2[1]) * (x1[1] - x2[1]) + (x1[2] - x2[2]) * (x1[2] - x2[2]));
}

// Compute cell diameter
double cellDiam(const Cell &c)
{
	ElementArray<Node> nodes = c.getNodes();
	unsigned m = static_cast<unsigned>(nodes.size());
	double diam = 0.;
	for (unsigned i = 0; i < m; i++)
	{
		for (unsigned j = 0; j < m; j++)
		{
			diam = max(diam, nodeDist(nodes[i], nodes[j]));
		}
	}
	return diam;
}

// Compute and print mesh diameter
void main_mesh_diam(Mesh &m)
{
	double diam = 0.0;
	for (Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++)
	{
		diam = max(diam, cellDiam(icell->self()));
	}
	cout << "Mesh diameter is " << diam << endl;
}


int main(int argc, char **argv)
{
	if (argc < 2)
	{
		printf("Usage: %s mesh_file\n", argv[0]);
		return -1;
	}

	Mesh m;
	m.Load(argv[1]);
	Problem P(m);
	P.initProblem();
	P.run();

	Mesh m1;
	m1.Load("res.vtk");
	main_mesh_diam(m1);
	cout << "|u - u_approx|_C = " << get_c_norm(m1) << endl;
	cout << "|u - u_approx|_L2 = " << P.get_L2_norm() << endl;
	printf("Success\n");
}