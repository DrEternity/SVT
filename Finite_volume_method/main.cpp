#include "inmost.h"
#include <iostream>
#include <cstdio>

using namespace INMOST;
using namespace std;

const double dx = 1.0;
const double dy = 1.0;
const double dxy = 0.0;
const double pi = 3.1415926535898;
const double a = 1;

double conc_an(double x, double y)
{
    return (sin(pi * x) * sinh(pi * (1 - y)) + sin(pi * y) * sinh(pi * (1 - x))) / sinh(pi);
}

double source(double x, double y)
{
    return 0;
}

enum BoundCondType
{
    BC_DIR = 1,
    BC_NEUM = 2
};

// Class including everything needed
class Problem
{
private:
    /// Mesh
    Mesh &m;
    // =========== Tags =============
    /// Solution tag: 1 real value per cell
    Tag tagConc;
    /// Diffusion tensor tag: 3 real values (Dx, Dy, Dxy) per cell
    Tag tagD;
    /// Boundary condition type tag: 1 integer value per face, sparse on faces
    Tag tagBCtype;
    /// Boundary condition value tag: 1 real value per face, sparse on faces
    Tag tagBCval;
    /// Right-hand side tag: 1 real value per cell
    Tag tagSource;
    /// Analytical solution tag: 1 real value per cell
    Tag tagConcAn;
    /// Global index tag: 1 integer value per cell
    Tag tagGlobInd;

    // =========== Tag names ===========
    const string tagNameConc = "Concentration";
    const string tagNameD = "Diffusion_tensor";
    const string tagNameBCtype = "BC_type";
    const string tagNameBCval = "BC_value";
    const string tagNameSource = "Source";
    const string tagNameConcAn = "Concentration_analytical";
    const string tagNameGlobInd = "Global_Index";

public:
    Problem(Mesh &m_);
    ~Problem();
    void initProblem();
    void assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs);
    void run();
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
    tagConc = m.CreateTag(tagNameConc, DATA_REAL, CELL, NONE, 1);
    tagD = m.CreateTag(tagNameD, DATA_REAL, CELL, NONE, 3);
    tagBCtype = m.CreateTag(tagNameBCtype, DATA_INTEGER, FACE, FACE, 1);
    tagBCval = m.CreateTag(tagNameBCval, DATA_REAL, FACE, FACE, 1);
    tagSource = m.CreateTag(tagNameSource, DATA_REAL, CELL, CELL, 1);
    tagConcAn = m.CreateTag(tagNameConcAn, DATA_REAL, CELL, CELL, 1);
    tagGlobInd = m.CreateTag(tagNameGlobInd, DATA_INTEGER, CELL, NONE, 1);

    // Cell loop
    // 1. Set diffusion tensor values
    // 2. Write analytical solution and source tags
    // 3. Assign global indices
    int glob_ind = 0;
    for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
        Cell c = icell->getAsCell();
        c.RealArray(tagD)[0] = dx; // Dx
        c.RealArray(tagD)[1] = dy; // Dy
        c.RealArray(tagD)[2] = dxy; // Dxy
        double xc[2];
        c.Barycenter(xc);
        c.Real(tagConcAn) = conc_an(xc[0], xc[1]);
        c.Real(tagSource) = source(xc[0], xc[1]);
        c.Integer(tagGlobInd) = glob_ind;
        glob_ind++;
    }

    // Face loop:
    // 1. Set BC
    for(Mesh::iteratorFace iface = m.BeginFace(); iface != m.EndFace(); iface++){
        Face f = iface->getAsFace();
        if(!f.Boundary())
            continue;
        f.Integer(tagBCtype) = BC_DIR;
        double xf[2];
        f.Barycenter(xf);
        f.Real(tagBCval) = conc_an(xf[0], xf[1]);
    }
}

void Problem::assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs)
{
    // Face loop
    // Calculate transmissibilities using
    // two-point flux approximation (TPFA)
    for(Mesh::iteratorFace iface = m.BeginFace(); iface != m.EndFace(); iface++){
        Face f = iface->getAsFace();
        double xf[2];
        rMatrix nf(2,1);
        f.UnitNormal(nf.data());
        f.Barycenter(xf);
        if(f.Boundary()){
            int BCtype = f.Integer(tagBCtype);
            //грань Дирихле
            if(BCtype == BC_DIR){
                Cell cA;
                cA = f.BackCell();
                double t = 0.0; // transmissibility
                double xA[2];
                cA.Barycenter(xA);
                // implement by yourself
                rMatrix D(2, 2);
                D(0, 0) = cA.RealArray(tagD)[0];
                D(1, 1) = cA.RealArray(tagD)[1];
                D(0, 1) = cA.RealArray(tagD)[2];
                D(1, 0) = cA.RealArray(tagD)[2];

                double dA[2];
                dA[0] = xf[0] - xA[0];
                dA[1] = xf[1] - xA[1];

                double cf;
                double tmp1 = D(0, 0) * dA[0] + D(0, 1) * dA[1];
                double tmp2 = D(1, 0) * dA[0] + D(1, 1) * dA[1];
                cf = (tmp1 * nf(0,0) + tmp2*nf(1,0)) / (tmp1 * tmp1 + tmp2*tmp2);
                int i = cA.Integer(tagGlobInd);
                A[i][i] += cf * f.Area();
                rhs[i] += cf * conc_an(xf[0], xf[1]) * f.Area();
            }
        }
        else{
            // Internal face
            Cell cA, cB;
            cA = f.BackCell();
            cB = f.FrontCell();
            if(!cB.isValid()){
                cout << "Invalid FrontCell!" << endl;
                exit(1);
            }
            //аналогично найти x_A но еще и x_B
            double xA[2], xB[2];
            cA.Barycenter(xA), cB.Barycenter(xB);
            //и еще d_a and d_B
            double dA[2], dB[2];
            dA[0] = xf[0] - xA[0];
            dA[1] = xf[1] - xA[1];
            dB[0] = xf[0] - xB[0];
            dB[1] = xf[1] - xB[1];

            rMatrix D(2, 2);
            D(0, 0) = cA.RealArray(tagD)[0];
            D(1, 1) = cA.RealArray(tagD)[1];
            D(0, 1) = cA.RealArray(tagD)[2];
            D(1, 0) = cA.RealArray(tagD)[2];
            rMatrix D1(2, 2);
            D1(0, 0) = cB.RealArray(tagD)[0];
            D1(1, 1) = cB.RealArray(tagD)[1];
            D1(0, 1) = cB.RealArray(tagD)[2];
            D1(1, 0) = cB.RealArray(tagD)[2];
            double cfA;
            double tmp1a = D(0, 0) * dA[0] + D(0, 1) * dA[1];
            double tmp2a = D(1, 0) * dA[0] + D(1, 1) * dA[1];
            cfA = (tmp1a * nf(0,0) + tmp2a*nf(1,0)) / (tmp1a * tmp1a + tmp2a*tmp2a);
            double cfB;
            double tmp1b = D1(0, 0) * dB[0] + D1(0, 1) * dB[1];
            double tmp2b = D1(1, 0) * dB[0] + D1(1, 1) * dB[1];
            cfB = (tmp1b * nf(0,0) + tmp2b*nf(1,0)) / (tmp1b * tmp1b + tmp2b*tmp2b);
            double cf = (cfA*cfB)/(cfA-cfB);
            int i = cA.Integer(tagGlobInd);
            int j = cB.Integer(tagGlobInd);
            A[i][i] -= cf * f.Area(), A[j][j] -= cf * f.Area();
            A[i][j] += cf * f.Area(), A[j][i] += cf * f.Area();
        }
    }
    for(auto icell = m.BeginCell(); icell != m.EndCell(); icell++){
        Cell c = icell->getAsCell();
        int tmp = c.Integer(tagGlobInd);
        double center[2];
        c.Barycenter(center);
        rhs[tmp] -= source(center[0], center[1]) * c.Volume();
    }

}

void checkOrtho(Mesh *m)
{
    Tag T = m->CreateTag("Ortho", DATA_REAL, CELL, NONE, 1);
    for(auto iface = m->BeginFace(); iface != m->EndFace(); iface++){
        Face f = iface->getAsFace();
        Cell cp, cm = f.BackCell();
        double xp[3], xm[3], xf[3], nf[3], l[3];
        f.Barycenter(xf);
        f.UnitNormal(nf);
        cm.Barycenter(xm);

        if(f.FrontCell().isValid()){
            cp = f.FrontCell();
            cp.Barycenter(xp);
            for(int i = 0; i < 3; i++)
                l[i] = xp[i] - xm[i];
        }
        else{
            for(int i = 0; i < 3; i++)
                l[i] = xf[i] - xm[i];
        }

        double val = fabs(l[1]*nf[2] - l[2]*nf[1]);
        f.BackCell().Real(T) += val;
        if(val > 1e-5){
            //printf("Non-orthogonal face!\n");
            continue;
        }
        val = fabs(l[0]*nf[2] - l[2]*nf[0]);
        f.BackCell().Real(T) += val;
        if(val > 1e-5){
            //printf("Non-orthogonal face!\n");
            continue;
        }
        val = fabs(l[0]*nf[1] - l[1]*nf[0]);
        f.BackCell().Real(T) += val;
        if(val > 1e-5){
            //printf("Non-orthogonal face!\n");
            continue;
        }
    }
}

void Problem::run()
{
    // Matrix size
    unsigned N = static_cast<unsigned>(m.NumberOfCells());
    // Global matrix called 'stiffness matrix'
    Sparse::Matrix A;
    // Solution vector
    Sparse::Vector sol;
    // Right-hand side vector
    Sparse::Vector rhs;

    checkOrtho(&m);

    A.SetInterval(0, N);
    sol.SetInterval(0, N);
    rhs.SetInterval(0, N);

    assembleGlobalSystem(A, rhs);

    A.Save("A.mtx");
    rhs.Save("rhs.mtx");

    string solver_name = "inner_mptiluc";
    Solver S(solver_name);
    S.SetParameter("drop_tolerance", "0");

    S.SetMatrix(A);
    bool solved = S.Solve(rhs, sol);
    printf("Number of iterations: %d\n", S.Iterations());
    printf("Residual:             %e\n", S.Residual());
    if(!solved){
        printf("Linear solver failed: %s\n", S.GetReason().c_str());
        exit(1);
    }

    double normC = 0.0, normL2 = 0.0;
    for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
        Cell c = icell->getAsCell();
        unsigned ind = static_cast<unsigned>(c.Integer(tagGlobInd));
        c.Real(tagConc) = sol[ind];
        double diff = fabs(c.Real(tagConc) - c.Real(tagConcAn));
        normL2 += diff * c.Volume();
        normC = max(normC, diff);
    }
    printf("Error C-norm:  %e\n", normC);
    printf("Error L2-norm: %e\n", normL2);

    m.Save("res.vtk");
}


// Compute distance between two nodes
double nodeDist(const Node &n1, const Node &n2)
{
    double x1[3], x2[3];
    n1.Centroid(x1);
    n2.Centroid(x2);
    return sqrt(
            (  x1[0] - x2[0])*(x1[0] - x2[0])
            + (x1[1] - x2[1])*(x1[1] - x2[1])
            + (x1[2] - x2[2])*(x1[2] - x2[2])
    );
}

// Compute cell diameter
double cellDiam(const Cell &c)
{
    ElementArray<Node> nodes = c.getNodes();
    unsigned m = static_cast<unsigned>(nodes.size());
    double diam = 0.;
    for(unsigned i = 0; i < m; i++){
        for(unsigned j = 0; j < m; j++){
            diam = max(diam, nodeDist(nodes[i], nodes[j]));
        }
    }
    return diam;
}


// Compute and print mesh diameter
void main_mesh_diam(Mesh &m)
{
    double diam = 0.0;
    for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
        diam = max(diam, cellDiam(icell->self()));
    }
    cout << "Mesh diameter is " << diam << endl;
}



int main(int argc, char ** argv)
{
    if( argc < 2 )
    {
        printf("Usage: %s mesh_file\n",argv[0]);
        return -1;
    }

    Mesh m;
    m.Load(argv[1]);
    Mesh::GeomParam table;
    table[BARYCENTER] = CELL | FACE;
    table[CENTROID] = CELL | NODE;
    table[MEASURE] = CELL | FACE;
    m.PrepareGeometricData(table);
    Problem P(m);
    P.initProblem();
    P.run();
    main_mesh_diam(m);
    printf("Success\n");
    return 0;
}
