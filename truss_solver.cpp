#include <iostream>
#include <vector>
#include <string>

struct Node {
    double x, y;            // Coordinates
    bool fixedX, fixedY;    // Boundary conditions (true if fixed)
    int id;                 // Unique identifier
};

struct Member {
    int node1, node2;       // Indices of connected nodes
    int id;                 // Unique identifier
    double length;          // Computed length (optional, for caching)
};

struct Load {
    int node;               // Node index where load is applied
    double fx, fy;          // Force components in x and y directions
};

class Truss {
private:
    std::vector<Node> nodes;
    std::vector<Member> members;
    std::vector<Load> loads;
    std::vector<double> displacements; // Solved displacements
    std::vector<double> forces;        // Solved member forces
    std::vector<double> reactions;     // Solved reaction forces
public:
    std::vector<double> getDisplacements() { return displacements; }
    std::vector<double> getForces() { return forces; }
    std::vector<double> getReactions() { return reactions; }
    size_t getNumNodes() { return nodes.size(); }
    void addNode(double x, double y, bool fixedX, bool fixedY);
    void addMember(int node1, int node2);
    void addLoad(int node, double fx, double fy);
    bool readFromFile(const std::string& filename);
    bool readFromJSON(const std::string& filename);
    void computeMemberLengths(); // Calculate and cache member lengths
    // void solve();                // Solve truss using direct stiffness method
    // using Solver class to change solver class separately
    void printResults() const;   // Output displacements, forces, reactions
};

class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows, cols;
public:
    Matrix(int r, int c);
    double& operator()(int i, int j);
    Matrix operator*(const Matrix& other) const;
    std::vector<double> solve(const std::vector<double>& b) const; // Gaussian elimination
};

class Solver {
private:
    const Truss& truss; // Reference to truss being solved
    Matrix globalStiffness; // Global stiffness matrix
public:
    Solver(Truss& t) : truss(t), globalStiffness(t.getNumNodes() * 2, t.getNumNodes() * 2) {}
    void assembleStiffness(); // Build global stiffness matrix
    void applyBoundaryConditions(); // Modify for fixed nodes
    void solveDisplacements(); // Solve for unknown displacements
    void computeForces(); // Calculate member forces
    void computeReactions(); // Calculate reaction forces
};

int main() {
    Truss truss;
    // Option 1: Console input
    int numNodes;
    std::cout << "Enter number of nodes: ";
    std::cin >> numNodes;
    for (int i = 0; i < numNodes; ++i) {
        double x, y;
        bool fx, fy;
        std::cout << "Enter x, y, fixedX, fixedY for node " << i + 1 << ": ";
        std::cin >> x >> y >> fx >> fy;
        truss.addNode(x, y, fx, fy);
    }
    Solver a(truss);
    truss.printResults();
    return 0;
}