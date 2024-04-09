// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    SparseMatrix<size_t> sparseMatrix(geometry->edgeIndices.size(), geometry->vertexIndices.size());

    for (Edge e : mesh->edges()) {
        auto eIndex = geometry->edgeIndices[e];
        auto beforPoint = e.halfedge().tipVertex();
        auto nextPoint = e.halfedge().tailVertex();

        auto beforIndex = geometry->vertexIndices[beforPoint];
        auto nextIndex = geometry->vertexIndices[nextPoint];

        sparseMatrix.insert(eIndex, beforIndex) = static_cast<size_t>(1);
        sparseMatrix.insert(eIndex, nextIndex) = static_cast<size_t>(1);
    }
    sparseMatrix.makeCompressed();
    return sparseMatrix;
    // return identityMatrix<size_t>(1); // placeholder
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    geometry->requireFaceIndices();
    geometry->requireEdgeIndices();

    SparseMatrix<size_t> spariseMatrix(geometry->faceIndices.size(), geometry->edgeIndices.size());

    for (Face f : mesh->faces()) {
        auto fIndex = geometry->faceIndices[f];
        auto adjEdges = f.adjacentEdges();
        for (Edge e : adjEdges) {
            auto eIndex = geometry->edgeIndices[e];
            spariseMatrix.insert(fIndex, eIndex) = static_cast<size_t>(1);
        }
    }
    spariseMatrix.makeCompressed();
    return spariseMatrix;
    // return identityMatrix<size_t>(1); // placeholder
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    if (subset.vertices.empty()) {
        return Vector<size_t>::Zero(1);
    }
    auto vertexSet = subset.vertices;
    Vector<size_t> ret(vertexSet.size());
    int i = 0;
    for (auto v : vertexSet) {
        // auto idx = geometry->vertexIndices[v];
        ret[i++] = v;
    }
    return ret;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    if (subset.edges.empty()) {
        return Vector<size_t>::Zero(1);
    }
    auto edgesSet = subset.edges;
    Vector<size_t> ret(edgesSet.size());
    int i = 0;
    for (auto e : edgesSet) {
        // auto idx = geometry->edgeIndices[e];
        ret[i++] = e;
    }
    return ret;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    if (subset.faces.empty()) {
        return Vector<size_t>::Zero(1);
    }
    auto faceSet = subset.faces;
    Vector<size_t> ret(faceSet.size());
    int i = 0;
    for (auto f : faceSet) {
        // auto idx = geometry->faceIndices[f];
        ret[i++] = f;
    }
    return ret;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO

    MeshSubset ret;
    auto ve_matrix = buildVertexEdgeAdjacencyMatrix();
    auto fe_matrix = buildFaceEdgeAdjacencyMatrix();
    for (auto v : subset.vertices) {
        ret.addVertex(v);
        for (Eigen::SparseMatrix<size_t>::InnerIterator it(ve_matrix, v); it; ++it) {
            auto edge = it.row();
            ret.addEdge(edge);

            for (Eigen::SparseMatrix<size_t>::InnerIterator it2(fe_matrix, edge); it2; ++it2) {
                auto face = it2.row();
                ret.addFace(face);
            }
        }
    }
    for (auto e : subset.edges) {
        ret.addEdge(e);
        for (Eigen::SparseMatrix<size_t>::InnerIterator it2(fe_matrix, e); it2; ++it2) {
            auto face = it2.row();
            ret.addFace(face);
        }
    }
    for (auto f : subset.faces) {
        ret.addFace(f);
    }

    return ret; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset ret;
    auto ve_matrix = buildVertexEdgeAdjacencyMatrix();
    auto fe_matrix = buildFaceEdgeAdjacencyMatrix();

    for (auto v : subset.vertices) {
        ret.addVertex(v);
    }
    for (auto e : subset.edges) {
        ret.addEdge(e);
        // std::cout << "add e: " << e << "\n";
        for (auto k = 0; k < ve_matrix.outerSize(); k++) {
            for (Eigen::SparseMatrix<size_t>::InnerIterator it2(ve_matrix, k); it2; ++it2) {
                if (it2.row() == e) {
                    auto vet = it2.col();
                    // std::cout << "add vet: " << vet << "\n";
                    ret.addVertex(vet);
                }
            }
        }
    }
    for (auto f : subset.faces) {
        ret.addFace(f);
        for (auto k = 0; k < fe_matrix.outerSize(); k++) {
            for (Eigen::SparseMatrix<size_t>::InnerIterator it(fe_matrix, k); it; ++it) {
                if (it.row() == f) {
                    auto edge = it.col();
                    ret.addEdge(edge);
                }
            }
        }
    }
    for (auto e : ret.edges) {
        ret.addEdge(e);
        // std::cout << "add e: " << e << "\n";
        for (auto k = 0; k < ve_matrix.outerSize(); k++) {
            for (Eigen::SparseMatrix<size_t>::InnerIterator it2(ve_matrix, k); it2; ++it2) {
                if (it2.row() == e) {
                    auto vet = it2.col();
                    // std::cout << "add vet: " << vet << "\n";
                    ret.addVertex(vet);
                }
            }
        }
    }
    return ret; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    auto cl_st = closure(star(subset));
    auto st_cl = star(closure(subset));

    cl_st.deleteSubset(st_cl);
    return cl_st; // placeholder
    // return subset; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}