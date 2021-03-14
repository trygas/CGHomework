#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

MinSurf::MinSurf(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>())
{
	Init(triMesh);
}

void MinSurf::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool MinSurf::Init(Ptr<TriMesh> triMesh) {
	Clear();

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}

	// init half-edge structure
	size_t nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;
	triangles.reserve(triMesh->GetTriangles().size());
	for (auto triangle : triMesh->GetTriangles())
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary()) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundaries\n");
		heMesh->Clear();
		return false;
	}

	// triangle mesh's positions ->  half-edge structure's positions
	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<vecf3>();
	}

	this->triMesh = triMesh;
	return true;
}

bool MinSurf::Run() {
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	Minimize();

	// half-edge structure -> triangle mesh
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	for (auto v : heMesh->Vertices())
		positions.push_back(v->pos.cast_to<pointf3>());
	for (auto f : heMesh->Polygons()) { // f is triangle
		for (auto v : f->BoundaryVertice()) // vertices of the triangle
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
	}

	triMesh->Init(indice, positions);

	return true;
}

void MinSurf::Minimize() {
	int nV = heMesh->NumVertices();
	vector<Eigen::Triplet<double> > triplets;

	for (int i = 0; i < nV; ++i) {
		V* vi = heMesh->Vertices()[i];
		triplets.push_back(Eigen::Triplet<double>(i, i, 1));

		if (!vi->IsBoundary()) {
			double adjVertexSize = vi->AdjVertices().size();
			for(int j = 0; j < adjVertexSize; ++j){
				triplets.push_back(Eigen::Triplet<double>(i, heMesh->Index(vi->AdjVertices()[j]), -1.0 / adjVertexSize));
			}
		}
	}

	Eigen::SparseMatrix<double> A(nV, nV);
	A.setZero();
	A.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	
	if (solver.info() != Eigen::Success) {
		cout << "compute is error" << endl;
		return;
	}

	Eigen::VectorXd resX(nV), resY(nV), resZ(nV), bX(nV), bY(nV), bZ(nV);
	bX.setZero(); bY.setZero(); bZ.setZero();
	for (int i = 0; i < nV; ++i) {
		V* vi = heMesh->Vertices()[i];
		if (vi->IsBoundary()) {
			bX(i) = vi->pos.at(0);
			bY(i) = vi->pos.at(1);
			bZ(i) = vi->pos.at(2);
		}
	}

	resX = solver.solve(bX);
	resY = solver.solve(bY);
	resZ = solver.solve(bZ);

	for (int i = 0; i < nV; ++i) {
		V* vi = heMesh->Vertices()[i];

		vi->pos.at(0) = resX(i);
		vi->pos.at(1) = resY(i);
		vi->pos.at(2) = resZ(i);
		//cout << resX(i) << " " << resY(i) << " " << resZ(i) << endl;
	}

	return;
}
