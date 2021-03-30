#include <Engine/MeshEdit/Paramaterize.h>

#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;

Paramaterize::Paramaterize(Ptr<TriMesh> triMesh) : heMesh(make_shared<HEMesh<V>>()) {
	Init(triMesh);
}

void Paramaterize::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool Paramaterize::Init(Ptr<TriMesh> triMesh) {
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

bool Paramaterize::Run() {
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	Paramaterization();

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

	if (show)
		triMesh->Update(texcoords);
	else
		triMesh->Update(positions);

	return true;
}

void Paramaterize::Paramaterization() {
	SetBoundaryPoints();

	int nV = heMesh->NumVertices();
	vector<Eigen::Triplet<double> > triplets;

	for (int i = 0; i < nV; ++i) {
		V* vi = heMesh->Vertices()[i];
		triplets.push_back(Eigen::Triplet<double>(i, i, 1));

		if (!vi->IsBoundary()) {
			double adjVertexSize = vi->AdjVertices().size();
			for (int j = 0; j < adjVertexSize; ++j) {
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

	Eigen::VectorXd resX(nV), resY(nV), bX(nV), bY(nV);
	bX.setZero(); bY.setZero();
	for (int i = 0; i < boundaryIndex.size(); ++i) {
		V* vi = heMesh->Vertices()[i];
		bX(boundaryIndex[i]) = boundaryFixedPoints[i][0];
		bY(boundaryIndex[i]) = boundaryFixedPoints[i][1];
	}

	resX = solver.solve(bX);
	resY = solver.solve(bY);

	for (int i = 0; i < nV; ++i) {
		V* vi = heMesh->Vertices()[i];

		vi->pos.at(0) = resX(i);
		vi->pos.at(1) = resY(i);
		vi->pos.at(2) = 0;

		texcoords.push_back(pointf2(resX(i), resY(i)));
	}

	return;
}

void Paramaterize::SetBoundaryPoints() {
	int nB = heMesh->Boundaries()[0].size();
	for (int i = 0; i < nB; ++i) {
		boundaryIndex.push_back(heMesh->Index(heMesh->Boundaries()[0][i]->Origin()));
	}

	int pointsPerEdge = std::ceil(nB / 4);
	double step = 1.0 / pointsPerEdge;
	for (int i = 0; i < nB; ++i) {
		if (i < pointsPerEdge) {
			boundaryFixedPoints.push_back(pointf2(0, step * i));
		}
		else if (i >= pointsPerEdge && i < 2 * pointsPerEdge) {
			boundaryFixedPoints.push_back(pointf2((i - pointsPerEdge) * step, 1));
		}
		else if (i >= pointsPerEdge * 2 && i < pointsPerEdge * 3) {
			boundaryFixedPoints.push_back(pointf2(1, (3 * pointsPerEdge - i) * step));
		}
		else
			boundaryFixedPoints.push_back(pointf2((nB - i) * step, 0));
	}

	return;
}
