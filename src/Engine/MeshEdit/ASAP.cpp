#include <Engine/MeshEdit/ASAP.h>

#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

using namespace Ubpa;

using namespace std;

ASAP::ASAP(Ptr<TriMesh> triMesh): heMesh(make_shared<HEMesh<V>>()) {
	Init(triMesh);
}

void ASAP::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool ASAP::Init(Ptr<TriMesh> triMesh) {
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

bool ASAP::Run() {
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
		this->triMesh->Update(texcoords);
	else
		this->triMesh->Update(positions);

	return true;
}

void ASAP::Paramaterization() {
	int nV = heMesh->NumVertices();
	int nP = heMesh->NumPolygons();
	A = Eigen::SparseMatrix<double>(2 * (nV + nP), 2 * (nV + nP));
	A.setZero();
	b = Eigen::VectorXd(2 * (nV + nP));
	b.setZero();

	std::cout << "LocalFlatern" << std::endl;
	LocalFlatern();
	std::cout << "InitLaplacianMatrix" << std::endl;
	InitLaplacianMatrix();

	std::cout << "solve matrix" << std::endl;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A.transpose() * A);
	if (solver.info() != Eigen::Success) {
		cout << "compute is error" << endl;
		return;
	}

	Eigen::VectorXd res = solver.solve(A.transpose() * b);
	for (int i = 0; i < nV; ++i) {
		heMesh->Vertices()[i]->pos.at(0) = res(i);
		heMesh->Vertices()[i]->pos.at(1) = res(nV + i);
		heMesh->Vertices()[i]->pos.at(2) = 0;
		texcoords.push_back(pointf2(res(i), res(nV + i)));
	}

	return;
}

// map each polygon to (0, 0) (u0u1, 0), (u0u2cos(theta), u0u2sin(theta))
void ASAP::LocalFlatern() {
	int nP = heMesh->NumPolygons();
	int nV = heMesh->NumVertices();
	auto polys = heMesh->Polygons();

	for (int i = 0; i < nP; ++i) {
		auto p = polys[i];
		if (p != nullptr) {
			int polyIdx = heMesh->Index(p);
			V* vec0 = p->BoundaryVertice()[0];
			V* vec1 = p->BoundaryVertice()[1];
			V* vec2 = p->BoundaryVertice()[2];

			double vec0vec1Dist = pointf3::distance(vec0->pos.cast_to<pointf3>(), vec1->pos.cast_to<pointf3>());
			double vec0vec2Dist = pointf3::distance(vec0->pos.cast_to<pointf3>(), vec2->pos.cast_to<pointf3>());
			double cos12 = vecf3::cos_theta((vec1->pos - vec0->pos), (vec2->pos - vec0->pos));
			double sin12 = sqrt(1 - cos12 * cos12);

			// map pos
			std::map<V*, pointf2> mappedPos;
			mappedPos[vec0] = { 0,0 };
			mappedPos[vec1] = { vec0vec1Dist, 0 };
			mappedPos[vec2] = { vec0vec2Dist * cos12, vec0vec2Dist * sin12 };
			verticesPerPoly[polyIdx] = mappedPos;

			std::map<V*, double> cot;
			for (size_t j = 0; j < 3; j++) {
				double cos = vecf3::cos_theta((p->BoundaryVertice()[j]->pos - p->BoundaryVertice()[(j + 2) % 3]->pos),
					(p->BoundaryVertice()[(j + 1) % 3]->pos - p->BoundaryVertice()[(j + 2) % 3]->pos));
				double sin = sqrt(1 - cos * cos);
				cot[p->BoundaryVertice()[(j + 2) % 3]] = (cos / sin); // store the angle of edge i, i+1. idx is V_{(i+2) %3}
			}
			cotPerPoly[polyIdx] = cot;
		}
	}
}

void ASAP::InitLaplacianMatrix() {
	std::cout << "Init u" << std::endl;
	Initu();
	std::cout << "Init Lt" << std::endl;
	InitLt();
	A.setFromTriplets(triplets.begin(), triplets.end());
	A.makeCompressed();
}

void ASAP::Initu() {
	int nV = heMesh->NumVertices();
	int nP = heMesh->NumPolygons();

	// fix two points;
	auto boundaries = heMesh->Boundaries()[0];
	auto v0 = boundaries[0]->Origin();
	auto v1 = boundaries[boundaries.size() / 2]->End();
	int Idx0 = heMesh->Index(v0);
	int Idx1 = heMesh->Index(v1);

	triplets.push_back(Eigen::Triplet<double>(Idx0, Idx0, 1));
	triplets.push_back(Eigen::Triplet<double>(nV + Idx0, nV + Idx0, 1));
	b(Idx0) = 0;
	b(nV + Idx0) = 0;

	triplets.push_back(Eigen::Triplet<double>(Idx1, Idx1, 1));
	triplets.push_back(Eigen::Triplet<double>(nV + Idx1, nV + Idx1, 1));
	b(Idx1) = 1;
	b(nV + Idx1) = 1;

	// others points
	for (int i = 0; i < nV; ++i) {
		auto vertex = heMesh->Vertices()[i];
		int vertexIdx = heMesh->Index(vertex);

		if (vertexIdx != Idx0 && vertexIdx != Idx1) {
			double cotSum = 0;

			for (auto adjVertex : vertex->AdjVertices()) {
				int adjIdx = heMesh->Index(adjVertex);
				auto edge = vertex->EdgeWith(adjVertex);
				auto he0 = edge->HalfEdge();
				auto he1 = edge->HalfEdge()->Pair();

				double cot0 = 0;
				double cot1 = 0;

				auto triangle0 = he0->Polygon();
				// 当这条边只与一个三角形相接
				if (triangle0 != nullptr) {
					auto tr0Ver2 = he0->Next()->End();
					int tri0Idx = heMesh->Index(triangle0);

					cot0 = cotPerPoly[tri0Idx][tr0Ver2];
					map<V*, pointf2> mappedPoints = verticesPerPoly[tri0Idx];
					// Lt
					triplets.push_back(Eigen::Triplet<double>(vertexIdx, 2 * nV + tri0Idx, -cot0 * (mappedPoints[vertex][0] - mappedPoints[adjVertex][0])));
					triplets.push_back(Eigen::Triplet<double>(vertexIdx, 2 * nV + nP + tri0Idx, -cot0 * (mappedPoints[vertex][1] - mappedPoints[adjVertex][1])));
					triplets.push_back(Eigen::Triplet<double>(nV + vertexIdx, 2 * nV + nP + tri0Idx, cot0 * (mappedPoints[vertex][0] - mappedPoints[adjVertex][0])));
					triplets.push_back(Eigen::Triplet<double>(nV + vertexIdx, 2 * nV + tri0Idx, -cot0 * (mappedPoints[vertex][1] - mappedPoints[adjVertex][1])));
				}

				auto triangle1 = he1->Polygon();
				if (triangle1 != nullptr) {
					auto tr1Ver2 = he1->Next()->End();
					int tri1Idx = heMesh->Index(triangle1);

					cot1 = cotPerPoly[tri1Idx][tr1Ver2];
					map<V*, pointf2> mappedPoints = verticesPerPoly[tri1Idx];
					// Lt
					triplets.push_back(Eigen::Triplet<double>(vertexIdx, 2 * nV + tri1Idx, -cot1 * (mappedPoints[vertex][0] - mappedPoints[adjVertex][0])));
					triplets.push_back(Eigen::Triplet<double>(vertexIdx, 2 * nV + nP + tri1Idx, -cot1 * (mappedPoints[vertex][1] - mappedPoints[adjVertex][1])));
					triplets.push_back(Eigen::Triplet<double>(nV + vertexIdx, 2 * nV + nP + tri1Idx, cot1 * (mappedPoints[vertex][0] - mappedPoints[adjVertex][0])));
					triplets.push_back(Eigen::Triplet<double>(nV + vertexIdx, 2 * nV + tri1Idx, -cot1 * (mappedPoints[vertex][1] - mappedPoints[adjVertex][1])));
				}

				cotSum += (cot0 + cot1);
				triplets.push_back(Eigen::Triplet<double>(vertexIdx, adjIdx, -(cot0 + cot1)));
				triplets.push_back(Eigen::Triplet<double>(nV + vertexIdx, nV + adjIdx, -(cot0 + cot1)));
			}

			triplets.push_back(Eigen::Triplet<double>(vertexIdx, vertexIdx, cotSum));
			triplets.push_back(Eigen::Triplet<double>(nV + vertexIdx, nV + vertexIdx, cotSum));
		}
	}
}

void ASAP::InitLt() {
	int nP = heMesh->NumPolygons();
	int nV = heMesh->NumVertices();
	auto polys = heMesh->Polygons();
	
	for (int i = 0; i < nP; ++i) {
		auto p = polys[i];
		auto vertice = p->BoundaryVertice();
		int pIdx = heMesh->Index(p);
		
		auto vertex0 = vertice[0];
		auto vertex1 = vertice[1];
		auto vertex2 = vertice[2];
		int vertex0Idx = heMesh->Index(vertex0);
		int vertex1Idx = heMesh->Index(vertex1);
		int vertex2Idx = heMesh->Index(vertex2);

		double cot0 = cotPerPoly[i][vertex0];
		double cot1 = cotPerPoly[i][vertex1];
		double cot2 = cotPerPoly[i][vertex2];

		pointf2 mappedVec0 = verticesPerPoly[i][vertex0];
		pointf2 mappedVec1 = verticesPerPoly[i][vertex1];
		pointf2 mappedVec2 = verticesPerPoly[i][vertex2];

		double deltaX01 = mappedVec0[0] - mappedVec1[0];
		double deltaX12 = mappedVec1[0] - mappedVec2[0];
		double deltaX20 = mappedVec2[0] - mappedVec0[0];

		double deltaY01 = mappedVec0[1] - mappedVec0[1];
		double deltaY12 = mappedVec1[1] - mappedVec2[1];
		double deltaY20 = mappedVec2[1] - mappedVec0[1];

		// a
		triplets.push_back(Eigen::Triplet<double>(2 * nV + pIdx, 2 * nV + pIdx, -(cot2 * (deltaX01 * deltaX01 + deltaY01 * deltaY01) +
																				  cot0 * (deltaX12 * deltaX12 + deltaY12 * deltaY12) +
																				  cot1 * (deltaX20 * deltaX20 + deltaY20 * deltaY20))));
		// ux
		triplets.push_back(Eigen::Triplet<double>(2 * nV + pIdx, vertex0Idx, cot2 * deltaX01 - cot1 * deltaX20));
		triplets.push_back(Eigen::Triplet<double>(2 * nV + pIdx, vertex1Idx, cot0 * deltaX12 - cot2 * deltaX01));
		triplets.push_back(Eigen::Triplet<double>(2 * nV + pIdx, vertex2Idx, cot1 * deltaX20 - cot0 * deltaX12));
		// uy
		triplets.push_back(Eigen::Triplet<double>(2 * nV + pIdx, nV + vertex0Idx, cot2 * deltaY01 - cot1 * deltaY20));
		triplets.push_back(Eigen::Triplet<double>(2 * nV + pIdx, nV + vertex1Idx, cot0 * deltaY12 - cot2 * deltaY01));
		triplets.push_back(Eigen::Triplet<double>(2 * nV + pIdx, nV + vertex2Idx, cot1 * deltaY20 - cot0 * deltaY12));
		
		// b
		triplets.push_back(Eigen::Triplet<double>(2 * nV + nP + pIdx, 2 * nV + nP + pIdx, -(cot2 * (deltaX01 * deltaX01 + deltaY01 * deltaY01) +
																							cot0 * (deltaX12 * deltaX12 + deltaY12 * deltaY12) +
																							cot1 * (deltaX20 * deltaX20 + deltaY20 * deltaY20))));
		// ux
		triplets.push_back(Eigen::Triplet<double>(2 * nV + nP + pIdx, vertex0Idx, cot2 * deltaY01 - cot1 * deltaY20));
		triplets.push_back(Eigen::Triplet<double>(2 * nV + nP + pIdx, vertex1Idx, cot0 * deltaY12 - cot2 * deltaY01));
		triplets.push_back(Eigen::Triplet<double>(2 * nV + nP + pIdx, vertex2Idx, cot1 * deltaY20 - cot0 * deltaY12));
		// uy
		triplets.push_back(Eigen::Triplet<double>(2 * nV + nP + pIdx, nV + vertex0Idx, - (cot2 * deltaX01 - cot1 * deltaX20)));
		triplets.push_back(Eigen::Triplet<double>(2 * nV + nP + pIdx, nV + vertex1Idx, - (cot0 * deltaX12 - cot2 * deltaX01)));
		triplets.push_back(Eigen::Triplet<double>(2 * nV + nP + pIdx, nV + vertex2Idx, - (cot1 * deltaX20 - cot0 * deltaX12)));
	}
}

