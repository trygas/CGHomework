#include <Engine/MeshEdit/HybridParamaterize.h>
#include <Engine/MeshEdit/ASAP.h>
#include <Engine/MeshEdit/Paramaterize.h>

#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseQR>
#include <Eigen/SVD>

#include <cmath>

using namespace Ubpa;

using namespace std;

HybridParamaterize::HybridParamaterize(Ptr<TriMesh> triMesh) : heMesh(make_shared<HEMesh<V>>()) {
	Init(triMesh);
}

void HybridParamaterize::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool HybridParamaterize::Init(Ptr<TriMesh> triMesh) {
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

bool HybridParamaterize::Run() {
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

void HybridParamaterize::Paramaterization() {
	int nV = heMesh->NumVertices();
	int nP = heMesh->NumPolygons();
	A = Eigen::SparseMatrix<double>(nV, nV);
	A.setZero();
	b = Eigen::MatrixXd(nV, 2);
	b.setZero();

	std::cout << "LocalFlatern" << std::endl;
	LocalFlatern();
	std::cout << "InitA" << std::endl;
	InitA();

	auto asap = ASAP::New(triMesh);
	asap->SetDisplay(true);
	asap->Run();
	for (int i = 0; i < heMesh->NumVertices(); i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = vecf3(triMesh->GetTexcoords()[i][0], triMesh->GetTexcoords()[i][1], 0);
	}

	for (int i = 0; i < iterTimes; ++i) {
		LocalUpdate();
		double maxError = GlobalUpdate();
		std::cout << "iter: " << i << ", maxError: " << maxError << std::endl;
		if (maxError < errorThreshold)
			break;
	}

	return;
}

// map each polygon to (0, 0) (u0u1, 0), (u0u2cos(theta), u0u2sin(theta))
void HybridParamaterize::LocalFlatern() {
	int nP = heMesh->NumPolygons();
	int nV = heMesh->NumVertices();
	auto polys = heMesh->Polygons();

	for (auto p : polys) {
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

// 和ASAP不一样，ARAP由于采用了迭代的思想，这里只构造拉普拉斯矩阵中不变的部分
void HybridParamaterize::InitA() {
	int nV = heMesh->NumVertices();

	// fix one points;
	auto triangle = heMesh->Polygons().back();
	auto v0 = triangle->BoundaryVertice()[0];
	int Idx0 = heMesh->Index(v0);

	triplets.push_back(Eigen::Triplet<double>(Idx0, Idx0, 1));

	// others points
	for (auto vertex : heMesh->Vertices()) {
		int vertexIdx = heMesh->Index(vertex);

		if (vertexIdx != Idx0) {
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
				}

				auto triangle1 = he1->Polygon();
				if (triangle1 != nullptr) {
					auto tr1Ver2 = he1->Next()->End();
					int tri1Idx = heMesh->Index(triangle1);

					cot1 = cotPerPoly[tri1Idx][tr1Ver2];
				}

				cotSum += (cot0 + cot1);
				triplets.push_back(Eigen::Triplet<double>(vertexIdx, adjIdx, -(cot0 + cot1)));
			}

			triplets.push_back(Eigen::Triplet<double>(vertexIdx, vertexIdx, cotSum));
		}
	}

	A.setFromTriplets(triplets.begin(), triplets.end());
	A.makeCompressed();
	solver.compute(A.transpose() * A);
}

double Pow13(double first) {
	if (first < 0) {
		return -pow(-first, 1.0 / 3.0);
	}
	else
		return pow(first, 1.0 / 3.0);
}

void HybridParamaterize::LocalUpdate() {
	LtMatrix.clear();
	for (auto triangle : heMesh->Polygons()) {
		double C1 = 0, C2 = 0, C3 = 0;
		int triIdx = heMesh->Index(triangle);
		auto trianglePoints = triangle->BoundaryVertice();
		auto mappedPoints = verticesPerPoly[triIdx];

		for (int i = 0; i < 3; ++i) {
			double cot = cotPerPoly[triIdx][trianglePoints[(i + 2) % 3]];
			V* u0 = trianglePoints[i];
			V* u1 = trianglePoints[(i + 1) % 3];

			pointf2 x0 = mappedPoints[u0];
			pointf2 x1 = mappedPoints[u1];

			double deltaU0 = u0->pos[0] - u1->pos[0], deltaU1 = u0->pos[1] - u1->pos[1];
			double deltaX0 = x0[0] - x1[0], deltaX1 = x0[1] - x1[1];

			C1 += cot * (deltaX0 * deltaX0 + deltaX1 * deltaX1);
			C2 += cot * (deltaX0 * deltaU0 + deltaX1 * deltaU1);
			C3 += cot * (deltaU0 * deltaX1 - deltaU1 * deltaX0);
		}

		double a = 0, b = 0;
		if (lambda == 0) {
			a = C2 / C1;
			b = C3 / C1;
		}
		else if (lambda >= 1e6) {
			a = C2 / sqrt(C2 * C2 + C3 * C3);
			b = C3 / sqrt(C2 * C2 + C3 * C3);
		}
		else {
			double thirdCoff = 2 * lambda * (1 + (C3 / C2) * (C3 / C2));
			double firstCoff = C1 - 2 * lambda;
			double zeroCoff = -C2;

			double p = firstCoff / thirdCoff;
			double q = zeroCoff / thirdCoff;
			double delta = pow(q / 2, 2) + pow(p / 3, 3);
			if (delta < 0) {
				double r = pow(-pow(p / 3, 3), 0.5);
				double theta = acos(-q / (2 * r)) / 3.0;
				a = 2 * pow(r, 1.0 / 3.0) * cos(theta + (3.141592653) / 3.0);
				b = (C3 / C2) * a;
				if (fabs(a) < 5e-1 || fabs(b) < 5e-1) {
					a = C2 / sqrt(C2 * C2 + C3 * C3);
					b = C3 / sqrt(C2 * C2 + C3 * C3);
				}
			}
			else {
				double midRes = pow(delta, 0.5);
				a = Pow13(-q / 2 + midRes) + Pow13(-q / 2 - midRes);
				b = (C3 / C2) * a;
			}
		}

		Eigen::Matrix2d Lt;
		Lt.setZero();
		Lt << a, b, -b, a;
		LtMatrix[triIdx] = Lt;
	}
}

double HybridParamaterize::GlobalUpdate() {
	Initb();
	double maxError = -1;
	Eigen::MatrixXd res = solver.solve(A.transpose() * b);

	texcoords.clear();
	int nV = heMesh->NumVertices();
	for (int i = 0; i < nV; ++i) {
		double dist = pointf3::distance(pointf3(res(i, 0), res(i, 1), 0.0), heMesh->Vertices()[i]->pos.cast_to<pointf3>());
		if (maxError < 0 || maxError < dist)
			maxError = dist;

		heMesh->Vertices()[i]->pos[0] = res(i, 0);
		heMesh->Vertices()[i]->pos[1] = res(i, 1);
		heMesh->Vertices()[i]->pos[2] = 0.0;
		texcoords.push_back(pointf2(res(i, 0), res(i, 1)));
	}

	return maxError;
}

void HybridParamaterize::Initb() {
	b.setZero();

	auto triangle = heMesh->Polygons().back();
	auto v0 = triangle->BoundaryVertice()[0];
	int Idx0 = heMesh->Index(v0);
	b(Idx0, 0) = 1;
	b(Idx0, 1) = 1;

	for (auto vertex : heMesh->Vertices()) {
		int idx = heMesh->Index(vertex);
		if (idx != Idx0) {
			Eigen::Matrix<double, 2, 1> bt;
			bt.setZero();
			for (auto adjVertex : vertex->AdjVertices()) {
				int adjIdx = heMesh->Index(adjVertex);
				auto edge = vertex->EdgeWith(adjVertex);
				auto he0 = edge->HalfEdge();
				auto he1 = edge->HalfEdge()->Pair();

				auto triangle0 = he0->Polygon();
				if (triangle0 != nullptr) {
					auto tr0Ver2 = he0->Next()->End();
					int triIdx = heMesh->Index(triangle0);
					double cot0 = cotPerPoly[triIdx][tr0Ver2];

					std::map<V*, pointf2> mappedPoints = verticesPerPoly[triIdx];
					Eigen::Matrix<double, 2, 2> Lt = LtMatrix[triIdx];
					Eigen::Matrix<double, 2, 1> deltaX;
					deltaX << mappedPoints[vertex][0] - mappedPoints[adjVertex][0], mappedPoints[vertex][1] - mappedPoints[adjVertex][1];

					bt += cot0 * Lt * deltaX;
				}

				auto triangle1 = he1->Polygon();
				if (triangle1 != nullptr) {
					auto tr1Ver2 = he1->Next()->End();
					int triIdx = heMesh->Index(triangle1);
					double cot1 = cotPerPoly[triIdx][tr1Ver2];

					std::map<V*, pointf2> mappedPoints = verticesPerPoly[triIdx];
					Eigen::Matrix<double, 2, 2> Lt = LtMatrix[triIdx];
					Eigen::Matrix<double, 2, 1> deltaX;
					deltaX << mappedPoints[vertex][0] - mappedPoints[adjVertex][0], mappedPoints[vertex][1] - mappedPoints[adjVertex][1];

					bt += cot1 * Lt * deltaX;
				}
			}

			b(idx, 0) = bt(0);
			b(idx, 1) = bt(1);
		}
	}
}
