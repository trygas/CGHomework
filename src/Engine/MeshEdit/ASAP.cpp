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

	triMesh->Init(indice, positions);
	triMesh->Update(texcoords);

	return true;
}

void ASAP::Paramaterization() {
	LocalFlatern();

	int nP = heMesh->NumPolygons();
	int nV = heMesh->NumVertices();
	auto polys = heMesh->Polygons();
	
	Eigen::MatrixXd A(6 * nP, 2 * nV + 2 * nP);
	for (int i = 0; i < nP; ++i) {
		auto p = polys[i];
		for (int j = 0; j < 3; ++j) {
			auto vtj = p->BoundaryVertice()[j];
			int index = heMesh->Index(vtj);
			if (index == 0 || index == nV - 1) {

			}
			
		}
	}
}

void ASAP::LocalFlatern() {
	int nP = heMesh->NumPolygons();
	int nV = heMesh->NumVertices();
	auto polys = heMesh->Polygons();

	for (int i = 0; i < nP; ++i) {
		auto p = polys[i];
		vecf3 vec0 = p->BoundaryVertice()[0]->pos;
		vecf3 vec1 = p->BoundaryVertice()[1]->pos;
		vecf3 vec2 = p->BoundaryVertice()[2]->pos;
		
		int index0 = heMesh->Index(p->BoundaryVertice()[0]);
		int index1 = heMesh->Index(p->BoundaryVertice()[1]);
		int index2 = heMesh->Index(p->BoundaryVertice()[2]);
		std::array<int, 3> tempIndex = {index0, index1, index2};
		indexPerPoly.push_back(tempIndex);

		double cosine[3] = { 0 };
		cosine[0] = vecf3::cos_theta(vec1 - vec0, vec2 - vec0);
		cosine[1] = vecf3::cos_theta(vec0 - vec1, vec2 - vec1);
		cosine[2] = vecf3::cos_theta(vec0 - vec2, vec1 - vec2);
		
		std::vector<pointf2> thisPoly;
		thisPoly.push_back(pointf2(0, 0));
		thisPoly.push_back(pointf2((vec1 - vec0).norm(), 0));
		thisPoly.push_back(pointf2((vec2 - vec0).norm() * cosine[0], (vec2 - vec0).norm() * sqrt(1 - cosine[0] * cosine[0])));
		verticesPerPoly.push_back(thisPoly);

		
	}
}

//void ASAP::LocalFlatern() {
//	random_set<V*> boundary_points;
//	random_set<V*> inner_points;
//
//	auto boundaries = this->heMesh->Boundaries();
//	if (boundaries.size() != 1) {
//		cout << "ERROR::Parameterize::DoPara:" << endl
//			<< "\t" << "got boundaries = " << boundaries.size()
//			<< " (expect 1)" << endl;
//		return;
//	}
//
//	for (auto v : boundaries[0]) {
//		boundary_points.insert(v->Origin());
//	}
//
//	for (auto v : heMesh->Vertices()) {
//		if (!boundary_points.contains(v)) {
//			inner_points.insert(v);
//		}
//	}
//
//	//const float boost_factor = 100;
//	const float boost_factor = 1;
//	// Fix our boundary
//
//	int points_total = boundary_points.size();
//	float step = 4.0f / points_total;
//
//	float curr = 0;
//	for (auto v : boundary_points) {
//		vecf3 new_pos;
//		if (curr >= 0 && curr < 1) {
//			new_pos[0] = curr;
//			new_pos[1] = new_pos[2] = 0;
//		}
//		else if (curr >= 1 && curr < 2) {
//			new_pos[0] = 1;
//			new_pos[1] = curr - 1;
//			new_pos[2] = 0;
//		}
//		else if (curr >= 2 && curr < 3) {
//			new_pos[0] = 1 - (curr - 2);
//			new_pos[1] = 1;
//			new_pos[2] = 0;
//		}
//		else { // curr >= 3; remember to cut off as fp precision is an issue
//			new_pos[0] = 0;
//			new_pos[1] = curr > 4 ? 4 : 4 - curr;
//			new_pos[2] = 0;
//		}
//		v->pos = new_pos * boost_factor;
//		curr += step;
//	}
//
//	// Build sparse matrix
//	size_t n = inner_points.size();
//	Eigen::SparseMatrix<double> coeff_mat(n, n);
//	coeff_mat.setZero();
//	Eigen::VectorXd b_vec_x = Eigen::VectorXd::Zero(n);
//	Eigen::VectorXd b_vec_y = Eigen::VectorXd::Zero(n);
//	Eigen::VectorXd b_vec_z = Eigen::VectorXd::Zero(n);
//
//	cout << "coeff mat build start" << endl;
//
//	int current_row = 0;
//	for (auto v : inner_points) {
//		// vidx CERTAINLY follows order (and it's redundant)
//		size_t vidx = inner_points.idx(v);
//		auto adj = v->AdjVertices();
//		size_t degree = v->Degree();
//		for (auto adjv : adj) {
//			// check type
//			if (boundary_points.contains(adjv)) { // this set is usually smaller
//				b_vec_x(current_row) += (1.0f / degree) * adjv->pos[0];
//				b_vec_y(current_row) += (1.0f / degree) * adjv->pos[1];
//				b_vec_z(current_row) += (1.0f / degree) * adjv->pos[2];
//			}
//			else { // inner
//				assert(inner_points.contains(adjv));
//				size_t adjidx = inner_points.idx(adjv);
//				// todo add assert = 0
//				coeff_mat.insert(current_row, adjidx) = -1.0f / degree;
//			}
//		}
//
//		// add itself
//		// todo add assert
//		coeff_mat.insert(current_row, vidx) = 1;
//		current_row++;
//	}
//
//
//	cout << "coeff mat build complete" << endl;
//
//	// Solve
//	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
//
//	cout << "begin makeCompressed()" << endl;
//	coeff_mat.makeCompressed();
//
//	cout << "begin compute()" << endl;
//	solver.compute(coeff_mat);
//	if (solver.info() != Eigen::Success) {
//		cout << "solver: decomposition was not successful." << endl;
//		return;
//	}
//
//	cout << "begin solve() for x" << endl;
//	Eigen::VectorXd res_x = solver.solve(b_vec_x);
//
//	cout << "begin solve() for y" << endl;
//	Eigen::VectorXd res_y = solver.solve(b_vec_y);
//
//	cout << "begin solve() for z" << endl;
//	Eigen::VectorXd res_z = solver.solve(b_vec_z);
//
//	// Update vertex coordinates
//	for (int i = 0; i < n; i++) {
//		// find the corresponding point
//		auto v = inner_points[i];
//		vecf3 new_pos = { res_x(i), res_y(i), res_z(i) }; // works?
//
//		//cout << new_pos << endl;
//		v->pos = new_pos;
//	}
//}

void ASAP::InitLaplacianMatrix() {

}

