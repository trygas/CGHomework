#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>

#include <map>
#include <Eigen/Sparse>

namespace Ubpa {
	class TriMesh;
	class MinSurf;

	// mesh boundary == 1
	class ARAP : public HeapObj {
	public:
		ARAP(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<ARAP> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<ARAP>(triMesh);
		}
	public:
		void Clear();
		bool Init(Ptr<TriMesh> triMesh);

		bool Run();

		void SetDisplay(bool isShow) {
			show = isShow;
		}
		void SetIter(int times) {
			iterTimes = times;
		}
		void SetErrorThreshold(double threshold) {
			errorThreshold = threshold;
		}

	private:
		class V;
		class E;
		class P;
		class V : public TVertex<V, E, P> {
		public:
			vecf3 pos;
		};
		class E : public TEdge<V, E, P> { };
		class P :public TPolygon<V, E, P> { };

	private:
		void Paramaterization();
		void LocalFlatern();
		void LocalUpdate();
		double GlobalUpdate();
		void InitA();
		void Initb();

	private:
		Ptr<TriMesh> triMesh;
		const Ptr<HEMesh<V>> heMesh; // vertice order is same with triMesh

		bool                 show;
		int                  iterTimes = 5;
		double               errorThreshold = 0.01;
		std::map<int, std::map<V*, double> > cotPerPoly;
		std::map<int, std::map<V*, pointf2> > verticesPerPoly;
		std::map<int, Eigen::Matrix<double, 2, 2> > LtMatrix;

		std::vector<pointf2> texcoords;
		Eigen::SparseMatrix<double> A;
		Eigen::MatrixXd b;
		std::vector<Eigen::Triplet<double> > triplets;
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	};
}
