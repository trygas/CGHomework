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
	class ASAP : public HeapObj {
	public:
		ASAP(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<ASAP> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<ASAP>(triMesh);
		}
	public:
		void Clear();
		bool Init(Ptr<TriMesh> triMesh);

		bool Run();

		void SetDisplay(bool isShow) {
			show = isShow;
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
		void InitLaplacianMatrix();
		void InitLt();
		void Initu();
		
	private:
		Ptr<TriMesh> triMesh;
		const Ptr<HEMesh<V>> heMesh; // vertice order is same with triMesh

		bool                 show;
		std::map<int, std::map<V*, double> > cotPerPoly;
		std::map<int, std::map<V*, pointf2> > verticesPerPoly;

		std::vector<pointf2> texcoords;
		Eigen::SparseMatrix<double> A;
		Eigen::VectorXd b;
		std::vector<Eigen::Triplet<double> > triplets;
	};
}
