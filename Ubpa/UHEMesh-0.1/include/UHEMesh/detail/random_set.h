#pragma once

#include <vector>
#include <unordered_map>

namespace Ubpa {
	template<typename T, class Hasher = std::hash<T>>
	class random_set {
	public:
		void insert(const T& e) {
			assert(Tmap.find(e) == Tmap.end());
			Tmap[e] = Tvec.size();
			Tvec.push_back(e);
		}

		void erase(const T& e) {
			auto target = Tmap.find(e);
			assert(target != Tmap.end());
			auto eIdx = target->second;
			if (eIdx != Tvec.size() - 1) {
				Tmap[Tvec.back()] = eIdx;
				Tvec[eIdx] = Tvec.back();
			}
			Tvec.pop_back();
			Tmap.erase(target);
		}

		auto begin() { return Tvec.begin(); }
		auto begin() const { return Tvec.begin(); }
		auto end() { return Tvec.end(); }
		auto end() const { return Tvec.end(); }

		T& operator[](size_t i) { return Tvec[i]; }
		const T& operator[](size_t i) const { return Tvec[i]; }

		size_t size() const { return Tvec.size(); }
		void reserve(size_t n) {
			Tvec.reserve(n);
			Tmap.reserve(n);
		}
		void clear() {
			Tvec.clear();
			Tmap.clear();
		}

		const std::vector<T>& vec() const { return Tvec; }

		size_t idx(const T& e) const {
			assert(Tmap.find(e) != Tmap.end());
			return Tmap.find(e)->second;
		}

		bool contains(const T& e) const { return Tmap.find(e) != Tmap.end(); }

		bool empty() const { return Tvec.empty(); }
	private:
		std::unordered_map<T, size_t, Hasher> Tmap;
		std::vector<T> Tvec;
	};
}
