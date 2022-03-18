#include "F4_classes_and_functions.h"

using namespace::shk_galoiscpp;
using namespace boost;

std::vector<polynomial> symbolic_preprocessingv2(std::vector<polynomial_pair> L, std::vector<polynomial> G) {
  
	std::set<monomial> terms;
	std::vector<polynomial> t;
	int size = L.size();
	int i = 0;
	for (i = 0; i < size; ++i) {
		t.push_back(L[i].flatten());
	}
	std::set<monomial> Done;
	for (i = 0; i < size; ++i) {
		Done.insert(t[i].lt());
	}
	for (i = 0; i < size; ++i) {
		terms.insert(t[i].begin(), t[i].end());
	}
	std::set<monomial> comp;
	set_difference(terms, Done, std::inserter(comp, comp.begin()));

	monomial m;
	monomial M;
	std::set<monomial> history;
	while (comp.size() != 0) {
		m = *(--comp.end());
		history.insert(m);
		comp.erase(--comp.end());
		for (auto j = G.begin(); j != G.end(); ++j) {
			if (is_reducible(m, (*j).lt())) {
				M = (m / (*j).lt());
				for (auto g = (*j).begin(); g != (*j).end(); ++g) {
					if ((history.find((*g) * M)) != history.end()) {
						goto breakpoint;
					}
					comp.insert((*g) * M);
				breakpoint:;
				}
				t.push_back((*j) * M);
				break;
			}
		}
	}
	int count = 0;
	for (auto i = t.begin(); i != t.end(); ++i) {
		for (auto j = i->begin(); j != i->end(); ++j) {
			++count;
		}
	}
	std::sort(t.begin(), t.end(), compare_polynomials);
	auto end = std::unique(t.begin(), t.end());
	return t;
}


std::vector<polynomial> rowechelonv3(std::vector<polynomial> f) {
	std::set<monomial, decltype(&compare)> temp_stores(&compare);

	for (auto i = f.begin(); i != f.end(); ++i) {
		for (auto j = ((*i).begin()); j != (*i).end(); ++j) {
			temp_stores.insert(*j);
		}
	}
	//conversion of key to contiguous memory
	std::vector<monomial> power;
	power.reserve(temp_stores.size());
	for (auto i = temp_stores.begin(); i != temp_stores.end(); ++i) {
		power.push_back(*i);
	}
	temp_stores.clear();
	std::vector<compressed_polynomial> sparse_matrix;
	sparse_matrix.reserve(f.size());
	for (auto i = f.begin(); i != f.end(); ++i) {
		auto previous = power.begin();
		compressed_polynomial temp_poly;
		for (auto j = (*i).begin(); j != (*i).end(); ++j) {
			previous = std::lower_bound(previous, power.end(), *j, compare);
			temp_poly.data.push_back(compressed_monomial(return_coefficients(*j), previous - power.begin()));
		}
		sparse_matrix.push_back(temp_poly);
	}
	std::vector < std::vector<compressed_polynomial>::iterator> leading_terms;
	leading_terms.push_back(sparse_matrix.begin());
	std::vector<compressed_monomial> leading_monomials;
	leading_monomials.push_back(*(*leading_terms.begin())->data.begin());
	for (auto i = sparse_matrix.begin() + 1; i != sparse_matrix.end(); ++i) {
		if (*(*i).data.begin() != (*(*(i - 1)).data.begin())) {
			leading_terms.push_back(i);
			leading_monomials.push_back(*(*i).data.begin());
		}
	}
	leading_terms.push_back(sparse_matrix.end());
	unsigned int thread_count;
	if (std::thread::hardware_concurrency() < leading_terms.size()) {
		thread_count = std::thread::hardware_concurrency() - 1;
	}
	else {
		thread_count = leading_terms.size() - 1;
	}
	std::vector<std::thread> threads;

	//reducing sparse matrix

	for (unsigned int i = 0; i < thread_count; ++i) {
		threads.push_back(std::thread(multi_threaded_reductionv3, std::ref(sparse_matrix), std::ref(leading_terms), std::ref(leading_monomials), i, thread_count));
	}
	for (auto i = threads.begin(); i != threads.end(); ++i) {
		(*i).join();
	}


	//removing leading term and empty polynomials
	std::vector<compressed_polynomial> J;
	for (auto i = leading_terms.begin(); i != leading_terms.end() - 1; ++i) {
		for (auto j = (*i) + 1; j != *(i + 1); ++j) {
			if ((*j).data.size() != 0) {
				J.push_back(*j);
			}
		}
	}


	//decompression of polynomial
	std::vector<polynomial> F;
	for (auto i = J.begin(); i != J.end(); ++i) {
		std::vector<monomial> temp_term;
		for (auto j = i->data.begin(); j != i->data.end(); ++j) {
			temp_term.push_back(monomial(power[j->position], j->value));
		}
		F.push_back(polynomial(temp_term));
	}
	std::set<monomial, decltype(&compare)> temp_store(&compare);

	for (auto i = F.begin(); i != F.end(); ++i) {
		for (auto j = ((*i).begin()); j != (*i).end(); ++j) {
			temp_store.insert(*j);
		}
	}
	//conversion of key to contiguous memory
	std::vector<monomial> powers;
	powers.reserve(temp_store.size());
	for (auto i = temp_store.begin(); i != temp_store.end(); ++i) {
		powers.push_back(*i);
	}

	//creation of Dense matrix
	using namespace Eigen;
	MatrixX<Fint> M = MatrixX<Fint>::Zero(F.size(), powers.size());
	for (auto i = F.begin(); i != F.end(); ++i) {
		auto previous = powers.begin();
		for (auto j = (*i).begin(); j != (*i).end(); ++j) {
			previous = std::lower_bound(previous, powers.end(), *j, compare);
			M(i - F.begin(), previous - powers.begin()) = return_coefficients(*j);
		}
	}


	Gauss_jordan_elimination(M, F.size(), powers.size());

	//rebuilding polynomial array
	int k = 0;
	bool found_entry;
	int j;
	int i;
	std::vector<polynomial> output;
	std::vector<monomial> temp;

	for (i = 0; i < F.size(); ++i) {
		found_entry = false;
		for (j = k; j < powers.size(); ++j) {
			if (M(i, j) != 0) {
				found_entry = true;
				k = j;
				goto foundvalue;
			}
		}
		if (!found_entry) {
			break;
		}
	mainloop:
		temp.clear();

	}
	std::cout << std::endl;



	return output;

foundvalue:
	for (; j < powers.size(); ++j) {
		if (M(i, j) != 0) {
			temp.push_back(monomial(powers[j], M(i, j)));
		}
	}
	output.push_back(polynomial(temp));
	goto mainloop;
}


std::set<Crit_pair> Create_pairs(const std::vector<polynomial> input) {
	std::set<Crit_pair> output;
	for (int i = 0; i < input.size(); ++i) {
		for (int j = 0; j < i; ++j) {
			output.emplace(input[i], input[j]);
		}
	}
	return output;
}

std::vector<polynomial> reduction(std::vector<polynomial_pair> L, std::vector<polynomial>  G) {
	return (rowechelonv3(symbolic_preprocessingv2(L, G)));
}


std::vector<polynomial> F4_2(std::vector<polynomial> input) {
	std::sort(input.begin(), input.end());
	std::vector<polynomial> G;
	std::set<Crit_pair> P;
	std::set<pair<monomial, int>> labeled_monomials;
	int label = 0;
	std::vector<Crit_pair> temp_store;
	//implementation of the M and F Gebauer Moller Criterion
	for (auto q = input.begin(); q != input.end(); ++q) {
		for (auto j = input.begin(); j != q; ++j) {
			bool test = true;
			monomial lcm = monomial((*j).lt(), (*q).lt());
			auto iterator = labeled_monomials.lower_bound({ lcm,label });
			for (auto w = labeled_monomials.begin(); w != iterator; ++w) {
				if (M_F_criterion(lcm, (*w).first, j - input.begin(), (*w).second, (*q).lt())) {
					test = false;
					break;
				}
			}
			if (test) {
				temp_store.push_back(Crit_pair(*j, *q));
			}
		}
		//erasing redundant pairs
		{auto U = P.begin(); while (U != P.end()) {
			if (B_k_criterion((*q).lt(), *U)) {
				++U;
				P.erase(std::prev(U));
				continue;
			}
			++U;
		}}

		for (auto i = temp_store.begin(); i != temp_store.end(); ++i) {
			P.insert(*i);
		}
		temp_store.clear();
		G.push_back(*q);
		labeled_monomials.insert({ (*q).lt(),label });
		++label;
	}

	std::vector<polynomial_pair> L_d;
	std::vector<polynomial> F_d;
	auto i = P.begin();
	int degree;
	while (P.size() != 0) {

		//selection criteria
		degree = (*P.begin()).pair_degree();
		for (i = P.begin(); i != P.end(); ++i) {
			if ((*i).pair_degree() == degree) {
				L_d.push_back((*i).left()); L_d.push_back((*i).right());
			}
			else {
				break;
			}
		}
		--i;
		if (i == P.begin()) {
			P.erase(P.begin());
		}
		else {
			P.erase(P.begin(), i);
		}


		F_d = reduction(L_d, G);

		L_d.clear();

		//implementation of the Gebauer Moller critereon
		for (auto q = F_d.begin(); q != F_d.end(); ++q) {
			for (auto j = G.begin(); j != G.end(); ++j) {
				bool test = true;
				monomial lcm = monomial((*j).lt(), (*q).lt());
				auto iterator = labeled_monomials.lower_bound({ lcm,label });
				for (auto w = labeled_monomials.begin(); w != iterator; ++w) {
					if (!M_F_criterion(lcm, (*w).first, j - G.begin(), (*w).second, (*q).lt())) {
						test = false;
						break;
					}
					if (test) {
						temp_store.push_back(Crit_pair(*j, *q));
					}
				}

			}

			//erasing reduntant pairs
			{auto U = P.begin(); while (U != P.end()) {
				if (B_k_criterion((*q).lt(), *U)) {
					++U;
					P.erase(std::prev(U));
					continue;
				}
				++U;
			}}
			for (auto U = temp_store.begin(); U != temp_store.end(); ++U) {
				P.insert(*U);
			}
			temp_store.clear();
			G.push_back(*q);
			labeled_monomials.insert({ (*q).lt(),label });
			++label;
		}
	}
	return G;
}
