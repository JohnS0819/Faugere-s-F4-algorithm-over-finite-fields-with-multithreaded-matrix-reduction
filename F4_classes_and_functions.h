#include <vector>
#include <set>
#include <iostream>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <algorithm>
#include <Eigen/Dense>
#include <thread>
#include "gfelement.h"



using namespace::shk_galoiscpp;
using namespace boost;

//size of finite field, only valid for prime powers less than 2e+16
const Fint Modulus = 31991;

// dimension of affine space
const int R = 9;
// variable names


//const char variable[R + 1] = { 'x','y','z','s','w','l','\0' };

//const char variable[R + 1] = { 'x','y','z','s','w','l','q','\0' };

//const char variable[R + 1] = { 'x','y','z','s','w','l','q', 'u', '\0' };

const char variable[R + 1] = { 'x','y','z','s','w','l','q', 'u', 'h', '\0' };

int Division_lookup_table[Modulus];

void initialize_division_table() {
	for (auto i = 1; i < Modulus; ++i) {
		Division_lookup_table[i] = inverseModular(i, Modulus);
	}
}

int DivideModular(Fint a, Fint b, Fint mod)
{
	return (a * Division_lookup_table[b]) % mod;
}




class monomial {

private:
	Fint coefficient;
	short exponents[R];
	int degree;
public:

	bool check_power_validity() {
		for (int i = 0; i < R; ++i) {
			if (exponents[i] < 0) {
				return false;
			}
		}
		if (degree > 1000000) {
			return false;
		}
		return true;
	}


	//constructor functions
	//default iniziallization
	monomial() {
		int i;
		for (i = 0; i < R; ++i) {
			exponents[i] = 0;
		}
		coefficient = 0;
		degree = 0;
	}
  
	monomial(const char* X) {
		int i = -1;
		bool isnumber = false;
		int counter;
		int track;
		int retrack;
		int power;
		int j;
		for (j = 0; j < R; j++) {
			exponents[j] = 0;
		}
		if (*X == '-') {
			coefficient = -1;
			++X;
		}
		else if (*X == '+') {
			coefficient = 1;
			++X;
		}
		else {
			coefficient = 1;
		}
		while (*X != '\0' and *X != '+' and *X != '-') {
			for (j = 0; j < R; j++) {
				if (variable[j] == *X) {
					if (isnumber) {
						exponents[i] = 1;
					}
					i = j;
					isnumber = true;
					break;
				}
			}
			if (*X <= 57 and *X >= 48) {
				track = 0;
				counter = 0;
				power = 1;
				while (*X <= 57 and *X >= 48) {
					++track;
					++X;
				}
				retrack = track;
				while (track != 0) {
					--X;
					counter = counter + ((*X - 48) * power);
					power = power * 10;
					--track;
				}
				X = X + retrack;
				if (i != -1) {
					exponents[i] = counter;
				}
				else {
					coefficient = coefficient * counter;
				}
				isnumber = false;
			}
			else {
				++X;
			}
		}
		if (isnumber) {
			exponents[i] = 1;
		}
		int deg = 0;
		for (j = 0; j < R; ++j) {
			deg = deg + exponents[j];
		}
		degree = deg;

	}
	//nonstring based initiallization
	monomial(Fint X, int Y[R]) {
		coefficient = X;
		int i;
		int deg = 0;
		for (i = 0; i < R; ++i) {
			exponents[i] = Y[i];
			deg = deg + exponents[i];
		}
		degree = deg;
	}
	//does not initialize coefficient;
	monomial(int Y[R]) {
		degree = 0;
		for (int i = 0; i < R; ++i) {
			exponents[i] = Y[i];
			degree += Y[i];
		}
		coefficient = 1;
	}

	//copy constructor
	monomial(const monomial& pnt) {
		coefficient = pnt.coefficient;
		degree = pnt.degree;
		int i;
		for (i = 0; i < R; ++i) {
			exponents[i] = pnt.exponents[i];

		}

	}

	//switches coefficients of monomial with another monomial
	monomial(const monomial& pnt, Fint a) {
		coefficient = a;
		degree = pnt.degree;
		int i;
		for (i = 0; i < R; ++i) {
			exponents[i] = pnt.exponents[i];

		}

	}

	//lcm constructor
	monomial(monomial a, monomial b) {
		coefficient = 1;
		int i;
		degree = 0;
		for (i = 0; i < R; ++i) {
			exponents[i] = std::max(a.exponents[i], b.exponents[i]);
			degree = degree + exponents[i];
		}
	}

	//operator definitions

	bool operator==(const monomial a) const {
		int i;
		for (i = 0; i < R; ++i) {
			if (a.exponents[i] != exponents[i]) {
				return false;
			}
		}
		return true;
	}

	monomial operator*(const monomial a) const {
		int i;
		monomial X;
		for (i = 0; i < R; ++i) {
			X.exponents[i] = a.exponents[i] + exponents[i];
		}
		X.degree = a.degree + degree;
		X.coefficient = multiplyModular(a.coefficient, coefficient, Modulus);
		return X;
	}


	monomial operator*(const Fint a) const {
		monomial b = *this;
		b.coefficient = multiplyModular(a, b.coefficient, Modulus);
		return b;
	}

	//does not check for negative exponents also does not change coefficients
	monomial operator/(const monomial a) {
		int diff[R];
		int i;
		for (i = 0; i < R; ++i) {
			diff[i] = exponents[i] - a.exponents[i];
		}
		return monomial(1, diff);
	}

	bool operator!=(const monomial a) const {
		return (1 - (*this == a));
	}

	friend std::ostream& operator<<(std::ostream& output, const monomial& a) {
		int i;
		output << "(" << a.coefficient << ")";
		for (i = 0; i < R; ++i) {
			switch (a.exponents[i]) {
			case 0:
				break;
			case 1:
				output << "*" << variable[i];
				break;
			default:
				output << "*" << variable[i] << "^" << a.exponents[i];
				break;
			}
		}
		return output;
	}

	// monomial ordering operations
	// only supports using degrevlex

	bool operator<(monomial a) const {
		if (degree < a.degree) {
			return true;
		}
		if (degree > a.degree) {
			return false;
		}
		int i;
		int j;
		for (i = R - 1; i > -1; --i) {
			j = exponents[i] - a.exponents[i];
			if (j != 0) {
				if (j > 0) {
					return true;
				}
				return false;
			}
		}
		return false;
	}



	bool operator>(monomial a) const {
		return (a < *this);
	}

	bool operator<=(const monomial a) const {
		return (1 - (*this > a));
	}

	bool operator>=(monomial a) {
		return (1 - (*this < a));

	}

	friend class Crit_pair;
	friend bool is_reducible(monomial a, monomial b);
	friend class polynomial;
	friend bool is_properly_divisible(monomial a, monomial b);
	friend Fint return_coefficients(monomial a);
	friend bool M_F_criterion(monomial a, monomial b, int label_of_a, int label_of_b, monomial c);
};

class compressed_monomial {
public:
	int position;
	Fint value;

	compressed_monomial() {

	}

	compressed_monomial(Fint a, int b) {
		position = b;
		value = a;
	}

	compressed_monomial(const compressed_monomial& a) {
		position = a.position;
		value = a.value;
	}
	compressed_monomial(int a[2]) {
		value = a[0];
		position = a[1];
	}


	bool operator<(compressed_monomial a) const {
		return (position < a.position);
	}
	bool operator>(compressed_monomial a) const {
		return (position > a.position);
	}
	bool operator==(compressed_monomial a) const {
		return (position == a.position);
	}

	bool operator!=(compressed_monomial a) const {
		return (position != a.position);
	}

	friend std::ostream& operator<<(std::ostream& output, const compressed_monomial& a) {
		int i;
		output << '(' << a.value << ',' << a.position << ')';
		return output;
	}

};
bool compare_compressed_monomial(compressed_monomial a, compressed_monomial b) {
	return b < a;
}

//used during sparse matrix operations
class compressed_polynomial {
public:
	std::vector<compressed_monomial> data;
	compressed_polynomial() {

	}

	compressed_polynomial(const compressed_polynomial& a) {
		data = a.data;
	}

	compressed_polynomial(std::vector<compressed_monomial> a) {
		data = a;
	}

	friend std::ostream& operator<<(std::ostream& output, const compressed_polynomial& a) {
		int i;
		for (auto i = a.data.begin(); i != a.data.end(); ++i) {
			output << *i << ' ';
		}
		return output;
	}


	compressed_polynomial operator*(Fint a) {
		std::vector<compressed_monomial> b = data;
		for (auto i = b.begin(); i != b.end(); ++i) {
			(*i).value = multiplyModular((*i).value, a, Modulus);
		}
		return compressed_polynomial(b);
	}

	void operator-=(const compressed_polynomial a) {
		for (auto i = a.data.begin(); i != a.data.end(); ++i) {
			auto j = lower_bound(data.begin(), data.end(), *i);
			if (j != data.end() && *j == *i) {
				if (j->value == i->value) {
					data.erase(j);
				}
				else {
					j->value = subtractModular(j->value, i->value, Modulus);
				}
			}
			else {
				data.insert(j, compressed_monomial(subtractModular(0, i->value, Modulus), i->position));
			}
		}
	}


};


//Part of the Gebauer moller criterion
bool M_F_criterion(monomial a, monomial b, int label_of_a, int label_of_b, monomial c) {
	for (int i = 0; i < R; ++i) {
		if (a.exponents[i] < b.exponents[i]) {
			return false;
		}
	}
	if (label_of_b < label_of_a) {
		return true;
	}
	for (int i = 0; i < R; ++i) {
		if (a.exponents[i] != max(b.exponents[i], c.exponents[i])) {
			return true;
		}
	}
	return false;
}




//sorts in descending rather than ascending order
bool compare(const monomial a1, const monomial a2) {
	return (a2 < a1);
}

//dummy functions to transfer friendship between functions that need monomial data but take polynomial/derivative data
Fint return_coefficients(monomial a) {
	return a.coefficient;
}


class polynomial {
private:
	int Length;
	// is exclusively used by addition and subtraction operations
	polynomial(const polynomial& a, int j) {
		Length = a.Length - 1;
		terms.reserve(Length);
		int i;
		for (i = 0; i < j; ++i) {
			terms.push_back(a.terms[i]);
		}
		for (i = j; i < Length; ++i) {
			terms.push_back(a.terms[i + 1]);
		}



	}

	std::vector<monomial> terms;

	void normalize() {
		if (terms[0].coefficient == 1) {
			return;
		}
		*this = *this * DivideModular(1, (terms[0].coefficient), Modulus);
		return;
	}


public:
	//returns pointer to start of terms
	auto begin() {
		return terms.begin();
	}
	//returns pointer to end of terms
	auto end() {
		return terms.end();
	}

	auto rbegin() {
		return terms.rbegin();
	}

	auto rend() {
		return terms.rbegin();
	}



	//automatically sorts polynomial but does not collect terms
	polynomial(const char* X) {
		int L = 1;
		const char* Y = X;
		if (*Y == '-') {
			L = 0;
		}
		while (*Y != '\0') {
			if (*Y == '+' or *Y == '-') {
				++L;
			}
			++Y;
		}
		Length = L;
		Y = X;
		terms.reserve(L);
		terms.push_back(monomial(Y));
		++Y;
		while (*Y != '\0') {
			if (*Y == '+' or *Y == '-') {
				terms.push_back(monomial(Y));
			}
			++Y;
			std::sort(terms.begin(), terms.end(), compare);
		}

	}

	polynomial() {
		Length = 0;
		terms = {};
	}

	//does not automatically sort polynomial by ordering
	//should not be used for user defined polynomials
	polynomial(std::vector<monomial> X) {
		Length = X.size();
		terms.reserve(Length);
		int i;
		for (i = 0; i < Length; ++i) {
			terms.push_back(X[i]);
		}


	}

	polynomial(monomial x) {
		Length = 1;
		terms = { x };
	}

	// assumes polynomial is currently sorted
	monomial lt() {
		return terms[0];
	}

	Fint lc() {
		return terms[0].coefficient;
	}


	friend std::ostream& operator<<(std::ostream& output, const polynomial& a) {
		int i;
		for (i = 0; i < a.Length; ++i) {
			if (i != (a.Length - 1)) {
				output << a.terms[i] << " + ";
			}
			else {
				output << a.terms[i];
			}
		}
		return output;
	}

	Fint has_term(monomial j) {
		auto i = std::lower_bound(terms.begin(), terms.end(), j, compare);
		if (i == terms.end()) {
			return 0;
		}
		if (*i == j) {
			return (*i).coefficient;
		}
		return 0;

	}

	polynomial(const polynomial& pnt) {
		terms = pnt.terms;
		Length = pnt.Length;
	}

	polynomial operator*(monomial a) {
		int i;
		std::vector<monomial> new_monomial;
		new_monomial.reserve(terms.size());
		for (i = 0; i < terms.size(); ++i) {
			new_monomial.push_back(terms[i] * a);
		}
		return polynomial(new_monomial);
	}
	polynomial operator*(Fint a) {
		int i;
		std::vector<monomial> new_monomial;
		new_monomial.reserve(terms.size());
		for (i = 0; i < terms.size(); ++i) {
			new_monomial.push_back(terms[i] * a);
		}
		return polynomial(new_monomial);
	}

	polynomial operator+(monomial a) {
		polynomial copy = *this;
		copy += a;
		return copy;
	}

	polynomial operator+(polynomial a) {
		int i;
		polynomial output = *this;
		for (i = 0; i < a.Length; ++i) {
			output = (output + a.terms[i]);
		}
		return output;


	}

	polynomial operator-(monomial a) {
		return(*this + (a * -1));
	}




	void operator+=(const monomial a) {
		auto i = lower_bound(terms.begin(), terms.end(), a, compare);
		if (i != terms.end() && *i == a) {
			if (i->coefficient == multiplyModular(-1, a.coefficient,Modulus)) {
				terms.erase(i);
				--Length;
				return;

			}
			else {
				i->coefficient = addModular(i->coefficient, a.coefficient, Modulus);
				return;
			}
		}
		else {
			terms.insert(i,a);
			++Length;
			return;
		}
	}

	

	void operator-=(const monomial a) {
		auto i = lower_bound(terms.begin(), terms.end(), a,compare);
		if (i != terms.end() && *i == a) {
			if (i->coefficient == a.coefficient) {
				terms.erase(i);
				--Length;
				return;

			}
			else {
				i->coefficient = subtractModular(i->coefficient, a.coefficient, Modulus);
				return;
			}
		}
		else {
			terms.insert(i,a * -1);
			++Length;
			return;
		}
	}

	void operator+=(const polynomial a) {
		for (int i = 0; i < a.Length; ++i) {
			(*this) += (a.terms[i]);
		}
		return;
	}

	void operator-=(const polynomial a) {
		for (int i = 0; i < a.Length; ++i) {
			(*this) -= a.terms[i];
		}
		return;
	}


	bool operator==(const polynomial a) const {
		if (a.Length != Length) {
			return false;
		}
		int i;
		Fint comp = DivideModular(terms[0].coefficient, a.terms[0].coefficient, Modulus);
		for (i = 0; i < Length; ++i) {
			if (terms[i] != a.terms[i]) {
				return false;
			}
			if (terms[i].coefficient != multiplyModular(a.terms[i].coefficient, comp, Modulus)) {
				return false;
			}
		}
		return true;

	}

	bool operator!=(const polynomial a) const {
		return (1 - (*this == a));
	}


	//Note this is not a total ordering, rather only a partial ordering on leading terms
	bool operator<(const polynomial a) const {
		return (terms[0] < a.terms[0]);
	}

	bool operator<(const monomial a) const {
		return (terms[0] < a);
	}
	bool operator>(const monomial a) const {
		return (terms[0] > a);
	}



	friend class Crit_pair;

	friend std::vector<polynomial> rowechelon(std::vector<polynomial> f);

	friend std::vector<polynomial> rowechelonv2(std::vector<polynomial> f);

	friend std::vector<polynomial> rowechelonv3(std::vector<polynomial> f);


	friend void single_threaded_reduction(std::vector<polynomial>& F);

	friend std::vector<polynomial> reduce_basis(std::vector<polynomial> F);

};


bool compare_polynomials(const polynomial a1, const polynomial a2) {
	return (a2 < a1);
}

//is use only for left and right operations
struct polynomial_pair {
	monomial A;
	polynomial B;
	polynomial_pair(monomial a, polynomial b) {
		A = a;
		B = b;
	}
	polynomial flatten() {
		return (B * A);
	}
	bool operator==(const polynomial_pair a) const {
		if (A == a.A && B == a.B) {
			return true;
		}
		return false;
	}
};


class Crit_pair {
private:
	monomial lcm;
	polynomial polys[2];
	monomial leading_terms[2];
public:
	Crit_pair(polynomial a, polynomial b) {
		lcm = monomial(a.lt(), b.lt());
		polys[0] = a; polys[1] = b;
		leading_terms[0] = lcm / a.lt(); leading_terms[1] = lcm / b.lt();
	}

	//creates a trivial pair with the goal of searching more efficiently
	Crit_pair(monomial a) {
		lcm = a;
		polys[0] = polynomial();
		polys[1] = polynomial();
	}

	int pair_degree() const {
		return lcm.degree;
	}

	polynomial_pair left() const {
		return polynomial_pair(leading_terms[0], polys[0]);
	}

	polynomial_pair right() const {
		return polynomial_pair(leading_terms[1], polys[1]);
	}

	bool operator==(const Crit_pair a) {
		if ((polys[0] == a.polys[0]) and (polys[1] == a.polys[1])) {
			return true;
		}
		if ((polys[1] == a.polys[0]) and (polys[0] == a.polys[1])) {
			return true;
		}
		return false;
	}
	bool operator!=(const Crit_pair a) {
		return (1 - (*this == a));
	}

	//will be removed once a topological sorting method has been implented and replaced with a partial ordering that compares
	//only degree
	bool operator<(const Crit_pair a) const {
		if ((*this).pair_degree() != a.pair_degree()) {
			return (*this).pair_degree() < a.pair_degree();
		}
		if (polys[0].Length != a.polys[0].Length) {
			return (polys[0].Length < a.polys[0].Length);
		}
		if (polys[1].Length != a.polys[1].Length) {
			return (polys[1].Length < a.polys[1].Length);
		}
		for (int i = 0; i < polys[0].Length; ++i) {
			if (polys[0].terms[i].coefficient != a.polys[0].terms[i].coefficient) {
				return (polys[0].terms[i].coefficient < a.polys[0].terms[i].coefficient);
			}
			if (polys[0].terms[i] != a.polys[0].terms[i]) {
				return (polys[0].terms[i] < a.polys[0].terms[i]);
			}
		}
		for (int i = 0; i < polys[1].Length; ++i) {
			if (polys[1].terms[i].coefficient != a.polys[1].terms[i].coefficient) {
				return (polys[1].terms[i].coefficient < a.polys[1].terms[i].coefficient);
			}
			if (polys[1].terms[i] != a.polys[1].terms[i]) {
				return (polys[1].terms[i] < a.polys[1].terms[i]);
			}
		}
		return false;

	}


	friend std::ostream& operator<<(std::ostream& output, const Crit_pair& a) {
		output << a.polys[0] << "\t" << a.leading_terms[0] << "\n"
			<< a.polys[1] << "\t" << a.leading_terms[1] << "\n" << a.lcm;
		return output;
	}

	friend bool B_k_criterion(monomial k, Crit_pair a);


};


bool B_k_criterion(monomial k, Crit_pair a) {
	if (!is_reducible(a.lcm, k)) {
		return false;
	}
	if (monomial(a.polys[0].lt(), k) == a.lcm) {
		return false;
	}
	if (monomial(a.polys[1].lt(), k) == a.lcm) {
		return false;
	}
	return true;

}



void Gauss_jordan_elimination(Eigen::MatrixX<Fint>& M, int size1, int size2);

bool is_reducible(monomial a, monomial b) {
	for (int i = 0; i < R; ++i) {
		if ((a.exponents[i] - b.exponents[i]) < 0) {
			return false;
		}
	}
	return true;
}

//produces a minimal representation of a groebner basis, note this representation is not unique
std::vector<polynomial> minimialize_basis(std::vector<polynomial> input) {
	std::vector<monomial> terms;
	terms.reserve(input.size());
	for (auto i = input.begin(); i != input.end(); ++i) {
		terms.push_back((*i).lt());
	}
	std::sort(terms.begin(), terms.end());
	std::sort(input.begin(), input.end());
	for (auto i = terms.begin(); i != terms.end(); ++i) {
		for (auto j = terms.begin(); j != i; ++j) {
			if (is_reducible(*i, *j)) {
				j = i - 1;
				input.erase(input.begin() + (i - terms.begin()));
				terms.erase(i);
				i = j;
				break;
			}
		}
	}
	return input;
}

//produces a unique reduced basis from an already minimal one, note this implementation is highly inefficient and shouldn't be used in affine spaces of dimension greater than 9
std::vector<polynomial> reduce_basis(std::vector<polynomial> F) {
	std::sort(F.begin(), F.end(), compare_polynomials);
	for (auto i = F.begin(); i != F.end(); ++i) {
		(*i).normalize();
	}
	std::vector<monomial> terms;
	terms.reserve(F.size());
	for (auto i = F.begin(); i != F.end(); ++i) {
		terms.push_back((*i).lt());
	}
	std::sort(terms.begin(), terms.end(), compare);

	for (auto i = F.begin(); i != F.end(); ++i) {
		for (auto j = (*i).begin() + 1; j != (*i).end(); ++j) {
			for (auto k = terms.begin() + (i - F.begin()); k != terms.end(); ++k) {
				if (is_reducible(*j, *k)) {
					/*std::cout << (*j) << "\tis reducible by \t" << (*k) << std::endl;
					std::cout << (*i).lt() <<std::endl;*/
					int position = (j - (*i).begin());
					(*i) -= (*(F.begin() + (k - terms.begin())) * ((*j) / (*k))) * return_coefficients(*j);
					j = (*i).begin() + position - 1;
				}
			}
		}
	}
	return F;



}



bool is_properly_divisible(monomial a, monomial b) {
	for (int i = 0; i < R; ++i) {
		if ((a.exponents[i] - b.exponents[i]) < 0) {
			return false;
		}
	}
	for (int i = 0; i < R; ++i) {
		if (a.exponents[i] > b.exponents[i]) {
			return true;
		}
	}
	return false;
}

//This implementation is very naive, however because the vast majority of reduction was already done inside the sparse reduction step,
//The performance gains of a more efficient implentation are negligible 
void Gauss_jordan_elimination(Eigen::MatrixX<Fint>& M, int size1, int size2) {
	using namespace Eigen;
	int r = 0;
	bool piv_found;
	for (int i = 0; i < size2; ++i) {
		piv_found = false;
		for (int j = r; j < size1; ++j) {
			if (M(j, i) != 0) {
				M.row(j).swap(M.row(r));
				if (M(r, i) != 1) {
					for (auto k = size2 - 1; k >= i; --k) {
						M(r, k) = DivideModular(M(r, k), M(r, i), Modulus);
					}
				}
				piv_found = true;
				break;
			}
		}
		if (piv_found) {
			for (int j = 0; j < r; ++j) {
				if (M(j, i) != 0) {
					for (auto k = size2 - 1; k >= i; --k) {
						M(j, k) = subtractModular(M(j, k), multiplyModular(M(r, k), M(j, i), Modulus), Modulus);
					}
				}
			}
			for (int j = r + 1; j < size1; ++j) {
				if (M(j, i) != 0) {
					for (auto k = size2 - 1; k >= i; --k) {
						M(j, k) = subtractModular(M(j, k), multiplyModular(M(r, k), M(j, i), Modulus), Modulus);
					}
				}
			}
			++r;
		}
	}
	return;
}



void reduceonrangev3(std::vector<compressed_polynomial>::iterator first, std::vector<compressed_polynomial>::iterator last, std::vector < std::vector<compressed_polynomial>::iterator> leading_terms, std::vector<compressed_monomial> leading_monomials) {
	for (auto i = first; i != last; ++i) {
		bool test = true;
		while (test) {
			test = false;
			for (auto j = (*i).data.begin(); j != (*i).data.end(); ++j) {

				//this line doesnt work
				auto k = std::lower_bound(leading_monomials.begin(), leading_monomials.end(), *j);

				if (k != leading_monomials.end() && *k == *j) {
					*i -= **(leading_terms.begin() + (k - leading_monomials.begin())) * (*j).value;
					test = true;
					break;
				}
			}
		}

	}
}



void multi_threaded_reductionv3(std::vector<compressed_polynomial>& F, std::vector<std::vector<compressed_polynomial>::iterator>& leading_term_polynomials, std::vector<compressed_monomial> leading_monomials, unsigned int current_thread, unsigned int thread_count) {
	unsigned int iterations;
	for (iterations = current_thread; (iterations + 1) < leading_term_polynomials.size(); iterations += thread_count) {
		reduceonrangev3(leading_term_polynomials[iterations] + 1, leading_term_polynomials[iterations + 1], leading_term_polynomials, leading_monomials);
	}
	return;
}


