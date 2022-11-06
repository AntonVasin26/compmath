#include <iostream>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include "polinom.hpp"


typedef long double ftype;
//constexpr std::pair<bool, int> max_index(long double n)
//{
//	std::pair<bool, int> answer(true, 0);
//	if (std::numeric_limits<long double>::is_integer == true)
//	{
//		int n_1 = static_cast<int>(n);
//		answer.second = 2*n_1;
//	}
//	else
//	{
//		answer.first = false;
//	}
//	return answer;
//}

int max_index(ftype n)
{
		const int n_1 = std::round(n);
	return 2*n_1;
}

void print_v(std::vector<ftype> vec)
{
	std::cout << "\nPrint index\n";

	for (auto elem : vec)
	{
		std::cout << elem << '\t';
	}
	std::cout << '\n';
}

std::pair<ftype, ftype> localization(std::vector<ftype> v)
{
	int n = v.size() - 1;
	auto B = abs(*std::max_element(v.begin()+1, v.end(), [](ftype a, ftype b) { return abs(a) < abs(b); }));
	auto A = abs(*std::max_element(v.begin(), v.end()-1, [](ftype a, ftype b) { return abs(a) < abs(b); }));
	auto a_0 = abs(v[n]);
	auto a_n = abs(v[0]);
	auto min = a_n / (a_n + B);
	auto max = 1 + A / a_0;
	return std::pair<ftype, ftype> (min, max);
}

ftype func(ftype x, std::vector<ftype> v)
{
	ftype y = 0.0;
	int n = 0;
	for (auto elm : v)
	{
		y += elm * powl(x, n++);
	}
	return y;
}

//typedef std::vector<ftype> polinom;
//typedef std::pair<polinom, polinom> pair_f;

//std::vector<ftype> minus(polinom p1, polinom p2)
//{
//	auto n1 = p1.size();
//	auto n2 = p2.size();
//	auto k = std::max(n1, n2);
//	for ()
//}
//
//pair_f polidiv(polinom p1, polinom p2)
//{
//	pair_f answer;
//	auto n0 = p1.size();
//	auto n1 = p2.size();
//	if (n1 = 0)
//	{
//		std::cout << "\n\npolidiv - invalid value\n\n";
//		return answer;
//	}
//	if (n1 > n0)
//	{
//		answer.second = p1;
//		return answer;
//	}
//	ftype c1 = p1[0]/p2[0]
//
//	
//}

setm::Polynomial<ftype> derivative(setm::Polynomial<ftype> p0)
{
	std::vector<ftype> der_p0;
	int m = p0.size();
	if (m == 0)
	{
		setm::Polynomial<ftype> answer;
		return answer;
	}
	for (int i = 1; i < m; i++)
	{
		der_p0.push_back(p0[i] * i);
	}
	setm::Polynomial<ftype> answer(der_p0.cbegin(), der_p0.cend());
	return answer;
}



int main()
{
	std::cout.setf(std::ios::fixed);
	std::cout << "I am Work\n\n\n";

	//Дано
	constexpr ftype gamma_0 = 5.0 / 3.0;
	constexpr ftype ro_0 = 1E-5;
	constexpr ftype V_0 = 0;
	constexpr ftype P_0 = 3.848E3;

	constexpr ftype gamma_3 = 5.0 / 3.0;
	constexpr ftype ro_3 = 11.36997486377774;
	constexpr ftype C_3 = 1.31478E-4;
	constexpr ftype V_3 = 5E4;
	constexpr ftype P_3 = 1.17928E9;


	//Замена
	constexpr ftype alpha_0 = (gamma_0 + 1.0) / (gamma_0 - 1.0);
	constexpr ftype n = (2 * gamma_0) / (gamma_3 - 1.0); //n = 5
	constexpr ftype X = P_3 / P_0;
	ftype mu = (V_3 - V_0) * std::sqrtl(((gamma_0 - 1.0) * ro_0) / (2 * P_0));
	ftype nu = 2/(gamma_3 - 1) * std::sqrtl((gamma_3 * (gamma_0 - 1) * P_3 * ro_0) / (2 * P_0 * ro_3));

	std::cout << "\nX = " << X;
	std::cout << "\nro_3 = " << ro_3;
	std::cout << "\nalpha_0 = " << alpha_0;
	std::cout << "\nn = " << n;
	std::cout << "\nmu = " << mu;
	std::cout << "\nnu = " << nu;
	std::cout << "\n";



	ftype Z;
	ftype P_1;


	//Рассчитаем индексы 
	ftype a_2n = X * X;
	ftype a_n2 = -(alpha_0 * nu * nu * X);
	ftype a_n1 = 2 * alpha_0 * nu * (mu + nu) * X;
	ftype a_n = -(2 + (mu + nu) * (mu + nu) * alpha_0) * X;
	ftype a_2 = -(nu * nu);
	ftype a_1 = 2 * nu * (mu + nu);
	ftype a_0 = 1 - (mu + nu) * (mu + nu);

	std::vector <ftype> indexs(max_index(n) + 1, 0.0);
	//print_v(indexs);
	indexs[2*n] = a_2n;
	indexs[n+2] = a_n2;
	indexs[n+1] = a_n1;
	indexs[n] = a_n;
	indexs[2] = a_2;
	indexs[1] = a_1;
	indexs[0] = a_0;

	//печатаем индексы
	//print_v(indexs);
	//std::cout << "\n";
	setm::Polynomial<ftype> p0(indexs.cbegin(), indexs.cend());
	std::cout << "\n\n" << p0;

	//локализация корней
	auto loc = localization(indexs);
	std::cout << "\n\nlocalization\n" << loc.first << " < x < " << loc.second << '\n';
	std::cout << "f(l) = " << func(loc.first, indexs) << "\t\tf(r) = " << func(loc.second, indexs);
	std::cout << "\n\n";

	//ряд штурма
	std::vector<setm::Polynomial<ftype>> row_sthurm;
	row_sthurm.push_back(p0);
	std::cout << p0 << "\n\n";
	setm::Polynomial<ftype> p1 = derivative(p0);
	row_sthurm.push_back(p1);
	std::cout << p1 << "\n\n";
	setm::Polynomial<ftype> pn(0);
	auto k = row_sthurm.size();
	auto pn_2 = row_sthurm[k - 2];
	auto pn_1 = row_sthurm[k - 1];
	pn = (pn_2 / pn_1);
	//setm::Polynomial<ftype> p2 = p0 - (p0 / p1)*p1;
	//std::cout << p2;
	do
	{
		k = row_sthurm.size();
		pn_2 = row_sthurm[k - 2];
		pn_1 = row_sthurm[k - 1];
		pn = (pn_2/ pn_1) * pn_1;
		std::cout << (pn_2 / pn_1) << "\n";
		std::cout << pn_1 << '\n';
		std::cout << pn << "\n";
		auto n = pn.size();
		auto a0 = pn_2[n - 1];
		auto b0 = pn[n - 1];
		if ((a0 - b0) / b0 < 0.0000001)
		{
			pn_2[n - 1] = 0;
			pn[n - 1] = 0;
			pn = pn_2 - pn;
		}
		row_sthurm.push_back(pn);
		std::cout << row_sthurm.back() << "\n\n";
	} while (k<9);
}



