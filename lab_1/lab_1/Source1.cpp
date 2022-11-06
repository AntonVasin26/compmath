#include <iostream>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include "polinom.hpp"
#include <iomanip>


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

using arr_t = std::vector <ftype>;
arr_t polynomial_remainder(const arr_t a, const arr_t b)
{
	auto m = std::size(b) - 1;
	auto n = std::size(a) - 1;

	if (m > n)
	{
		throw std::logic_error("n should be greater than m");
	}
	if (b[m] == 0.0f)
	{
		throw std::logic_error("b_m should be greater than 0");
	}

	arr_t c(n + 1);
	std::copy(std::begin(a), std::end(a), std::begin(c));

	float k_i;
	for (auto i = 0U; i < n - m + 1U; i++)
	{
		k_i = c[n - i] / b[m];
		for (auto j = i; j <= m + i; j++)
		{
			c[n - j] -= k_i * b[m - j + i];
		}
		c.pop_back();
	}
	while (c[std::size(c) - 1] == 0.0)
	{
		c.pop_back();
	}
	return c;
}

int max_index(ftype n)
{
	const int n_1 = std::round(n);
	return 2 * n_1;
}

void print_v(std::vector<ftype> vec)
{
	std::cout << "\nPrint index\n";

	for (int i = 0; i < vec.size(); i++)
	{
		std::cout << vec[i] << "x^" << i << '\t';
	}
	std::cout << '\n';
}

std::pair<ftype, ftype> localization(std::vector<ftype> v)
{
	int n = v.size() - 1;
	auto B = abs(*std::max_element(v.begin() + 1, v.end(), [](ftype a, ftype b) { return abs(a) < abs(b); }));
	auto A = abs(*std::max_element(v.begin(), v.end() - 1, [](ftype a, ftype b) { return abs(a) < abs(b); }));
	auto a_0 = abs(v[n]);
	auto a_n = abs(v[0]);
	auto min = a_n / (a_n + B);
	auto max = 1 + A / a_0;
	return std::pair<ftype, ftype>(min, max);
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

std::vector<ftype> derivative_v(std::vector<ftype> p0)
{
	std::vector<ftype> der_p0;
	int m = p0.size();
	if (m == 0)
	{
		return der_p0;
	}
	for (int i = 1; i < m; i++)
	{
		der_p0.push_back(p0[i] * i);
	}
	return der_p0;
}

bool is_positiv(ftype value)
{
	if (value > 0)
		return true;
	if (value <= 0)
		return false;
}

ftype count_sign_change(std::vector<std::vector<ftype>> row, ftype value)
{
	int count = 0;
	bool plus = is_positiv(func(value, row[0]));
	for (auto elm : row)
	{
		ftype sum = func(value, elm);
		if (sum == 0)
			continue;
		bool sign = is_positiv(sum);
		if (sign == plus)
			continue;
		plus = sign;
		count++;
	}
	return count;
}

std::vector< std::pair<ftype, ftype>> sthurm_grid(
	std::vector<std::vector<ftype>>row_sthurm, ftype first, ftype second, ftype eps)
{
	int k = ceil(abs(second - first) / eps);
	std::vector<std::pair<ftype, ftype>> grid(k+2);
	ftype step = eps;
	for (int i = 0; i < k; i++)
	{
		grid[i].first = first + step * i;
		grid[i].second = count_sign_change(row_sthurm, grid[i].first);
	}
	grid[k+1].first = second;
	grid[k+1].second = count_sign_change(row_sthurm, grid[k+1].first);
	std::vector< std::pair<ftype, ftype>> answer;

	for (int i = 1; i < grid.size(); i++)
	{
		if (grid[i - 1].second != grid[i].second)
		{
			auto x1 = grid[i - 1].first;
			auto x2 = grid[i].first;
			std::pair<ftype, ftype> interwal(x1, x2);
			answer.push_back(interwal);
		}
	}
	return answer;
}
	


int main()
{
	//std::cout.setf(std::ios::fixed);
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
	ftype nu = 2 / (gamma_3 - 1) * std::sqrtl((gamma_3 * (gamma_0 - 1) * P_3 * ro_0) / (2 * P_0 * ro_3));

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
	indexs[2 * n] = a_2n;
	indexs[n + 2] = a_n2;
	indexs[n + 1] = a_n1;
	indexs[n] = a_n;
	indexs[2] = a_2;
	indexs[1] = a_1;
	indexs[0] = a_0;

	//печатаем индексы
	//print_v(indexs);
	//std::cout << "\n";
	/*setm::Polynomial<ftype> p0(indexs.cbegin(), indexs.cend());*/
	auto p0 = indexs;
	//std::cout << "\n\n" << p0;
	print_v(p0);
	//локализация корней
	auto loc = localization(p0);
	std::cout << "\n\nlocalization\n" << loc.first << " < x < " << loc.second << '\n';
	std::cout << "f(1) = " << func(-1*loc.second, p0) << "\t\tf(2) = " << func(-1*loc.first, p0)
		<< "\t\tf(3) = " << func(loc.first, p0) << "\t\tf(4) = " << func(loc.second, p0);
	std::cout << "\n\n";

	//ряд штурма
	//std::vector<setm::Polynomial<ftype>> row_sthurm;
	//row_sthurm.push_back(p0);
	//std::cout << p0 << "\n\n";
	//setm::Polynomial<ftype> p1 = derivative(p0);
	//row_sthurm.push_back(p1);
	//std::cout << p1 << "\n\n";
	//setm::Polynomial<ftype> pn(0);
	//auto k = row_sthurm.size();
	//setm::Polynomial<ftype> p2 = p0 - (p0 / p1)*p1;
	//std::cout << p2;
	//do
	//{
	//	k = row_sthurm.size();
	//	auto pn_2 = row_sthurm[k - 2];
	//	auto pn_1 = row_sthurm[k - 1];
	//	pn = pn_2 - (pn_2 / pn_1) * pn_1;
	//	row_sthurm.push_back(pn);
	//	std::cout << pn << "\n\n";
	//} while (k < 9);

	
	
	//ряд штурма
	std::vector<std::vector<ftype>> row_sthurm;
	row_sthurm.push_back(p0);
	print_v(p0);
	std::vector<ftype> p1 = derivative_v(p0);
	row_sthurm.push_back(p1);
	print_v(p1);
	std::vector<ftype> pn;
	auto k = row_sthurm.size();
	do
	{
		k = row_sthurm.size();
		auto pn_2 = row_sthurm[k - 2];
		auto pn_1 = row_sthurm[k - 1];
		pn = polynomial_remainder(pn_2, pn_1);
		row_sthurm.push_back(pn);
		print_v(pn);
	} while (pn.size()>1);

	//Составлен ряд штурма.
	//проверим число корней на концах интервалов локализации.
	int i1 = count_sign_change(row_sthurm, -loc.second);
	int i2 = count_sign_change(row_sthurm, -loc.first);
	int i3 = count_sign_change(row_sthurm, +loc.first);
	int i4 = count_sign_change(row_sthurm, -loc.second);
	std::cout << "\n\n (" << i1 << " ; " << i2 << ") \t\t" << "(" << i3 << ";" << i4 << ")\n\n";

	//видно что по однуму корню на каждом участке
	auto eps = 0.00001;
	
	
	
	//поиск корней
	std::cout << "f(1) = " << func(-1.0 * loc.second, p0) << "\t\tf(2) = " << func(-1.0 * loc.first, p0)
		<< "\t\tf(3) = " << func(loc.first, p0) << "\t\tf(4) = " << func(loc.second, p0);
	std::cout << "\n\n";

	//Видно что можно воспользоватся методом методом половинного деления.
	ftype x1 = -1.0 * loc.second;
	auto f1 = func(x1, p0);
	ftype x2 = -1.0 * loc.first;
	auto f2 = func(x2, p0);
	int l = 100;
	ftype x3;
	for (int i = 0; i < l; i++)
	{
		x3 = (x1 + x2)/2;
		auto f3 = func(x3, p0);
		if (f1 * f3 == 0)
			break;
		if (f1 * f3 < 0)
			x2 = x3;
		else
			x1 = x3;
	}

	//оценим точность вычисления
	auto delta_X = (loc.second - loc.first) / powl(2, l);
	std::cout << "\n\nx1 = " << std::setprecision(15) << x3 << " +|- " << delta_X;
	std::cout << "\n f(x1) = " << func(x3, p0);

	//Найдём второй корень
	x1 = 1.0 * loc.first;
	f1 = func(x1, p0);
	x2 = 1.0 * loc.second;
	f2 = func(x2, p0);
	l = 100;
	//ftype x3;
	for (int i = 0; i < l; i++)
	{
		x3 = (x1 + x2) / 2;
		auto f3 = func(x3, p0);
		if (f1 * f3 == 0)
			break;
		if (f1 * f3 < 0)
			x2 = x3;
		else
			x1 = x3;
	}
	std::cout << "\n\nx2 = " << std::setprecision(15) << x3 << " +|- " << delta_X;
	std::cout << "\n f(x2) = " << func(x3, p0);
}