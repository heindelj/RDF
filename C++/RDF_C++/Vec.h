#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <assert.h>

template<typename T>
class Vec
{
	std::vector<T> m_Vec;

public:
	// Constructors 
	template<typename ... Ps>
	Vec(Ps... elems)
	{
		m_Vec = get_data(elems...);
	}

	template<typename It>
	Vec(It begin, It end)
	{
		for (begin; begin != end; begin++)
			m_Vec.push_back(*begin);
	}

	// Operator overloads
	T operator [](int i)
	{
		return m_Vec[i];
	}

	friend std::ostream& operator<< (std::ostream& os, const Vec& V)
	{
		std::copy(V.m_Vec.begin(), V.m_Vec.end(),
			std::ostream_iterator<T>(os, " "));
		return os;
	}

	Vec<T> operator *(const Vec<T>& vec)
	{
		std::vector<T> mult;
		for (int i = 0; i < m_Vec.size(); i++)
			mult.push_back(m_Vec[i] * vec.m_Vec[i]);
		return Vec<T>(mult.begin(), mult.end());
	}

	Vec<T> operator +(const Vec<T>& vec)
	{
		std::vector<T> add;
		for (int i = 0; i < m_Vec.size(); i++)
			add.push_back(m_Vec[i] + vec.m_Vec[i]);
		return Vec<T>(add.begin(), add.end());
	}

	Vec<T> operator -(const Vec<T>& vec)
	{
		std::vector<T> minus;
		for (int i = 0; i < m_Vec.size(); i++)
			minus.push_back(m_Vec[i] - vec.m_Vec[i]);
		return Vec<T>(minus.begin(), minus.end());
	}

	T dot(const Vec<T>& vec)
	{
		T dotProd{};
		for (int i = 0; i < m_Vec.size(); i++)
			dotProd += m_Vec[i] * vec.m_Vec[i];
		return dotProd;
	}

	T elementProduct()
	{
		T product = m_Vec[0];
		for (int i = 1; i < m_Vec.size(); i++)
			product *= m_Vec[i];
		return product;
	}

	// returns the cross product of two 3-D vectors
	// cross product is only defined for 3-D (sort of)
	Vec<T> cross(const Vec<T>& vec)
	{
		assert(m_Vec.size() == 3);
		assert(vec.m_Vec.size() == 3);
		return Vec<T>(m_Vec[1] * vec.m_Vec[2] - m_Vec[2] * vec.m_Vec[1],
			m_Vec[2] * vec.m_Vec[0] - m_Vec[0] * vec.m_Vec[2],
			m_Vec[0] * vec.m_Vec[1] - m_Vec[1] * vec.m_Vec[0]);
	}

	// returns norm of the vector as tyep T
	// there is no simple way to normalize vector of ints without
	// loss of data, so be careful
	T norm()
	{
		T norm{};
		for (auto it = m_Vec.begin(); it != m_Vec.end(); it++)
			norm += ((*it) * (*it));
		return sqrt(norm);
	}

	// returns normalized version of vector
	Vec<T> normalize()
	{
		std::vector<T> normed;
		normed.reserve(m_Vec.size());
		auto length = norm();
		for (auto it = m_Vec.begin(); it != m_Vec.end(); it++)
			normed.emplace_back((*it) / length);
		return Vec(normed.begin(), normed.end());
	}

private:
	template<typename ... Ps>
	std::vector<T> get_data(const Ps& ... elems)
	{
		return { elems... };
	}
};