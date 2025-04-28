// Vec2D.h : Vec2D struct
//

#pragma once

#include <cmath>


struct Vec2D
{
	double x, y;

	Vec2D()
	{
	}

	Vec2D(const Vec2D &v) : x(v.x), y(v.y)
	{
	}

	Vec2D(double xx, double yy) : x(xx), y(yy)
	{
	}

	Vec2D &operator = (const Vec2D &v)
	{
		x = v.x;  y = v.y;  return *this;
	}

	Vec2D operator + (const Vec2D &v) const
	{
		return Vec2D(x + v.x, y + v.y);
	}

	Vec2D &operator += (const Vec2D &v)
	{
		x += v.x;  y += v.y;  return *this;
	}

	Vec2D operator - (const Vec2D &v) const
	{
		return Vec2D(x - v.x, y - v.y);
	}

	Vec2D &operator -= (const Vec2D &v)
	{
		x -= v.x;  y -= v.y;  return *this;
	}

	Vec2D operator - () const
	{
		return Vec2D(-x, -y);
	}

	Vec2D operator * (double a) const
	{
		return Vec2D(a * x, a * y);
	}

	Vec2D &operator *= (double a)
	{
		x *= a;  y *= a;  return *this;
	}

	double operator * (const Vec2D &v) const
	{
		return x * v.x + y * v.y;
	}

	Vec2D operator / (double a) const
	{
		return (*this) * (1 / a);
	}

	Vec2D operator /= (double a)
	{
		return (*this) *= (1 / a);
	}

	double operator ^ (const Vec2D &v) const
	{
		return x * v.y - y * v.x;
	}

	Vec2D operator ~() const
	{
		return Vec2D(y, -x);
	}

	double sqr() const
	{
		return x * x + y * y;
	}

	double len() const
	{
		return sqrt(x * x + y * y);
	}
};

inline Vec2D operator * (double a, const Vec2D &v)
{
	return v * a;
}
