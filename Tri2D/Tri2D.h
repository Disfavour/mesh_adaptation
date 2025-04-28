// Tri2D.h : TriGrid class definition
//

#pragma once

#include <new>
#include <limits>
#include <vector>
#include <windows.h>
#include "Vec2D.h"
#include "Array.h"
#include "Tree.h"

using namespace std;



extern const double Pi;

void init_random();



template<typename T> struct NumRef
{
	T *ptr;
	size_t i;

	NumRef()
	{
	}

	template<typename T1> NumRef(T1 *p, size_t index) : ptr(reinterpret_cast<T *>(p)), i(index)
	{
	}

	NumRef(T *p, size_t index) : ptr(p), i(index)
	{
	}

	operator T *()
	{
		return ptr;
	}

	T *operator -> ()
	{
		return ptr;
	}

	operator bool() const
	{
		return ptr != 0;
	}

	bool operator ! () const
	{
		return ptr == 0;
	}
};


struct Density;
struct BorderPoint;
struct BorderEdge;

class TriGrid
{
public:
	struct BadGeometry : public exception
	{
		BadGeometry() : exception("bad geometry")
		{
		}
	};

	struct Point : public Vec2D
	{
		Vec2D RoR;
		double Ro;
		bool border;

		Point()
		{
		}

		Point(double x, double y) : Vec2D(x, y), RoR(x, y), Ro(1), border(false)
		{
		}

		Point(const Vec2D &pt) : Vec2D(pt), RoR(pt), Ro(1), border(false)
		{
		}

		Vec2D &operator = (const Vec2D &pt)
		{
			return *static_cast<Vec2D *>(this) = pt;
		}

		bool operator < (const Point &pt) const
		{
			return x < pt.x;
		}
	};

	struct EdgeCross : public Vec2D
	{
		double t;
		EdgeCross *next[2];
		Density *area[2];

		Vec2D &operator = (const Vec2D &pt)
		{
			return *static_cast<Vec2D *>(this) = pt;
		}
	};

	struct Triangle : public Vec2D
	{
		double w;
		NumRef<Triangle> next[3];
		Point *pt[3];

		Density *area;
		EdgeCross *cross[3];
		NumRef<Triangle> queue[3];
		double skew[3];

		Vec2D &operator = (const Vec2D &pt)
		{
			return *static_cast<Vec2D *>(this) = pt;
		}
	};

	typedef NumRef<Triangle> TriRef;


private:
	struct BeachPoint;
	struct QueueEvent;

	typedef BlockList<EdgeCross, 1024> CrossArray;
	typedef BlockListAllocator<BeachPoint, 1024> BeachAlloc;
	typedef BlockListAllocator<QueueEvent, 1024> EventAlloc;

	struct BeachPoint : public TreeNode<BeachPoint>
	{
		Point *pt[2];
		Triangle *tr;
		QueueEvent *ev;

		BeachPoint(Point *pt1, Point *pt2, Triangle *t) : tr(t), ev(0)
		{
			pt[0] = pt1;  pt[1] = pt2;
		}

		struct Cmp : public Vec2D
		{
			Cmp(Vec2D &v) : Vec2D(v)
			{
			}

			bool operator () (BeachPoint *bp)
			{
				Vec2D D = *bp->pt[1] - *bp->pt[0];
				if(D.y == 0)
				{
					if(D.x <= 0)return false;
					if(y <= bp->pt[0]->y)return x > bp->pt[1]->x;
				}
				double y1 = y - bp->pt[0]->y, y2 = y - bp->pt[1]->y, w = 4 * y1 * y2;
				if(D.x > 0)
				{
					double t = (D.x * D.x - w) / (D.x * (y1 + y2) + sqrt(w * D.sqr()));
					return 2 * x > bp->pt[1]->x + bp->pt[0]->x + t * D.y;
				}
				else
				{
					double t = D.x * (y1 + y2) - sqrt(w * D.sqr());
					return 2 * x > bp->pt[1]->x + bp->pt[0]->x + t / D.y;
				}
			}
		};
	};

	struct QueueEvent : public Vec2D, public TreeNode<QueueEvent>
	{
		double yy;
		BeachPoint *seg;
		Point *pt;

		QueueEvent(Point *p) : Vec2D(*p), yy(p->y), seg(0), pt(p)
		{
		}

		QueueEvent(BeachPoint *s, Vec2D &v, double delta) : Vec2D(v), yy(v.y + delta), seg(s), pt(0)
		{
		}

		static QueueEvent *Triangle(EventAlloc &alloc, BeachPoint *s, Point *pt2)
		{
			if(s->pt[0] == pt2)return 0;
			Vec2D R0 = (*s->pt[0] + *pt2) / 2, R = *s->pt[1] - R0;
			Vec2D D = *pt2 - *s->pt[0];  double S = R ^ D;  if(S <= 0)return 0;
			double dd4 = D.sqr() / 4, h = (R.sqr() - dd4) / S;
			return alloc.Alloc(s, R0 + (~D) * (h / 2), sqrt(dd4 * (1 + h * h)));
		}

		bool operator < (QueueEvent &qe)
		{
			if(yy < qe.yy)return true;
			if(yy > qe.yy)return false;
			return x < qe.x;
		}
	};


	vector<Point> point;
	vector<Triangle> triangle;
	CrossArray cross;
	BorderEdge *border;
	size_t nBorder;

	double lambda, skew;
	size_t nTr, nSt;
	bool rebuild;



	static void Connect(Triangle *tr1, size_t ed1, Triangle *tr2, size_t ed2)
	{
		tr1->next[ed1] = TriRef(tr2, ed2);  tr2->next[ed2] = TriRef(tr1, ed1);
	}

	static void Connect(Triangle *tr1, size_t ed1, TriRef tr2)
	{
		tr1->next[ed1] = tr2;  tr2->next[tr2.i] = TriRef(tr1, ed1);
	}


	void AddPoint(Point *pt);
	void Triangulate();

	bool CheckEdge(TriRef &first);
	bool Retriangulate();

	void ProcessBorder(BorderEdge *brd, BorderPoint *&last);
	void ProcessEdge(TriRef &first);
	void Integrate();

	void SubdivideEdge(Triangle *tr, size_t i);

	void DrawEdge(HDC hdc, int x0, int y0, double h, Triangle *tr, size_t i, HPEN pen1, HPEN pen2);


public:
	TriGrid(BorderEdge *brd, size_t nBrd) : border(brd), nBorder(nBrd),
		lambda(1), skew(0), nTr(0), nSt(0), rebuild(true)
	{
	}

	~TriGrid()
	{
	}


	void CreateRandom(double R, size_t n);
	void CreateRegular(double R, size_t n);

	void ClearCounter()
	{
		nTr = nSt = 0;
	}

	void Step();
	void Subdivide();

	void Draw(HDC hdc, int x0, int y0, double h);
	void Dump(const char *file);

	void IncLambda()
	{
		lambda *= 2;
	}

	void DecLambda()
	{
		lambda /= 2;
	}
};


struct Density
{
	virtual void Integrate(TriGrid::Point *pt, const Vec2D *pt1, const Vec2D *pt2) = 0;
	virtual void IntegrateBorder(TriGrid::Point *pt, const Vec2D *pt1, const Vec2D *pt2, double lambda) = 0;
};

struct ConstDensity : public Density
{
	double Ro;

	ConstDensity(double ro) : Ro(ro)
	{
	}

	void Integrate(TriGrid::Point *pt, const Vec2D *pt1, const Vec2D *pt2)
	{
		double S = ((*pt1 - *pt) ^ (*pt2 - *pt)) * (Ro / 2);
		pt->RoR += (*pt + *pt1 + *pt2) * (S / 3);  pt->Ro += S;
	}

	void IntegrateBorder(TriGrid::Point *pt, const Vec2D *pt1, const Vec2D *pt2, double lambda)
	{
		double l = lambda * (*pt1 - *pt2).len() * Ro;
		pt->RoR += (*pt1 + *pt2) * (l / 2);  pt->Ro += l;  pt->border = true;
	}
};

template<typename F> struct FunctionalDensity : public Density
{
	F Ro;

	FunctionalDensity(F ro) : Ro(ro)
	{
	}

	void Integrate(TriGrid::Point *pt, const Vec2D *pt1, const Vec2D *pt2)
	{
		Vec2D R0 = (*pt * 4 + *pt1 + *pt2) / 6;
		Vec2D R1 = (*pt + *pt1 * 4 + *pt2) / 6;
		Vec2D R2 = (*pt + *pt1 + *pt2 * 4) / 6;
		double Ro0 = Ro(R0), Ro1 = Ro(R1), Ro2 = Ro(R2);

		double S3 = ((*pt1 - *pt) ^ (*pt2 - *pt)) / 6;
		pt->RoR += (R0 * Ro0 + R1 * Ro1 + R2 * Ro2) * S3;
		pt->Ro += (Ro0 + Ro1 + Ro2) * S3;
	}

	void IntegrateBorder(TriGrid::Point *pt, const Vec2D *pt1, const Vec2D *pt2, double lambda)
	{
		const double t1 = 0.78867513459481288225457439025098;  // (3 + sqrt(3)) / 6
		const double t2 = 0.21132486540518711774542560974902;  // (3 - sqrt(3)) / 6
		Vec2D R1 = *pt1 * t1 + *pt2 * t2;
		Vec2D R2 = *pt1 * t2 + *pt2 * t1;
		double Ro1 = sqrt(Ro(R1)), Ro2 = sqrt(Ro(R2));

		double l2 = lambda * (*pt1 - *pt2).len() / 2;
		pt->RoR += (R1 * Ro1 + R2 * Ro2) * l2;  pt->Ro += (Ro1 + Ro2) * l2;  pt->border = true;
	}
};

struct BorderPoint : public Vec2D
{
	double weight;
	TriGrid::Point *pt;
	TriGrid::TriRef tr;
};

struct BorderEdge
{
	double lambda;
	BorderPoint *pt[2];
	Density *area[2];
};

