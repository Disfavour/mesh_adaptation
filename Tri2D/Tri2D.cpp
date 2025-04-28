// Tri2D.cpp : TriGrid class implementation
//

#include <cmath>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "Tri2D.h"

using namespace std;


#define DRAW_BORDER
#define DRAW_GRID
//#define DRAW_VORONOY
//#define DRAW_POINT
//#define DRAW_CROSS



const double Pi = 3.14159265358979323846264338327950288;
const double large_val = 1E10;


unsigned long long rand_val = 0x576D050FCC16835;

void init_random()
{
	rand_val = GetTickCount();
}

double random()
{
	rand_val = rand_val * 0x46D920EE0D9B85A3 + 0x54B639321933838D;
	return ((rand_val & 0xFFFFFFFF) + 0.5) / 4294967296.0;
}


// TriGrid class

void TriGrid::CreateRandom(double R, size_t n)
{
	point.resize(n);
	for(vector<Point>::iterator pt = point.begin(); pt != point.end(); pt++)
	{
		double fi = 2 * Pi * random(), r = R * sqrt(random());
		*pt = Point(r * cos(fi), r * sin(fi));
	}
	rebuild = true;
}

void TriGrid::CreateRegular(double R, size_t n)
{
	point.resize(3 * n * (n + 1) + 1);
	vector<Point>::iterator pt = point.begin();

	const double h = 0.86602540378443864676372317075294;  // sqrt(3) / 2

	double Rn = R / n, R2 = 0.5 * Rn, Rh = h * Rn;
	Vec2D V1(Rn, 0), V2(-R2, Rh), V3(-R2, -Rh), Vi, Vj;  size_t i, j;
	for(i = 0, Vi = V1; i < n; i++, Vi += V1)for(j = 0, Vj = Vi; j <= n; j++, Vj += V2)*pt++ = Point(Vj);
	for(i = 0, Vi = V2; i < n; i++, Vi += V2)for(j = 0, Vj = Vi; j <= n; j++, Vj += V3)*pt++ = Point(Vj);
	for(i = 0, Vi = V3; i < n; i++, Vi += V3)for(j = 0, Vj = Vi; j <= n; j++, Vj += V1)*pt++ = Point(Vj);
	*pt = Point(0, 0);  rebuild = true;
}


inline size_t inc(size_t i)
{
	return i >= 2 ? 0 : i + 1;
}

inline size_t dec(size_t i)
{
	return !i ? 2 : i - 1;
}

inline size_t shift(size_t i, bool add)
{
	return add ? inc(i) : dec(i);
}


void TriGrid::Triangulate()  // nPoints >= 2
{
	nTr++;

	EventAlloc queue;  Tree<QueueEvent, EventAlloc::Freer> event(queue.GetFreer());
	for(vector<Point>::iterator pt = point.begin(); pt != point.end(); pt++)
		event.Insert(queue.Alloc(&*pt));

	BeachAlloc beach;  Tree<BeachPoint, BeachAlloc::Freer> line(beach.GetFreer());
	QueueEvent *pt1 = event.First();  pt1->Remove();
	QueueEvent *pt2 = event.First();  pt2->Remove();

	Triangle *tr = &triangle[0];
	tr[0].pt[0] = tr[1].pt[0] = 0;
	tr[0].pt[1] = tr[1].pt[2] = pt1->pt;
	tr[0].pt[2] = tr[1].pt[1] = pt2->pt;

	Connect(tr, 0, tr + 1, 0);  Connect(tr, 1, tr + 1, 2);  Connect(tr, 2, tr + 1, 1);
	tr[0] = -(tr[1] = ~(*pt2->pt - *pt1->pt));  tr[0].w = tr[1].w = 0;

	BeachPoint *bp1 = beach.Alloc(pt1->pt, pt2->pt, tr++);
	line.Insert(bp1, Tree<BeachPoint>::Pos(0, false));
	BeachPoint *bp2 = beach.Alloc(pt2->pt, pt1->pt, tr++);
	line.Insert(bp2, Tree<BeachPoint>::Pos(bp1, true));
	queue.Free(pt1);  queue.Free(pt2);

	while(event)
	{
		QueueEvent *ev = event.First();  ev->Remove();
		if(ev->pt)
		{
			BeachPoint::Cmp cmp(*ev->pt);
			OwningTree<BeachPoint>::Pos pos = line.Find(cmp);
			BeachPoint *prev = pos.Node(), *next = pos.Node();
			if(pos.After())
			{
				next = next->Next();  if(!next)next = line.First();
			}
			else
			{
				prev = prev->Prev();  if(!prev)prev = line.Last();
			}

			tr[0].pt[0] = tr[1].pt[0] = 0;
			tr[0].pt[1] = tr[1].pt[2] = prev->pt[1];
			tr[0].pt[2] = tr[1].pt[1] = ev->pt;

			Connect(tr, 0, tr + 1, 0);  Connect(tr, 1, tr + 1, 2);
			Connect(tr, 2, prev->tr, 1);  Connect(next->tr, 2, tr + 1, 1);
			tr[0] = -(tr[1] = ~(*ev->pt - *prev->pt[1]));  tr[0].w = tr[1].w = 0;

			BeachPoint *bp1 = beach.Alloc(prev->pt[1], ev->pt, tr++);  line.Insert(bp1, pos);
			BeachPoint *bp2 = beach.Alloc(ev->pt, prev->pt[1], tr++);
			line.Insert(bp2, Tree<BeachPoint>::Pos(bp1, true));

			queue.Free(prev->ev);
			if(prev->ev = QueueEvent::Triangle(queue, prev, ev->pt))event.Insert(prev->ev);
			if(bp2->ev = QueueEvent::Triangle(queue, bp2, next->pt[1]))event.Insert(bp2->ev);
		}
		else
		{
			BeachPoint *prev = ev->seg, *next = ev->seg->Next();

			prev->tr->pt[0] = next->pt[1];  next->tr->pt[1] = prev->pt[0];
			Connect(prev->tr, 1, next->tr->next[0]);  Connect(next->tr, 2, prev->tr->next[2]);
			Connect(prev->tr, 2, next->tr, 0);  *prev->tr = *ev;  prev->tr->w = 1;
			*next->tr = ~(*prev->pt[0] - *next->pt[1]);

			next->pt[0] = prev->pt[0];  beach.Free(prev);
			if(BeachPoint *ptr = next->Prev())
			{
				queue.Free(ptr->ev);
				if(ptr->ev = QueueEvent::Triangle(queue, ptr, next->pt[1]))event.Insert(ptr->ev);
			}
			if(BeachPoint *ptr = next->Next())
			{
				queue.Free(next->ev);
				if(next->ev = QueueEvent::Triangle(queue, next, ptr->pt[1]))event.Insert(next->ev);
			}
		}
		queue.Free(ev);
	}
}


inline bool Circumcenter(TriGrid::Triangle *tr)
{
	Vec2D D = *tr->pt[1] - *tr->pt[2];
	if(tr->pt[0])
	{
		Vec2D R0 = (*tr->pt[1] + *tr->pt[2]) / 2, R = *tr->pt[0] - R0;  double S = R ^ D;
		if(S < 0)return false;  double dd4 = D.sqr() / 4, h = (R.sqr() - dd4) / S;
		*tr = R0 + (~D) * (h / 2);
	}
	else *tr = ~D;  return true;
}

inline void QueueEdge(TriGrid::TriRef &first, TriGrid::Triangle *ptr, size_t i)
{
	TriGrid::TriRef tr = (ptr < ptr->next[i].ptr ? TriGrid::TriRef(ptr, i) : ptr->next[i]);
	if(tr->queue[i])return;  tr->queue[i] = first;  first = TriGrid::TriRef(tr, i);
}

bool TriGrid::CheckEdge(TriRef &first)  // TODO: bags...
{
	TriRef tr = first, next = tr->next[tr.i];  first = first->queue[first.i];

	if(!tr->pt[0])
	{
		if(next->pt[0])return true;
		if((tr.i == 2) == ((*tr ^ *next) >= 0))return true;
	}
	else
	{
		if(!next->pt[0])return true;
		if(((*tr - *next) ^ (*tr->pt[inc(tr.i)] - *tr->pt[dec(tr.i)])) >= 0)return true;
	}

	// flip

	tr->pt[dec(tr.i)] = next->pt[next.i];  next->pt[dec(next.i)] = tr->pt[tr.i];
	Connect(tr, tr.i, next->next[inc(next.i)]);  Connect(next, next.i, tr->next[inc(tr.i)]);
	Connect(tr, inc(tr.i), next, inc(next.i));
	
	if(!Circumcenter(tr) || !Circumcenter(next))return false;
	QueueEdge(first, tr, inc(tr.i));  QueueEdge(first, next, inc(next.i));
	QueueEdge(first, tr, dec(tr.i));  QueueEdge(first, next, dec(next.i));
	return true;
}

bool TriGrid::Retriangulate()
{
	TriRef first = TriRef(&triangle[0], 3);
	for(vector<Triangle>::iterator tr = triangle.begin(); tr != triangle.end(); tr++)
	{
		if(!Circumcenter(&*tr))return false;
		if(&*tr < tr->next[0])
		{
			tr->queue[0] = first;  first = TriRef(&*tr, 0);
		}
		if(&*tr < tr->next[1])
		{
			tr->queue[1] = first;  first = TriRef(&*tr, 1);
		}
		if(&*tr < tr->next[2])
		{
			tr->queue[2] = first;  first = TriRef(&*tr, 2);
		}
	}
	while(first.i < 3)if(!CheckEdge(first))return false;  return true;
}


void TriGrid::ProcessBorder(BorderEdge *brd, BorderPoint *&last)
{
	if(brd->pt[0] > last)throw BadGeometry();

	Vec2D D, P, R;  double h, w;
	Point *pt;  TriRef cur;  bool cw;
	if(!brd->pt[0])
	{
		if(brd->lambda > 0)throw BadGeometry();
		D = Vec2D(1, 0);  cur = TriRef(&triangle[0], 0);  while(cur->pt[0])cur.ptr++;
		while(cur->y >= 0)cur.ptr = cur->next[1];  while(cur->y < 0)cur.ptr = cur->next[1];
		R = *cur;  h = cur->y;  w = 0;  pt = cur->pt[1];
		cur = cur->next[0];  cw = true;
	}
	else
	{
		pt = brd->pt[0]->pt;  if(brd->pt[1] <= last && pt == brd->pt[1]->pt)return;
		D = *brd->pt[1] - *brd->pt[0];  P = *brd->pt[0];  cur = brd->pt[0]->tr;
		w = cur->w;  R = *cur - w * (*brd->pt[1]);  h = (D ^ R);  cw = (h > 0);
		cur = cw ? cur->next[cur.i] : cur->next[dec(cur.i)];
	}

	for(;;)
	{
		if(brd->pt[1] <= last && pt == brd->pt[1]->pt)break;

		Vec2D RR;  double hh;
		for(;;)
		{
			RR = *cur - cur->w * (*brd->pt[1]);  hh = (D ^ RR);  if(cw ? hh <= 0 : hh >= 0)break;
			h = hh;  R = RR;  w = cur->w;  cur = cur->next[shift(cur.i, cw)];
		}

		TriRef next = cw ? cur->next[cur.i] : cur;
		double f = (w * R + cur->w * RR) ^ (cur->w * R - w * RR);
		if((cw ? f >= 0 : f <= 0) && (brd->pt[1] > last || brd->pt[1]->pt != next->pt[dec(next.i)]))
		{
			if(brd->pt[1] <= last)throw BadGeometry();
			brd->pt[1]->pt = pt;  last = brd->pt[1];  brd->pt[1]->tr = next;  break;
		}

		double t = (w * h + cur->w * hh) / (cur->w * h - w * hh);
		if(cur.ptr > cur->next[cur.i].ptr)t = -t;
		if(!cw)
		{
			R = RR;  h = hh;  w = cur->w;  cur = cur->next[cur.i];
		}

		TriRef curnew = next->next[inc(next.i)];
		Point *ptnew = next->pt[dec(next.i)];  size_t sw = 0;
		if(cur.ptr > next.ptr)
		{
			swap(cur, next);  sw = 1;
		}
		double c1 = (cur->w + next->w * t) / (cur->w + next->w);
		double c2 = (next->w - cur->w * t) / (cur->w + next->w);
		Vec2D C = (*cur) * c1 + (*next) * c2;

		if(brd->area[0] != brd->area[1])
		{
			if(brd->area[0])
			{
				brd->area[0]->Integrate(pt, &P, &C);
				if(brd->lambda > 0)brd->area[0]->IntegrateBorder(pt, &P, &C, lambda * brd->lambda);
			}
			if(brd->area[1])
			{
				brd->area[1]->Integrate(pt, &C, &P);
				if(brd->lambda > 0)brd->area[1]->IntegrateBorder(pt, &C, &P, lambda * brd->lambda);
			}

			EdgeCross *crs = cross.append();  *crs = C;  crs->t = t;  
			crs->area[0] = brd->area[sw];  crs->area[1] = brd->area[sw ^ 1];

			EdgeCross *old = 0, **ptr = &cur->cross[cur.i];
			while(*ptr && (*ptr)->t < t)
			{
				old = *ptr;  ptr = &(*ptr)->next[0];
			}
			(*ptr ? (*ptr)->next[1] : next->cross[next.i]) = crs;
			crs->next[0] = *ptr;  crs->next[1] = old;  *ptr = crs;
		}
		else if(brd->lambda > 0 && brd->area[0])
			brd->area[0]->IntegrateBorder(pt, &P, &C, 2 * lambda * brd->lambda);

		P = C;  pt = ptnew;  cur = curnew;  cw = true;
	}

	if(brd->area[0] != brd->area[1])
	{
		if(brd->area[0])
		{
			brd->area[0]->Integrate(pt, &P, brd->pt[1]);
			if(brd->lambda > 0)brd->area[0]->IntegrateBorder(pt, &P, brd->pt[1], lambda * brd->lambda);
		}
		if(brd->area[1])
		{
			brd->area[1]->Integrate(pt, brd->pt[1], &P);
			if(brd->lambda > 0)brd->area[1]->IntegrateBorder(pt, brd->pt[1], &P, lambda * brd->lambda);
		}
	}
	else if(brd->lambda > 0 && brd->area[0])
		brd->area[0]->IntegrateBorder(pt, &P, brd->pt[1], 2 * lambda * brd->lambda);

	if(brd->pt[1]->weight > 0)
	{
		double w = lambda * brd->pt[1]->weight;  pt->RoR += *brd->pt[1] * w;  pt->Ro += w;
	}
}

void TriGrid::ProcessEdge(TriRef &first)
{
	TriRef tr = first;  first = first->queue[first.i];  // remove

	TriRef next = tr->next[tr.i];
	Point *pt1 = tr->pt[dec(tr.i)], *pt2 = tr->pt[inc(tr.i)];
	size_t dir = (tr.ptr < next.ptr ? 0 : 1);

	Density *area = tr->area;  Vec2D *pt = tr;
	for(;;)
	{
		EdgeCross **ptr = &tr->cross[tr.i], *cr = *ptr;
		for(;; ptr = &cr->next[dir], cr = *ptr)
		{
			if(!cr)
			{
				if(tr->cross[tr.i])throw BadGeometry();
				if(area)
				{
					area->Integrate(pt1, pt, next);  area->Integrate(pt2, next, pt);
				}
				if(!next->queue[next.i])
				{
					size_t j1 = inc(next.i);
					TriRef tr1 = next->next[j1];  next->queue[j1] = first;
					if(!tr1->queue[tr1.i])first = TriRef(next, j1);

					size_t j2 = dec(next.i);
					TriRef tr2 = next->next[j2];  next->queue[j2] = first;
					if(!tr2->queue[tr2.i])first = TriRef(next, j2);

					next->area = area;
				}
				else if(area != next->area)throw BadGeometry();

				double s = 0;
				if(area && tr->area)
				{
					Vec2D L = *pt1 - *pt2;  L /= 2 * L.sqr();
					s = fabs((*tr + *next - *pt1 - *pt2) ^ L);  if(skew < s)skew = s;
				}
				tr->skew[tr.i] = next->skew[next.i] = s;  return;
			}
			if(cr->area[dir ^ 1] != area)continue;
			if(area)
			{
				area->Integrate(pt1, pt, cr);  area->Integrate(pt2, cr, pt);
			}
			pt = cr;  area = cr->area[dir];  *ptr = cr->next[dir];  break;
		}
	}
}

void TriGrid::Integrate()  // nPoints >= 2
{
	TriRef start(0, 0);
	for(vector<Triangle>::iterator tr = triangle.begin(); tr != triangle.end(); tr++)
	{
		tr->cross[0] = tr->cross[1] = tr->cross[2] = 0;
		tr->queue[0] = tr->queue[1] = tr->queue[2] = TriRef(0, 0);
		if(!tr->pt[0])start.ptr = &*tr;
	}

	TriRef first(start, 3);
	for(Triangle *cur = start;;)
	{
		cur->queue[0] = first;  first = TriRef(cur, 0);
		cur->area = 0;  if((cur = cur->next[1]) == start)break;
	}
	
	BorderPoint *last = 0;
	for(BorderEdge *brd = border; brd < border + nBorder; brd++)ProcessBorder(brd, last);
	while(first.i < 3)ProcessEdge(first);
}

void TriGrid::Step()
{
	cross.resize(cross.begin());  vector<Point>::iterator cur = point.begin();
	for(vector<Point>::iterator pt = point.begin(); pt != point.end(); pt++)
		if(pt->Ro > 0)
		{
			*cur = pt->RoR / pt->Ro;  cur->RoR = Vec2D(0, 0);  cur->Ro = 0;  cur++;
		}
		else rebuild = true;

	if(rebuild)
	{
		point.resize(cur - point.begin());  
		if(point.size() < 2)
		{
			triangle.resize(0);  return;
		}
		triangle.resize(2 * (point.size() - 1));
	}

	skew = 0;  if(rebuild || !Retriangulate())Triangulate();
	rebuild = false;  Integrate();  nSt++;
}


inline void TriGrid::SubdivideEdge(Triangle *tr, size_t i)
{
	Point *pt1 = tr->pt[inc(i)], *pt2 = tr->pt[dec(i)];  if(!pt1 || !pt2)return;
	if(tr->area || tr->next[i]->area)point.push_back(Point((*pt1 + *pt2) / 2));
}

void TriGrid::Subdivide()
{
	Point *beg = &point[0];  point.reserve(4 * point.size());
	for(vector<Triangle>::iterator tr = triangle.begin(); tr != triangle.end(); tr++)
	{
		if(tr->pt[0])tr->pt[0] = &point[0] + (tr->pt[0] - beg);
		tr->pt[1] = &point[0] + (tr->pt[1] - beg);
		tr->pt[2] = &point[0] + (tr->pt[2] - beg);

		if(&*tr < tr->next[0])SubdivideEdge(&*tr, 0);
		if(&*tr < tr->next[1])SubdivideEdge(&*tr, 1);
		if(&*tr < tr->next[2])SubdivideEdge(&*tr, 2);
	}
	rebuild = true;
}


inline int round(double x)
{
	return (int)floor(x + 0.5);
}

void TriGrid::DrawEdge(HDC hdc, int x0, int y0, double h, Triangle *tr, size_t i, HPEN pen1, HPEN pen2)
{
	Point *pt1 = tr->pt[inc(i)], *pt2 = tr->pt[dec(i)];  if(!pt1 || !pt2)return;
	SelectObject(hdc, tr->skew[i] >= 0.9 * skew ? pen2 : pen1);

#ifdef DRAW_GRID
	if(tr->area || tr->next[i]->area)
	{
		MoveToEx(hdc, x0 + round(h * pt1->x), y0 - round(h * pt1->y), 0);
		LineTo(hdc, x0 + round(h * pt2->x), y0 - round(h * pt2->y));
	}
#endif

#ifdef DRAW_VORONOY
	if(tr->pt[0])MoveToEx(hdc, x0 + round(h * tr->x), y0 - round(h * tr->y), 0);
	else
	{
		Vec2D R = (*pt1 + *pt2) / 2 + *tr * (10 / tr->len());
		MoveToEx(hdc, x0 + round(h * R.x), y0 - round(h * R.y), 0);
	}

	tr = tr->next[i];
	if(tr->pt[0])LineTo(hdc, x0 + round(h * tr->x), y0 - round(h * tr->y));
	else
	{
		Vec2D R = (*pt1 + *pt2) / 2 + *tr * (10 / tr->len());
		LineTo(hdc, x0 + round(h * R.x), y0 - round(h * R.y));
	}
#endif
}

void TriGrid::Draw(HDC hdc, int x0, int y0, double h)
{
	HPEN hIn   = CreatePen(PS_SOLID, 0, 0x0000FF);
	HPEN hOut  = CreatePen(PS_SOLID, 0, 0xFF0000);
	HPEN hBrd  = CreatePen(PS_SOLID, 2, 0xCCFFCC);
	HPEN hGray = CreatePen(PS_SOLID, 0, 0xCCCCCC);
	HPEN hBold = CreatePen(PS_SOLID, 2, 0x000000);

#ifdef DRAW_BORDER
	SelectObject(hdc, hBrd);
	for(BorderEdge *brd = border; brd < border + nBorder; brd++)
	{
		if(brd->lambda <= 0 && brd->area[0] == brd->area[1])continue;
		MoveToEx(hdc, x0 + round(h * brd->pt[0]->x), y0 - round(h * brd->pt[0]->y), 0);
		LineTo(hdc, x0 + round(h * brd->pt[1]->x), y0 - round(h * brd->pt[1]->y));
	}
#endif

	for(vector<Triangle>::iterator tr = triangle.begin(); tr != triangle.end(); tr++)
	{
		if(&*tr < tr->next[0])DrawEdge(hdc, x0, y0, h, &*tr, 0, hGray, hBold);
		if(&*tr < tr->next[1])DrawEdge(hdc, x0, y0, h, &*tr, 1, hGray, hBold);
		if(&*tr < tr->next[2])DrawEdge(hdc, x0, y0, h, &*tr, 2, hGray, hBold);

#ifdef DRAW_VORONOY
		if(tr->w <= 0)continue;
		int x = x0 + round(h * tr->x), y = y0 - round(h * tr->y);
		SelectObject(hdc, tr->area ? hIn : hOut);
		Ellipse(hdc, x - 1, y - 1, x + 2, y + 2);
#endif
	}

#ifdef DRAW_POINT
	SelectObject(hdc, GetStockObject(BLACK_PEN));
	for(vector<Point>::iterator pt = point.begin(); pt != point.end(); pt++)
	{
		int x = x0 + round(h * pt->x), y = y0 - round(h * pt->y);
		MoveToEx(hdc, x, y, 0);  Vec2D R = pt->RoR / pt->Ro;
		LineTo(hdc, x0 + round(h * R.x), y0 - round(h * R.y));
		Ellipse(hdc, x - 1, y - 1, x + 2, y + 2);
	}
#endif


	SelectObject(hdc, GetStockObject(BLACK_PEN));
#ifdef DRAW_CROSS
	for(CrossArray::iterator crs = cross.begin(); crs != cross.end(); crs++)
	{
		int x = x0 + round(h * crs->x), y = y0 - round(h * crs->y);
		Ellipse(hdc, x - 2, y - 2, x + 3, y + 3);
	}
#endif

	ostringstream out;  out << "N = " << point.size();
	out << ", Lambda = " << lambda << ", Skew = " << skew;
	out << ", Triangulate " << round(100.0 * nTr / nSt) << '%';
	TextOut(hdc, 16, 16, out.str().c_str(), out.str().length());

	DeleteObject(hIn);  DeleteObject(hOut);  DeleteObject(hBrd);  DeleteObject(hGray);
}

void TriGrid::Dump(const char *file)
{
	size_t nTr = 0;
	for(vector<Triangle>::iterator tr = triangle.begin(); tr != triangle.end(); tr++)
		if(tr->area)nTr++;

	ofstream out(file);
	out << "VARIABLES = \"X\", \"Y\"\nZONE N=" << point.size();
	out << ", E=" << nTr << ", DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n\n";

	for(vector<Point>::iterator pt = point.begin(); pt != point.end(); pt++)
		out << pt->x << ' ' << pt->y << '\n';  out << '\n';

	for(vector<Triangle>::iterator tr = triangle.begin(); tr != triangle.end(); tr++)if(tr->area)
	{
		out << (tr->pt[0] - &point[0]) + 1 << ' ';
		out << (tr->pt[1] - &point[0]) + 1 << ' ';
		out << (tr->pt[2] - &point[0]) + 1 << '\n';
	}
}
