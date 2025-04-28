// Main.cpp : defines the entry point for the application
//

#include <ctime>
#include <crtdbg.h>
#include "Resource.h"
#include "Tri2D.h"
#include "Tree.h"


const size_t buf_size = 256;

HINSTANCE hInstance;
TriGrid *pGrid = 0;
bool play = false;


int CALLBACK About(HWND hDlg, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch(msg)
	{
	case WM_INITDIALOG:  return true;

	case WM_COMMAND:
		if(LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));  return true;
		}
		return false;
	}
	return false;
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch(msg)
	{
	case WM_COMMAND:
		switch(LOWORD(wParam))
		{
		case IDM_SAVE:
			{
				char buf[1024] = "";  OPENFILENAME ofn;
				memset(&ofn, 0, sizeof(ofn));  ofn.lStructSize = sizeof(ofn);
				ofn.hwndOwner = hWnd;  ofn.lpstrFilter = "All Files (*.*)\0*.*\0";
				ofn.lpstrFile = buf;  ofn.nMaxFile = sizeof(buf);
				ofn.Flags = OFN_OVERWRITEPROMPT;  ofn.lpstrDefExt = "dat";
				if(GetSaveFileName(&ofn))pGrid->Dump(buf);  return 0;
			}
		case IDM_ABOUT:  DialogBox(hInstance, (LPCTSTR)IDD_ABOUTBOX, hWnd, (DLGPROC)About);  return 0;
		case IDM_EXIT:  DestroyWindow(hWnd);  return 0;
		}
		break;

	case WM_KEYDOWN:
		switch(wParam)
		{
		case VK_ESCAPE:  DestroyWindow(hWnd);  return 0;

		case VK_SPACE:
			if(play)KillTimer(hWnd, 0);
			else SetTimer(hWnd, 0, 100, 0);
			play = !play;  return 0;
		}
		break;

	case WM_LBUTTONDOWN:
		pGrid->Step();  InvalidateRect(hWnd, 0, FALSE);  return 0;

	case WM_RBUTTONDOWN:
		pGrid->Subdivide();  InvalidateRect(hWnd, 0, FALSE);  return 0;

	case WM_MOUSEWHEEL:
		{
			int delta = GET_WHEEL_DELTA_WPARAM(wParam);
			if(delta > 0)pGrid->IncLambda();
			else pGrid->DecLambda();
			InvalidateRect(hWnd, 0, FALSE);  return 0;
		}

	case WM_PAINT:
		{
			PAINTSTRUCT ps;  HDC hdc = BeginPaint(hWnd, &ps);
			RECT rc;  GetClientRect(hWnd, &rc);
			FillRect(hdc, &rc, (HBRUSH)GetStockObject(WHITE_BRUSH));
			double h = min(rc.right - rc.left, rc.bottom - rc.top) / 2.1;
			pGrid->Draw(hdc, (rc.left + rc.right) / 2, (rc.top + rc.bottom) / 2, h);
			EndPaint(hWnd, &ps);  return 0;
		}

	case WM_TIMER:
		pGrid->ClearCounter();  for(size_t i = 0; i < 10; i++)pGrid->Step();
		InvalidateRect(hWnd, 0, FALSE);  return 0;

	case WM_DESTROY:
		if(play)KillTimer(hWnd, 0);  PostQuitMessage(0);  return 0;
	}
	return DefWindowProc(hWnd, msg, wParam, lParam);
}


struct Ro1
{
	double operator () (Vec2D &R)
	{
		double r = R.sqr() - 0.25;  return 0.01 + exp(-10 * r);
	}
};

struct Ro2
{
	double operator () (Vec2D &R)
	{
		double r = 0.25 - R.sqr();  return 0.01 + exp(-10 * r);
	}
};

int APIENTRY WinMain(HINSTANCE hInst, HINSTANCE, LPTSTR, int show)
{
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_CHECK_ALWAYS_DF |
		_CRTDBG_CHECK_CRT_DF | _CRTDBG_DELAY_FREE_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	hInstance = hInst;  init_random();

	const int N = 1024;
	Ro1 ro1;  FunctionalDensity<Ro1> area1(ro1);
	Ro2 ro2;  FunctionalDensity<Ro2> area2(ro2);
	BorderPoint bpt[2 * N];
	BorderEdge beg[2 * N + 2];
	for(int i = 0; i < N; i++)
	{
		double fi = 2 * Pi * i / N, r = 0.9 + 0.1 * cos(8 * fi);
		bpt[i].x = -r * cos(fi);  bpt[i].y = -r * sin(fi);  bpt[i].weight = 0;
		bpt[i + N].x = -0.5 * cos(fi);  bpt[i + N].y = -0.5 * sin(fi);  bpt[i + N].weight = 0;
		beg[i].pt[1] = beg[i + 1].pt[0] = &bpt[i];  beg[i + 1].lambda = 1;
		beg[i + 1].area[0] = &area1;  beg[i + 1].area[1] = 0;
		beg[i + N + 1].pt[1] = beg[i + N + 2].pt[0] = &bpt[i + N];  beg[i + N + 2].lambda = 1;
		beg[i + N + 2].area[0] = &area2;  beg[i + N + 2].area[1] = &area1;
	}
	beg[0].pt[0] = 0;  beg[N].pt[1] = &bpt[0];  beg[0].lambda = 0;
	beg[0].area[0] = beg[0].area[1] = 0;
	beg[N + 1].pt[0] = &bpt[0];  beg[2 * N + 1].pt[1] = &bpt[N];  beg[N + 1].lambda = 0;
	beg[N + 1].area[0] = beg[N + 1].area[1] = 0;

	TriGrid grid(beg, 2 * N + 2);  pGrid = &grid;
	//grid.CreateRandom(1, 1000);
	grid.CreateRegular(1.2, 10);


	char classname[] = "Tri2DMWC";

	WNDCLASSEX wcex;
	wcex.cbSize = sizeof(wcex);
	wcex.style = CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc = WndProc;
	wcex.cbClsExtra = wcex.cbWndExtra = 0;
	wcex.hInstance = hInst;
	wcex.hIcon = LoadIcon(hInst, (LPCTSTR)IDI_LARGE);
	wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
	wcex.hbrBackground = NULL;
	wcex.lpszMenuName = (LPCTSTR)IDM_MAIN;
	wcex.lpszClassName = classname;
	wcex.hIconSm = LoadIcon(hInst, (LPCTSTR)IDI_SMALL);
	RegisterClassEx(&wcex);


	char title[buf_size];
	LoadString(hInst, IDS_TITLE, title, buf_size);
	HWND hWnd = CreateWindowEx(0, classname, title, WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, 0, 0, hInst, 0);
	if(!hWnd)return -1;

	ShowWindow(hWnd, show);
	UpdateWindow(hWnd);


	MSG msg;
	while(GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return (int)msg.wParam;
}