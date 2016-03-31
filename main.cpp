#include <iostream>
#include "BSplineSurfIntp.h"

using namespace std;

/* parameter function */
//void paramFunc(double* ret, int sz)
//{
//	double de = 1.0 / double(sz);
//	for (int i = 0; i < sz; ++i)
//		ret[i] = double(i)*de;
//}

/* knot sequence generating function */
//void knotFunc(double* ret, int degree, int n_dat, int sz)
//{
//	for (int i = 0; i <= degree; ++i)
//		ret[i] = 0;
//	for (int i = n_dat; i < sz; ++i)
//		ret[i] = 1;
//	for (int i = 1; i < n_dat - degree; ++i)
//		ret[i + degree] = i / double(n_dat - degree);
//}

int main(int argc, char* argv[])
{
	int deg = 2;
	BSplineSurfIntp bsp(3,3,4,4);
	//bsp.setParamFunc(paramFunc);
	//bsp.setKnotSeqFunc(knotFunc);
	ChordParams chordParams;
	CentrParams centrParams;
	//bsp.setParamFunc(&chordParams);
	bsp.setParamFunc(&centrParams);
	AvgKnotSeq avgKnotSeq;
	bsp.setKnotSeqFunc(&avgKnotSeq);
	
	vector<Vector3d> data;
	data.push_back(Vector3d(0, 6, 0));
	data.push_back(Vector3d(2, 9, 0));
	data.push_back(Vector3d(4, 6, 0));
	data.push_back(Vector3d(6, 9, 0));
	data.push_back(Vector3d(0, 4, 0));
	data.push_back(Vector3d(2, 7, 0));
	data.push_back(Vector3d(4, 4, 0));
	data.push_back(Vector3d(6, 7, 0));
	data.push_back(Vector3d(0, 2, 0));
	data.push_back(Vector3d(2, 5, 0));
	data.push_back(Vector3d(4, 2, 0));
	data.push_back(Vector3d(6, 5, 0));
	data.push_back(Vector3d(0, 0, 0));
	data.push_back(Vector3d(2, 3, 0));
	data.push_back(Vector3d(4, 0, 0));
	data.push_back(Vector3d(6, 3, 0));
	bsp.setData(data);
	bsp.print();
	if (bsp.compIntp())
	{
		vector<Vector3d> cont = bsp.getControlPoints();
		for (Vector3d c : cont)
			cout << c << endl;
	}
	vector<double> u_knot = bsp.getUKnot();
	for (double k : u_knot)
		cout << k << " ";
	cout << endl;

	vector<double> v_knot = bsp.getVKnot();
	for (double k : v_knot)
		cout << k << " ";
	cout << endl;
	
	return 0;
}
