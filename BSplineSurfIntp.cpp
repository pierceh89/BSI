#include "BSplineSurfIntp.h"

/* Compute chord length parameters for global surface interpolation */
bool ChordParams::computeParams(const BSplineSurfIntp& surf, vector<double>& ret_u, vector<double>& ret_v)
{
	const vector<Vector3d>& dat = surf.getDataPoints();
	// compute parameters in u direction
	int sz_v = surf.getVData();
	int sz_u = surf.getUData();
	vector<double> _ret_u(ret_u);
	vector<double> _ret_v(ret_v);
	_ret_u[0] = 0.0;
	_ret_v[0] = 0.0;
	double total;
	vector<double> dist;
	dist.resize(sz_v);
	for (int l = 0; l < sz_v; ++l)
	{
		total = 0.0;
		for (int k = 1; k < sz_u; ++k)
		{
			dist[k] = sqrtf((dat[l*sz_u + k] - dat[l*sz_u + k - 1]).squaredNorm());
			//cout << "dist : " << dist[k] << endl;
			total = total + dist[k];
		}
		if (total == 0.0)
			return false;
		for (int k = 1; k < sz_u; ++k)
		{
			_ret_u[k] = _ret_u[k - 1] + dist[k] / total;
			ret_u[k] = ret_u[k] + _ret_u[k];
		}
	}
	double div_sz_v = 1 / double(sz_v);
	for (int i = 0; i < sz_u; ++i)
	{
		ret_u[i] = ret_u[i] * div_sz_v;
		//cout << ret_u[i] << " ";
	}
	//cout << endl;

	// compute parameters in v direction
	dist.resize(sz_u);
	for (int l = 0; l < sz_u; ++l)
	{
		total = 0.0;
		for (int k = 1; k < sz_v; ++k)
		{
			dist[k] = sqrtf((dat[k*sz_u + l] - dat[(k-1)*sz_u + l]).squaredNorm());
			//cout << "dist : " << dist[k] << endl;
			total = total + dist[k];
		}
		if (total == 0.0)
			return false;
		for (int k = 1; k < sz_v; ++k)
		{
			_ret_v[k] = _ret_v[k - 1] + dist[k] / total;
			ret_v[k] = ret_v[k] + _ret_v[k];
		}
	}
	double div_sz_u = 1 / double(sz_u);
	for (int i = 0; i < sz_v; ++i)
	{
		ret_v[i] = ret_v[i] * div_sz_u;
		//cout << ret_v[i] << " ";
	}
	//cout << endl;
	return true;
}

/* Compute centripetal parameters for global surface interpolation */
bool CentrParams::computeParams(const BSplineSurfIntp& surf, vector<double>& ret_u, vector<double>& ret_v)
{
	const vector<Vector3d>& dat = surf.getDataPoints();
	// compute parameters in u direction
	int sz_v = surf.getVData();
	int sz_u = surf.getUData();
	vector<double> _ret_u(ret_u);
	vector<double> _ret_v(ret_v);
	_ret_u[0] = 0.0;
	_ret_v[0] = 0.0;
	double total;
	vector<double> dist;
	dist.resize(sz_v);
	for (int l = 0; l < sz_v; ++l)
	{
		total = 0.0;
		for (int k = 1; k < sz_u; ++k)
		{
			dist[k] = sqrtf(sqrtf((dat[l*sz_u + k] - dat[l*sz_u + k - 1]).squaredNorm()));
			//cout << "dist : " << dist[k] << endl;
			total = total + dist[k];
		}
		if (total == 0.0)
			return false;
		for (int k = 1; k < sz_u; ++k)
		{
			_ret_u[k] = _ret_u[k - 1] + dist[k] / total;
			ret_u[k] = ret_u[k] + _ret_u[k];
		}
	}
	double div_sz_v = 1 / double(sz_v);
	for (int i = 0; i < sz_u; ++i)
	{
		ret_u[i] = ret_u[i] * div_sz_v;
		//cout << ret_u[i] << " ";
	}
	//cout << endl;

	// compute parameters in v direction
	dist.resize(sz_u);
	for (int l = 0; l < sz_u; ++l)
	{
		total = 0.0;
		for (int k = 1; k < sz_v; ++k)
		{
			dist[k] = sqrtf(sqrtf((dat[k*sz_u + l] - dat[(k - 1)*sz_u + l]).squaredNorm()));
			//cout << "dist : " << dist[k] << endl;
			total = total + dist[k];
		}
		if (total == 0.0)
			return false;
		for (int k = 1; k < sz_v; ++k)
		{
			_ret_v[k] = _ret_v[k - 1] + dist[k] / total;
			ret_v[k] = ret_v[k] + _ret_v[k];
		}
	}
	double div_sz_u = 1 / double(sz_u);
	for (int i = 0; i < sz_v; ++i)
	{
		ret_v[i] = ret_v[i] * div_sz_u;
		//cout << ret_v[i] << " ";
	}
	//cout << endl;
	return true;
}

bool AvgKnotSeq::computeSeq(const BSplineSurfIntp& surf, vector<double>& ret_u_knot, vector<double>& ret_v_knot)
{
	int u_deg = surf.getUOrder() - 1;
	int v_deg = surf.getVOrder() - 1;
	if (u_deg < 1 || v_deg < 1)
		return false;

	const vector<double>& u_param = surf.getUParam();
	const vector<double>& v_param = surf.getVParam();
	int sz_u = u_param.size();
	int sz_v = v_param.size();
	//for (double v : u_param)
	//	cout << v << " ";
	//cout << endl;

	for (int i = 0; i <= u_deg; ++i)
		ret_u_knot[i] = 0;
	for (int i = sz_u; i < ret_u_knot.size(); ++i)
		ret_u_knot[i] = 1;
	double total;
	for (int i = 1; i < sz_u - u_deg; ++i)
	{
		total = 0.0;
		for (int j = i; j < i + u_deg; ++j)
			total = total + u_param[j];
		ret_u_knot[i + u_deg] = total / u_deg;
		//cout << total << " , " << ret_u_knot[i + u_deg] << endl;
	}
	
	for (int i = 0; i <= v_deg; ++i)
		ret_v_knot[i] = 0;
	for (int i = sz_v; i < ret_v_knot.size(); ++i)
		ret_v_knot[i] = 1;
	for (int i = 1; i < sz_u - u_deg; ++i)
	{
		total = 0.0;
		for (int j = i; j < i + v_deg; ++j)
			total = total + v_param[j];
		ret_v_knot[i + v_deg] = total / v_deg;
		//cout << total << " , " << ret_u_knot[i + v_deg] << endl;
	}
	return true;
}

void BSplineSurfIntp::print()
{
	cout << "order : ( " << u_order << " , " << v_order << " )" << endl;
	cout << "Data  : ( " << u_dat_num << " , " << v_dat_num << " )" << endl;
	cout << "Size of knot vectors : ( " << u_knot.size() << " , " << v_knot.size() << " )" << endl;
	cout << "Size of parameters   : ( " << u_param.size() << " , " << v_param.size() << " )" << endl;
	cout << "Size of data and control points : ( " << data.size() << " , " << con.size() << " )" << endl;
}

/* Initialize variables and vectors
 u,v orders (degree + 1)
 u,v data length */
BSplineSurfIntp::BSplineSurfIntp(int u_o, int v_o, int u_dat, int v_dat):u_order(u_o), v_order(v_o), u_dat_num(u_dat), v_dat_num(v_dat)
{
	size_t sz_dat = u_dat_num * v_dat_num;
	data.resize(sz_dat);
	con.resize(sz_dat);
	u_param.resize(u_dat_num);
	v_param.resize(v_dat_num);
	size_t sz_u_knot = u_order + u_dat_num;
	size_t sz_v_knot = v_order + v_dat_num;
	u_knot.resize(sz_u_knot);
	v_knot.resize(sz_v_knot);
}

/* Compute control points */
bool BSplineSurfIntp::compIntp()
{
	int sz = u_dat_num * v_dat_num;
	MatrixXd mat;
	mat.resize(sz, sz);

	double *_v_ret = nullptr, *_u_ret = nullptr;
	_u_ret = new double[u_dat_num];
	_v_ret = new double[v_dat_num];
	if (!_v_ret || !_u_ret)
		return false;

	// compute parameters
	if (!paramFunc->computeParams(*this, u_param, v_param))
		return false;

	// compute knot sequence
	if (!knotSeqFunc->computeSeq(*this, u_knot, v_knot))
		return false;

	int cnt = 0;
	for (int j = 0; j < v_dat_num; ++j)
	{
		if (!basis(_v_ret, v_param[j], v_order, v_knot))
			return false;
		for (int i = 0; i < u_dat_num; ++i)
		{
			if (!basis(_u_ret, u_param[i], u_order, u_knot))
				return false;
			for (int m = 0; m < v_dat_num; ++m)
			{
				for (int n = 0; n < u_dat_num; ++n)
				{
					mat(cnt, m*u_dat_num + n) = _v_ret[m] * _u_ret[n];
				}
			}
			cnt++;
		}
	}
	
	// compute control points
	VectorXd _x, _y, _z;
	_x.resize(sz);
	_y.resize(sz);
	_z.resize(sz);
	for (int i = 0; i < sz; ++i)
	{
		_x(i) = data[i](0);
		_y(i) = data[i](1);
		_z(i) = data[i](2);
	}
	VectorXd _con_x, _con_y, _con_z;
	_con_x = mat.jacobiSvd(ComputeThinU | ComputeThinV).solve(_x);
	_con_y = mat.jacobiSvd(ComputeThinU | ComputeThinV).solve(_y);
	_con_z = mat.jacobiSvd(ComputeThinU | ComputeThinV).solve(_z);

	// copy the control points
	for (int i = 0; i < sz; ++i)
	{
		con[i](0) = _con_x(i);
		con[i](1) = _con_y(i);
		con[i](2) = _con_z(i);
	}

	delete[] _v_ret;
	delete[] _u_ret;
	return true;
}

/*
* ret is the array of N(i,d)(t) values 
* t is the parameter
* order is B-Spline's order
* knot is the knot vector array  */
bool BSplineSurfIntp::basis(double* ret, double t, int order, const vector<double>& knot) const
{
	// d == 0
	int degree = order - 1;
	int m = knot.size() - 1;
	double **lowb = new double*[degree];	// degree = order-1
	if (!lowb)
		return false;
	lowb[0] = new double[m];
	if (!lowb[0])
		return false;
	for (int i = 0; i < m; i++)
	{
		if (t >= knot[i] && t < knot[i + 1])
			lowb[0][i] = 1;
		else
			lowb[0][i] = 0;
	}
	for (int k = 1; k < order-1; k++)
	{
		lowb[k] = new double[m - k];
		if (!lowb[k])
			return false;
		int len = m - k;
		for (int i = 0; i < len; i++)
		{
			double op1 = 0, op2 = 0;
			double den = knot[i + k] - knot[i];
			if (den != 0.0)
				op1 = (t - knot[i]) / den * lowb[k - 1][i];
			den = knot[i + k + 1] - knot[i + 1];
			if (den != 0.0)
				op2 = (knot[i + k + 1] - t) / den * lowb[k - 1][i + 1];
			lowb[k][i] = op1 + op2;
		}
	}

	// finished computing lower degree basis
	int n = m - degree; // # of control points
	for (int i = 0; i < n; i++)
	{
		double op1 = 0, op2 = 0;
		double den = knot[i + degree] - knot[i];
		if (den != 0.0)
			op1 = (t - knot[i]) / den * lowb[degree - 1][i];
		den = knot[i + degree + 1] - knot[i + 1];
		if (den != 0)
			op2 = (knot[i + degree + 1] - t) / den * lowb[degree - 1][i + 1];
		ret[i] = op1 + op2;
	}

	for (int k = 0; k < degree; k++)
		delete[] lowb[k];
	delete[] lowb;
	return true;
}

/* By using user-defined knot sequence generating function, set knot sequence for B-spline surface
The function signature should be: void knotFunc(double* ret, int degree, int n_data, int n_knot)  */
bool BSplineSurfIntp::setKnotSeqFunc(void(*knotFunc)(double*, int, int, int))
{
	double *_u_knot, *_v_knot;
	int sz_u = u_order + u_dat_num, sz_v = v_order + v_dat_num;
	// In case variables are not initialized;
	if (sz_u == 1 || sz_v == 1)
		return false;

	_u_knot = new double[sz_u];
	_v_knot = new double[sz_v];

	if (!_u_knot || !_v_knot)
		return false;
	knotFunc(_u_knot, u_order-1, u_dat_num, sz_u);
	knotFunc(_v_knot, v_order-1, v_dat_num, sz_v);

	for (int i = 0; i < sz_u; ++i)
		u_knot[i] = _u_knot[i];
	for (int i = 0; i < sz_v; ++i)
		v_knot[i] = _v_knot[i];

	delete[] _u_knot;
	delete[] _v_knot;

	return true;
}

/* By using user-defined parametrization function, set parameter for each data point
 The function signature should be: void paramFunc(double* ret, int knot_size) */
bool BSplineSurfIntp::setParamFunc(void(*paramFunc)(double*, int))
{
	double *_u_arr, *_v_arr;
	_u_arr = new double[u_dat_num];
	_v_arr = new double[v_dat_num];
	if (!_u_arr || !_v_arr)
	{	
		cerr << "Error: dynamic allocation error" << endl;
		return false;
	}
	paramFunc(_u_arr, u_dat_num);
	paramFunc(_v_arr, v_dat_num);

	// copy
	for (size_t i = 0; i < u_dat_num; ++i)
		u_param[i] = _u_arr[i];
	for (size_t i = 0; i < v_dat_num; ++i)
		v_param[i] = _v_arr[i];
	
	delete[] _u_arr;
	delete[] _v_arr;
	return true;
}

void BSplineSurfIntp::setData(const vector<Vector3d>& _data)
{
	size_t sz = u_dat_num * v_dat_num;
	if (sz != _data.size())
	{
		cerr << "Error: data size does not match" << endl;
		return;
	}
	for (size_t i = 0; i < _data.size(); ++i)
		data[i] = _data[i];
}

void BSplineSurfIntp::setOrder(int _u, int _v)
{
	u_order = _u; v_order = _v;
}

void BSplineSurfIntp::setDataSize(int _u, int _v)
{
	if (_u != u_dat_num || _v != v_dat_num)
	{
		u_dat_num = _u;
		v_dat_num = _v;
		data.resize(_u*_v);
		con.resize(_u*_v);
		u_knot.resize(_u + u_order + 1);
		v_knot.resize(_v + v_order + 1);
	}
}
