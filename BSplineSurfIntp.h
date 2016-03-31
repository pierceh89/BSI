/* B-Spline Surface Interpolation Class
 * Given M x N data points, M x N parameter values and 
 * order of curves in both direction,
 * Obtain M x N control points and knot vectors in both direction
 * for B-Spline surface which passes given data points
 * 
 * Author: Hyounggap An
 * Email : hyounggap.an@gmail.com
 * Date  : 2016. Mar. 29th
 */
#ifndef H_BSPSURFINTP
#define H_BSPSURFINTP
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <utility>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using namespace Eigen;

class BSplineSurfIntp;	// implicit declaration of the class
class MeshParams
{
public:
	virtual bool computeParams(const BSplineSurfIntp& surf, vector<double>& ret_u, vector<double>& ret_v) = 0;
};

class ChordParams : public MeshParams
{
public:
	bool computeParams(const BSplineSurfIntp& surf, vector<double>& ret_u, vector<double>& ret_v);
};
class CentrParams : public MeshParams
{
public:
	bool computeParams(const BSplineSurfIntp& surf, vector<double>& ret_u, vector<double>& ret_v);
};

class KnotSeq
{
public:
	virtual bool computeSeq(const BSplineSurfIntp& surf, vector<double>& ret_u_knot, vector<double>& ret_v_knot) = 0;
};

class AvgKnotSeq : public KnotSeq
{
public:
	bool computeSeq(const BSplineSurfIntp& surf, vector<double>& ret_u_knot, vector<double>& ret_v_knot);
};

class BSplineSurfIntp
{
public:
	BSplineSurfIntp():u_order(0), v_order(0), u_dat_num(0), v_dat_num(0){};
	BSplineSurfIntp(int u_o, int v_o, int u_dat, int v_dat);
	BSplineSurfIntp(const BSplineSurfIntp& rhs);
	BSplineSurfIntp& operator=(const BSplineSurfIntp& rhs);
	~BSplineSurfIntp(){};
	friend class MeshParams;

	bool compIntp();

	// set
	void setData(const vector<Vector3d>& _data);
	void setOrder(int _u, int _v);
	void setDataSize(int _u, int _v);
	void setParamFunc(MeshParams* _paramFunc){ paramFunc = _paramFunc; }
	void setKnotSeqFunc(KnotSeq* _knotSeqFunc){ knotSeqFunc = _knotSeqFunc; }
	bool setParamFunc(void (*paramFunc)(double*,int));
	bool setKnotSeqFunc(void (*knotFunc)(double*,int,int,int));
	
	// get
	int getUOrder() const { return u_order; }
	int getVOrder() const { return v_order; }
	int getUData() const { return u_dat_num; }
	int getVData() const { return v_dat_num; }
	const vector<Vector3d>& getControlPoints() const { return con; }
	const vector<Vector3d>& getDataPoints() const { return data; }
	const vector<double>& getUParam() const { return u_param; }
	const vector<double>& getVParam() const { return v_param; }
	const vector<double>& getUKnot() const { return u_knot; }
	const vector<double>& getVKnot() const { return v_knot; }

	// for debug
	void print();
	
private:
	int u_order, v_order, u_dat_num, v_dat_num;
	vector<Vector3d> data, con;		// data points are given, con will be obtained
	vector<double> u_param, v_param;
	vector<double> u_knot, v_knot;
	MeshParams* paramFunc;
	KnotSeq* knotSeqFunc;

	bool basis(double* ret, double t, int order, const vector<double>& knot) const;
};

#endif

