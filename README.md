# BSI
B-Spline Surface Interpolation

BSI is a simple tool of interpolation for integral (not rational) B-Spline Surface.

Given data points in the form of _(m \* n)_ and degree of curves _(k,l)_ for each direction, control points of B-Spline surface lying on given points can be obtained. Additionally, parametrization for given data points and generating knot sequence is needed and it can be done by using a provided function or user-defined function.

To implement new method for parametrization or generating knot sequence, you can use the following:
```c++
class YourParameterizationClass : public MeshParams
{
	public:
		// implement this function
		bool computeParams(const BSplineSurfIntp& surf, vector<double>& ret_u, vector<double>& ret_v);
}

class YourKnotSequenceClass : public KnotSeq
{
	public:
		// implement this function
		bool computeSeq(const BSplineSurfIntp& surf, vector<double>& ret_u_knot, vector<double>& ret_v_knot);
}
```

You can find the example code is __main.cpp__. Please feel free to use, modify, or whatever!

## Eigen
This code uses [Eigen](http://eigen.tuxfamily.org/) to solve linear system.
