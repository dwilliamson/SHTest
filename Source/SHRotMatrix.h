// SHRotMatrix.h
// Spherical Harmonic Rotation Matrix
// rjgreen june 2002

#ifndef _SHROTMATRIX_H
#define _SHROTMATRIX_H

typedef double* double_ptr;

class SHRotMatrix
{
public:
	SHRotMatrix(int i);
	~SHRotMatrix();

	void rotate_sh( const double rot[],	// 4x4 rotation matrix
						 const int order,		// number of bands in the SH vector
						 const double sh[],	// SH vector with order^2 values
						 double result[] );	// SH vector to hold the result

private:
	SHRotMatrix();												// no unparametarised constructor
	SHRotMatrix(const SHRotMatrix& m);					// no copy constructor
	SHRotMatrix &operator=(const SHRotMatrix &m);	// no asignment

	void alloc_mat(double_ptr* &mat);
   void free_mat(double_ptr* &mat);
	void dump_mat(const double_ptr* mat);

	void store_d(double value, int l,int m,int n);
	double get_d(int l,int m,int n);
	void calc_d();
	double sigma(int m,double theta);
	double eqn_64(int l,int m,int n);
	double eqn_65(int l,int m,int n);
	double eqn_66(int l,int m,int n);
	double eqn_67(int l,int m,int n);
	double eqn_68(int l,int m,int n);
	double delta(int l,int m,int n);
	
	int highest_band;
	double **d_matrix;	// array of d matrices
	double **sh_matrix;	// array of sh rotaton matrices
};

#endif