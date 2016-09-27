// SHRotMatrix.cpp
// Spherical Harmonic Rotation Matrix
// from "Evaluation of the rotation matrices in the basis of real spherical harmonics"
// M.A.Blanco et al, Journal of Molecular Structure (Theochem), 419, 1997, p.19-27
// rjgreen june 2002

#include <cmath>
#include <cstdio>

#include "SHRotMatrix.h"

// -----------------------------------------------

inline double sign(double x) { return x<0.0?-1.0:1.0; }
inline double negative_one_power(int x) { return (abs(x)&1)?-1.0:1.0; }
inline int row_offset(int l) { return 2*l+1; }

static double sinalpha = 0.0;
static double cosalpha = 0.0;
static double sinbeta = 0.0;
static double cosbeta = 0.0;
static double singamma = 0.0;
static double cosgamma = 0.0;
static double tanbetaover2 = 0.0;

// -----------------------------------------------

double clamp(double x, const double lo, const double hi)
{
	if (x < lo) x = lo;
	if (x > hi) x = hi;
	return (x);
}

void SHRotMatrix::alloc_mat(double_ptr* &m)
{
	// allocate and clear memory for matrices
	m = new double_ptr[highest_band];
	// *** TODO: add error checking
	for(int i=0; i<highest_band; ++i) {
		int n_terms = 2*i+1;
		m[i] = new double[n_terms * n_terms];
		for(int j=0; j<n_terms*n_terms; ++j) m[i][j] = 0.0;
	}
}

void SHRotMatrix::free_mat(double_ptr* &m)
{
	// free memory for matrices.
	for(int i=highest_band-1; i>=0; --i) {
		delete[] m[i];
	}
	delete[] m;
	m = 0;
}

void SHRotMatrix::dump_mat(const double_ptr* m)
{
	for(int i=0; i<highest_band; ++i) {
		for(int y=0; y<(2*i+1); ++y) {
			for(int x=0; x<(2*i+1); ++x) {
				printf("%7.3f, ", m[i][row_offset(i)*y+x]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void SHRotMatrix::store_d(double value, int l,int m,int n) 
{
	// *** TODO: sanity check indices here
	double *ptr = d_matrix[l];
	*(ptr + row_offset(l) * (m+l) + (n+l) ) = value;
	*(ptr + row_offset(l) * (l-n) + (l-m) ) = value;
	*(ptr + row_offset(l) * (l-m) + (l-n) ) = negative_one_power(m+n) * value;
	*(ptr + row_offset(l) * (n+l) + (m+l) ) = negative_one_power(m+n) * value;

// from Maple worksheet
//   d[ a+l,  b+l] := val;
//   d[-b+l, -a+l] := val;
//   d[-a+l, -b+l] := (-1)^(a+b) * val;
//   d[ b+l,  a+l] := (-1)^(a+b) * val;

}

double SHRotMatrix::get_d(int l,int m,int n)
{
	int x = m+l;
	int y = n+l;
	return d_matrix[l][row_offset(l)*x+y];
}

// -----------------------------------------------

double SHRotMatrix::eqn_64(int l,int m,int n)
{
	double ll = l*l;
	double lm1 = l-1;
	double twolm1 = 2*l-1;
	double lm12 = lm1*lm1;
	double mm = m*m;
	double nn = n*n;

	double a = ( l*twolm1 ) / sqrt( (ll-mm)*(ll-nn) );
	double b = ( get_d(1,0,0) - (m*n)/(l*lm1) ) * get_d(l-1, m, n);
	double c = ( sqrt((lm12-mm)*(lm12-nn))/(lm1*twolm1) ) * get_d(l-2,m,n);
	return a * ( b - c );
}

double SHRotMatrix::eqn_65(int l,int m,int n)
{
	return get_d(1,1,1) * get_d(l-1, l-1, l-1);
}

double SHRotMatrix::eqn_66(int l,int m,int n)
{
	return (l * get_d(1,0,0) - l + 1.0) * get_d(l-1, l-1, l-1);
}

double SHRotMatrix::eqn_67(int l,int m,int n)
{
	double lpm = l+m;
	double lmm = l-m;
	return -sqrt(double(lpm)/(lmm+1.0)) * tanbetaover2 * get_d(l,l,m);
}

double SHRotMatrix::eqn_68(int l,int m,int n)
{
	double lpm = l+m;
	double lmm = l-m;
	double a = l * cosbeta - m;
	return ((a+1.0)/a) * sqrt((lpm)/(lmm+1.0)) * tanbetaover2 * get_d(l, l-1, m);
}

// ---------------------------------------------

void SHRotMatrix::calc_d()
{
	// order 0 matrix
	store_d(1.0, 0,0,0);

	// order 1 matrix
	store_d(cosbeta , 1, 0, 0);									// Rzz
	store_d(-sqrt(0.5*(1.0-cosbeta*cosbeta)) , 1, 1, 0);	// -sqrt((1-Rzz^2)/2)
	store_d(0.5*(1.0+cosbeta) , 1, 1, 1);						// (1+Rzz)/2
	store_d(0.5*(1.0-cosbeta) , 1, 1,-1);						// (1-Rzz)/2
	tanbetaover2 = sqrt(get_d(1,1,-1) / get_d(1,1,1));		// used in eqn67 and eqn68

	// calc higher order matrices 2..highest_band
	for(int l=2; l<highest_band; ++l) {
		// ----------------
		for(int m=0; m<=(l-2); ++m) {
			for(int n=-m; n<=m; ++n) {
				store_d(eqn_64(l,m,n), l, m, n);
			}
		}
		store_d(eqn_65(l,l,l), l,l,l);
		store_d(eqn_66(l,l-1,l-1), l, l-1, l-1);
		for(int n=l-1; n>=-l; --n) {
			store_d(eqn_67(l,l,n), l,l,n);
		}
		for(int n=l-2; n>=1-l; --n) {
			store_d(eqn_68(l,l-1,n), l,l-1,n);
		}
		// ----------------
	}
}

double sinmx(int m, double sinx, double cosx)
{
	double prev_sinmx = 0.0;
	double prev_cosmx = 1.0;
	for(int i=0; i<m; ++i) {
		double new_sinmx = sinx*prev_cosmx + cosx*prev_sinmx;
		double new_cosmx = cosx*prev_cosmx - sinx*prev_sinmx;
		prev_sinmx = new_sinmx;
		prev_cosmx = new_cosmx;
	}
	return prev_sinmx;
}

double cosmx(int m, double sinx, double cosx)
{
	double prev_sinmx = 0.0;
	double prev_cosmx = 1.0;
	for(int i=0; i<m; ++i) {
		double new_sinmx = sinx*prev_cosmx + cosx*prev_sinmx;
		double new_cosmx = cosx*prev_cosmx - sinx*prev_sinmx;
		prev_sinmx = new_sinmx;
		prev_cosmx = new_cosmx;
	}
	return prev_cosmx;
}

double sigma_alpha(int m)
{
	const double sqrt2 = sqrt(2.0f);
	if(m>0) return sqrt2 * cosmx(m, sinalpha, cosalpha);
	else if(m==0) return 1.0;
	else return sqrt2 * sinmx(-m, sinalpha, cosalpha);
}

double sigma_gamma(int m)
{
	const double sqrt2 = sqrt(2.0f);
	if(m>0) return sqrt2 * cosmx(m, singamma, cosgamma);
	else if(m==0) return 1.0;
	else return sqrt2 * sinmx(-m, singamma, cosgamma);
}

double SHRotMatrix::delta(int l,int m,int n)
{
	double a = sign(n) * sigma_alpha(m) * sigma_gamma(n);
	double b = 0.5 * ( get_d(l,abs(n),abs(m)) + negative_one_power(m) * get_d(l,abs(m),-abs(n)) );
	double c = sign(m) * sigma_alpha(-m) * sigma_gamma(-n);
	double d = 0.5 * ( get_d(l,abs(n),abs(m)) - negative_one_power(m) * get_d(l,abs(m),-abs(n)) );
	return a*b - c*d;
}

// -----------------------------------------------
// 5=1,1=y,y
// 4=1,0=y,x
// 8=2,0=z,x
void SHRotMatrix::rotate_sh(
		const double rot[],
		const int order,
		const double sh[],
		double result[] )
{
	const double epsilon = 1e-8;

	highest_band = order;

	// calc trig values from rotation matrix
	// NOTE: assumes 4x4 input
	cosbeta = rot[10];							// Rzz
	sinbeta = sqrt(1.0-rot[10]*rot[10]);	// sqrt(1-Rzz^2);
	if(sinbeta<epsilon && sinbeta>-epsilon) {
		// sinbeta == 0.0 so fake the results.
		cosbeta = 1.0;
		sinbeta = 0.0;
		cosgamma = rot[5];	//  Ryy
		singamma = -rot[4];	// -Ryx
		cosalpha = 1.0;
		sinalpha = 0.0;
	} else {
		cosgamma =  rot[ 8]/sinbeta;			//  Rzx/sqrt(1-Rzz^2)
		singamma =  rot[ 9]/sinbeta;			//  Rzy/sqrt(1-Rzz^2)
		cosalpha = -rot[ 2]/sinbeta;			// -Rxz/sqrt(1-Rzz^2)
		sinalpha =  rot[ 6]/sinbeta;			//  Ryz/sqrt(1-Rzz^2)
	}

	// account for small sinbeta giving screwy results
	cosalpha = clamp(cosalpha, -1.0, 1.0);
	sinalpha = clamp(sinalpha, -1.0, 1.0);
	cosgamma = clamp(cosgamma, -1.0, 1.0);
	singamma = clamp(singamma, -1.0, 1.0);

	// calculate d-matrix
	calc_d();

	// calculate sh rotation matrices
	for(int l=0; l<highest_band; ++l) {
		for(int m=-l; m<=l; ++m) {
			for(int n=-l; n<=l; ++n) {
				sh_matrix[l][row_offset(l)*(m+l) + (n+l)] = delta(l,m,n);
			}
		}
	}

	//dump_mat(sh_matrix);

	// transform the input sh vector
	result[0] = sh[0]; 	// no need to transform order 0
	for(int i=1; i<order; ++i) {
		// matrix-vector multiply
		for(int m=0; m<(2*i+1); ++m) {
		   double value = 0.0;
			for(int n=0; n<(2*i+1); ++n) {
				value += sh_matrix[i][row_offset(i)*m + n] * sh[n+(i*i)];
			}
			result[i*i + m] = value;
		}
	}

}

SHRotMatrix::SHRotMatrix(int order)
{
	// allocate and clear space for matrices.
	highest_band = order;
	
	alloc_mat(d_matrix);
	alloc_mat(sh_matrix);
}


SHRotMatrix::~SHRotMatrix()
{
	free_mat(sh_matrix);
	free_mat(d_matrix);
}

// -----------------------------------------------