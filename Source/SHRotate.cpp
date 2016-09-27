
//
// Spherical Harmonic Rotation
//
// Implementation Reference:
//    Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion
//    Joseph Ivanic and Klaus Ruedenberg
//    J. Phys. Chem. 1996, 100, 6342-5347
//

#include "SHRotate.h"
#include <sdla/Matrix.h>
#include <sdla/Exception.h>
#include <cmath>
#include <cstring>

using namespace sdla;


namespace
{
	inline real delta(const int m, const int n)
	{
		// Kronecker Delta
		return (m == n ? 1 : 0);
	}


	void uvw(const int l, const int m, const int n, real& u, real& v, real& w)
	{
		// Pre-calculate simple reusable terms
		real d = delta(m, 0);
		int abs_m = abs(m);

		// Only calculate the required denominator once
		real denom;
		if (abs(n) == l)
			denom = (2 * l) * (2 * l - 1);

		else
			denom = (l + n) * (l - n);

		// Now just calculate the scalars
		u = sqrt((l + m) * (l - m) / denom);
		v = 0.5f * sqrt((1 + d) * (l + abs_m - 1) * (l + abs_m) / denom) * (1 - 2 * d);
		w = -0.5f * sqrt((l - abs_m - 1) * (l - abs_m) / denom) * (1 - d);
	}


	struct PermutedMatrix
	{
		PermutedMatrix(const cMatrix& m) : matrix(m) { }

		static int permute(const int v)
		{
			if (v == 1) return (0);
			if (v == -1) return (1);
			if (v == 0) return (2);
			ASSERT(false);
			return (0);
		}

		// m=row, n=column
		real operator () (const int m, const int n) const
		{
			return (matrix.e[permute(n)][permute(m)]);
		}

		const cMatrix& matrix;
	};


	real P(const int i, const int l, const int a, const int b, const PermutedMatrix& R, const SHRotateMatrix& M)
	{
		if (b == -l)
			return (R(i, 1) * M(l - 1, a, -l + 1) + R(i, -1) * M(l - 1, a, l - 1));

		else if (b == l)
			return (R(i, 1) * M(l - 1, a, l - 1) - R(i, -1) * M(l - 1, a, -l + 1));

		else // |b|<l
			return (R(i, 0) * M(l - 1, a, b));
	}


	real U(const int l, const int m, const int n, const PermutedMatrix& R, const SHRotateMatrix& M)
	{
		if (m == 0)
			return (P(0, l, 0, n, R, M));

		// For both m>0, m<0
		return (P(0, l, m, n, R, M));
	}


	real V(const int l, const int m, const int n, const PermutedMatrix& R, const SHRotateMatrix& M)
	{
		if (m == 0)
		{
			real p0 = P(1, l, 1, n, R, M);
			real p1 = P(-1, l, -1, n, R, M);
			return (p0 + p1);
		}

		else if (m > 0)
		{
			real d = delta(m, 1);
			real p0 = P(1, l, m - 1, n, R, M);
			real p1 = P(-1, l, -m + 1, n, R, M);
			return (p0 * sqrt(1 + d) - p1 * (1 - d));
		}

		else // m < 0
		{
			real d = delta(m, -1);
			real p0 = P(1, l, m + 1, n, R, M);
			real p1 = P(-1, l, -m - 1, n, R, M);
			return (p0 * (1 - d) + p1 * sqrt(1 - d));
		}
	}


	real W(const int l, const int m, const int n, const PermutedMatrix& R, const SHRotateMatrix& M)
	{
		if (m == 0)
		{
			// Never gets called as kd=0
			ASSERT(false);
			return (0);
		}

		else if (m > 0)
		{
			real p0 = P(1, l, m + 1, n, R, M);
			real p1 = P(-1, l, -m - 1, n, R, M);
			return (p0 + p1);
		}

		else // m < 0
		{
			real p0 = P(1, l, m - 1, n, R, M);
			real p1 = P(-1, l, -m + 1, n, R, M);
			return (p0 - p1);
		}
	}


	real M(const int l, const int m, const int n, const PermutedMatrix& R, const SHRotateMatrix& M)
	{
		// First get the scalars
		real u, v, w;
		uvw(l, m, n, u, v, w);

		// Scale by their functions
		if (u)
			u *= U(l, m, n, R, M);
		if (v)
			v *= V(l, m, n, R, M);
		if (w)
			w *= W(l, m, n, R, M);

		return (u + v + w);
	}
}


void SHRotate(SHRotateMatrix& shrm, const cMatrix& rotation)
{
	// Start with already known 1x1 rotation matrix for band zero
	memset(shrm.e, 0, shrm.DIM * shrm.DIM * sizeof(real));
	shrm.e[0][0] = 1;

	// Create matrix index modifiers
	PermutedMatrix pm(rotation);

	for (int m = -1; m <= 1; m++)
		for (int n = -1; n <= 1; n++)
			shrm(1, m, n) = pm(m, n);

	// Calculate each block of the rotation matrix for each band
	for (int band = 2; band < SampleSphere::NB_BANDS; band++)
	{
		for (int m = -band; m <= band; m++)
			for (int n = -band; n <= band; n++)
				shrm(band, m, n) = M(band, m, n, pm, shrm);
	}
}