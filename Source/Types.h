
#pragma once


#include <cmath>


typedef double	real;
//typedef float	real;


struct Sample2D
{
	real	x, y;
};


struct SampleSphere
{
	static const int	NB_BANDS = 4;
	static const int	NB_BASES = NB_BANDS * NB_BANDS;

	SampleSphere(void) { }

	SampleSphere(const real u, const real v) :
		theta(u),
		phi(v),
		index(-1)
	{
		x = sin(theta) * cos(phi);
		y = sin(theta) * sin(phi);
		z = cos(theta);
	}

	// Spherical co-ordinates
	real	theta, phi;

	// Cartesian co-ordinates
	real	x, y, z;

	// y(l,m,theta,phi) for each l,m
	real	coeffs[NB_BASES];

	// Sample index in a set
	int		index;
};


struct RGBCoeff
{
	union
	{
		real	e[3];

		struct
		{
			real	r, g, b;
		};
	};

	RGBCoeff(void) { }
	RGBCoeff(const real _r, const real _g, const real _b) : r(_r), g(_g), b(_b) { }
	RGBCoeff(const RGBCoeff& rhs) : r(rhs.r), g(rhs.g), b(rhs.b) { }

	RGBCoeff& operator += (const RGBCoeff& rhs)
	{
		r += rhs.r;
		g += rhs.g;
		b += rhs.b;
		return (*this);
	}
};


inline RGBCoeff operator * (const RGBCoeff& lhs, const real rhs)
{
	return (RGBCoeff(lhs.r * rhs, lhs.g * rhs, lhs.b * rhs));
}


inline RGBCoeff operator + (const RGBCoeff& lhs, const RGBCoeff& rhs)
{
	return (RGBCoeff(lhs.r + rhs.r, lhs.g + rhs.g, lhs.b + rhs.b));
}


struct RGBCoeffList
{
	RGBCoeff	e[SampleSphere::NB_BASES];

	RGBCoeff& operator () (const int l, const int m)
	{
		return (e[l * (l + 1) + m]);
	}
	const RGBCoeff& operator () (const int l, const int m) const
	{
		return (e[l * (l + 1) + m]);
	}
};


struct SHRotateMatrix
{
	// The size of the matrix is the square of the number of bands
	static const int DIM = SampleSphere::NB_BANDS * SampleSphere::NB_BANDS;

	// Matrix array
	real e[DIM][DIM];

	// SH rotation  matrices are block diagonal sparse.
	// Given a band index this will offset the lookup into the block within the matrix
	// that is responsible for rotating that band and shift it by (a,b).
	real& operator () (const int l, const int a, const int b)
	{
		int centre = (l + 1) * l;
		return (e[centre + a][centre + b]);
	}
	real operator () (const int l, const int a, const int b) const
	{
		int centre = (l + 1) * l;
		return (e[centre + a][centre + b]);
	}

	void Transform(const RGBCoeffList& source, RGBCoeffList& dest) const
	{
		// Loop through each band
		for (int l = 0; l < SampleSphere::NB_BANDS; l++)
		{
			for (int mo = -l; mo <= l; mo++)
			{
				RGBCoeff& rgb = dest(l, mo);
				rgb.r = rgb.g = rgb.b = 0;

				for (int mi = -l; mi <= l; mi++)
					rgb += source(l, mi) * (*this)(l, mo, mi);
			}

		}

		/*for (int i = 0; i < SampleSphere::NB_BASES; i++)
		{
			dest.e[i].r = 0;
			dest.e[i].g = 0;
			dest.e[i].b = 0;

			for (int j = 0; j < DIM; j++)
				dest.e[i] += source.e[j] * e[j][i];
		}*/
	}
};