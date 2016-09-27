
#pragma once


#include "Types.h"


// Template meta-program to generate offsets into a linear array of SH sub-matrices
template <int D> struct SHRM_LookupEntry
{
	enum { val = D * D + SHRM_LookupEntry<D - 2>::val };
};
template <> struct SHRM_LookupEntry<1>
{
	enum { val = 0 };
};

// Sub-matrix offsets for a maximum of 12 bands minus the first
static const int SHRM_MATRIX_LOOKUP[] =
{
	0,							// 2 bands
	SHRM_LookupEntry<3>::val,	// 3 bands
	SHRM_LookupEntry<5>::val,	// 4 bands
	SHRM_LookupEntry<7>::val,	// ...
	SHRM_LookupEntry<9>::val,
	SHRM_LookupEntry<11>::val,
	SHRM_LookupEntry<13>::val,
	SHRM_LookupEntry<15>::val,
	SHRM_LookupEntry<17>::val,
	SHRM_LookupEntry<19>::val,
	SHRM_LookupEntry<21>::val,
	SHRM_LookupEntry<23>::val
};


class cSHRotateMatrix
{
public:
	class SubMatrix
	{
	public:
		// Constructor that needs the width of the sub-matrix
		SubMatrix(const int width) : m_Elements(0), m_Width(width), m_Size(width * width), m_Shift((width - 1) / 2)
		{
			m_Elements = new double[m_Size];
		}

		// Cleanup destructor
		~SubMatrix(void)
		{
			if (m_Elements)
				delete [] m_Elements;
		}

		// Element accessors with (-,+) lookup
		double& operator () (const int m, const int n)
		{
			return (m_Elements[(m + m_Shift) * m_Width + (n + m_Shift)]);
		}
		double operator () (const int m, const int n) const
		{
			return (m_Elements[(m + m_Shift) * m_Width + (n + m_Shift)]);
		}

	private:

		// Element array
		double*	m_Elements;

		// Width of the matrix
		int		m_Width;

		// Total size of the matrix in multiples of doubles
		int		m_Size;

		// Value to shift incoming matrix indices by
		int		m_Shift;
	};


	cSHRotateMatrix(const int order) :

		m_Order(order),
		m_Matrices(0)

	{
		m_Matrices = new SubMatrix*[m_Order - 1];
		for (int i = 1; i < m_Order; i++)
			m_Matrices[i - 1] = new SubMatrix(i * 2 + 1);
	}

	~cSHRotateMatrix(void)
	{
		if (m_Matrices)
			delete [] m_Matrices;
	}


	double& operator () (const int l, const int m, const int n)
	{
		return ((*m_Matrices[l - 1])(m, n));
	}
	double operator () (const int l, const int m, const int n) const
	{
		return ((*m_Matrices[l - 1])(m, n));
	}

	void Transform(const RGBCoeffList& source, RGBCoeffList& dest) const
	{
		// ... Need to check here that input coeff-list size matches the matrix ...

		dest.e[0] = source.e[0];

		// Loop through each band
		for (int l = 1; l < m_Order; l++)
		{
			// Now through each argument of the destination coeff-list
			for (int mo = -l; mo <= l; mo++)
			{
				RGBCoeff& rgb = dest(l, mo);
				rgb.r = rgb.g = rgb.b = 0;

				// Multiply-add with each argument off the source coeff-list
				for (int mi = -l; mi <= l; mi++)
					rgb += source(l, mi) * (*m_Matrices[l - 1])(mo, mi);
			}
		}
	}

	int GetOrder(void) const
	{
		return (m_Order);
	}

private:
	// Number of SH bands 
	int		m_Order;

	// Matrices for each band
	SubMatrix**	m_Matrices;
};