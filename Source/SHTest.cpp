
#include "SHTest.h"
#include "BitArray.h"
#include "SHRotate.h"
#include "SHRotMatrix.h"
#include "SHRotateMatrix.h"
#include <sdla/Maths.h>
#include <sdla/Font.h>
#include <sdla/ComponentManager.h>
#include <sdla/MeshLoaderFactory.h>
#include <sdla/ImageLoaderFactory.h>
#include <sdla/PolyMesh.h>
#include <sdla/DrawPolyMesh.h>
#include <sdla/VoxelRaytracer.h>
#include <cstdlib>
#include <ctime>

using namespace sdla;


#define USE_SSS_HACK


// http://mathworld.wolfram.com/LegendrePolynomial.html
// http://mathworld.wolfram.com/SphericalHarmonic.html
// http://mathworld.wolfram.com/DoubleFactorial.html
// http://www.gnu.org/software/gsl/
// http://wwwinfo.cern.ch/asd/lhc++/


namespace
{
	const int		SQRT_NB_SAMPLES = 20;
	const int		MAX_NB_SAMPLES = SQRT_NB_SAMPLES * SQRT_NB_SAMPLES;

	const int		NB_INDIRECT_PASSES = 3;


	real y(const int l, const int m, const real theta, const real phi);


	// Generate random number between 0 and 1
	real Random(void)
	{
		return ((real)rand() / RAND_MAX);
	}

	
	// Generate simple random 2D samples with high variance
	void GenerateSamples2D(Sample2D* dest, const int nb_samples)
	{
		for (int i = 0; i < nb_samples; i++)
		{
			dest[i].x = Random();
			dest[i].y = Random();
		}
	}


	// Reduce the variance of the generated samples by breaking the sample space
	// into an NxN grid, iterating through each grid cell, and generating a random
	// sample within that
	void GenerateSamples2D_Stratified(Sample2D* dest, const int sqrt_nb_samples)
	{
		for (int i = 0, k = 0; i < sqrt_nb_samples; i++)
		{
			for (int j = 0; j < sqrt_nb_samples; j++, k++)
			{
				dest[k].x = (i + Random()) / (real)sqrt_nb_samples;
				dest[k].y = (j + Random()) / (real)sqrt_nb_samples;
			}
		}
	}


	void DrawSamples2D(Sample2D* samples, const int nb_samples, const cMatrix& m)
	{
		glMatrixMode(GL_MODELVIEW);
		glLoadMatrix(cMatrix::Translate(cVector(-0.5f, -0.5f, 0)) * m * cMatrix::Translate(cVector(0, 0, 2)));

		glBegin(GL_POINTS);
		for (int i = 0; i < nb_samples; i++)
			glVertex3f((float)samples[i].x, (float)samples[i].y, 0);
		glEnd();
	}

	
	// Map a set of 2D samples into spherical co-ordinates
	// Also generate their position on the sphere in 3D
	void ProjectSpherical(SampleSphere* dest, const Sample2D* source, const int nb_samples)
	{
		for (int i = 0; i < nb_samples; i++)
		{
			dest[i] = SampleSphere(2 * acos(sqrt(1 - source[i].x)), 2 * PI * source[i].y);

			// Pre-calc the SH basis function values at this sample
			for (int l = 0, c = 0; l < SampleSphere::NB_BANDS; l++)
				for (int m = -l; m <= l; m++, c++)
					dest[i].coeffs[c] = y(l, m, dest[i].theta, dest[i].phi);

			dest[i].index = i;
		}
	}


	void DrawSamplesSphere(SampleSphere* samples, const int nb_samples, const cMatrix& m)
	{
		glMatrixMode(GL_MODELVIEW);
		glLoadMatrix(m * cMatrix::Translate(cVector(0, 0, 2)));

		glBegin(GL_POINTS);
		for (int i = 0; i < nb_samples; i++)
			glVertex3f((float)samples[i].x, (float)samples[i].y, (float)samples[i].z);
		glEnd();
	}


	void DrawSamples(Sample2D* samples_2d, SampleSphere* samples_sphere, const cMatrix& m)
	{
		glDisable(GL_LIGHTING);

		glPointSize(2);
		glColor3f(1, 1, 1);

		DrawSamples2D(samples_2d, MAX_NB_SAMPLES, m);
		DrawSamplesSphere(samples_sphere, MAX_NB_SAMPLES, m);
	}


	// Defined as: n!! = n.(n - 2).(n - 4)..., n!!(0,-1)=1
	int DoubleFactorial(int x)
	{
		if (x == 0 || x == -1)
			return (1);

		int result = x;
		while ((x -= 2) > 0)
			result *= x;
		return (result);
	}


	class LegendreNormal;
	class LegendreOptimised;

	typedef LegendreNormal	Legendre;


	// Calculate Associated Legendre Polynomial
	template <typename T> real P(const int l, const int m, const real x)
	{
	}


	// Standard method using recursion, no optimisations
	template <> real P<LegendreNormal>(const int l, const int m, const real x)
	{
		// Rule 2 needs no previous results
		if (l == m)
			return (pow(-1.0f, m) * DoubleFactorial(2 * m - 1) * pow(sqrt(1 - x * x), m));

		// Rule 3 requires the result for the same argument of the previous band
		if (l == m + 1)
			return (x * (2 * m + 1) * P<LegendreNormal>(m, m, x));

		// Main reccurence used by rule 1 that uses result of the same argument from
		// the previous two bands
		return ((x * (2 * l - 1) * P<LegendreNormal>(l - 1, m, x) - (l + m - 1) * P<LegendreNormal>(l - 2, m, x)) / (l - m));
	}


	// Optimised, non-recursive method
	template <> real P<LegendreOptimised>(const int l, const int m, const real x)
	{
		// Start with P(0,0) at 1
		real pmm = 1;

		// First calculate P(m,m) since that is the only rule that requires results
		// from previous bands
		// No need to check for m>0 since SH function always gives positive m

		// Precalculate (1 - x^2)^0.5
		real somx2 = sqrt(1 - x * x);

		// This calculates P(m,m). There are three terms in rule 2 that are being iteratively multiplied:
		//
		// 0: -1^m
		// 1: (2m-1)!!
		// 2: (1-x^2)^(m/2)
		//
		// Term 2 has been partly precalculated and the iterative multiplication by itself m times
		// completes the term.
		// The result of 2m-1 is always odd so the double factorial calculation multiplies every odd
		// number below 2m-1 together. So, term 3 is calculated using the 'fact' variable.
		real fact = 1;
		for (int i = 1; i <= m; i++)
		{
			pmm *= -1 * fact * somx2;
			fact += 2;
		}

		// No need to go any further, rule 2 is satisfied
		if (l == m)
			return (pmm);

		// Since m<l in all remaining cases, all that is left is to raise the band until the required
		// l is found

		// Rule 3, use result of P(m,m) to calculate P(m,m+1)
		real pmmp1 = x * (2 * m + 1) * pmm;

		// Is rule 3 satisfied?
		if (l == m + 1)
			return (pmmp1);

		// Finally, use rule 1 to calculate any remaining cases
		real pll = 0;
		for (int ll = m + 2; ll <= l; ll++)
		{
			// Use result of two previous bands
			pll = (x * (2.0 * ll - 1.0) * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);

			// Shift the previous two bands up
			pmm = pmmp1;
			pmmp1 = pll;
		}

		return (pll);
	}


	void DrawGrid(const float r, const float g, const float b)
	{
		struct V
		{
			V(void) : swap(false) { }
			bool swap;

			void operator () (const float x, const float y)
			{
				if (swap)
					glVertex3f(y, x, 0);
				else
					glVertex3f(x, y, 0);
			}
		} vertex;

		// Position in world space
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glBegin(GL_LINES);
		const float STEP = 0.1f;
		for (int d = 0; d < 2; d++)
		{
			for (float i = -1.5f; i < 1.5f; i += STEP)
			{
				const float GRID_SIZE = 1.5f;
				const float MEDIUM_SIZE = 0.05f;
				const float MINI_SIZE = 0.01f;

				// The main grid
				glColor3f(r, g, b);
				vertex(i, GRID_SIZE);
				vertex(i, -GRID_SIZE);

				// The grid cell sized marks
				glColor3f(1, 1, 1);
				vertex(i, MEDIUM_SIZE);
				vertex(i, -MEDIUM_SIZE);

				// The mini-marks inbetween the grid cells
				const float SMALL_STEP = 0.025f;
				for (float i_small = i + SMALL_STEP; i_small < i + STEP; i_small += SMALL_STEP)
				{
					vertex(i_small, MINI_SIZE);
					vertex(i_small, -MINI_SIZE);
				}
			}

			glColor3f(1, 1, 1);
			vertex(0, 1.5f);
			vertex(0, -1.5f);

			vertex.swap = true;
		}
		glEnd();
	}


	cColour g_LColours[] =
	{
		#define C 1
		cColour(C, 0, 0),
		cColour(0, C, 0),
		cColour(0, 0, C),
		cColour(C, C, 0),
		cColour(C, 0, C),
		cColour(0, 1, 1),
		#undef	C
		#define C 0.5f
		cColour(C, 0, 0),
		cColour(0, C, 0),
		cColour(0, 0, C),
		cColour(C, C, 0),
		cColour(C, 0, C),
		cColour(0, 1, 1),
		#undef	C
	};	


	void DrawLegendre(const int l, const int m, cFont* font_ptr, const int colour)
	{
		const float STEP = 0.01f;

		// Position in world space
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		real x, y;

		cColour& c = g_LColours[colour];
		glColor3f(c.r, c.g, c.b);

		glBegin(GL_LINES);
		for (real i = -1; i < 1 - STEP; i += STEP)
		{
			x = i;
			y = P<Legendre>(l, m, x) / 10;
			glVertex3d(x, y, 0);

			x += STEP;
			y = P<Legendre>(l, m, x) / 10;
			glVertex3d(x, y, 0);
		}
		glEnd();

		font_ptr->WriteText(0.8f, 0.8f - font_ptr->GetLineHeight() * colour, "P(%d,%d)", l, m);
	}


	// Basic integer factorial
	int Factorial(int v)
	{
		if (v == 0)
			return (1);

		int result = v;
		while (--v > 0)
			result *= v;
		return (result);
	}


	// Re-normalisation constant for SH function
	// Wonder what all this would look like in C++ TMP :)
	real K(const int l, const int m)
	{
		// Note that |m| is not used here as the SH function always passes positive m
		return (sqrt(((2 * l + 1) * Factorial(l - m)) / (4 * PI * Factorial(l + m))));
	}


	// Return a point sample of an SH basis function
	// l is the band
	// m is the argument, in the range [-l...l]
	// theta is in the range [0...PI] (altitude)
	// phi is in the range [0...2PI] (azimuth)
	real y(const int l, const int m, const real theta, const real phi)
	{	
		if (m == 0)
			return (K(l, 0) * P<Legendre>(l, 0, cos(theta)));

		if (m > 0)
			return (sqrt(2.0f) * K(l, m) * cos(m * phi) * P<Legendre>(l, m, cos(theta)));

		// m < 0, m is negated in call to K
		return (sqrt(2.0f) * K(l, -m) * sin(-m * phi) * P<Legendre>(l, -m, cos(theta)));
	}


	cVector CalculateNormal(const cVector& a, const cVector& b, const cVector&c )
	{
		// Get edges
		cVector e0 = b - a;
		cVector e1 = c - a;

		// Get perpendicular
		return ((e0 ^ e1).Normalise());
	}


	struct SphericalFunction
	{
		virtual real ToReal(const SampleSphere& s) const = 0;
		virtual cColour ToColour(const SampleSphere& s) const = 0;

		// Convert radius to vector
		// Map from spherical co-ordinates onto the sphere using the unsigned sample as the radius
		cVector ToVector(const real u, const real v)
		{
			SampleSphere s(u, v);

			// Sample the function
			real r = fabs(ToReal(s));

			// Map from spherical co-ordinates onto the sphere using the unsigned sample as the radius
			cVector vec;
			vec.x = (float)(r * s.x);
			vec.y = (float)(r * s.y);
			vec.z = (float)(r * s.z);

			return (vec);
		}
	};


	// For drawing an SH basis
	struct SHSample : public SphericalFunction
	{
		SHSample(const int _l, const int _m) : l(_l), m(_m) { }

		real ToReal(const SampleSphere& s) const
		{
			return (y(l, m, s.theta, s.phi));
		}

		cColour ToColour(const SampleSphere& s) const
		{
			return (cColour((float)ToReal(s)));
		}

		int l, m;
	};


	void DrawSphericalFunction(SphericalFunction& func, const int name)
	{
		const int	NB_SAMPLES = 50;

		real du = PI / NB_SAMPLES;
		real dv = 2 * PI / NB_SAMPLES;

		glNewList(name, GL_COMPILE);

		glBegin(GL_QUADS);
		for (real u = 0; u < PI; u += du)
		{
			for (real v = 0; v < 2 * PI; v += dv)
			{
				// Sample 4 points in the neighbourhood
				cVector p[4] =
				{
					func.ToVector(u, v),
					func.ToVector(u + du, v),
					func.ToVector(u + du, v + dv),
					func.ToVector(u, v + dv)
				};

				real ddu = du / 10;
				real ddv = dv / 10;

				// Sample the normals as even smaller patches in the neighbourhood
				cVector n[4] =
				{
					CalculateNormal(p[0], func.ToVector(u + ddu, v), func.ToVector(u, v + ddv)),
					CalculateNormal(p[1], func.ToVector(u + du + ddu, v), func.ToVector(u + du, v + ddv)),
					CalculateNormal(p[2], func.ToVector(u + du + ddu, v + dv), func.ToVector(u + du, v + dv + ddv)),
					CalculateNormal(p[3], func.ToVector(u + ddu, v + dv), func.ToVector(u, v + dv + ddv))
				};

				// Positive SH samples are green, negative are red
				float c[3] = { 0, 0, 0 };
				if (func.ToReal(SampleSphere(u, v)) < 0)
					c[0] = 1;
				else
					c[1] = 1;
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);

				// Draw the quad into the display list
				glNormal3fv(n[0]);
				glVertex3fv(p[0]);
				glNormal3fv(n[1]);
				glVertex3fv(p[1]);
				glNormal3fv(n[2]);
				glVertex3fv(p[2]);
				glNormal3fv(n[3]);
				glVertex3fv(p[3]);
			}
		}
		glEnd();

		glEndList();
	}


	template <typename T> T min(const T& a, const T& b)
	{
		return (a > b ? b : a);
	}


	template <typename T> T max(const T& a, const T& b)
	{
		return (a < b ? b : a);
	}


	struct ExampleLight : public SphericalFunction
	{
		real ToReal(const SampleSphere& s) const
		{
			// Example lighting function provided in documentation
			//double theta = s.theta + 3.141f;
			//return (max(0.0, 5 * cos(s.theta) - 4));
			return (max(0.0, 5 * cos(s.theta) - 4) + max(0.0, -4 * sin(s.theta - PI) * cos(s.phi - 2.5) - 3));
			//return (max(0.0, 10 * cos(theta) - 4) + max(0.0, -4 * sin(theta - PI) * cos(s.phi - 2.5) - 3));
		}

		cColour ToColour(const SampleSphere& s) const
		{
			return (cColour((float)ToReal(s)));
		}
	};


	struct ReconstructSH : public SphericalFunction
	{
		ReconstructSH(const real* _coeffs)
		{
			// Copy the SH basis weights
			for (int i = 0; i < SampleSphere::NB_BASES; i++)
				coeffs[i] = _coeffs[i];
		}

		real ToReal(const SampleSphere& s) const
		{
			real value = 0;

			// Calculate a weighted sum of all the SH bases
			for (int l = 0, i = 0; l < SampleSphere::NB_BANDS; l++)
				for (int m = -l; m <= l; m++, i++)
					value += coeffs[i] * y(l, m, s.theta, s.phi);

			return (value);
		}

		cColour ToColour(const SampleSphere& s) const
		{
			return (cColour((float)ToReal(s)));
		}

		real coeffs[SampleSphere::NB_BASES];
	};


	// Project a spherical function onto the SH bases
	void ProjectSphericalFunction(SphericalFunction& func, const SampleSphere* samples, real* coeffs)
	{
		// Clear out sums
		for (int i = 0; i < SampleSphere::NB_BASES; i++)
			coeffs[i] = 0;

		for (int i = 0; i < MAX_NB_SAMPLES; i++)
		{
			// Take the sample at this point on the sphere
			const SampleSphere& s = samples[i];
			real sample = func.ToReal(s);

			// Project the sample onto each SH basis
			if (sample)
				for (int j = 0; j < SampleSphere::NB_BASES; j++)
					coeffs[j] += sample * s.coeffs[j];
		}

		// Divide each coefficient by number of samples and multiply by weights
		for (int i = 0; i < SampleSphere::NB_BASES; i++)
			coeffs[i] = coeffs[i] * ((4 * PI) / MAX_NB_SAMPLES);
	}


	void RGBProjectSphericalFunction(SphericalFunction& func, const SampleSphere* samples, RGBCoeffList& coeffs)
	{
		// Clear out sums
		for (int i = 0; i < SampleSphere::NB_BASES; i++)
			coeffs.e[i].r = coeffs.e[i].g = coeffs.e[i].b = 0;

		// For each sample sum the product of the lighting function and each SH basis
		for (int i = 0; i < MAX_NB_SAMPLES; i++)
		{
			const SampleSphere& s = samples[i];
			cColour sample = func.ToColour(s);

			if (sample.r || sample.g || sample.b)
			{
				for (int j = 0; j < SampleSphere::NB_BASES; j++)
				{
					// Project onto band/arg
					cColour colour = sample * (float)s.coeffs[j];

					// Sum
					coeffs.e[j].r += colour.r;
					coeffs.e[j].g += colour.g;
					coeffs.e[j].b += colour.b;
				}
			}
		}

		// Divide each coefficient by number of samples and multiply by weights
		for (int i = 0; i < SampleSphere::NB_BASES; i++)
			coeffs.e[i] = coeffs.e[i] * ((4 * PI) / MAX_NB_SAMPLES);
	}


	struct DiffuseUnshadowedTransfer : public SphericalFunction
	{
		DiffuseUnshadowedTransfer(const cPolyMesh::Vertex& vertex, const cColour& _colour) :
			normal(vertex.normal),
			colour(_colour)
		{
		}

		real ToReal(const SampleSphere& s) const
		{
			// Calculate cosine term for this sample
			float H = cVector((float)s.x, (float)s.y, (float)s.z) * normal;

			// Only valid if in the upper hemisphere
			if (H > 0)
				return (H);

			return (0);
		}

		cColour ToColour(const SampleSphere& s) const
		{
			return (colour * (float)ToReal(s));
		}

		// Normal at the sample point
		cVector normal;

		// Colour at the sample point
		cColour colour;
	};


	struct DiffuseShadowedTransfer : public SphericalFunction
	{
		DiffuseShadowedTransfer(const cPolyMesh::Vertex& v, const cVoxelRaytracer& rt, const cColour& c) :
			vertex(v),
			raytracer(rt),
			colour(c)
		{
		}

		real ToReal(const SampleSphere& s) const
		{
			// Calculate cosine term for this sample
			cVector dir((float)s.x, (float)s.y, (float)s.z);
			float H = dir * vertex.normal;

			float total = 0;

#ifdef	USE_SSS_HACK
			// Only valid if in the lower hemisphere
			if (H < 0)
			{
				// Trace back into the model looking for the smallest depth (won't work correctly if the ray reenters the model)
				float t, u, v;
				raytracer.GetClosestHitTriangle(vertex.position, dir, false, t, u, v);

				// Convert the depth into contributed radiance
				t = std::min(std::max(0.0f, t), 1.0f);
				t = 1 - t * 4;
				t = std::min(std::max(0.0f, t), 1.0f);
				t /= 4;
				total += t;
			}
#endif

			// Only valid if in the upper hemisphere
			if (H > 0)
			{
				// If not self-shadowing use the cosine term
				if (!raytracer.RayHitsTriangle(vertex.position, dir, false))
					total += H;
			}

			return (total);
		}

		cColour ToColour(const SampleSphere& s) const
		{
			return (colour * (float)ToReal(s));
		}

		// Vertex being processed
		const cPolyMesh::Vertex& vertex;

		// Raytracer used for visibility term
		const cVoxelRaytracer& raytracer;

		// Colour at the sample point
		cColour colour;
	};


	struct DiffuseInterreflectedTransfer
	{
		// Direct light in the form of diffuse shadowed transfer
		struct DirectPass : public SphericalFunction
		{
			DirectPass(const cPolyMesh::Vertex& v, const cVoxelRaytracer& rt, const cColour& c, cBitArray& rh, const int bi) :
				vertex(v),
				raytracer(rt),
				colour(c),
				rayhits(rh),
				base_index(bi)
			{
			}

			real ToReal(const SampleSphere& s) const
			{
				// Calculate cosine term for this sample
				cVector dir(float(s.x), float(s.y), float(s.z));
				float H = dir * vertex.normal;

				float total = 0;

#ifdef	USE_SSS_HACK
				// Only valid if in the lower hemisphere
				if (H < 0)
				{
					// Trace back into the model looking for the smallest depth (won't work correctly if the ray reenters the model)
					float t, u, v;
					raytracer.GetClosestHitTriangle(vertex.position, dir, false, t, u, v);

					// Convert the depth into contributed radiance
					t = std::min(std::max(0.0f, t), 1.0f);
					t = 1 - t * 4;
					t = std::min(std::max(0.0f, t), 1.0f);
					t /= 4;
					total += t;
				}
#endif

				// Only valid if in the upper hemisphere
				if (H > 0)
				{
					// Get closest intersection with self
					float t, u, v;
					int hit = raytracer.GetClosestHitTriangle(vertex.position, dir, false, t, u, v);

					// Only mark those triangles facing the vertex as hit as those are the only causes
					// of inter-reflectance. This reduces the hit count by a large amount.
					if (hit != -1 && raytracer.GetTriangle(hit).normal * dir < 0)
						rayhits[base_index + s.index] = 1;

					// If not self-shadowing use the cosine term
					if (hit == -1)
						total += H;
				}

				return (total);
			}

			cColour ToColour(const SampleSphere& s) const
			{
				return (colour * (float)ToReal(s));
			}

			// Vertex being processed
			const cPolyMesh::Vertex& vertex;

			// Raytracer used for indirect lighting
			const cVoxelRaytracer& raytracer;

			// Colour at the sample point
			cColour colour;

			// List of ray hit information per sphere sample point
			cBitArray& rayhits;
			int		base_index;
		};
	};


	struct HDRLightSource : public SphericalFunction
	{
		HDRLightSource(const char* filename, cComponentManager* cm_ptr)
		{
			iImageLoader* image_loader = cImageLoaderFactory::Instance().LoadImage(cm_ptr, filename);

			// Allocate the image buffer
			image_size = image_loader->GetWidth();
			ASSERT(image_size == image_loader->GetHeight());
			image_ptr = new cColour[image_size * image_size];

			// Read the entire image
			for (int i = 0; i < image_size; i++)
				image_loader->ReadScanline((iImageLoader::Pixel*)image_ptr + image_size * i);

			cImageLoaderFactory::Instance().Destroy(image_loader);
		}

		~HDRLightSource(void)
		{
			delete [] image_ptr;
		}

		real ToReal(const SampleSphere& s) const
		{
			cColour c = ToColour(s);
			return ((c.r + c.g + c.b) / 3);
		}

		cColour ToColour(const SampleSphere& s) const
		{
			static const real one_over_pi = 1 / PI;
			real invl = 1 / sqrt(s.x * s.x + s.y * s.y);
			real r = one_over_pi * acos(s.z) * invl;
			real u = s.x * r;
			real v = -s.y * r;

			int x = int(u * image_size + image_size) >> 1;
			int y = int(v * image_size + image_size) >> 1;

			cColour col = image_ptr[y * image_size + x];
			col = col * 0.35f;

			return (col);
		}

		// Width & Height
		int		image_size;

		// Image colour array
		cColour*	image_ptr;
	};


	// Light with a monochromatic light source
	cColour SHLight(const RGBCoeffList& source_coeffs, real* light_coeffs)
	{
		cColour colour(0);

		for (int k = 0; k < SampleSphere::NB_BASES; k++)
		{
			const RGBCoeff& coeff = source_coeffs.e[k];

			colour.r += (float)(coeff.r * light_coeffs[k]);
			colour.g += (float)(coeff.g * light_coeffs[k]);
			colour.b += (float)(coeff.b * light_coeffs[k]);
		}

		return (colour);
	}


	// Light with an RGB light source
	cColour SHLight(const RGBCoeffList& source_coeffs, const RGBCoeffList& light_coeffs)
	{
		cColour colour(0);

		for (int k = 0; k < SampleSphere::NB_BASES; k++)
		{
			const RGBCoeff& scoeff = source_coeffs.e[k];
			const RGBCoeff& lcoeff = light_coeffs.e[k];

			colour.r += (float)(scoeff.r * lcoeff.r);
			colour.g += (float)(scoeff.g * lcoeff.g);
			colour.b += (float)(scoeff.b * lcoeff.b);
		}

		return (colour);
	}


	void BlancoSHRotate(const RGBCoeffList& source, const cMatrix& rotation, RGBCoeffList& dest)
	{
		SHRotMatrix shrm(SampleSphere::NB_BANDS);

		double rotm[16];
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				rotm[i * 4 + j] = rotation.e[i][j];

		for (int i = 0; i < 3; i++)
		{
			double shsrc[SampleSphere::NB_BASES], shdest[SampleSphere::NB_BASES];
			for (int j = 0; j < SampleSphere::NB_BASES; j++)
				shsrc[j] = source.e[j].e[i];

			shrm.rotate_sh(rotm, SampleSphere::NB_BANDS, shsrc, shdest);

			for (int j = 0; j < SampleSphere::NB_BASES; j++)
				dest.e[j].e[i] = shdest[j];
		}
	}


	void BlancoSHRotate(const real* source, const cMatrix& rotation, real* dest)
	{
		SHRotMatrix shrm(SampleSphere::NB_BANDS);

		double rotm[16];
		memset(rotm, 0, sizeof(rotm));
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				rotm[i * 4 + j] = rotation.e[i][j];

		double shsrc[SampleSphere::NB_BASES], shdest[SampleSphere::NB_BASES];
		for (int j = 0; j < SampleSphere::NB_BASES; j++)
			shsrc[j] = source[j];

		shrm.rotate_sh(rotm, SampleSphere::NB_BANDS, shsrc, shdest);

		for (int j = 0; j < SampleSphere::NB_BASES; j++)
			dest[j] = shdest[j];
	}
}


cSHTest::cSHTest(const int width, const int height) :

	cApplication(width, height, false),
	m_TrackBall(width, height),
	m_Font(0),
	m_ViewMode(0),
	m_Model(0)

{
	printf("Initialisation start: %d\n", clock());

	// Ensure needed components are active
	m_ComponentManager->LoadComponent("LWOLoader");
	m_ComponentManager->LoadComponent("HDRLoader");

	m_Samples = new Sample2D[MAX_NB_SAMPLES];
	m_SphereSamples = new SampleSphere[MAX_NB_SAMPLES];

	printf("Generating %d random samples on the sphere.\n", MAX_NB_SAMPLES);
	//GenerateSamples2D(m_Samples, MAX_NB_SAMPLES);
	GenerateSamples2D_Stratified(m_Samples, SQRT_NB_SAMPLES);
	ProjectSpherical(m_SphereSamples, m_Samples, MAX_NB_SAMPLES);

	// Generate a list of spherical harmonic visualisations
	printf("Generating 3D models of %d Spherical Harmonic bases.\n", SampleSphere::NB_BASES);
	m_SHLists = glGenLists(SampleSphere::NB_BASES);
	Enable3D();
	for (int l = 0, i = 0; l < SampleSphere::NB_BANDS; l++)
	{
		for (int m = -l; m <= l; m++, i++)
		{
			SHSample sample(l, m);
			DrawSphericalFunction(sample, m_SHLists + i);
		}
	}

	// Project the example lighting function onto the SH bases
	printf("Projecting example light onto SH bases.\n");
	ExampleLight el;
	ProjectSphericalFunction(el, m_SphereSamples, m_ExampleCoeffs);

	// Create a 3D model of the lighting function
	printf("Generating a 3D model of the example lighting function.\n");
	m_ExampleLists = glGenLists(2);
	DrawSphericalFunction(el, m_ExampleLists);

	// Reconstruct an approximation to the example lighting function from its SH co-efficients
	printf("Reconstructing example lighting model from SH co-efficients.\n");
	ReconstructSH rsh(m_ExampleCoeffs);
	DrawSphericalFunction(rsh, m_ExampleLists + 1);

	// Project the HDR image onto SH bases
	static const char* HDR_NAME = "rnl_probe.hdr";
	printf("Projecting HDR image %s onto SH bases.\n", HDR_NAME);
	HDRLightSource hdr(HDR_NAME, m_ComponentManager);
	RGBProjectSphericalFunction(hdr, m_SphereSamples, m_HDRCoeffs);

	//printf("\n");
	for (int i = 0; i < SampleSphere::NB_BASES; i++)
		printf("%.5f\n", m_HDRCoeffs.e[i].r);

	// Load the mesh
	static const char* MODEL_NAME = "beethoven.lwo";
	printf("Loading model %s.\n", MODEL_NAME);
	iMeshLoader* mesh_loader = cMeshLoaderFactory::Instance().LoadMesh(m_ComponentManager, MODEL_NAME);
	m_Model = new cPolyMesh(mesh_loader);

	// Get the polygon colours
	m_PolygonColours = new cColour[mesh_loader->GetNbPolygons()];
	for (int i = 0; i < mesh_loader->GetNbPolygons(); i++)
		mesh_loader->GetMaterialDiffuse(mesh_loader->GetPolygonMaterial(i), m_PolygonColours[i]);

	// Colour the vertices using the polygon material colours
	// This relies upon a vertex not sharing any material properties so it's pretty crap
	m_VertexColours = new cColour[mesh_loader->GetNbVertices()];
	for (int i = 0; i < mesh_loader->GetNbPolygons(); i++)
	{
		for (int j = 0; j < mesh_loader->GetPolygonNbIndices(i); j++)
			m_VertexColours[mesh_loader->GetPolygonVertexIndex(i, j)] = m_PolygonColours[i];
	}

	cMeshLoaderFactory::Instance().Destroy(mesh_loader);

	// Need a closed mesh to continue
	if (m_Model->IsClosed())
		printf("   model is closed.\n");
	else
		printf("   model is NOT closed.\n");

	// Allocate lists of SH co-efficients
	const std::vector<cPolyMesh::Vertex*>& vertices = m_Model->GetVertices();
	const std::vector<cPolyMesh::Face*>& faces = m_Model->GetFaces();
	printf("   contains %d vertices.\n", vertices.size());
	m_CoeffsDU = new RGBCoeffList[vertices.size()];
	m_CoeffsDS = new RGBCoeffList[vertices.size()];
	m_CoeffsDI = new RGBCoeffList[vertices.size()];

	if (FILE* fp = fopen("diffuse_unshadowed.shc", "rb"))
	{
		// Read in cached co-efficients
		printf("Reading diffuse unshadowed transfer co-efficients from file.\n");
		fread(m_CoeffsDU, 1, vertices.size() * sizeof(RGBCoeffList), fp);
		fclose(fp);
	}

	else
	{
		// First generate them
		printf("Generating co-efficients for diffuse unshadowed transfer...\n");
		for (size_t i = 0; i < vertices.size(); i++)
		{
			printf("\r   %.2f%%", ((float)i / (float)(vertices.size() - 1)) * 100);
			DiffuseUnshadowedTransfer dut(*vertices[i], m_VertexColours[i] * (1.0f / float(PI)));
			RGBProjectSphericalFunction(dut, m_SphereSamples, m_CoeffsDU[i]);
		}

		// Cache to file
		FILE* fp = fopen("diffuse_unshadowed.shc", "wb");
		fwrite(m_CoeffsDU, 1, vertices.size() * sizeof(RGBCoeffList), fp);
		fclose(fp);

		printf("\n");
	}

	if (FILE* fp = fopen("diffuse_shadowed.shc", "rb"))
	{
		// Read in cached co-efficients
		printf("Reading diffuse shadowed transfer co-efficients from file.\n");
		fread(m_CoeffsDS, 1, vertices.size() * sizeof(RGBCoeffList), fp);
		fclose(fp);
	}

	else
	{
		// Build the mesh vertex list
		std::vector<cVector> rt_vertices;
		rt_vertices.resize(vertices.size());
		for (size_t i = 0; i < vertices.size(); i++)
			rt_vertices[i] = vertices[i]->position;

		// Build the mesh triangle list
		std::vector<cVoxelRaytracer::Triangle> rt_triangles;
		rt_triangles.resize(faces.size());
		for (size_t i = 0; i < faces.size(); i++)
		{
			cVoxelRaytracer::Triangle& tri = rt_triangles[i];
			cPolyMesh::Face& face = *faces[i];

			tri.normal = face.normal;
			
			// Copy triangle vertex indices
			int c = 0;
			for (cPolyMesh::VertexList::iterator i = face.vertices.begin(); i != face.vertices.end(); ++i)
			{
				ASSERT(c < 3);
				tri.v[c++] = (*i)->index;
			}
		}

		// Need a raytracer for hit tests
		cVoxelRaytracer voxel_rt(rt_vertices, rt_triangles, 80);

		// Now generate diffuse shadowed co-efficients
		printf("Generating co-efficients for diffuse shadowed transfer...\n");
		std::vector<int> ignore_list;
		for (size_t i = 0; i < vertices.size(); i++)
		{
			cPolyMesh::Vertex& v = *vertices[i];

			// Ignore all faces that share this vertex
			ignore_list.clear();
			for (cPolyMesh::FaceList::iterator j = v.faces.begin(); j != v.faces.end(); ++j)
				ignore_list.push_back((*j)->index);
			voxel_rt.SetIgnoreList(ignore_list);

			printf("\r   %.2f%%", ((float)i / (float)(vertices.size() - 1)) * 100);
			DiffuseShadowedTransfer dst(*vertices[i], voxel_rt, m_VertexColours[i] * (1.0f / float(PI)));
			RGBProjectSphericalFunction(dst, m_SphereSamples, m_CoeffsDS[i]);
		}

		// Cache to file
		FILE* fp = fopen("diffuse_shadowed.shc", "wb");
		fwrite(m_CoeffsDS, 1, vertices.size() * sizeof(RGBCoeffList), fp);
		fclose(fp);

		printf("\n");
	}

	if (FILE* fp = fopen("diffuse_interreflected.shc", "rb"))
	{
		// Read in cached co-efficients
		printf("Reading diffuse interreflected transfer co-efficients from file.\n");
		fread(m_CoeffsDI, 1, vertices.size() * sizeof(RGBCoeffList), fp);
		fclose(fp);
	}

	else
	{
		// Build the mesh vertex list
		std::vector<cVector> rt_vertices;
		rt_vertices.resize(vertices.size());
		for (size_t i = 0; i < vertices.size(); i++)
			rt_vertices[i] = vertices[i]->position;

		// Build the mesh triangle list
		std::vector<cVoxelRaytracer::Triangle> rt_triangles;
		rt_triangles.resize(faces.size());
		for (size_t i = 0; i < faces.size(); i++)
		{
			cVoxelRaytracer::Triangle& tri = rt_triangles[i];
			cPolyMesh::Face& face = *faces[i];

			tri.normal = face.normal;
			
			// Copy triangle vertex indices
			int c = 0;
			for (cPolyMesh::VertexList::iterator i = face.vertices.begin(); i != face.vertices.end(); ++i)
			{
				ASSERT(c < 3);
				tri.v[c++] = (*i)->index;
			}
		}

		// Need a raytracer for hit tests
		cVoxelRaytracer voxel_rt(rt_vertices, rt_triangles, 80);

		// A bit for each sample of each vertex telling whether that sample is involved in self-intersection or not
		cBitArray rayhits(MAX_NB_SAMPLES * vertices.size());

		// Now generate diffuse interreflected co-efficients
		printf("Generating co-efficients for diffuse interreflected transfer...\n");
		std::vector<int> ignore_list;
		for (size_t i = 0; i < vertices.size(); i++)
		{
			cPolyMesh::Vertex& v = *vertices[i];

			// Ignore all faces that share this vertex
			ignore_list.clear();
			for (cPolyMesh::FaceList::iterator j = v.faces.begin(); j != v.faces.end(); ++j)
				ignore_list.push_back((*j)->index);
			voxel_rt.SetIgnoreList(ignore_list);

			printf("\r   Direct Illumination:              %.2f%%", ((float)i / (float)(vertices.size() - 1)) * 100);
			DiffuseInterreflectedTransfer::DirectPass dit_dp(*vertices[i], voxel_rt, m_VertexColours[i] * (1.0f / float(PI)), rayhits, i * MAX_NB_SAMPLES);
			RGBProjectSphericalFunction(dit_dp, m_SphereSamples, m_CoeffsDI[i]);
		}

		printf("\n");

		// Setup the coefficient lists
		// Vertices first of all take incoming light from the diffuse shadowed co-efficients
		// Then each subsequent pass takes incoming REFLECTED light only
		RGBCoeffList* coeff_db[NB_INDIRECT_PASSES + 1];
		coeff_db[0] = m_CoeffsDI;
		for (int i = 1; i <= NB_INDIRECT_PASSES; i++)
		{
			coeff_db[i] = new RGBCoeffList[vertices.size()];
			memset(coeff_db[i], 0, vertices.size() * sizeof(RGBCoeffList));
		}

		for (int pass = 1; pass <= NB_INDIRECT_PASSES; pass++)
		{
			// Figure buffer pointers
			RGBCoeffList* dbuf = coeff_db[pass];
			RGBCoeffList* sbuf = coeff_db[pass - 1];

			for (size_t i = 0; i < vertices.size(); i++)
			{
				cColour material_albedo = m_VertexColours[i] * (1.0f / float(PI));

				cPolyMesh::Vertex& vtx = *vertices[i];

				RGBCoeffList& d = dbuf[i];

				// Ignore all faces that share this vertex
				ignore_list.clear();
				for (cPolyMesh::FaceList::iterator j = vtx.faces.begin(); j != vtx.faces.end(); ++j)
					ignore_list.push_back((*j)->index);
				voxel_rt.SetIgnoreList(ignore_list);

				printf("\r   Indirect Illumination (pass %d/%d): %.2f%%", pass, NB_INDIRECT_PASSES, ((float)i / (float)(vertices.size() - 1)) * 100);

				for (int j = 0; j < MAX_NB_SAMPLES; j++)
				{
					// Ignore cases where there's no inter-reflectance
					if (rayhits[i * MAX_NB_SAMPLES + j] == false)
						continue;

					SampleSphere& s = m_SphereSamples[j];

					// Calculate cosine term for this sample
					cVector dir(float(s.x), float(s.y), float(s.z));
					float H = dir * vtx.normal;
					ASSERT(H > 0);

					// Do a hitscan through the triangle list again finding the nearest intersection
					float t, u, v;
					int hit = voxel_rt.GetClosestHitTriangle(vtx.position, dir, false, t, u, v);
					ASSERT(hit != -1);
					const cVoxelRaytracer::Triangle& tri = voxel_rt.GetTriangle(hit);

					// Complete barycentric co-ordinates
					float w = 1 - (u + v);

					// Ignore erroneous
					if (fabs(u + v + w - 1) > 1e-3)
						continue;

					// Clamp coords (some come out as tiny negative)
					u = min(1.0f, max(0.0f, u));
					v = min(1.0f, max(0.0f, v));
					w = min(1.0f, max(0.0f, w));

					// Get SH co-efficients for each vertex in the source buffer
					RGBCoeffList& a = sbuf[tri.v[0]];
					RGBCoeffList& b = sbuf[tri.v[1]];
					RGBCoeffList& c = sbuf[tri.v[2]];

					// Sum SH co-efficients at the triangle intersection point
					// (as the graphics hardware would interpolate them)
					for (int k = 0; k < SampleSphere::NB_BASES; k++)
					{
						RGBCoeff c((a.e[k] * u + b.e[k] * v + c.e[k] * w) * H);
						c.r *= material_albedo.r * 0.5f;
						c.g *= material_albedo.g * 0.5f;
						c.b *= material_albedo.b * 0.5f;
						d.e[k] += c;
					}
				}
			}

			// Divide each coefficient by number of samples and multiply by weights
			for (size_t i = 0; i < vertices.size(); i++)
				for (int j = 0; j < SampleSphere::NB_BASES; j++)
					dbuf[i].e[j] = dbuf[i].e[j] * ((4 * PI) / MAX_NB_SAMPLES);
		}

		// Finally, sum all indirect light passes with the direct light pass
		for (int i = 1; i <= NB_INDIRECT_PASSES; i++)
		{
			for (size_t j = 0; j < vertices.size(); j++)
				for (int k = 0; k < SampleSphere::NB_BASES; k++)
					m_CoeffsDI[j].e[k] += coeff_db[i][j].e[k];
		}

		// Cache to file
		FILE* fp = fopen("diffuse_interreflected.shc", "wb");
		fwrite(m_CoeffsDI, 1, vertices.size() * sizeof(RGBCoeffList), fp);
		fclose(fp);

		printf("\n");
	}

	printf("Initialisation end: %d\n", clock());
}


cSHTest::~cSHTest(void)
{
	delete [] m_SphereSamples;
	delete [] m_Samples;
}


bool cSHTest::ProcessFrame(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// List of viewing modes - just add to this list to have it usable
	ViewFunc vfs[] =
	{
		&cSHTest::DrawSamples,
		&cSHTest::DrawLegendre,
		&cSHTest::DrawSHBases,
		&cSHTest::DrawExampleLight,
		&cSHTest::DrawModelOpenGL,
		&cSHTest::DrawModelDiffuseUnshadowed,
		&cSHTest::DrawModelDiffuseShadowed,
		&cSHTest::DrawModelDiffuseInterreflected
	};

	// Move onto the next mode on space bar
	if (m_KeysReleased[SDLK_SPACE])
		m_ViewMode = (m_ViewMode + 1) % (sizeof(vfs) / sizeof(vfs[0]));

	// Execute current view mode
	(this->*vfs[m_ViewMode])();

	return (true);
}


void cSHTest::BeforeSwitch(void)
{
	delete m_Font;
}


void cSHTest::AfterSwitch(void)
{
	m_Font = new cFont("Courier New", 0, 12);
}


void cSHTest::DrawSamples(void)
{
	const cMatrix& rot = m_TrackBall.UpdateMat(m_MouseX, m_MouseY, m_MouseButtons & SDL_BUTTON(SDL_BUTTON_LEFT));

	Enable3D();
	::DrawSamples(m_Samples, m_SphereSamples, rot);
}


void cSHTest::DrawLegendre(void)
{
	Enable2D();
	DrawGrid(0.5f, 0.5f, 0.5f);

	for (int i = 0, k = 0; i <= 3; i++)
		for (int j = 0; j <= i; j++, k++)
			::DrawLegendre(i, j, m_Font, k);
}


void cSHTest::DrawSHBases(void)
{
	const cMatrix& rot = m_TrackBall.UpdateMat(m_MouseX, m_MouseY, m_MouseButtons & SDL_BUTTON(SDL_BUTTON_LEFT));

	Enable3D();

	for (int l = 0, i = 0; l < SampleSphere::NB_BANDS; l++)
	{
		for (int m = -l; m <= l; m++, i++)
		{
			glMatrixMode(GL_MODELVIEW);
			glLoadMatrix(rot * cMatrix::Translate(cVector((float)m, 2.5f - l * 1.2f, 6)));
			glCallList(m_SHLists + i);
		}
	}
}


void cSHTest::DrawExampleLight(void)
{
	const cMatrix& rot = m_TrackBall.UpdateMat(m_MouseX, m_MouseY, m_MouseButtons & SDL_BUTTON(SDL_BUTTON_LEFT));

	Enable3D();
	glMatrixMode(GL_MODELVIEW);

	glLoadMatrix(rot * cMatrix::Translate(cVector(-1, 0, 3)));
	glCallList(m_ExampleLists);

	glLoadMatrix(rot * cMatrix::Translate(cVector(1, 0, 3)));
	glCallList(m_ExampleLists + 1);
}


void cSHTest::DrawModelOpenGL(void)
{
	const cMatrix& rot = m_TrackBall.UpdateMat(m_MouseX, m_MouseY, m_MouseButtons & SDL_BUTTON(SDL_BUTTON_LEFT));

	Enable3D();

	// Load the model-view matrix for the matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrix(rot * cMatrix::Translate(cVector(0, 0, 1)));

	// Retrieve the mesh elements
	const std::vector<sdla::cPolyMesh::Vertex*>& vertices = m_Model->GetVertices();
	const std::vector<sdla::cPolyMesh::Face*>& faces = m_Model->GetFaces();

	// For each face
	for (size_t i = 0; i < faces.size(); i++)
	{
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, m_PolygonColours[i]);

		// Draw a polygon as a triangle fan
		glBegin(GL_TRIANGLE_FAN);
		for (cPolyMesh::VertexList::iterator j = faces[i]->vertices.begin(); j != faces[i]->vertices.end(); ++j)
		{
			const sdla::cPolyMesh::Vertex& v = *(*j);
			glNormal3fv(v.normal);
			glVertex3fv(v.position);
		}
		glEnd();
	}
}


void cSHTest::DrawModel(RGBCoeffList* coeffs)
{
	const cMatrix& rot = m_TrackBall.UpdateMat(m_MouseX, m_MouseY, m_MouseButtons & SDL_BUTTON(SDL_BUTTON_LEFT));

	Enable3D();

	// Load the model-view matrix for the orientation
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrix(rot * cMatrix::Translate(cVector(0, 0, 1)));

	// Retrieve the mesh elements
	const std::vector<sdla::cPolyMesh::Vertex*>& vertices = m_Model->GetVertices();
	const std::vector<sdla::cPolyMesh::Face*>& faces = m_Model->GetFaces();

	glDisable(GL_LIGHTING);

	RGBCoeffList transformed;
	static float t = 0;
	cMatrix rotation(cMatrix::Identity());
	cSHRotateMatrix shrm(SampleSphere::NB_BANDS);
	SHRotateMatrix newone;
	SHRotate(newone, cMatrix::Inverse(rot));
	//for (int i = 0; i < shrm.DIM; i++)
	//	for (int j = 0; j < shrm.DIM; j++)
	//		shrm.e[i][j] = (i == j);
	newone.Transform(m_HDRCoeffs, transformed);
	//newone.Transform(m_ExampleCoeffs, transformed);

	static float exp = 1;

	if (m_Keys[SDLK_UP])
		exp += 0.01f;
	if (m_Keys[SDLK_DOWN])
		exp -= 0.01f;

	//BlancoSHRotate(m_HDRCoeffs, cMatrix::Inverse(rot.ToMatrix()), transformed);

	//printf("%f - %f\n", m_HDRCoeffs.e[3].r, transformed.e[3].r);

	// For each face
	for (size_t i = 0; i < faces.size(); i++)
	{
		// Draw a polygon as a triangle fan
		glBegin(GL_TRIANGLE_FAN);
		for (cPolyMesh::VertexList::iterator j = faces[i]->vertices.begin(); j != faces[i]->vertices.end(); ++j)
		{
			const sdla::cPolyMesh::Vertex& v = *(*j);

			glColor3fv(SHLight(coeffs[v.index], transformed));
			glVertex3fv(v.position);
		}
		glEnd();
	}

	glEnable(GL_LIGHTING);
}


void cSHTest::DrawModelDiffuseUnshadowed(void)
{
	// Could do with functors for this...
	DrawModel(m_CoeffsDU);
}


void cSHTest::DrawModelDiffuseShadowed(void)
{
	DrawModel(m_CoeffsDS);
}


void cSHTest::DrawModelDiffuseInterreflected(void)
{
	DrawModel(m_CoeffsDI);
}
