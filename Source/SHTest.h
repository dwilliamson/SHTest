
#pragma once


#include <sdla/Application.h>
#include <sdla/TrackBall.h>
#include "Types.h"


namespace sdla
{
	class cFont;
	class cPolyMesh;
	struct cColour;
};


class cSHTest : public sdla::cApplication
{
public:
	// Constructor/destructor
	cSHTest(const int width, const int height);
	~cSHTest(void);

private:
	// cApplication implementations
	bool	ProcessFrame(void);
	void	BeforeSwitch(void);
	void	AfterSwitch(void);

	// List of viewing modes
	typedef void (cSHTest::*ViewFunc)(void);
	void	DrawSamples(void);
	void	DrawLegendre(void);
	void	DrawSHBases(void);
	void	DrawExampleLight(void);
	void	DrawModelOpenGL(void);
	void	DrawModel(RGBCoeffList* coeffs);
	void	DrawModelDiffuseUnshadowed(void);
	void	DrawModelDiffuseShadowed(void);
	void	DrawModelDiffuseInterreflected(void);

	// Randomly generated samples
	Sample2D*		m_Samples;
	SampleSphere*	m_SphereSamples;

	// Current view mode
	int		m_ViewMode;

	// Allow the user to rotate the objects
	sdla::cTrackBall	m_TrackBall;

	// Font to write text with
	sdla::cFont*	m_Font;

	// Polygonal generated SH basis functions
	int		m_SHLists;

	// Coefficient per base for example lighting function
	real	m_ExampleCoeffs[SampleSphere::NB_BASES];

	// Coefficient per base for the loaded HDR light source
	RGBCoeffList	m_HDRCoeffs;

	// 3D representations of the example lighting functions
	int		m_ExampleLists;

	// Model to light
	sdla::cPolyMesh*	m_Model;

	// Colours taken from materials
	sdla::cColour*	m_PolygonColours;
	sdla::cColour*	m_VertexColours;

	// Diffuse unshadowed SH co-efficients for each vertex in the model
	RGBCoeffList*	m_CoeffsDU;

	// Diffuse shadowed SH co-efficients for each vertex in the model
	RGBCoeffList*	m_CoeffsDS;

	// Diffuse interreflected SH co-efficients for each vertex in the model
	RGBCoeffList*	m_CoeffsDI;
};