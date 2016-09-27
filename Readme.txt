
Spherical Harmonic Lighting
Don Williamson, June 2003
don@donw.co.uk
---------------------------

This is an implementation of Spherical Harmonic Lighting to partially simulate Global Illumination
on a per-vertex level in real-time. The structure of the program follows directly from the material
written in Robin Green's excellent GDC2003 presentation. It basically functions as an object
viewer with many modes. You can rotate all objects you can see by clicking the left mouse button on
the window and dragging in any direction. You can change modes by pressing the space bar. The modes
are:

	1. Monte-carlo estimator samples on the sphere.
	2. Associated Legendre plot.
	3. 3D representation of spherical harmonic bases.
	4. Comparison of an original spherical lighting function (left) and it's SH projected equivalent (right, 4 bands).
	5. Simple OpenGL-lit beethoven bust.
	6. Same model diffusely lit by HDR lightsource (rnl_probe.hdr).
	7. Diffusely lit with self-shadowing.
	8. Diffusely lit with self-shadowing and inter-reflectance.

This is a Win32 OpenGL program that any recent (last 3 years) video card should be able to run.

The simulation starts by generating SH co-efficients for each vertex in 3 types of lighting model:
diffuse unshadowed, diffuse shadowed, and diffuse inter-reflected. For each vertex, integration over
the sphere is performed using a Monte-Carlo integrator with 10,000 samples, and then projected onto
4 bands of spherical harmonics. For the first time the program starts you're going to have

SH rotation is performed using my verbose implementation of the method of Ivanic/Ruedenberg and is
unoptimised. I've also implemented the SH rotation method of Choi et. al. but it is not present in this
release. I am currently co-ordinating an SH rotation specific release with Robin Green for later.
I have also removed Robin's implementation of Blanco SH rotation as it's not mine so the code won't
compile as-is: you will just have to uncomment any errors that occur (everything is included for you
to be able to compile the demo itself).

I apologise for the messy state of the code but this is purely and simply a spare time test which I
hope you all can benefit from. In order to compile this you need to download my SDLApp framework which
is available on my website.


Acknowledgements
----------------

Thanks to Robin Green for his excellent GDC2003 presentation and accompany documentation, and for
help and discussion with regard to SH lighting. Thanks to Joseph Ivanic and Klaus Ruedenberg for
help with SH rotation.

The HDR image rnl_probe.hdr is taken directly from Paul Debevec's HDR gallery at http://www.debevec.org.

References
----------

Precomputed Radiance Transfer for Real-Time Rendering in Dynamic, Low-Frequency Lighting Environments
Peter-Pike Sloan, Jan Kautz, and John Snyder
SIGGRAPH 2002, July 2002
http://research.microsoft.com/~ppsloan/

Spherical Harmonic Lighting: The Gritty Details
Robin Green
GDC 2003, Presentation
http://www.research.scea.com/gdc2003/spherical-harmonic-lighting.html

Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion
Joseph Ivanic and Klaus Ruedenberg
J. Phys. Chem. 1996, 100, 6342-5347

Additions and Corrections (to previous paper)
Joseph Ivanic and Klaus Ruedenberg
J. Phys. Chem. A, 1998, Vol. 102, No. 45, 9099

Rapid and stable determination of rotation matrices between spherical harmonics by direct recursion
Cheol Ho Choi, Joseph Ivanic, Mark S. Gordon, and Klaus Ruedenberg
J. Chem. Phys, 1999, Vol. 111, No. 19
