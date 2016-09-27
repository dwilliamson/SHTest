
// The thing you need is the 3x3 real roation matrix, which is equation 5.1
// in our paper.(Note that there is a small thypo.) The way we do this is


      /*subroutine rotr(A,B,R)
//
//    This routine returns the 3x3 coordinate rotation matrix R which
//     rotate A to B, where A and B are row vectors.
//    B=A*R
//
//     C. H. Choi March 1999
//
      implicit double precision(a-h,o-z)
      parameter (zero=0.0d+00,one=1.0d+00,pi=3.141592653589793238d+00,Eps=1.0d-14)
      dimension A(3), B(3), R(3,3)
//     Calculate the angle between A and B.
      an=sqrt(A(1)*A(1)+A(2)*A(2)+A(3)*A(3))
      bn=sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))
      dot=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      Ang=acos(dot/(an*bn))
      if (ang.eq.zero) then
         do i=1,3
            do j=1,i
               R(i,j)=zero
               R(j,i)=zero
               if (i.eq.j) then
                  R(i,j)=one
               endif
            enddo
         enddo
      elseif (abs(abs(ang)-pi).le.Eps) then
         do i=1,3
            do j=1,i
               R(i,j)=zero
               R(j,i)=zero
               if (i.eq.j) then
                  R(i,j)=-one
               endif
            enddo
         enddo
      else
         Cw=cos(Ang)
         Sw=sin(Ang)
         Cwn=1.0d+00-Cw
//        Calculate the cross product between A and B.
         x=B(2)*A(3)-B(3)*A(2)
         y=B(3)*A(1)-B(1)*A(3)
         z=B(1)*A(2)-B(2)*A(1)
//        Normalize the cross product.
         cn=sqrt(x*x+y*y+z*z)
         x=x/cn
         y=y/cn
         z=z/cn
//        Calculate R.
         R(1,1)=Cwn*x*x+Cw
         R(2,1)=Cwn*x*y+z*Sw
         R(3,1)=Cwn*x*z-y*Sw
         R(1,2)=Cwn*x*y-z*Sw
         R(2,2)=Cwn*y*y+Cw
         R(3,2)=Cwn*y*z+x*Sw
         R(1,3)=Cwn*x*z+y*Sw
         R(2,3)=Cwn*y*z-x*Sw
         R(3,3)=Cwn*z*z+Cw
      endif

      return
      end*/

// Ok, now once you have R, you can construct D, the first order Wigner D
// matrix. (eq. 5.4 and 5.5)
// We do this by

      subroutine rotc(F,G,R)
//
//     This routine returns the Wigner rotation matrix (D) of order 1
//     on the basis of 3x3 orthogonal coordinate transformation matrix
//     R.
//     F and G are real and imaginary parts of it.
//     Consequently, D = F + i*G.
//    And F and G are determined entirely by R.
//
//     Ref : Choi, C.H.; Ivanic, J.; Gordon, M.S.; Ruedenberg, K,
//           J. Chem. Phys.
//
//     C. H. Choi March 1999
//
      implicit double precision(a-h,o-z)
      parameter (fs=0.70710678118655d+00,f2=0.5d+00,zero=0.0d+00)
      dimension R(3,3),F(-1:1,-1:1),G(-1:1,-1:1)

      F(-1,-1)=f2*(R(2,2)+R(1,1))
      F( 0,-1)= fs*R(3,1)
      F( 1,-1)=f2*(R(2,2)-R(1,1))
      F(-1, 0)= fs*R(1,3)
      F( 0, 0)=    R(3,3)
      F( 1, 0)=-fs*R(1,3)
      F(-1, 1)=f2*(R(2,2)-R(1,1))
      F( 0, 1)=-fs*R(3,1)
      F( 1, 1)=f2*(R(2,2)+R(1,1))
      G(-1,-1)=f2*(R(2,1)-R(1,2))
      G( 0,-1)=-fs*R(3,2)
      G( 1,-1)=f2*(R(2,1)+R(1,2))
      G(-1, 0)= fs*R(2,3)
      G( 0, 0)= zero
      G( 1, 0)= fs*R(2,3)
      G(-1, 1)=-f2*(R(2,1)+R(1,2))
      G( 0, 1)=-fs*R(3,2)
      G( 1, 1)=f2*(R(1,2)-R(2,1))

      return
      end

// So far so good, now the real thing!
// We use the ROTL to calculate D of order L on the basis of D of order L-1
// and first order D.

// Note that before we do this, we initially calculate and save
// coefficients CLM by

      subroutine getclm(CLM)

//     This routine gets the CLM data.
//
//     C. H. Choi April 1999
//
      implicit double precision (a-h,o-z)
      logical QOPS,qfmm
      common /qmfm  /Size,Eps,dpgd,qfmm,Np,Ns,iws,NpGP,MPMTHD,NUMRD,ITERMS,QOPS,ISCUT
      dimension CLM(-NP:NP)
      do LLM=-Np,Np
         CLM(LLM)=sqrt(dfloat(Np+LLM))
      enddo

      return
      end

// You don't have to worry about the common block. All you need to know is
// parameter Np, which is the highest
// order of D you want you calculate. CLM does not depend on the particular
// rotation. In oder words, CLM is composed of contant values. So You don't
// have to calculate this everytime you use it. You can define them as
// constants.

// We also use function eo. I think you can find better way to do this.

      function eo(ix)

//     This function returns either 1 or -1.
//     eo = 1, if ix is even or zero
//     eo =-1, if ix is odd.
//
//     C. H. Choi April 1999

      implicit double precision (a-h,o-z)
      parameter (one=1.0d+00,zero=0.0d+00)
      eo=-mod(abs(ix),2)
      if (eo.eq.zero) then
         eo=one
      endif

      return
      end



// Once you have first order D, CLM, then you are ready to calculate higher
// order D by



      subroutine rotl(F,G,Fo,Go,Fn,Gn,CLM,L)

//     This routine calculates the rotation matrix (D=Fn+i*Gn)
//     of order L which rotates complex spherical harmonics of order L
//     on the basis of order 1 (F,G) and order L-1 (Fo,Go).
//     There is exact correspondence between the rows(i) and columns(j)
//     of the matrices, and m and m' indices.
//
//        Fn/Gn(m,m') : (m,m') element of the Fn/Gn
//
//     C.H. Choi, March 1999
//
//     Ref : Choi, C.H.; Ivanic, J.; Gordon, M.S.; Ruedenberg, K.
//           J. Chem. Phys. 1999, 111, 8825.
//
      implicit double precision(a-h,o-z)
      parameter (sq=1.414213562373095d+00)
      dimension F(-1:1,-1:1),G(-1:1,-1:1),Fo(-L+1:L-1,-L+1:L-1),Go(-L+1:L-1,-L+1:L-1),Fn(-L:L,-L:L),Gn(-L:L,-L:L),CLM(-L:L)
//     Let's define frequent indices.
      Ip  =  L-1
      In  = -L+1
      I2p =  L-2
      I2n = -L+2
//     ---------------------------------------------
//            The cases of abs(m') = L
//     ---------------------------------------------
//
//     Case 1 : m=L,-L    m'=L,-L
      d=CLM(L)*CLM(Ip)/(CLM(L)*CLM(Ip))
      Fn( L, L)=d*(F( 1,1)*Fo(Ip,Ip)-G( 1,1)*Go(Ip,Ip))
      Gn( L, L)=d*(F( 1,1)*Go(Ip,Ip)+G( 1,1)*Fo(Ip,Ip))
      Fn(-L, L)=d*(F(-1,1)*Fo(In,Ip)-G(-1,1)*Go(In,Ip))
      Gn(-L, L)=d*(F(-1,1)*Go(In,Ip)+G(-1,1)*Fo(In,Ip))
//     Due to the parity relations.
      Fn( L,-L)= Fn(-L,L)
      Gn( L,-L)=-Gn(-L,L)
      Fn(-L,-L)= Fn( L,L)
      Gn(-L,-L)=-Gn( L,L)
//
//     Case 2 : m= L-1,  m'=L,-L
//     Case 3 : m=-L+1,  m'=L,-L
//
      d=CLM(Ip)*CLM(I2p)/(CLM(L)*CLM(Ip))
      c=sq*CLM(Ip)*CLM(In)/(CLM(L)*CLM(Ip))
      Fn(Ip, L)=d*(F(1, 1)*Fo(I2p,Ip)-G(1, 1)*Go(I2p,Ip))+c*(F(0, 1)*Fo(Ip, Ip)-G(0, 1)*Go(Ip, Ip))
      Gn(Ip, L)=d*(F(1, 1)*Go(I2p,Ip)+G(1, 1)*Fo(I2p,Ip))+c*(F(0, 1)*Go(Ip, Ip)+G(0, 1)*Fo(Ip, Ip))
      Fn(Ip,-L)=d*(F(1,-1)*Fo(I2p,In)-G(1,-1)*Go(I2p,In))+c*(F(0,-1)*Fo(Ip, In)-G(0,-1)*Go(Ip, In))
      Gn(Ip,-L)=d*(F(1,-1)*Go(I2p,In)+G(1,-1)*Fo(I2p,In))+c*(F(0,-1)*Go(Ip, In)+G(0,-1)*Fo(Ip, In))
      Fn(In,-L)= -Fn(Ip, L)
      Gn(In,-L)=  Gn(Ip, L)
      Fn(In, L)= -Fn(Ip,-L)
      Gn(In, L)=  Gn(Ip,-L)
//     Case 4 : abs(m) < L-1, m'=L, m'=-L
      sig=-eo(L)
      do LMm=0,I2p
         LMmp=LMm+1
         LMmm=LMm-1
         d1=CLM(-LMm)*CLM(-LMm-1)/(CLM(L)*CLM(Ip))
         d2=CLM(LMm)*CLM(LMmm)/(CLM(L)*CLM(Ip))
         c=sq*CLM(LMm)*CLM(-LMm)/(CLM(L)*CLM(Ip))
         Fn(LMm, L)=d1*(F(-1, 1)*Fo(LMmp,Ip)-G(-1, 1)*Go(LMmp,Ip))+d2*(F( 1, 1)*Fo(LMmm,Ip)-G( 1, 1)*Go(LMmm,Ip))+c* (F( 0, 1)*Fo(LMm, Ip)-G( 0, 1)*Go(LMm, Ip))
         Gn(LMm, L)=d1*(F(-1, 1)*Go(LMmp,Ip)+G(-1, 1)*Fo(LMmp,Ip))+d2*(F( 1, 1)*Go(LMmm,Ip)+G( 1, 1)*Fo(LMmm,Ip))+c* (F( 0, 1)*Go(LMm, Ip)+G( 0, 1)*Fo(LMm, Ip))
         Fn(LMm,-L)=d1*(F(-1,-1)*Fo(LMmp,In)-G(-1,-1)*Go(LMmp,In))+d2*(F( 1,-1)*Fo(LMmm,In)-G( 1,-1)*Go(LMmm,In))+c* (F( 0,-1)*Fo(LMm, In)-G( 0,-1)*Go(LMm, In))
         Gn(LMm,-L)=d1*(F(-1,-1)*Go(LMmp,In)+G(-1,-1)*Fo(LMmp,In))+d2*(F( 1,-1)*Go(LMmm,In)+G( 1,-1)*Fo(LMmm,In))+c* (F( 0,-1)*Go(LMm, In)+G( 0,-1)*Fo(LMm, In))
         sig=-sig
         Fn(-LMm,-L)= sig*Fn(LMm, L)
         Gn(-LMm,-L)=-sig*Gn(LMm, L)
         Fn(-LMm, L)= sig*Fn(LMm,-L)
         Gn(-LMm, L)=-sig*Gn(LMm,-L)
      enddo
//     ---------------------------------------------
//            Now the cases of abs(m') < L
//     ---------------------------------------------
//
//     Case 5 : m=L,-L    abs(m') < L
      sig=-eo(L)
      do LMp=0,Ip
        b=CLM(L)*CLM(Ip)/(sq*CLM(LMp)*CLM(-LMp))
        Fn( L,LMp)=b*(F( 1,0)*Fo(Ip,LMp)-G( 1,0)*Go(Ip,LMp))
        Gn( L,LMp)=b*(F( 1,0)*Go(Ip,LMp)+G( 1,0)*Fo(Ip,LMp))
        Fn(-L,LMp)=b*(F(-1,0)*Fo(In,LMp)-G(-1,0)*Go(In,LMp))
        Gn(-L,LMp)=b*(F(-1,0)*Go(In,LMp)+G(-1,0)*Fo(In,LMp))
        sig=-sig
        Fn(-L,-LMp)= sig*Fn( L,LMp)
        Gn(-L,-LMp)=-sig*Gn( L,LMp)
        Fn( L,-LMp)= sig*Fn(-L,LMp)
        Gn( L,-LMp)=-sig*Gn(-L,LMp)
      enddo
//
//     Case 6 : m=L-1,-L+1    abs(m') < L
//
      sig=-eo(Ip)
      do LMp=0,Ip
        b=CLM(Ip)*CLM(I2p)/(sq*CLM(LMp)*CLM(-LMp))
        a=CLM(Ip)*CLM(In)/(CLM(LMp)*CLM(-LMp))
        Fn(Ip,LMp)=b*(F( 1,0)*Fo(I2p,LMp)-G( 1,0)*Go(I2p,LMp))+a*F( 0,0)*Fo(Ip, LMp)
        Gn(Ip,LMp)=b*(F( 1,0)*Go(I2p,LMp)+G( 1,0)*Fo(I2p,LMp))+a*F( 0,0)*Go(Ip, LMp)
        Fn(In,LMp)=b*(F(-1,0)*Fo(I2n,LMp)-G(-1,0)*Go(I2n,LMp))+a*F( 0,0)*Fo(In, LMp)
        Gn(In,LMp)=b*(F(-1,0)*Go(I2n,LMp)+G(-1,0)*Fo(I2n,LMp))+a*F( 0,0)*Go(In, LMp)
        sig=-sig
        Fn(In,-LMp) = sig*Fn(Ip,LMp)
        Gn(In,-LMp) =-sig*Gn(Ip,LMp)
        Fn(Ip,-LMp) = sig*Fn(In,LMp)
        Gn(Ip,-LMp) =-sig*Gn(In,LMp)
      enddo
//
//     Case 7 : abs(m) < L-1,   abs(m') < L
//
      do LMp=In,Ip
         sig=-eo(LMp)
         do LMm=0,I2p
            LMmm=LMm-1
            LMmp=LMm+1
            b1=CLM(-LMm)*CLM(-LMm-1)/(sq*CLM(LMp)*CLM(-LMp))
            b2=CLM(LMm)*CLM(LMmm)/(sq*CLM(LMp)*CLM(-LMp))
            a =CLM(LMm)*CLM(-LMm)/(CLM(LMp)*CLM(-LMp))
            Fn(LMm,LMp)=b1*(F(-1,0)*Fo(LMmp,LMp)-G(-1,0)*Go(LMmp,LMp))+b2*(F( 1,0)*Fo(LMmm,LMp)-G( 1,0)*Go(LMmm,LMp))+a*  F( 0,0)*Fo(LMm, LMp)
            Gn(LMm,LMp)=b1*(F(-1,0)*Go(LMmp,LMp)+G(-1,0)*Fo(LMmp,LMp))+b2*(F( 1,0)*Go(LMmm,LMp)+G( 1,0)*Fo(LMmm,LMp))+a*  F( 0,0)*Go(LMm, LMp)
            sig=-sig
            Fn(-LMm,-LMp)= sig*Fn(LMm,LMp)
            Gn(-LMm,-LMp)=-sig*Gn(LMm,LMp)
         enddo
      enddo

      return
      end
