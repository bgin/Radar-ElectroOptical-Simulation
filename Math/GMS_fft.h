

#ifndef __GMS_FFT_H__
#define __GMS_FFT_H__

/********************************************************************
 *                                                                  *
 * THIS FILE IS PART OF THE OggSQUISH SOFTWARE CODEC SOURCE CODE.   *
 *                                                                  *
 ********************************************************************



> From xiphmont@MIT.EDU Tue Jun  3 14:50:15 1997
> Subject: Re: fft.c 
> Date: Tue, 03 Jun 1997 14:50:11 EDT
> From: Monty <xiphmont@MIT.EDU>
> 
> The include files are not really necessary; they define a few macros
> depending on the compiler in use (eg, STIN becomes one of "static" or
> "static inline") and the prototypes for the functions in fft.c.  You
> should be able to compile and use the code without the headers with
> minimal modification.  The C code behaves nearly identically to the
> fortran code; inspection of the code should clarify the differences
> (eg, the 'ifac' argument).  Let me know if you need any other help
> with the code.
> 
> Monty




  file: fft.c
  function: Fast discrete Fourier and cosine transforms and inverses
  author: Monty <xiphmont@mit.edu>
  modifications by: Monty
  last modification date: Jul 1 1996

 ********************************************************************/

/* These Fourier routines were originally based on the Fourier
   routines of the same names from the NETLIB bihar and fftpack
   fortran libraries developed by Paul N. Swarztrauber at the National
   Center for Atmospheric Research in Boulder, CO USA.  They have been
   reimplemented in C and optimized in a few ways for OggSquish. */

/* As the original fortran libraries are public domain, the C Fourier
   routines in this file are hereby released to the public domain as
   well.  The C routines here produce output exactly equivalent to the
   original fortran routines.  Of particular interest are the facts
   that (like the original fortran), these routines can work on
   arbitrary length vectors that need not be powers of two in
   length. */

/* Notes from RFB: 

Looks like the user-level routines are:

Real FFT

void __ogg_fdrffti(int n, float *wsave, int *ifac)
void __ogg_fdrfftf(int n, float *r, float *wsave, int *ifac)
void __ogg_fdrfftb(int n, float *r, float *wsave, int *ifac)

__ogg_fdrffti == initialization
__ogg_fdrfftf == forward transform
__ogg_fdrfftb == backward transform

Parameters are
n == length of sequence
r == sequence to be transformed (input)
== transformed sequence (output)
wsave == work array of length 2n (allocated by caller)
ifac == work array of length 15 (allocated by caller)

Cosine quarter-wave FFT

void __ogg_fdcosqi(int n, float *wsave, int *ifac)
void __ogg_fdcosqf(int n, float *x, float *wsave, int *ifac)
void __ogg_fdcosqb(int n, float *x, float *wsave, int *ifac)
*/


void __ogg_fdrffti(int n, float *wsave, int *ifac)
                                 __attribute__((aligned(32)))
								         __attribute__((hot));
void __ogg_fdrfftf(int n, float *r, float *wsave, int *ifac)
                                 __attribute__((aligned(32)))
								 __attribute__((hot));
void __ogg_fdrfftb(int n, float *r, float *wsave, int *ifac)
                                 __attribute__((aligned(32)))
								 __attribute__((hot));

void __ogg_fdcosqi(int n, float *wsave, int *ifac)
                                  __attribute__((aligned(32)))
								  __attribute__((hot));
void __ogg_fdcosqf(int n, float *x, float *wsave, int *ifac)
                                   __attribute__((aligned(32)))
								   __attribute__((hot));
void __ogg_fdcosqb(int n, float *x, float *wsave, int *ifac)
                                   __attribute__((aligned(32)))
								   __attribute__((hot));











#endif /*__GMS_FFT_H__*/