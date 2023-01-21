/*
 * fft.c: Code for generic fft routines.
 * [Part of simple-fft-1.5.tar.Z (ver 1.5)]
 *
 * Copyright (C) 1991 Steve Haehnichen.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 1, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 * 
 * 
 *  fft.c was derived from US Navy radar code written by Steve Sampson
 *  in June, 1988 for the MS-DOS platform, and placed in the Public Domain.
 *  The fft() function was copied directly from fft.c in his FFT26.ZIP.
 *    Steve Sampson
 *    Box 45668
 *    Tinker AFB, OK, 73145
 *
 *  Steve Haehnichen  (shaehnic@ucsd.edu)      March, 1992.
 *
 *  $Id: fft.c,v 1.8 1992/04/02 06:56:34 steve Exp steve $
 */

#include "fft.h"


# define		SINE(x)		SineTable[(x)]
# define		COSINE(x)       CosineTable[(x)]
# define		WINDOW(x)	WinTable[(x)]
static FLOAT	 	*SineTable, *CosineTable, *WinTable;

/*
 *  This PERMUTE could be defined as a table lookup (Sampson did this),
 *  but I found that to be significantly slower than just computing it
 *  inline.  Of course, this is w/ GCC, and not MS-DOS Turbo C.
 *  I'll leave it as a macro; you can tweak with it.
 */
#define PERMUTE(x, y)	reverse((x), (y))

/* Number of samples in one "frame" */
#define SAMPLES 	(1 << bits)
#define REAL(x)		wave[(x)].re
#define IMAG(x)		wave[(x)].im
#define ALPHA		0.54	/* ala hanning */
  
/*
 *  Bit reverser for unsigned ints
 *  Reverses 'bits' bits.
 */
inline const unsigned int
reverse (unsigned int val, int bits)
{
  unsigned int retn = 0;
  
  while (bits--)
    {
      retn <<= 1;
      retn |= (val & 1);
      val >>= 1;
    }
  return (retn);
}

/*
 *  Here is the real work-horse.
 *  It's a generic FFT, so nothing is lost or approximated.
 *  The samples in wave[] should be in order, and they
 *  will be decimated when fft() returns.
 */
void fft (complx wave[], int bits)
{
  register int	loop, loop1, loop2;
  unsigned	i1;		/* going to right shift this */
  int		i2, i3, i4, y;
  FLOAT		a1, a2, b1, b2, z1, z2;

  i1 = SAMPLES / 2;
  i2 = 1;

  /* perform the butterfly's */

  for (loop = 0; loop < bits; loop++)
    {
      i3 = 0;
      i4 = i1;

      for (loop1 = 0; loop1 < i2; loop1++)
	{
	  y  = PERMUTE(i3 / (int)i1, bits);
	  z1 = COSINE(y);
	  /* z2 = -SINE(y);             VGZ: fix the sign */
	  z2 = SINE(y);                

	  for (loop2 = i3; loop2 < i4; loop2++)
	    {
	      a1 = REAL(loop2);
	      a2 = IMAG(loop2);

	      b1 = z1 * REAL(loop2+i1) - z2 * IMAG(loop2+i1);
	      b2 = z2 * REAL(loop2+i1) + z1 * IMAG(loop2+i1);

	      REAL(loop2) = a1 + b1;
	      IMAG(loop2) = a2 + b2;

	      REAL(loop2+i1) = a1 - b1;
	      IMAG(loop2+i1) = a2 - b2;
	    }

	  i3 += (i1 << 1);
	  i4 += (i1 << 1);
	}

      i1 >>= 1;
      i2 <<= 1;
    }
}
  
/*
 *  Put the samples back in order after the FFT scrambles them.
 *  (Decimation-in-frequecy)
 *  Untested and unused because I haven't done IFFT yet.
 */
void
reorder (complx wave[], int bits)
{
  int i, j;
  complx temp;
  
  for (i = 0; i < SAMPLES; i++)
    {
      j = PERMUTE (i, bits);
      if (i > j)		/* So we only do each pair once */
	{
	  temp = wave[i];
	  wave[i] = wave[j];
	  wave[j] = temp;
	}
    }
}

/*
 *  Scale sampled values.
 *  Do this *before* the fft.
 */
void scale (complx wave[], int bits)
{
  int i;

  for (i = 0; i < SAMPLES; i++)  {
    wave[i].re /= (SAMPLES/2);       /* VGZ: add devide by 2 */
    wave[i].im /= (SAMPLES/2);
  }
  //wave[0].re /=2;                    /* fix 0th element */
}


/*
 *  Calculate amplitude of component n in the decimated wave[] array.
 */
FLOAT amp (int n, complx wave[], int bits)
{
  n = PERMUTE (n, bits);
  return (hypot (REAL(n), IMAG(n)));
}


/*
 *  Calculate phase of component n in the decimated wave[] array.
 */
FLOAT phase (int n, complx wave[], int bits)
{
  n = PERMUTE (n, bits);
  if (REAL(n) != 0.0)
    return (atan (IMAG(n) / REAL(n)));
  else
    return (0.0);
}


/*
 *  Initializer for FFT routines.  Currently only sets up tables.
 *  - Generates scaled lookup tables for sin() and cos()
 *  - Fills a table for the Hamming window function
 */
void fft_init (int bits)
{
  int i;
  const FLOAT  	TWOPIoN   = (atan(1.0) * 8.0) / (FLOAT)SAMPLES;
  const FLOAT 	TWOPIoNm1 = (atan(1.0) * 8.0) / (FLOAT)(SAMPLES - 1);

  SineTable   = (FLOAT *) malloc (sizeof(FLOAT) * SAMPLES);
  CosineTable = (FLOAT *) malloc (sizeof(FLOAT) * SAMPLES);
  WinTable    = (FLOAT *) malloc (sizeof(FLOAT) * SAMPLES);
  for (i=0; i < SAMPLES; i++)
    {
      SineTable[i]   = sin((FLOAT) i * TWOPIoN);
      CosineTable[i] = cos((FLOAT) i * TWOPIoN);
      /*
       * Generalized Hamming window function.
       * Set ALPHA to 0.54 for a hanning window. (Good idea)
       */
      WinTable[i] = ALPHA + ((1.0 - ALPHA)
				* cos (TWOPIoNm1 * (i - SAMPLES/2))); 
    }
}


/*
 *  Apply some windowing function to the samples.
 */
void fft_window (complx wave[], int bits)
{
  int i;

  for (i = 0; i < SAMPLES; i++)
    {
      REAL(i) *= WINDOW(i);
      IMAG(i) *= WINDOW(i);
    }
}

/*
 *  Undo the window function to restore full magnitude.
 *  Only works if there are NO zeros in WinTable!  (like hanning's)
 */
void fft_unwindow (complx wave[], int bits)
{
  int i;

  for (i = 0; i < SAMPLES; i++)
    {
      REAL(i) /= WINDOW(i);
      IMAG(i) /= WINDOW(i);
    }
}
