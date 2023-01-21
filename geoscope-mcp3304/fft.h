/*
 * fft.h: Defs for fft routines
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
 * $Id: fft.h,v 1.6 1992/04/02 06:56:44 steve Exp steve $
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
 * If you want double-precision floats used everywhere,
 * define this to be "double"
 */
#define FLOAT 	double

/*
 * Since every math.h seems to have their own idea of what the
 * elements in struct complex are called, we define our own.
 */
struct _complx
{
  FLOAT re;
  FLOAT im;
};
typedef struct _complx complx;

/*
 * Functions provided by fft.c
 */
void		fft (complx wavetrum[], int bits);
void		ifft (complx wavetrum[], int bits);
void		scale (complx wavetrum[], int bits);
void		reorder (complx wavetrum[], int bits);
void		fft_window (complx wavetrum[], int bits);
void		fft_unwindow (complx wavetrum[], int bits);
void		fft_init (int bits);
FLOAT		phase (int n, complx wavetrum[], int bits);
FLOAT		amp (int n, complx wavetrum[], int bits);
inline const unsigned int reverse (unsigned int val, int bits);
