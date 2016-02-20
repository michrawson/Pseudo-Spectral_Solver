/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Tue Mar  4 13:50:24 EST 2014 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_r2cb.native -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 32 -name r2cb_32 -include r2cb.h */

/*
 * This function contains 156 FP additions, 84 FP multiplications,
 * (or, 72 additions, 0 multiplications, 84 fused multiply/add),
 * 82 stack variables, 9 constants, and 64 memory accesses
 */
#include "r2cb.h"

static void r2cb_32(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP1_662939224, +1.662939224605090474157576755235811513477121624);
     DK(KP668178637, +0.668178637919298919997757686523080761552472251);
     DK(KP1_961570560, +1.961570560806460898252364472268478073947867462);
     DK(KP198912367, +0.198912367379658006911597622644676228597850501);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(128, rs), MAKE_VOLATILE_STRIDE(128, csr), MAKE_VOLATILE_STRIDE(128, csi)) {
	       E T1F, T1C, T1H, T1z, T1G, T1I;
	       {
		    E T8, T1t, Tz, T1R, T5, T1S, T1u, TE, T1w, TP, T1U, Tg, T2m, T1X, T1x;
		    E TK, T1D, T1d, T20, To, T2p, T28, T1A, TW, T11, T1e, Tv, T25, T23, T2q;
		    E T16, T1f, TA, TD;
		    {
			 E T4, Ty, T1, T2, T6, T7;
			 T4 = Cr[WS(csr, 8)];
			 Ty = Ci[WS(csi, 8)];
			 T1 = Cr[0];
			 T2 = Cr[WS(csr, 16)];
			 T6 = Cr[WS(csr, 4)];
			 T7 = Cr[WS(csr, 12)];
			 {
			      E TB, Tx, T3, TC;
			      TB = Ci[WS(csi, 4)];
			      Tx = T1 - T2;
			      T3 = T1 + T2;
			      TA = T6 - T7;
			      T8 = T6 + T7;
			      TC = Ci[WS(csi, 12)];
			      T1t = FMA(KP2_000000000, Ty, Tx);
			      Tz = FNMS(KP2_000000000, Ty, Tx);
			      T1R = FNMS(KP2_000000000, T4, T3);
			      T5 = FMA(KP2_000000000, T4, T3);
			      TD = TB + TC;
			      T1S = TB - TC;
			 }
		    }
		    {
			 E Td, TG, Tc, T1V, TO, Te, TH, TI;
			 {
			      E Ta, Tb, TM, TN;
			      Ta = Cr[WS(csr, 2)];
			      T1u = TA + TD;
			      TE = TA - TD;
			      Tb = Cr[WS(csr, 14)];
			      TM = Ci[WS(csi, 2)];
			      TN = Ci[WS(csi, 14)];
			      Td = Cr[WS(csr, 10)];
			      TG = Ta - Tb;
			      Tc = Ta + Tb;
			      T1V = TM - TN;
			      TO = TM + TN;
			      Te = Cr[WS(csr, 6)];
			      TH = Ci[WS(csi, 10)];
			      TI = Ci[WS(csi, 6)];
			 }
			 {
			      E Tl, TS, Tk, T26, T1c, Tm, TT, TU;
			      {
				   E Ti, Tj, T1a, T1b;
				   Ti = Cr[WS(csr, 1)];
				   {
					E TL, Tf, T1W, TJ;
					TL = Td - Te;
					Tf = Td + Te;
					T1W = TH - TI;
					TJ = TH + TI;
					T1w = TO - TL;
					TP = TL + TO;
					T1U = Tc - Tf;
					Tg = Tc + Tf;
					T2m = T1W + T1V;
					T1X = T1V - T1W;
					T1x = TG + TJ;
					TK = TG - TJ;
					Tj = Cr[WS(csr, 15)];
				   }
				   T1a = Ci[WS(csi, 1)];
				   T1b = Ci[WS(csi, 15)];
				   Tl = Cr[WS(csr, 9)];
				   TS = Ti - Tj;
				   Tk = Ti + Tj;
				   T26 = T1a - T1b;
				   T1c = T1a + T1b;
				   Tm = Cr[WS(csr, 7)];
				   TT = Ci[WS(csi, 9)];
				   TU = Ci[WS(csi, 7)];
			      }
			      {
				   E Ts, TX, Tr, T22, T10, Tt, T13, T14;
				   {
					E Tp, Tq, TY, TZ;
					Tp = Cr[WS(csr, 5)];
					{
					     E T19, Tn, T27, TV;
					     T19 = Tl - Tm;
					     Tn = Tl + Tm;
					     T27 = TT - TU;
					     TV = TT + TU;
					     T1D = T1c - T19;
					     T1d = T19 + T1c;
					     T20 = Tk - Tn;
					     To = Tk + Tn;
					     T2p = T27 + T26;
					     T28 = T26 - T27;
					     T1A = TS + TV;
					     TW = TS - TV;
					     Tq = Cr[WS(csr, 11)];
					}
					TY = Ci[WS(csi, 5)];
					TZ = Ci[WS(csi, 11)];
					Ts = Cr[WS(csr, 3)];
					TX = Tp - Tq;
					Tr = Tp + Tq;
					T22 = TY - TZ;
					T10 = TY + TZ;
					Tt = Cr[WS(csr, 13)];
					T13 = Ci[WS(csi, 3)];
					T14 = Ci[WS(csi, 13)];
				   }
				   {
					E T12, Tu, T21, T15;
					T11 = TX - T10;
					T1e = TX + T10;
					T12 = Ts - Tt;
					Tu = Ts + Tt;
					T21 = T14 - T13;
					T15 = T13 + T14;
					Tv = Tr + Tu;
					T25 = Tr - Tu;
					T23 = T21 - T22;
					T2q = T22 + T21;
					T16 = T12 - T15;
					T1f = T12 + T15;
				   }
			      }
			 }
		    }
		    {
			 E T1B, T1E, T1l, T1m, T1p, T1o, T1T, T1Y, T29, T2g, T2j, T2f, T2h, T24;
			 {
			      E T1g, T17, T2n, T2t, T2u, T2s;
			      {
				   E T2o, Tw, T2w, T2r, T2l, T9, Th, T2v;
				   T2o = To - Tv;
				   Tw = To + Tv;
				   T2w = T2q + T2p;
				   T2r = T2p - T2q;
				   T1g = T1e - T1f;
				   T1B = T1e + T1f;
				   T17 = T11 + T16;
				   T1E = T16 - T11;
				   T2l = FNMS(KP2_000000000, T8, T5);
				   T9 = FMA(KP2_000000000, T8, T5);
				   Th = FMA(KP2_000000000, Tg, T9);
				   T2v = FNMS(KP2_000000000, Tg, T9);
				   T2n = FNMS(KP2_000000000, T2m, T2l);
				   T2t = FMA(KP2_000000000, T2m, T2l);
				   R0[WS(rs, 4)] = FNMS(KP2_000000000, T2w, T2v);
				   R0[WS(rs, 12)] = FMA(KP2_000000000, T2w, T2v);
				   R0[0] = FMA(KP2_000000000, Tw, Th);
				   R0[WS(rs, 8)] = FNMS(KP2_000000000, Tw, Th);
				   T2u = T2o + T2r;
				   T2s = T2o - T2r;
			      }
			      {
				   E T1j, TR, T18, T1h, TF, TQ;
				   T1l = FNMS(KP1_414213562, TE, Tz);
				   TF = FMA(KP1_414213562, TE, Tz);
				   TQ = FNMS(KP414213562, TP, TK);
				   T1m = FMA(KP414213562, TK, TP);
				   R0[WS(rs, 2)] = FMA(KP1_414213562, T2s, T2n);
				   R0[WS(rs, 10)] = FNMS(KP1_414213562, T2s, T2n);
				   R0[WS(rs, 6)] = FNMS(KP1_414213562, T2u, T2t);
				   R0[WS(rs, 14)] = FMA(KP1_414213562, T2u, T2t);
				   T1j = FNMS(KP1_847759065, TQ, TF);
				   TR = FMA(KP1_847759065, TQ, TF);
				   T1p = FNMS(KP707106781, T17, TW);
				   T18 = FMA(KP707106781, T17, TW);
				   T1h = FMA(KP707106781, T1g, T1d);
				   T1o = FNMS(KP707106781, T1g, T1d);
				   {
					E T2d, T2e, T1k, T1i;
					T1T = FNMS(KP2_000000000, T1S, T1R);
					T2d = FMA(KP2_000000000, T1S, T1R);
					T2e = T1U + T1X;
					T1Y = T1U - T1X;
					T29 = T25 + T28;
					T2g = T28 - T25;
					T1k = FMA(KP198912367, T18, T1h);
					T1i = FNMS(KP198912367, T1h, T18);
					T2j = FMA(KP1_414213562, T2e, T2d);
					T2f = FNMS(KP1_414213562, T2e, T2d);
					R1[WS(rs, 4)] = FNMS(KP1_961570560, T1k, T1j);
					R1[WS(rs, 12)] = FMA(KP1_961570560, T1k, T1j);
					R1[0] = FMA(KP1_961570560, T1i, TR);
					R1[WS(rs, 8)] = FNMS(KP1_961570560, T1i, TR);
					T2h = T20 - T23;
					T24 = T20 + T23;
				   }
			      }
			 }
			 {
			      E T1v, T1y, T1M, T1P, T1L, T1N;
			      {
				   E T1r, T1n, T2k, T2i;
				   T2k = FMA(KP414213562, T2g, T2h);
				   T2i = FNMS(KP414213562, T2h, T2g);
				   T1r = FMA(KP1_847759065, T1m, T1l);
				   T1n = FNMS(KP1_847759065, T1m, T1l);
				   R0[WS(rs, 7)] = FNMS(KP1_847759065, T2k, T2j);
				   R0[WS(rs, 15)] = FMA(KP1_847759065, T2k, T2j);
				   R0[WS(rs, 11)] = FMA(KP1_847759065, T2i, T2f);
				   R0[WS(rs, 3)] = FNMS(KP1_847759065, T2i, T2f);
				   {
					E T1J, T1K, T1s, T1q;
					T1v = FNMS(KP1_414213562, T1u, T1t);
					T1J = FMA(KP1_414213562, T1u, T1t);
					T1K = FMA(KP414213562, T1w, T1x);
					T1y = FNMS(KP414213562, T1x, T1w);
					T1F = FNMS(KP707106781, T1E, T1D);
					T1M = FMA(KP707106781, T1E, T1D);
					T1s = FMA(KP668178637, T1o, T1p);
					T1q = FNMS(KP668178637, T1p, T1o);
					T1P = FMA(KP1_847759065, T1K, T1J);
					T1L = FNMS(KP1_847759065, T1K, T1J);
					R1[WS(rs, 6)] = FNMS(KP1_662939224, T1s, T1r);
					R1[WS(rs, 14)] = FMA(KP1_662939224, T1s, T1r);
					R1[WS(rs, 10)] = FMA(KP1_662939224, T1q, T1n);
					R1[WS(rs, 2)] = FNMS(KP1_662939224, T1q, T1n);
					T1N = FMA(KP707106781, T1B, T1A);
					T1C = FNMS(KP707106781, T1B, T1A);
				   }
			      }
			      {
				   E T2b, T1Z, T1Q, T1O, T2c, T2a;
				   T1Q = FMA(KP198912367, T1M, T1N);
				   T1O = FNMS(KP198912367, T1N, T1M);
				   T2b = FNMS(KP1_414213562, T1Y, T1T);
				   T1Z = FMA(KP1_414213562, T1Y, T1T);
				   R1[WS(rs, 7)] = FNMS(KP1_961570560, T1Q, T1P);
				   R1[WS(rs, 15)] = FMA(KP1_961570560, T1Q, T1P);
				   R1[WS(rs, 11)] = FMA(KP1_961570560, T1O, T1L);
				   R1[WS(rs, 3)] = FNMS(KP1_961570560, T1O, T1L);
				   T2c = FMA(KP414213562, T24, T29);
				   T2a = FNMS(KP414213562, T29, T24);
				   T1H = FMA(KP1_847759065, T1y, T1v);
				   T1z = FNMS(KP1_847759065, T1y, T1v);
				   R0[WS(rs, 5)] = FNMS(KP1_847759065, T2c, T2b);
				   R0[WS(rs, 13)] = FMA(KP1_847759065, T2c, T2b);
				   R0[WS(rs, 1)] = FMA(KP1_847759065, T2a, T1Z);
				   R0[WS(rs, 9)] = FNMS(KP1_847759065, T2a, T1Z);
			      }
			 }
		    }
	       }
	       T1G = FNMS(KP668178637, T1F, T1C);
	       T1I = FMA(KP668178637, T1C, T1F);
	       R1[WS(rs, 5)] = FNMS(KP1_662939224, T1I, T1H);
	       R1[WS(rs, 13)] = FMA(KP1_662939224, T1I, T1H);
	       R1[WS(rs, 1)] = FMA(KP1_662939224, T1G, T1z);
	       R1[WS(rs, 9)] = FNMS(KP1_662939224, T1G, T1z);
	  }
     }
}

static const kr2c_desc desc = { 32, "r2cb_32", {72, 0, 84, 0}, &GENUS };

void X(codelet_r2cb_32) (planner *p) {
     X(kr2c_register) (p, r2cb_32, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_r2cb.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 32 -name r2cb_32 -include r2cb.h */

/*
 * This function contains 156 FP additions, 50 FP multiplications,
 * (or, 140 additions, 34 multiplications, 16 fused multiply/add),
 * 54 stack variables, 9 constants, and 64 memory accesses
 */
#include "r2cb.h"

static void r2cb_32(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP1_662939224, +1.662939224605090474157576755235811513477121624);
     DK(KP1_111140466, +1.111140466039204449485661627897065748749874382);
     DK(KP1_961570560, +1.961570560806460898252364472268478073947867462);
     DK(KP390180644, +0.390180644032256535696569736954044481855383236);
     DK(KP765366864, +0.765366864730179543456919968060797733522689125);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(128, rs), MAKE_VOLATILE_STRIDE(128, csr), MAKE_VOLATILE_STRIDE(128, csi)) {
	       E T9, T2c, TB, T1y, T6, T2b, Ty, T1v, Th, T2e, T2f, TD, TK, T1C, T1F;
	       E T1h, Tp, T2i, T2m, TN, T13, T1K, T1Y, T1k, Tw, TU, T1l, TW, T1V, T2j;
	       E T1R, T2l;
	       {
		    E T7, T8, T1w, Tz, TA, T1x;
		    T7 = Cr[WS(csr, 4)];
		    T8 = Cr[WS(csr, 12)];
		    T1w = T7 - T8;
		    Tz = Ci[WS(csi, 4)];
		    TA = Ci[WS(csi, 12)];
		    T1x = Tz + TA;
		    T9 = KP2_000000000 * (T7 + T8);
		    T2c = KP1_414213562 * (T1w + T1x);
		    TB = KP2_000000000 * (Tz - TA);
		    T1y = KP1_414213562 * (T1w - T1x);
	       }
	       {
		    E T5, T1u, T3, T1s;
		    {
			 E T4, T1t, T1, T2;
			 T4 = Cr[WS(csr, 8)];
			 T5 = KP2_000000000 * T4;
			 T1t = Ci[WS(csi, 8)];
			 T1u = KP2_000000000 * T1t;
			 T1 = Cr[0];
			 T2 = Cr[WS(csr, 16)];
			 T3 = T1 + T2;
			 T1s = T1 - T2;
		    }
		    T6 = T3 + T5;
		    T2b = T1s + T1u;
		    Ty = T3 - T5;
		    T1v = T1s - T1u;
	       }
	       {
		    E Td, T1A, TG, T1E, Tg, T1D, TJ, T1B;
		    {
			 E Tb, Tc, TE, TF;
			 Tb = Cr[WS(csr, 2)];
			 Tc = Cr[WS(csr, 14)];
			 Td = Tb + Tc;
			 T1A = Tb - Tc;
			 TE = Ci[WS(csi, 2)];
			 TF = Ci[WS(csi, 14)];
			 TG = TE - TF;
			 T1E = TE + TF;
		    }
		    {
			 E Te, Tf, TH, TI;
			 Te = Cr[WS(csr, 10)];
			 Tf = Cr[WS(csr, 6)];
			 Tg = Te + Tf;
			 T1D = Te - Tf;
			 TH = Ci[WS(csi, 10)];
			 TI = Ci[WS(csi, 6)];
			 TJ = TH - TI;
			 T1B = TH + TI;
		    }
		    Th = KP2_000000000 * (Td + Tg);
		    T2e = T1A + T1B;
		    T2f = T1E - T1D;
		    TD = Td - Tg;
		    TK = TG - TJ;
		    T1C = T1A - T1B;
		    T1F = T1D + T1E;
		    T1h = KP2_000000000 * (TJ + TG);
	       }
	       {
		    E Tl, T1I, TZ, T1X, To, T1W, T12, T1J;
		    {
			 E Tj, Tk, TX, TY;
			 Tj = Cr[WS(csr, 1)];
			 Tk = Cr[WS(csr, 15)];
			 Tl = Tj + Tk;
			 T1I = Tj - Tk;
			 TX = Ci[WS(csi, 1)];
			 TY = Ci[WS(csi, 15)];
			 TZ = TX - TY;
			 T1X = TX + TY;
		    }
		    {
			 E Tm, Tn, T10, T11;
			 Tm = Cr[WS(csr, 9)];
			 Tn = Cr[WS(csr, 7)];
			 To = Tm + Tn;
			 T1W = Tm - Tn;
			 T10 = Ci[WS(csi, 9)];
			 T11 = Ci[WS(csi, 7)];
			 T12 = T10 - T11;
			 T1J = T10 + T11;
		    }
		    Tp = Tl + To;
		    T2i = T1I + T1J;
		    T2m = T1X - T1W;
		    TN = Tl - To;
		    T13 = TZ - T12;
		    T1K = T1I - T1J;
		    T1Y = T1W + T1X;
		    T1k = T12 + TZ;
	       }
	       {
		    E Ts, T1L, TT, T1M, Tv, T1O, TQ, T1P;
		    {
			 E Tq, Tr, TR, TS;
			 Tq = Cr[WS(csr, 5)];
			 Tr = Cr[WS(csr, 11)];
			 Ts = Tq + Tr;
			 T1L = Tq - Tr;
			 TR = Ci[WS(csi, 5)];
			 TS = Ci[WS(csi, 11)];
			 TT = TR - TS;
			 T1M = TR + TS;
		    }
		    {
			 E Tt, Tu, TO, TP;
			 Tt = Cr[WS(csr, 3)];
			 Tu = Cr[WS(csr, 13)];
			 Tv = Tt + Tu;
			 T1O = Tt - Tu;
			 TO = Ci[WS(csi, 13)];
			 TP = Ci[WS(csi, 3)];
			 TQ = TO - TP;
			 T1P = TP + TO;
		    }
		    Tw = Ts + Tv;
		    TU = TQ - TT;
		    T1l = TT + TQ;
		    TW = Ts - Tv;
		    {
			 E T1T, T1U, T1N, T1Q;
			 T1T = T1L + T1M;
			 T1U = T1O + T1P;
			 T1V = KP707106781 * (T1T - T1U);
			 T2j = KP707106781 * (T1T + T1U);
			 T1N = T1L - T1M;
			 T1Q = T1O - T1P;
			 T1R = KP707106781 * (T1N + T1Q);
			 T2l = KP707106781 * (T1N - T1Q);
		    }
	       }
	       {
		    E Tx, T1r, Ti, T1q, Ta;
		    Tx = KP2_000000000 * (Tp + Tw);
		    T1r = KP2_000000000 * (T1l + T1k);
		    Ta = T6 + T9;
		    Ti = Ta + Th;
		    T1q = Ta - Th;
		    R0[WS(rs, 8)] = Ti - Tx;
		    R0[WS(rs, 12)] = T1q + T1r;
		    R0[0] = Ti + Tx;
		    R0[WS(rs, 4)] = T1q - T1r;
	       }
	       {
		    E T1i, T1o, T1n, T1p, T1g, T1j, T1m;
		    T1g = T6 - T9;
		    T1i = T1g - T1h;
		    T1o = T1g + T1h;
		    T1j = Tp - Tw;
		    T1m = T1k - T1l;
		    T1n = KP1_414213562 * (T1j - T1m);
		    T1p = KP1_414213562 * (T1j + T1m);
		    R0[WS(rs, 10)] = T1i - T1n;
		    R0[WS(rs, 14)] = T1o + T1p;
		    R0[WS(rs, 2)] = T1i + T1n;
		    R0[WS(rs, 6)] = T1o - T1p;
	       }
	       {
		    E TM, T16, T15, T17;
		    {
			 E TC, TL, TV, T14;
			 TC = Ty - TB;
			 TL = KP1_414213562 * (TD - TK);
			 TM = TC + TL;
			 T16 = TC - TL;
			 TV = TN + TU;
			 T14 = TW + T13;
			 T15 = FNMS(KP765366864, T14, KP1_847759065 * TV);
			 T17 = FMA(KP765366864, TV, KP1_847759065 * T14);
		    }
		    R0[WS(rs, 9)] = TM - T15;
		    R0[WS(rs, 13)] = T16 + T17;
		    R0[WS(rs, 1)] = TM + T15;
		    R0[WS(rs, 5)] = T16 - T17;
	       }
	       {
		    E T2t, T2x, T2w, T2y;
		    {
			 E T2r, T2s, T2u, T2v;
			 T2r = T2b + T2c;
			 T2s = FMA(KP1_847759065, T2e, KP765366864 * T2f);
			 T2t = T2r - T2s;
			 T2x = T2r + T2s;
			 T2u = T2i + T2j;
			 T2v = T2m - T2l;
			 T2w = FNMS(KP1_961570560, T2v, KP390180644 * T2u);
			 T2y = FMA(KP1_961570560, T2u, KP390180644 * T2v);
		    }
		    R1[WS(rs, 11)] = T2t - T2w;
		    R1[WS(rs, 15)] = T2x + T2y;
		    R1[WS(rs, 3)] = T2t + T2w;
		    R1[WS(rs, 7)] = T2x - T2y;
	       }
	       {
		    E T1a, T1e, T1d, T1f;
		    {
			 E T18, T19, T1b, T1c;
			 T18 = Ty + TB;
			 T19 = KP1_414213562 * (TD + TK);
			 T1a = T18 - T19;
			 T1e = T18 + T19;
			 T1b = TN - TU;
			 T1c = T13 - TW;
			 T1d = FNMS(KP1_847759065, T1c, KP765366864 * T1b);
			 T1f = FMA(KP1_847759065, T1b, KP765366864 * T1c);
		    }
		    R0[WS(rs, 11)] = T1a - T1d;
		    R0[WS(rs, 15)] = T1e + T1f;
		    R0[WS(rs, 3)] = T1a + T1d;
		    R0[WS(rs, 7)] = T1e - T1f;
	       }
	       {
		    E T25, T29, T28, T2a;
		    {
			 E T23, T24, T26, T27;
			 T23 = T1v - T1y;
			 T24 = FMA(KP765366864, T1C, KP1_847759065 * T1F);
			 T25 = T23 - T24;
			 T29 = T23 + T24;
			 T26 = T1K - T1R;
			 T27 = T1Y - T1V;
			 T28 = FNMS(KP1_662939224, T27, KP1_111140466 * T26);
			 T2a = FMA(KP1_662939224, T26, KP1_111140466 * T27);
		    }
		    R1[WS(rs, 10)] = T25 - T28;
		    R1[WS(rs, 14)] = T29 + T2a;
		    R1[WS(rs, 2)] = T25 + T28;
		    R1[WS(rs, 6)] = T29 - T2a;
	       }
	       {
		    E T2h, T2p, T2o, T2q;
		    {
			 E T2d, T2g, T2k, T2n;
			 T2d = T2b - T2c;
			 T2g = FNMS(KP1_847759065, T2f, KP765366864 * T2e);
			 T2h = T2d + T2g;
			 T2p = T2d - T2g;
			 T2k = T2i - T2j;
			 T2n = T2l + T2m;
			 T2o = FNMS(KP1_111140466, T2n, KP1_662939224 * T2k);
			 T2q = FMA(KP1_111140466, T2k, KP1_662939224 * T2n);
		    }
		    R1[WS(rs, 9)] = T2h - T2o;
		    R1[WS(rs, 13)] = T2p + T2q;
		    R1[WS(rs, 1)] = T2h + T2o;
		    R1[WS(rs, 5)] = T2p - T2q;
	       }
	       {
		    E T1H, T21, T20, T22;
		    {
			 E T1z, T1G, T1S, T1Z;
			 T1z = T1v + T1y;
			 T1G = FNMS(KP765366864, T1F, KP1_847759065 * T1C);
			 T1H = T1z + T1G;
			 T21 = T1z - T1G;
			 T1S = T1K + T1R;
			 T1Z = T1V + T1Y;
			 T20 = FNMS(KP390180644, T1Z, KP1_961570560 * T1S);
			 T22 = FMA(KP390180644, T1S, KP1_961570560 * T1Z);
		    }
		    R1[WS(rs, 8)] = T1H - T20;
		    R1[WS(rs, 12)] = T21 + T22;
		    R1[0] = T1H + T20;
		    R1[WS(rs, 4)] = T21 - T22;
	       }
	  }
     }
}

static const kr2c_desc desc = { 32, "r2cb_32", {140, 34, 16, 0}, &GENUS };

void X(codelet_r2cb_32) (planner *p) {
     X(kr2c_register) (p, r2cb_32, &desc);
}

#endif				/* HAVE_FMA */
