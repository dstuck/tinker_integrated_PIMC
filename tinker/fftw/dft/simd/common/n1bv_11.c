/*
 * Copyright (c) 2003, 2007-11 Matteo Frigo
 * Copyright (c) 2003, 2007-11 Massachusetts Institute of Technology
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Wed Jul 27 06:13:38 EDT 2011 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw_c.native -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -sign 1 -n 11 -name n1bv_11 -include n1b.h */

/*
 * This function contains 70 FP additions, 60 FP multiplications,
 * (or, 15 additions, 5 multiplications, 55 fused multiply/add),
 * 67 stack variables, 11 constants, and 22 memory accesses
 */
#include "n1b.h"

static void n1bv_11(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP959492973, +0.959492973614497389890368057066327699062454848);
     DVK(KP876768831, +0.876768831002589333891339807079336796764054852);
     DVK(KP918985947, +0.918985947228994779780736114132655398124909697);
     DVK(KP989821441, +0.989821441880932732376092037776718787376519372);
     DVK(KP778434453, +0.778434453334651800608337670740821884709317477);
     DVK(KP830830026, +0.830830026003772851058548298459246407048009821);
     DVK(KP372785597, +0.372785597771792209609773152906148328659002598);
     DVK(KP634356270, +0.634356270682424498893150776899916060542806975);
     DVK(KP715370323, +0.715370323453429719112414662767260662417897278);
     DVK(KP342584725, +0.342584725681637509502641509861112333758894680);
     DVK(KP521108558, +0.521108558113202722944698153526659300680427422);
     {
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ii;
	  xo = io;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	       V T1, Tb, T4, Tq, Tg, Tm, T7, Tp, Ta, To, Tc, T11;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       {
		    V T2, T3, Te, Tf;
		    T2 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T3 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
		    Te = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    Tf = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    {
			 V T5, T6, T8, T9;
			 T5 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 T6 = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 T8 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 T9 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 Tb = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 T4 = VADD(T2, T3);
			 Tq = VSUB(T2, T3);
			 Tg = VADD(Te, Tf);
			 Tm = VSUB(Te, Tf);
			 T7 = VADD(T5, T6);
			 Tp = VSUB(T5, T6);
			 Ta = VADD(T8, T9);
			 To = VSUB(T8, T9);
			 Tc = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
		    }
	       }
	       T11 = VFMA(LDK(KP521108558), Tm, Tq);
	       {
		    V TA, TS, TE, TW, Td, Tn, Ts, Tw, Tr, Tv, TT, TF;
		    Tr = VFNMS(LDK(KP521108558), Tq, Tp);
		    Tv = VFNMS(LDK(KP342584725), T7, Tg);
		    TA = VFMA(LDK(KP715370323), To, Tq);
		    TS = VFMA(LDK(KP521108558), To, Tm);
		    TE = VFNMS(LDK(KP342584725), T4, Ta);
		    TW = VFNMS(LDK(KP342584725), Ta, T7);
		    Td = VADD(Tb, Tc);
		    Tn = VSUB(Tb, Tc);
		    Ts = VFNMS(LDK(KP715370323), Tr, To);
		    Tw = VFNMS(LDK(KP634356270), Tv, T4);
		    TT = VFNMS(LDK(KP715370323), TS, Tp);
		    TF = VFNMS(LDK(KP634356270), TE, Tg);
		    {
			 V Tu, TV, TD, TL, T14, TP, TZ, Tj, Tz, TI, TB, TJ, TM;
			 TB = VFMA(LDK(KP372785597), Tn, TA);
			 TJ = VFNMS(LDK(KP521108558), Tp, Tn);
			 {
			      V T12, TN, TX, Th;
			      T12 = VFMA(LDK(KP715370323), T11, Tn);
			      ST(&(xo[0]), VADD(Tg, VADD(Td, VADD(Ta, VADD(T7, VADD(T4, T1))))), ovs, &(xo[0]));
			      TN = VFNMS(LDK(KP342584725), Td, T4);
			      TX = VFNMS(LDK(KP634356270), TW, Td);
			      Th = VFNMS(LDK(KP342584725), Tg, Td);
			      {
				   V Tt, Tx, TU, TG;
				   Tt = VFNMS(LDK(KP830830026), Ts, Tn);
				   Tx = VFNMS(LDK(KP778434453), Tw, Ta);
				   TU = VFMA(LDK(KP830830026), TT, Tq);
				   TG = VFNMS(LDK(KP778434453), TF, Td);
				   {
					V TC, TK, T13, TO;
					TC = VFNMS(LDK(KP830830026), TB, Tm);
					TK = VFMA(LDK(KP715370323), TJ, Tm);
					T13 = VFMA(LDK(KP830830026), T12, Tp);
					TO = VFNMS(LDK(KP634356270), TN, T7);
					{
					     V TY, Ti, Ty, TH;
					     TY = VFNMS(LDK(KP778434453), TX, T4);
					     Ti = VFNMS(LDK(KP634356270), Th, Ta);
					     Tu = VMUL(LDK(KP989821441), VFNMS(LDK(KP918985947), Tt, Tm));
					     Ty = VFNMS(LDK(KP876768831), Tx, Td);
					     TV = VMUL(LDK(KP989821441), VFNMS(LDK(KP918985947), TU, Tn));
					     TH = VFNMS(LDK(KP876768831), TG, T7);
					     TD = VMUL(LDK(KP989821441), VFMA(LDK(KP918985947), TC, Tp));
					     TL = VFNMS(LDK(KP830830026), TK, To);
					     T14 = VMUL(LDK(KP989821441), VFMA(LDK(KP918985947), T13, To));
					     TP = VFNMS(LDK(KP778434453), TO, Tg);
					     TZ = VFNMS(LDK(KP876768831), TY, Tg);
					     Tj = VFNMS(LDK(KP778434453), Ti, T7);
					     Tz = VFNMS(LDK(KP959492973), Ty, T1);
					     TI = VFNMS(LDK(KP959492973), TH, T1);
					}
				   }
			      }
			 }
			 TM = VMUL(LDK(KP989821441), VFNMS(LDK(KP918985947), TL, Tq));
			 {
			      V TQ, T10, Tk, TR, Tl;
			      TQ = VFNMS(LDK(KP876768831), TP, Ta);
			      T10 = VFNMS(LDK(KP959492973), TZ, T1);
			      Tk = VFNMS(LDK(KP876768831), Tj, T4);
			      ST(&(xo[WS(os, 7)]), VFMAI(TD, Tz), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 4)]), VFNMSI(TD, Tz), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 8)]), VFNMSI(TM, TI), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 3)]), VFMAI(TM, TI), ovs, &(xo[WS(os, 1)]));
			      TR = VFNMS(LDK(KP959492973), TQ, T1);
			      ST(&(xo[WS(os, 10)]), VFNMSI(T14, T10), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 1)]), VFMAI(T14, T10), ovs, &(xo[WS(os, 1)]));
			      Tl = VFNMS(LDK(KP959492973), Tk, T1);
			      ST(&(xo[WS(os, 9)]), VFMAI(TV, TR), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 2)]), VFNMSI(TV, TR), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 6)]), VFNMSI(Tu, Tl), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 5)]), VFMAI(Tu, Tl), ovs, &(xo[WS(os, 1)]));
			 }
		    }
	       }
	  }
     }
     VLEAVE();
}

static const kdft_desc desc = { 11, XSIMD_STRING("n1bv_11"), {15, 5, 55, 0}, &GENUS, 0, 0, 0, 0 };

void XSIMD(codelet_n1bv_11) (planner *p) {
     X(kdft_register) (p, n1bv_11, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw_c.native -simd -compact -variables 4 -pipeline-latency 8 -sign 1 -n 11 -name n1bv_11 -include n1b.h */

/*
 * This function contains 70 FP additions, 50 FP multiplications,
 * (or, 30 additions, 10 multiplications, 40 fused multiply/add),
 * 32 stack variables, 10 constants, and 22 memory accesses
 */
#include "n1b.h"

static void n1bv_11(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP959492973, +0.959492973614497389890368057066327699062454848);
     DVK(KP654860733, +0.654860733945285064056925072466293553183791199);
     DVK(KP142314838, +0.142314838273285140443792668616369668791051361);
     DVK(KP415415013, +0.415415013001886425529274149229623203524004910);
     DVK(KP841253532, +0.841253532831181168861811648919367717513292498);
     DVK(KP540640817, +0.540640817455597582107635954318691695431770608);
     DVK(KP909631995, +0.909631995354518371411715383079028460060241051);
     DVK(KP989821441, +0.989821441880932732376092037776718787376519372);
     DVK(KP755749574, +0.755749574354258283774035843972344420179717445);
     DVK(KP281732556, +0.281732556841429697711417915346616899035777899);
     {
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ii;
	  xo = io;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	       V Th, T3, Tm, Tf, Ti, Tc, Tj, T9, Tk, T6, Tl, Ta, Tb, Ts, Tt;
	       Th = LD(&(xi[0]), ivs, &(xi[0]));
	       {
		    V T1, T2, Td, Te;
		    T1 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T2 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
		    T3 = VSUB(T1, T2);
		    Tm = VADD(T1, T2);
		    Td = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    Te = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
		    Tf = VSUB(Td, Te);
		    Ti = VADD(Td, Te);
	       }
	       Ta = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
	       Tb = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
	       Tc = VSUB(Ta, Tb);
	       Tj = VADD(Ta, Tb);
	       {
		    V T7, T8, T4, T5;
		    T7 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    T8 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    T9 = VSUB(T7, T8);
		    Tk = VADD(T7, T8);
		    T4 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    T5 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
		    T6 = VSUB(T4, T5);
		    Tl = VADD(T4, T5);
	       }
	       ST(&(xo[0]), VADD(Th, VADD(Tm, VADD(Ti, VADD(Tl, VADD(Tj, Tk))))), ovs, &(xo[0]));
	       {
		    V Tg, Tn, Tu, Tv;
		    Tg = VBYI(VFMA(LDK(KP281732556), T3, VFMA(LDK(KP755749574), T6, VFNMS(LDK(KP909631995), Tc, VFNMS(LDK(KP540640817), Tf, VMUL(LDK(KP989821441), T9))))));
		    Tn = VFMA(LDK(KP841253532), Ti, VFMA(LDK(KP415415013), Tj, VFNMS(LDK(KP142314838), Tk, VFNMS(LDK(KP654860733), Tl, VFNMS(LDK(KP959492973), Tm, Th)))));
		    ST(&(xo[WS(os, 5)]), VADD(Tg, Tn), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 6)]), VSUB(Tn, Tg), ovs, &(xo[0]));
		    Tu = VBYI(VFMA(LDK(KP755749574), T3, VFMA(LDK(KP540640817), T6, VFNMS(LDK(KP909631995), T9, VFNMS(LDK(KP989821441), Tf, VMUL(LDK(KP281732556), Tc))))));
		    Tv = VFMA(LDK(KP841253532), Tl, VFMA(LDK(KP415415013), Tk, VFNMS(LDK(KP959492973), Tj, VFNMS(LDK(KP142314838), Ti, VFNMS(LDK(KP654860733), Tm, Th)))));
		    ST(&(xo[WS(os, 4)]), VADD(Tu, Tv), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 7)]), VSUB(Tv, Tu), ovs, &(xo[WS(os, 1)]));
	       }
	       Ts = VBYI(VFMA(LDK(KP909631995), T3, VFNMS(LDK(KP540640817), T9, VFNMS(LDK(KP989821441), Tc, VFNMS(LDK(KP281732556), T6, VMUL(LDK(KP755749574), Tf))))));
	       Tt = VFMA(LDK(KP415415013), Tm, VFMA(LDK(KP841253532), Tk, VFNMS(LDK(KP142314838), Tj, VFNMS(LDK(KP959492973), Tl, VFNMS(LDK(KP654860733), Ti, Th)))));
	       ST(&(xo[WS(os, 2)]), VADD(Ts, Tt), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 9)]), VSUB(Tt, Ts), ovs, &(xo[WS(os, 1)]));
	       {
		    V Tq, Tr, To, Tp;
		    Tq = VBYI(VFMA(LDK(KP540640817), T3, VFMA(LDK(KP909631995), Tf, VFMA(LDK(KP989821441), T6, VFMA(LDK(KP755749574), Tc, VMUL(LDK(KP281732556), T9))))));
		    Tr = VFMA(LDK(KP841253532), Tm, VFMA(LDK(KP415415013), Ti, VFNMS(LDK(KP959492973), Tk, VFNMS(LDK(KP654860733), Tj, VFNMS(LDK(KP142314838), Tl, Th)))));
		    ST(&(xo[WS(os, 1)]), VADD(Tq, Tr), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 10)]), VSUB(Tr, Tq), ovs, &(xo[0]));
		    To = VBYI(VFMA(LDK(KP989821441), T3, VFMA(LDK(KP540640817), Tc, VFNMS(LDK(KP909631995), T6, VFNMS(LDK(KP281732556), Tf, VMUL(LDK(KP755749574), T9))))));
		    Tp = VFMA(LDK(KP415415013), Tl, VFMA(LDK(KP841253532), Tj, VFNMS(LDK(KP654860733), Tk, VFNMS(LDK(KP959492973), Ti, VFNMS(LDK(KP142314838), Tm, Th)))));
		    ST(&(xo[WS(os, 3)]), VADD(To, Tp), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 8)]), VSUB(Tp, To), ovs, &(xo[0]));
	       }
	  }
     }
     VLEAVE();
}

static const kdft_desc desc = { 11, XSIMD_STRING("n1bv_11"), {30, 10, 40, 0}, &GENUS, 0, 0, 0, 0 };

void XSIMD(codelet_n1bv_11) (planner *p) {
     X(kdft_register) (p, n1bv_11, &desc);
}

#endif				/* HAVE_FMA */
