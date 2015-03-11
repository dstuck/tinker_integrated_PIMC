/*
 * CallTinker.h
 *
 *  Created on: Oct 30, 2012
 *      Author: dstuck
 */

#ifndef CALLTINKER_H_
#define CALLTINKER_H_

extern"C" {
//extern void dstuckenergy_(const int&, int*, double*, double*, const char*, const int&, double&, bool&);
extern void dstuckenergy_(const int&, int*, double*, double*, const char*, const int&, double&, int&);
extern void dstuckvibrate_(int*, double*, double*, const int&, const int&, double*, double*, double*);
extern void dstuckoptimize_(int*, double*, double*, const int&, const int&);
//extern void dstuckmechanic_(const char*, const int&);
}

#endif /* CALLTINKER_H_ */
