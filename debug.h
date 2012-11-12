/*
 * debug.h
 *
 *  Created on: May 4, 2012
 *      Author: dstuck
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#ifdef DEBUG
#define LINE std::cout << "file:" << __FILE__ << " line:" << __LINE__ << std::endl;
#else
#define LINE
#endif


#endif /* DEBUG_H_ */
