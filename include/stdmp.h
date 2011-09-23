/*
 *  stdmp.h
 *  
 *
 *  Created by Andy Eschbacher on 2010-10-12.
 *  
 *
 */
#ifndef _STDMP_H_
#define _STDMP_H_

//void stdmpInit(const double, double*, double*, const int);
void stdmpInit(double*, double*);
//void stdmp(const double, double*, double*, const int, const int);
void stdmp(double*, double*);
//double stdmpLifted(const double, double, double, const int);
double stdmpLifted(double, double, const int);
/*
void noise(double*, double*, const int);
void noise(double*, double*, const int, const int);

*/
#endif // _STDMP_H