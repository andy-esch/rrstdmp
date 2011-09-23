
#ifndef _USAGE_H
#define _USAGE_H

#include "gmp.h"
#include "gmpxx.h"

void version(void);
void usage_rrstdmp(double, double, mpz_t, unsigned int, unsigned int, double, const int);
void usage_rrmapgen(double, double, mpz_t, unsigned int, unsigned int, double, const int);
//void changeLog(void);

#endif