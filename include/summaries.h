
#ifndef _SUMMARIES_H_
#define _SUMMARIES_H_

#include <ctime>

//void datasummary(double, double, unsigned long int, unsigned long int, mpz_t, mpz_t, double, double*, double*, time_t);
void timesummary(time_t, time_t&);
void errorMsg(char*);
void errorMsg(char*, char*);

#endif