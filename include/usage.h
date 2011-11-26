
#ifndef USAGE_H
#define USAGE_H

#include <iostream>
#include <string>
#include <cstring>

using std::cerr;
using std::endl;
using std::string;

void version(void);
void usage_rrstdmp(double, double, long, unsigned int, unsigned int, double, const int);
void usage_rrmapgen(double, double, long, unsigned int, unsigned int, double, const int);

#endif // USAGE_H