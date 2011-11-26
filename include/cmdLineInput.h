/*
 *  cmdLineInput.cpp
 *  
 *  Description:
 *
 *
 *  Created by Peter Eschbacher on 11/10/11.
 *
 */

#ifndef CMDLINEINPUT_H
#define CMDLINEINPUT_H

#include <iostream> // cin, cout, endl
#include <cstdlib>  // atoi, atof, exit(),   
#include <cstring> // strcpy

#include "usage.h"

using std::cin;
using std::cerr;
using std::cout;
using std::endl;

extern const double TWOPI;
extern double k, ge;
extern int globalWindow, globalOverlap;

void cmdLineInput(int, char**, char *, char*, long&, double&, int&, \
                  bool&, bool&, bool&, bool&);

#endif // CMDLINEINPUT_H