# makefile for rrstdmp -- a program to calculate the recurrence rate of the 
#  standard map written by andy eschbacher, Feb. 2011
vpath %.cpp src
vpath %.h include

CXX = g++-mp-4.4
CXXFLAGS = -L/opt/local/lib -lgmpxx -lgmp -lgsl -I/opt/local/include -I include

OBJECTS = usage.o stdmp.o rr.o rrmean.o rrtester.o summaries.o misc.o
HEADERS = usage.h stdmp.h rr.h rrmean.h rrtester.h summaries.h misc.h
BLT = "--->  Making"
SEPR = "\n\t\c"

ifeq "$(MAKECMDGOALS)" "rrmapgen"
	OBJECTS += stdmpout.o
	HEADERS += stdmpout.h
endif

# =-=-=-=-=-=-=-=- make everything -=-=-=-=-=-=-=-=-= #
.PHONY: all
all: rrstdmp rrmapgen rrmeangen

# -=-=-=-=-=-=-=-=-=-= rrstdmp -=-=-=-=-=-=-=-=-=-=-= #
# rrstdmp is the main program for generating and analyzing
#   standard map recurrence rate.  Its output is a histogram
#   of sticky events, and, if chosen, the fits to that distribution
#

rrstdmp: rrstdmp.o $(OBJECTS)
	@echo $(BLT) $@ $(SEPR)
	$(CXX) $(CXXFLAGS) $^ -o $@

rrstdmp.o: rrstdmp.cpp
	@echo $(BLT) $@ $(SEPR)
	$(CXX) -I include -c $< -o $@
#	@echo $(BLT) "rrstdmp.o: \n\t\c"
#	$(CXX) $(CXXFLAGS) -c rrstdmp.cpp

# -=-=-=-=-=-=-=-=-=-= rrmapgen =-=-=-=-=-=-=-=-=-=-= #
#  rrmapgen shares much of the source code with rrstdmp but 
#    outputs trajectories of sticking regions to file for later
#    analysis, such as with txt2png.sh
#

rrmapgen: rrmapgen.o $(OBJECTS)
	@echo $(BLT) $@ $(SEPR)
	$(CXX) $(CXXFLAGS) $^ -o $@

rrmapgen.o: rrmapgen.cpp
	@echo $(BLT) $@ $(SEPR)
	$(CXX) -I include -c $< -o $@

# -=-=-=-=-=-=-=-=-=-= rrmeangen -=-=-=-=-=-=-=-=-=-=-= #
# rrmeangen is the main program for generating and analyzing
#   standard map recurrence rate according to mean recurrence
#   and its relation to its sticking event
#

rrmeangen: rrmeangen.o $(OBJECTS)
	@echo $(BLT) $@ $(SEPR)
	$(CXX) $(CXXFLAGS) $^ -o $@

rrmeangen.o: rrmeangen.cpp
	$(CXX) -I include -c $< -o $@

# -=-=-=-=-=-=-=-=-=-= dependents -=-=-=-=-=-=-=-=-=-=-= #
#
$(OBJECTS): %.o: %.cpp
	@echo $(BLT) $@ $(SEPR)
	$(CXX) -I include -c $< -o $@

# -=-=-=-=-=-=-=-=-= other files -=-=-=-=-=-=-=-=-=-= #
#  other files...
#
#  standard non-twist map module
snm.o: snm.h

# -=-=-=-=-=-=-=-=- Utilities, etc. =-=-=-=-=-=-=-=-= #
# Installs rrstdmp into home bin
.PHONY: install
install:
	@echo "--->  Staging rrstdmp into ~/bin"
	mv rrstdmp ~/bin/

#.PHONY: zip
#zip:
#	@echo "---> Zipping file to backup."
#	cd ~/Desktop
#	mkdir backup_rrstdmp
#	cd backup_rrstdmp
#	mkdir src/
#	mkdir include/
#	cd ~/Project/rrstdmp/
#	cp -r src/ ~/Desktop/backup-rrstdmp/src/
#	cp -r include/ ~/Desktop/backup-rrstdmp/include/

# Cleans object files
.PHONY: clean
clean:
ifneq ("$(shell ls | grep '.*\.o')","")
	@echo "--->  Cleaning object files \n\t\c"
	rm $(shell ls | grep ".*\.o")
else
	@echo "--->  No object files to remove"
endif

