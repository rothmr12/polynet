# Makefile for PolyNet, which performs mesoscale
# simulations of a polymer network.
#
# John L. Barber
# 2012 February
# jlbarber@lanl.gov
# 505-664-0605
#
# 2014 April 25: Revised to change names of files
#
# 2018 Feb 23: Added "ReadInputDeck" and "Point" and divided PolyNet into multiple files
#
# 2018 August 14: Replaced Point with faster Point3D

OBJECT_FILES = Point3D.o ReadInputDeck.o PolyNet_Setup.o PolyNet_Stepping.o PolyNet_Diagnostics.o PolyNet_Output.o PolyNet_Misc.o
HEADER_FILES = Point3D.h ReadInputDeck.h PolyNet.h

#CFLAGS = -Og
#CFLAGS = -Og -march=native
#CFLAGS = -O3
CFLAGS = -Ofast
#CFLAGS = -O3 -march=native
#CFLAGS = -march=native
#CFLAGS = -Og -pg
#CFLAGS = -Ofast -pg

CPP_COMPILER = g++

PolyNet : PolyNet.cpp $(OBJECT_FILES) $(HEADER_FILES)
	$(CPP_COMPILER) $(CFLAGS) -o PolyNet PolyNet.cpp $(OBJECT_FILES)

PolyNet_Setup.o : PolyNet_Setup.cpp $(HEADER_FILES)
	$(CPP_COMPILER) $(CFLAGS) -c PolyNet_Setup.cpp

PolyNet_Stepping.o : PolyNet_Stepping.cpp PolyNet.h Point3D.h
	$(CPP_COMPILER) $(CFLAGS) -c PolyNet_Stepping.cpp

PolyNet_Diagnostics.o : PolyNet_Diagnostics.cpp PolyNet.h Point3D.h
	$(CPP_COMPILER) $(CFLAGS) -c PolyNet_Diagnostics.cpp

PolyNet_Output.o : PolyNet_Output.cpp PolyNet.h Point3D.h ReadInputDeck.h
	$(CPP_COMPILER) $(CFLAGS) -c PolyNet_Output.cpp

PolyNet_Misc.o : PolyNet_Misc.cpp
	$(CPP_COMPILER) $(CFLAGS) -c PolyNet_Misc.cpp

ReadInputDeck.o : ReadInputDeck.cpp ReadInputDeck.h
	$(CPP_COMPILER) $(CFLAGS) -c ReadInputDeck.cpp

Point3D.o : Point3D.cpp Point3D.h
	$(CPP_COMPILER) $(CFLAGS) -c Point3D.cpp

clean :
	-rm PolyNet.exe $(OBJECT_FILES)
