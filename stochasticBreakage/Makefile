COMPILER       = g++
LINKER         = g++

COMPILE_FLAGS  = -O2 -Wall -fPIC\

LINK_FLAGS     = 
LIBS           = 

EXEC = testStoB
SRC  = $(wildcard *.cpp)
OBJ  = $(SRC:.cpp=.o)


exec : $(OBJ)
	$(LINKER) $(LINK_FLAGS) $(OBJ) $(LIBS) -o $(EXEC)
	@size $(EXEC)

#Meaning: for every .o file, the corresponding .cpp file is a prerequisite
#Recipe for the cpp file is to compile it with the compile flags, and produce a .o file
%.o : %.cpp
	$(COMPILER) $(COMPILE_FLAGS) -c -MMD $<

clean : 
	rm *.o *.d $(EXEC)

-include *.d
