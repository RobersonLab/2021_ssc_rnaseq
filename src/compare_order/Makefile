CXX = g++

# Any custom compiler flags here ( -g, -march, -mtune, -m64, etc)
# For 64-bit compilation add -m64 to CXXFLAGS
CXXFLAGS = 
# --------------------------------------------------------------------
# Do not edit below this line
# --------------------------------------------------------------------

CXXFLAGS += -O3 -I. -Wall
OUTPUT = compare_order
SRC = $(addsuffix .cpp, $(OUTPUT))
OBJ = $(addsuffix .o, $(OUTPUT))

all : $(OUTPUT) 

$(OUTPUT) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

.PHONEY : clean install uninstall

clean:
	rm -f *.o $(OUTPUT)

install:
	install -m 0755 ./$(OUTPUT) /usr/local/bin

uninstall:
	rm /usr/local/bin/$(OUTPUT)
