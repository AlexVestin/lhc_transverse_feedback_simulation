
BUILD_PATH=./dependencies/build/

# FFTW and sciplot
INC=-I$(BUILD_PATH)include
LIB=-L$(BUILD_PATH)lib
FLAGS=-lhdf5 -lfftw3 -lm 
# Skia

SK_INC=-I/home/alex/skia/
SK_LIB=-L/home/alex/skia/out/Static
SK_FLAGS=-lskia -ldl -lpthread -ljpeg -lfreetype -lz -lpng -lglfw -lfontconfig -lwebp -lwebpmux -lwebpdemux -lGL

DB_FLAGS=-Wall -fsanitize=address -g


NSRC=build/NAFFlib/source
NAFFLIB=$(NSRC)/fft.o $(NSRC)/brent.o $(NSRC)/frequency.o $(NSRC)/signal_processing.o $(NSRC)/windows.o

skia: 
	g++ -std=c++17 -DSKIA main.cpp $(NAFFLIB)  $(DB_FLAGS)  $(INC) $(SK_INC) $(LIB) $(SK_LIB) $(SK_FLAGS) $(FLAGS) -o main

sciplot:
	g++ -std=c++17 -DSCIPLOT main.cpp $(DB_FLAGS) $(INC) $(LIB) $(FLAGS) -o main

