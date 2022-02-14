OPT = -O3 
CFLAGS = -std=c++11

all:
	g++ $(CFLAGS) $(OPT) main.cpp -o main.exe

clean:
	rm -f *.exe
