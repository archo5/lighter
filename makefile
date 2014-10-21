
CFLAGS = -Wall -static-libgcc -static-libstdc++ -DLTRBUILD

bin/lighter.dll: lighter.cpp lighter.h lighter_int.h
	$(CXX) -o $@ lighter.cpp -shared $(CFLAGS)
