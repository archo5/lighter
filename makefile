
CFLAGS = -Wall -static-libgcc -static-libstdc++ -DLTRBUILD

lighter.dll: lighter.cpp lighter.h lighter_int.h
	$(CXX) -o $@ lighter.cpp -shared $(CFLAGS)
