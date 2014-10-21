
CFLAGS = -Wall -static-libgcc -static-libstdc++ -DLTRBUILD

bin/lighter.dll: lighter.cpp lighter_math.cpp lighter.h lighter_int.hpp
	$(CXX) -o $@ lighter.cpp lighter_math.cpp -shared $(CFLAGS)
