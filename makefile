
CFLAGS = -Wall -fno-exceptions -static-libgcc -static-libstdc++ -DLTRBUILD -g

bin/lighter_test.exe: lighter_test.cpp lighter.h bin/lighter.dll
	$(CXX) -o $@ lighter_test.cpp $(CFLAGS) bin/lighter.dll

bin/lighter.dll: lighter.cpp lighter_math.cpp lighter.h lighter_int.hpp
	$(CXX) -o $@ lighter.cpp lighter_math.cpp -shared $(CFLAGS)

.PHONY: test
test: bin/lighter_test.exe
	bin/lighter_test.exe
