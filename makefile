
CFLAGS = -Wall -fno-exceptions -fno-rtti -static-libgcc -static-libstdc++ -DLTRBUILD -g -DLTRDEBUG -m32

bin/lighter_test.exe: lighter_test.cpp lighter.h bin/lighter.dll
	$(CXX) -o $@ lighter_test.cpp $(CFLAGS) bin/lighter.dll

bin/lighter.dll: lighter.cpp lighter_math.cpp lighter.h lighter_int.hpp
	$(CXX) -o $@ lighter.cpp lighter_math.cpp -shared $(CFLAGS)

.PHONY: test clean
test: bin/lighter_test.exe
	bin/lighter_test.exe
clean:
	del bin\lighter_test.exe bin\lighter.dll
