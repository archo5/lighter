
CFLAGS = -Wall -fno-exceptions -fno-rtti -static-libgcc -static-libstdc++ -DLTRBUILD -g -DLTRDEBUG -m32 -O1
ifeq ($(OS),Windows_NT)
	fnREMOVE_ALL = del /F /S /Q
	fnFIX_PATH = $(subst /,\,$1)
	TARGET_TEST = bin/lighter_test.exe
	TARGET_LIB = bin/lighter.dll
else
	fnREMOVE_ALL = rm -rf
	fnFIX_PATH = $1
	TARGET_TEST = bin/lighter_test
	TARGET_LIB = bin/liblighter.so
	CFLAGS += -lpthread
endif

$(TARGET_TEST): lighter_test.cpp lighter.h $(TARGET_LIB)
	$(CXX) -o $@ lighter_test.cpp $(CFLAGS) $(TARGET_LIB)

$(TARGET_LIB): lighter.cpp lighter_math.cpp lighter.h lighter_int.hpp
	$(CXX) -o $@ lighter.cpp lighter_math.cpp -shared $(CFLAGS)

.PHONY: test clean
test: $(TARGET_TEST)
	$(TARGET_TEST)
clean:
	-$(fnREMOVE_ALL) $(call fnFIX_PATH,$(TARGET_TEST) $(TARGET_LIB))

