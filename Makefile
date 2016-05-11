# This makefile will build the C++ components of the semi-global long read aligner.

# Example commands:
#   make (build in release mode)
#   make debug (build in debug mode)
#   make clean (deletes *.o files, which aren't required to run the aligner)
#   make distclean (deletes *.o files and the *.so file, which is required to run the aligner)
#   make CXX=g++-5 (build with a particular compiler)
#   make CXXFLAGS="-Werror -g3" (build with particular compiler flags)


# CXX and CXXFLAGS can be overriden by the user.
CXX         ?= g++-5
CXXFLAGS    ?= -Wall -Wextra -pedantic -march=native

# These flags are required for the build to work.
FLAGS        = -std=c++11 -Iinclude -fPIC
LDFLAGS      = -shared

# Different debug/optimisation levels for debug/release builds.
DEBUGFLAGS   = -DSEQAN_ENABLE_DEBUG=1 -g
RELEASEFLAGS = -O3 -D NDEBUG

TARGET       = seqan_align.so
SHELL        = /bin/sh
SOURCES      = $(shell echo src/*.cpp)
HEADERS      = $(shell echo include/*.h)
OBJECTS      = $(SOURCES:.cpp=.o)

.PHONY: release
release: FLAGS+=$(RELEASEFLAGS)
release: $(TARGET)

.PHONY: debug
debug: FLAGS+=$(DEBUGFLAGS)
debug: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(LDFLAGS) -Wl,-install_name,$(TARGET) -o $(TARGET) $(OBJECTS)

clean:
	$(RM) $(OBJECTS)

distclean: clean
	$(RM) $(TARGET)

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) $(CXXFLAGS) -c -o $@ $<
