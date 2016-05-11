SHELL        = /bin/sh
CXX          = g++-5

FLAGS        = -std=c++11 -Iinclude -fPIC
CXXFLAGS     = -Wall -Wextra -pedantic -march=native
LDFLAGS      = -shared
DEBUGFLAGS   = -DSEQAN_ENABLE_DEBUG=1
RELEASEFLAGS = -O3 -D NDEBUG

TARGET       = seqan_align.so
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
