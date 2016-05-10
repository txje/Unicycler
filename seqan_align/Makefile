CXX       = g++-5
CXXFLAGS  = -fPIC
CXXFLAGS += -O3
CXXFLAGS += -Wall -Wextra -pedantic
CXXFLAGS += -lpthread -DSEQAN_ENABLE_DEBUG=1
CXXFLAGS += -std=c++11
CXXFLAGS += -Iinclude

TARGET    = seqan_align.so
SHELL     = /bin/sh
SOURCES   = $(shell echo src/*.cpp)
HEADERS   = $(shell echo include/*.h)
OBJECTS   = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -shared -Wl,-install_name,$(TARGET) -o $(TARGET) $(OBJECTS)

clean:
	$(RM) $(TARGET) $(OBJECTS)
