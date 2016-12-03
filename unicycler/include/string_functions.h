#ifndef STRING_CONVERSION_H
#define STRING_CONVERSION_H

#include <string>

char * cppStringToCString(std::string cpp_string);

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {
    void freeCString(char * p);
}


#endif // STRING_CONVERSION_H
