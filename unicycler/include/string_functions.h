#ifndef STRING_CONVERSION_H
#define STRING_CONVERSION_H

#include <string>
#include <vector>

char * cppStringToCString(std::string cpp_string);

std::vector<std::string> splitString(char * inString, char delimiter);
std::vector<std::string> splitString(std::string inString, char delimiter);

std::string getReverseComplement(std::string sequence);


// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {
    void freeCString(char * p);
}


#endif // STRING_CONVERSION_H
