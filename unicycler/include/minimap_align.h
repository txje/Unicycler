
#ifndef MINIMAP_ALIGN_H
#define MINIMAP_ALIGN_H

#include "minimap/minimap.h"
#include "minimap/kseq.h"
#include <string>

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    char * minimapAlignReads(char * referenceFasta, char * readsFastq, int n_threads,
                             int verbosity, int sensitivityLevel);
}

#endif // MINIMAP_ALIGN_H

