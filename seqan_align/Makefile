CXX=g++-5

seqan_align: alignmentline.cpp kmers.cpp semiglobalalignment.cpp seqan_align.cpp
	$(CXX) -I/Users/Ryan/Applications/seqan/include -std=c++11 -shared -Wl,-install_name,seqan_align.so -o seqan_align.so -fPIC alignmentline.cpp kmers.cpp semiglobalalignment.cpp seqan_align.cpp -O3
