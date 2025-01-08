#include <seqan3/core/debug_stream.hpp>
#include <htslib/hts.h> 
int main()
{
    seqan3::debug_stream << "Hello World!\n";
    return 0;
}
