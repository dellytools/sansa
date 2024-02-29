#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "util.h"
#include "version.h"
#include "annotate.h"
#include "compvcf.h"
#include "markdup.h"

using namespace sansa;

inline void
displayUsage() {
  std::cerr << "Usage: sansa <command> <arguments>" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Commands:" << std::endl;
  std::cerr << std::endl;
  std::cerr << "    annotate     annotate VCF file" << std::endl;
  std::cerr << "    markdup      mark duplicate SV sites based on SV allele and GT concordance" << std::endl;
  std::cerr << "    compvcf      compare multi-sample VCF to a ground truth VCF" << std::endl;
  std::cerr << std::endl;
  std::cerr << std::endl;
}

int main(int argc, char **argv) {
  if (argc < 2) { 
    printTitle("Sansa");
    displayUsage();
    return 0;
  }
  
  if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
    std::cerr << "Sansa version: v" << sansaVersionNumber << std::endl;
    std::cerr << " using Boost: v" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
    std::cerr << " using HTSlib: v" << hts_version() << std::endl;
    return 0;
  }
  else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
    printTitle("Sansa");
    displayUsage();
    return 0;
  }
  else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
    displayWarranty();
    return 0;
  }
  else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
    bsd();
    return 0;
  }
  else if ((std::string(argv[1]) == "annotate")) {
    return annotate(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "markdup")) {
    return markdup(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "compvcf")) {
    return compvcf(argc-1,argv+1);
  }
  std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
  return 1;
}

