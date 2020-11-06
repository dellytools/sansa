#ifndef ANNOTATE_H
#define ANNOTATE_H

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include "parsedb.h"
#include "query.h"

namespace sansa
{

  struct AnnotateConfig {
    bool hasCT;
    bool matchSvType;
    int32_t bpwindow;
    float sizediff;
    boost::filesystem::path annofile;
    boost::filesystem::path db;
    boost::filesystem::path matchfile;
    boost::filesystem::path infile;
  };


  template<typename TConfig>
  inline int32_t
  runAnnotate(TConfig& c) {
    
#ifdef PROFILE
    ProfilerStart("sansa.prof");
#endif

    // Chromosome map
    typedef std::map<std::string, int32_t> TChrMap;
    TChrMap chrMap;

    // Structural variants
    typedef std::vector<SV> TSV;
    TSV svs;
    
    // Parse DB
    if (!parseDB(c, svs, chrMap)) {
      std::cerr << "Sansa couldn't parse database!" << std::endl;
      return 1;
    }

    // Fix chr names
    fixChrNames(chrMap);

    // Debug chr names
    //for(TChrMap::const_iterator itcm = chrMap.begin(); itcm != chrMap.end(); ++itcm) std::cerr << itcm->first << '=' << itcm->second << std::endl;

    // Query SV
    query(c, svs, chrMap);
    
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int annotate(int argc, char** argv) {
    AnnotateConfig c;
    c.hasCT = false;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("db,d", boost::program_options::value<boost::filesystem::path>(&c.db)->default_value("database.bcf"), "database VCF/BCF file")
      ("bpoffset,b", boost::program_options::value<int32_t>(&c.bpwindow)->default_value(50), "max. breakpoint offset")
      ("sizediff,s", boost::program_options::value<float>(&c.sizediff)->default_value(0.8), "min. size ratio smaller SV to larger SV")
      ("anno,a", boost::program_options::value<boost::filesystem::path>(&c.annofile)->default_value("anno.bcf"), "output annotation VCF/BCF file")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.matchfile)->default_value("query.txt.gz"), "gzipped output file for query SVs")
      ("notype,n", "Do not require matching SV types")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "query VCF/BCF file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: sansa " << argv[0] << " [OPTIONS] input.bcf" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Match SV types
    if (vm.count("notype")) c.matchSvType = false;
    else c.matchSvType = true;

    // Check size ratio
    if (c.sizediff < 0) c.sizediff = 0;
    else if (c.sizediff > 1) c.sizediff = 1;
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "sansa ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runAnnotate(c);
  }

}

#endif
