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


namespace sansa
{

  struct AnnotateConfig {
    bool hasDumpFile;
    bool hasCT;
    boost::filesystem::path dumpfile;
    boost::filesystem::path db;
    boost::filesystem::path outfile;
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
    for(TChrMap::const_iterator itcm = chrMap.begin(); itcm != chrMap.end(); ++itcm) std::cerr << itcm->first << '=' << itcm->second << std::endl;
	
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
      ("anno,a", boost::program_options::value<boost::filesystem::path>(&c.db)->default_value("annotation.bcf"), "annotation database VCF/BCF file")
      ("dump,d", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped output file for DB SVs")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("test.tsv"), "output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input VCF/BCF file")
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

    // Dump DB
    if (vm.count("dump")) c.hasDumpFile = true;
    else c.hasDumpFile = false;
    
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
