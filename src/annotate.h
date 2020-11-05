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
    
    // Parse DB
    if (!parseDB(c, chrMap)) {
      std::cerr << "Sansa couldn't parse database!" << std::endl;
      return 1;
    }

    // Take care of 1 vs. chr1, X vs chrX, ...
    std::vector<std::string> altName(chrMap.size(), "NA");
    for(TChrMap::iterator itcm = chrMap.begin(); itcm != chrMap.end(); ++itcm) {
      std::string cn = itcm->first;
      if (cn.size() == 1) {
	if (cn == "1") altName[itcm->second] = "chr1";
	else if (cn == "2") altName[itcm->second] = "chr2";
	else if (cn == "3") altName[itcm->second] = "chr3";
	else if (cn == "4") altName[itcm->second] = "chr4";
	else if (cn == "5") altName[itcm->second] = "chr5";
	else if (cn == "6") altName[itcm->second] = "chr6";
	else if (cn == "7") altName[itcm->second] = "chr7";
	else if (cn == "8") altName[itcm->second] = "chr8";
	else if (cn == "9") altName[itcm->second] = "chr9";
	else if (cn == "X") altName[itcm->second] = "chrX";
	else if (cn == "Y") altName[itcm->second] = "chrY";
	else if (cn == "M") altName[itcm->second] = "chrM";
      } else if (cn.size() == 2) {
	if (cn == "10") altName[itcm->second] = "chr10";
	else if (cn == "11") altName[itcm->second] = "chr11";
	else if (cn == "12") altName[itcm->second] = "chr12";
	else if (cn == "13") altName[itcm->second] = "chr13";
	else if (cn == "14") altName[itcm->second] = "chr14";
	else if (cn == "15") altName[itcm->second] = "chr15";
	else if (cn == "16") altName[itcm->second] = "chr16";
	else if (cn == "17") altName[itcm->second] = "chr17";
	else if (cn == "18") altName[itcm->second] = "chr18";
	else if (cn == "19") altName[itcm->second] = "chr19";
	else if (cn == "20") altName[itcm->second] = "chr20";
	else if (cn == "21") altName[itcm->second] = "chr21";
	else if (cn == "22") altName[itcm->second] = "chr22";
	else if (cn == "MT") altName[itcm->second] = "chrM";
      }
    }
    for(uint32_t i = 0; i < altName.size(); ++i) {
      if (altName[i] != "NA") chrMap.insert(std::make_pair(altName[i], i));
    }
    for(TChrMap::const_iterator itcm = chrMap.begin(); itcm != chrMap.end(); ++itcm) {
      std::cerr << itcm->first << '=' << itcm->second << std::endl;
    }
	
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int annotate(int argc, char** argv) {
    AnnotateConfig c;
  
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
