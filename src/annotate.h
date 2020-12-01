#ifndef ANNOTATE_H
#define ANNOTATE_H

#include <fstream>
#include <iomanip>

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
#include <boost/icl/split_interval_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include "parsedb.h"
#include "query.h"
#include "bed.h"
#include "gtf.h"
#include "gff3.h"

namespace sansa
{

  struct AnnotateConfig {
    typedef std::map<std::string, int32_t> TChrMap;
    bool hasCT;
    bool matchSvType;
    bool bestMatch;
    bool reportNoMatch;
    int32_t gtfFileFormat;   // 0 = gtf, 1 = bed, 2 = gff3
    int32_t bpwindow;
    int32_t maxDistance;
    float sizediff;
    std::string idname;
    std::string feature;
    TChrMap nchr;
    boost::filesystem::path gtfFile;
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

    // Unify sequence dictionaries
    int32_t maxRID;
    uint32_t numseq = 0;
    for(uint32_t k = 0; k < 2; ++k) {
      htsFile* ifile = NULL;
      if (k == 0) ifile = bcf_open(c.db.string().c_str(), "r");
      else ifile = bcf_open(c.infile.string().c_str(), "r");
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      int32_t nseq=0;
      const char** seqnames = bcf_hdr_seqnames(hdr, &nseq);
      for(int32_t i = 0; i<nseq;++i) {
	std::string chrName(bcf_hdr_id2name(hdr, i));
	if (c.nchr.find(chrName) == c.nchr.end()) c.nchr[chrName] = numseq++;
      }
      if (seqnames!=NULL) free(seqnames);
      bcf_hdr_destroy(hdr);
      bcf_close(ifile);

      // Fix chrX vs X naming inconsistencies
      maxRID = fixChrNames(c);

      // Debug
      //typedef typename TConfig::TChrMap TChrMap;
      //for(typename TChrMap::const_iterator itcm = c.nchr.begin(); itcm != c.nchr.end(); ++itcm) std::cerr << itcm->first << ',' << itcm->second << std::endl;
    }

    // Structural variants
    typedef std::vector<SV> TSV;
    TSV svs;
    
    // Parse DB
    if (!parseDB(c, svs)) {
      std::cerr << "Sansa couldn't parse database!" << std::endl;
      return 1;
    }

    // Optionally parse GFF/GTF/BED file
    typedef std::vector<IntervalLabel> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    TGenomicRegions gRegions;
    gRegions.resize(maxRID, TChromosomeRegions());
    typedef std::vector<std::string> TGeneIds;
    TGeneIds geneIds;
    if (c.gtfFileFormat != -1) {
      int32_t tf = 0;
      if (c.gtfFileFormat == 0) tf = parseGTF(c, gRegions, geneIds);
      else if (c.gtfFileFormat == 1) tf = parseBED(c, gRegions, geneIds);
      else if (c.gtfFileFormat == 2) tf = parseGFF3(c, gRegions, geneIds);
      if (tf == 0) {
	std::cerr << "Error parsing GTF/GFF3/BED file!" << std::endl;
	return 1;
      }
      for(int32_t refIndex = 0; refIndex < maxRID; ++refIndex) {
	// Sort by position
	std::sort(gRegions[refIndex].begin(), gRegions[refIndex].end(), SortIntervalStart<IntervalLabel>());
      }
    }

    // Query SV
    query(c, svs, gRegions, geneIds);
    
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int annotate(int argc, char** argv) {
    AnnotateConfig c;
    c.hasCT = false;
    std::string strategy = "best";
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("anno,a", boost::program_options::value<boost::filesystem::path>(&c.annofile)->default_value("anno.bcf"), "output annotation VCF/BCF file")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.matchfile)->default_value("query.tsv.gz"), "gzipped output file for query SVs")
      ;


    boost::program_options::options_description svopt("SV annotation file options");
    svopt.add_options()
      ("db,d", boost::program_options::value<boost::filesystem::path>(&c.db), "database VCF/BCF file")
      ("bpoffset,b", boost::program_options::value<int32_t>(&c.bpwindow)->default_value(50), "max. breakpoint offset")
      ("ratio,r", boost::program_options::value<float>(&c.sizediff)->default_value(0.8), "min. size ratio smaller SV to larger SV")
      ("strategy,s", boost::program_options::value<std::string>(&strategy)->default_value("best"), "matching strategy [best|all]")
      ("notype,n", "do not require matching SV types")
      ("nomatch,m", "report SVs without match in database (ANNOID=None)")
      ;
      
    boost::program_options::options_description gtfopt("BED/GTF/GFF3 annotation file options");
    gtfopt.add_options()
      ("gtf,g", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf/gff3/bed file")
      ("id,i", boost::program_options::value<std::string>(&c.idname)->default_value("gene_name"), "gtf/gff3 attribute")
      ("feature,f", boost::program_options::value<std::string>(&c.feature)->default_value("gene"), "gtf/gff3 feature")
      ("distance,t", boost::program_options::value<int32_t>(&c.maxDistance)->default_value(1000), "max. distance (0: overlapping features only)")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "query VCF/BCF file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(svopt).add(gtfopt).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(svopt).add(gtfopt);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: sansa " << argv[0] << " [OPTIONS] input.bcf" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // SV database
    if (!vm.count("db")) {
      // Set input SV file as DB to fill chr array
      c.db = c.infile;
    }
    
    // Match SV types
    if (vm.count("notype")) c.matchSvType = false;
    else c.matchSvType = true;

    // Report no matches
    if (vm.count("nomatch")) c.reportNoMatch = true;
    else c.reportNoMatch = false;
    
    // Check size ratio
    if (c.sizediff < 0) c.sizediff = 0;
    else if (c.sizediff > 1) c.sizediff = 1;

    // Check strategy
    if (strategy == "all") c.bestMatch = false;
    else c.bestMatch = true;

    // Check output directory
    if (!_outfileValid(c.matchfile)) return 1;
    if (!_outfileValid(c.annofile)) return 1;

    // GTF/GFF3/BED
    if (vm.count("gtf")) {
      if (is_gff3(c.gtfFile)) c.gtfFileFormat = 2; // GFF3
      else if (is_gtf(c.gtfFile)) c.gtfFileFormat = 0; // GTF/GFF2
      else c.gtfFileFormat = 1;  // BED
    } else c.gtfFileFormat = -1;

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
