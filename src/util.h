#ifndef UTIL_H
#define UTIL_H

#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>

namespace sansa
{

  #ifndef DELLY_SVT_TRANS
  #define DELLY_SVT_TRANS 5
  #endif


  // Structural variant record
  struct StructuralVariantRecord {
    int32_t chr;
    int32_t svStart;
    int32_t chr2;
    int32_t svEnd;
    int32_t svt;
    int32_t id;
    
    StructuralVariantRecord() : chr(0), svStart(0), chr2(0), svEnd(0), svt(-1), id(-1) {}
  };

  
  template<typename TSV>
  struct SortSVs : public std::binary_function<TSV, TSV, bool>
  {
    inline bool operator()(TSV const& sv1, TSV const& sv2) {
      return ((sv1.chr<sv2.chr) || ((sv1.chr==sv2.chr) && (sv1.svStart<sv2.svStart)) || ((sv1.chr==sv2.chr) && (sv1.svStart==sv2.svStart) && (sv1.chr2<sv2.chr2)) || ((sv1.chr==sv2.chr) && (sv1.svStart==sv2.svStart) && (sv1.chr2==sv2.chr2) && (sv1.svEnd<sv2.svEnd)) || ((sv1.chr==sv2.chr) && (sv1.svStart==sv2.svStart) && (sv1.chr2==sv2.chr2) && (sv1.svEnd==sv2.svEnd) && (sv1.id < sv2.id)));
    }
  };

  
  inline int32_t
    _decodeOrientation(std::string const& value) {
    if (value=="3to3") return 0;
    else if (value=="5to5") return 1;
    else if (value=="3to5") return 2;
    else if (value=="5to3") return 3;
    else return 4;
  }
  
  // Decode Orientation
  inline int32_t
  _decodeOrientation(std::string const& value, std::string const& svt) {
    if (svt == "BND") {
      if (value=="3to3") return DELLY_SVT_TRANS + 0;
      else if (value=="5to5") return DELLY_SVT_TRANS + 1;
      else if (value=="3to5") return DELLY_SVT_TRANS + 2;
      else if (value=="5to3") return DELLY_SVT_TRANS + 3;
      else return -1;
    } else {
      if (value=="3to3") return 0;
      else if (value=="5to5") return 1;
      else if (value=="3to5") return 2;
      else if (value=="5to3") return 3;
      else return 4;
    }
  }
  
  // Output directory/file checks
  inline bool
  _outfileValid(boost::filesystem::path& outfile) {
    try {
      boost::filesystem::path outdir;
      if (outfile.has_parent_path()) outdir = outfile.parent_path();
      else outdir = boost::filesystem::current_path();
      if (!boost::filesystem::exists(outdir)) {
	std::cerr << "Output directory does not exist: " << outdir << std::endl;
	return false;
      } else {
	boost::filesystem::file_status s = boost::filesystem::status(outdir);
	boost::filesystem::ofstream file(outfile.string());
	file.close();
	if (!(boost::filesystem::exists(outfile) && boost::filesystem::is_regular_file(outfile))) {
	  std::cerr << "Fail to open output file " << outfile.string() << std::endl;
	  std::cerr << "Output directory permissions: " << s.permissions() << std::endl;
	  return false;
	} else {
	  boost::filesystem::remove(outfile.string());
	}
      }
    } catch (boost::filesystem::filesystem_error const& e) {
      std::cerr << e.what() << std::endl;
      return false;
    }
    return true;
  }

}

#endif
