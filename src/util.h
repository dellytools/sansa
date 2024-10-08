#ifndef UTIL_H
#define UTIL_H

#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>

namespace sansa
{

  #ifndef DELLY_SVT_TRANS
  #define DELLY_SVT_TRANS 5
  #endif

  struct IntervalLabel {
    int32_t start;
    int32_t end;
    char strand;
    int32_t lid;

    explicit IntervalLabel(int32_t s) : start(s), end(s+1), strand('*'), lid(-1) {}
    IntervalLabel(int32_t s, int32_t e, char t, int32_t l) : start(s), end(e), strand(t), lid(l) {}

    bool operator<(const IntervalLabel& s2) const {
      return start < s2.start;
    }
  };

  
  template<typename TRecord>
  struct SortIntervalLabel {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return s1.lid < s2.lid;
    }
  };

  inline void
  _insertInterval(std::vector<IntervalLabel>& cr, int32_t s, int32_t e, char strand, int32_t lid, int32_t) {
    // Uniqueness not necessary because we flatten the interval map
    cr.push_back(IntervalLabel(s, e, strand, lid));
  }
  
  // Structural variant record
  struct SV {
    int32_t chr;
    int32_t svStart;
    int32_t chr2;
    int32_t svEnd;
    int32_t id;
    int32_t qual;
    int32_t svt;
    int32_t svlen;
    
    SV() : chr(0), svStart(0), chr2(0), svEnd(0), id(-1), qual(0), svt(-1), svlen(-1) {}

    SV(int32_t const c1, int32_t const pos1, int32_t const c2, int32_t pos2) : chr(c1), svStart(pos1), chr2(c2), svEnd(pos2), id(-1), qual(0), svt(-1), svlen(-1) {}

    SV(int32_t const c1, int32_t const pos1, int32_t const c2, int32_t pos2, int32_t const ival, int32_t const qval, int32_t const svtval, int32_t const svl) : chr(c1), svStart(pos1), chr2(c2), svEnd(pos2), id(ival), qual(qval), svt(svtval), svlen(svl) {}

    bool operator<(const SV& sv2) const {
      return ((chr<sv2.chr) || ((chr==sv2.chr) && (svStart<sv2.svStart)) || ((chr==sv2.chr) && (svStart==sv2.svStart) && (chr2<sv2.chr2)) || ((chr==sv2.chr) && (svStart==sv2.svStart) && (chr2==sv2.chr2) && (svEnd<sv2.svEnd)) || ((chr==sv2.chr) && (svStart==sv2.svStart) && (chr2==sv2.chr2) && (svEnd==sv2.svEnd) && (id < sv2.id)));
    }
  };

  
  inline bool
  _translocation(int32_t const svt) {
    return (DELLY_SVT_TRANS <= svt) && (svt < 9);
  }
  
  inline bool
  _translocation(bam1_t* rec) {
    return (rec->core.tid != rec->core.mtid);
  }

  inline uint8_t
  _getSpanOrientation(int32_t const svt) {
    if (_translocation(svt)) {
      return svt - DELLY_SVT_TRANS;
    } else {
      return svt;
    }
  }
  
  inline std::string
  _addID(int32_t const svt) {
    if (svt == 0) return "INV";
    else if (svt == 1) return "INV";
    else if (svt == 2) return "DEL";
    else if (svt == 3) return "DUP";
    else if (svt == 4) return "INS";
    else if (svt == 9) return "CNV";
    else return "BND";
  }

  // Add Orientation
  inline std::string
  _addOrientation(int32_t const svt) {
    uint8_t ct = _getSpanOrientation(svt);
    if (ct==0) return "3to3";
    else if (ct==1) return "5to5";
    else if (ct==2) return "3to5";
    else if (ct==3) return "5to3";
    else return "NtoN";
  }

  inline void
  _remove_info_tag(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& tag) {
    bcf_update_info(hdr, rec, tag.c_str(), NULL, 0, BCF_HT_INT);  // Type does not matter for n = 0
  }
  
  inline void
  _remove_format_tag(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& tag) {
    bcf_update_format(hdr, rec, tag.c_str(), NULL, 0, BCF_HT_INT);  // Type does not matter for n = 0
  }

  inline bool
  _isKeyPresent(bcf_hdr_t const* hdr, std::string const& key) {
    return (bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str())>=0);
  }
  
  inline int
  _getInfoType(bcf_hdr_t const* hdr, std::string const& key) {
    return bcf_hdr_id2type(hdr, BCF_HL_INFO, bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str()));
  }

  inline int
  _getFormatType(bcf_hdr_t const* hdr, std::string const& key) {
    return bcf_hdr_id2type(hdr, BCF_HL_FMT, bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str()));
  }
  
  inline bool _missing(bool const value) {
    return !value;
  }
  
  inline bool _missing(float const value) {
    return bcf_float_is_missing(value);
  }
  
  inline bool _missing(int8_t const value) {
    return (value == bcf_int8_missing);
  }
  
  inline bool _missing(int16_t const value) {
    return (value == bcf_int16_missing);
  }
  
  inline bool _missing(int32_t const value) {
    return (value == bcf_int32_missing);
  }
  
  inline bool _missing(std::string const& value) {
    return ((value.empty()) || (value == "."));
  }
  
  inline bool
  _parse_bcf_string(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& key, std::string& strval) {
    if (_isKeyPresent(hdr, key)) {
      char* chr2 = NULL;
      int32_t nchr2 = 0;
      if (bcf_get_info_string(hdr, rec, key.c_str(), &chr2, &nchr2) > 0) strval = std::string(chr2);
      else return false;
      if (chr2 != NULL) free(chr2);
    } else return false;
    return true;
  }

  inline bool
  _parseSVTYPE(bcf_hdr_t* hdr, bcf1_t* rec, std::string& svtval) {
    if (_parse_bcf_string(hdr, rec, "SVTYPE", svtval)) return true;
    else {
      std::string altAllele = rec->d.allele[1];
      if ((altAllele.size() > 2) && (altAllele[0] == '<') && (altAllele[altAllele.size() - 1] == '>')) {
	svtval = altAllele.substr(1,altAllele.size()-2);
      } else return false;
    }
    return true;
  }

  inline bool
  _parse_bcf_int32(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& key, int32_t& intval) {
    if (_isKeyPresent(hdr, key)) {
      int32_t nsvlen = 0;
      int32_t* svlen = NULL;
      if (bcf_get_info_int32(hdr, rec, key.c_str(), &svlen, &nsvlen) > 0) intval = *svlen;
      else return false;
      if (svlen != NULL) free(svlen);
    } else return false;
    return true;
  }


  inline int32_t
  deriveEndPos(bcf1_t* rec, std::string const& svtval, int32_t const pos2val, int32_t const endval) {
    if ((pos2val != -1) && (endval != -1)) {
      if ((svtval == "BND") || (svtval == "TRA")) return pos2val;
      else return endval;
    }
    else if (pos2val != -1) return pos2val;
    else if (endval != -1) return endval;
    else {
      if (svtval == "INS") return rec->pos + 2;
      if (svtval == "DEL") {
	std::string refAllele = rec->d.allele[0];
	std::string altAllele = rec->d.allele[1];
	if ((altAllele.size()) && (refAllele.size() > altAllele.size()) && (altAllele[0] != '<')) {
	  return rec->pos + 1 + (refAllele.size() - altAllele.size());
	}
      }
    }
    return -1;
  }

  inline bool
  parseAltBnd(bcf_hdr_t* hdr, bcf1_t* rec, std::string& svtval, std::string& ctval, std::string& chr2Name, int32_t& endval) {
    if (endval != -1) return true;  // Done.
    if (svtval != "BND") return false;
    std::string altAllele = rec->d.allele[1];
    std::size_t found1 = altAllele.find_first_of("[");
    std::size_t found2 = altAllele.find_last_of("[");
    if ((found1 == std::string::npos) && (found2 == std::string::npos)) {
      found1 = altAllele.find_first_of("]");
      found2 = altAllele.find_last_of("]");
    }
    if ((found1 == std::string::npos) || (found2 == std::string::npos)) return false;
    if (found1 + 1 >= found2) return false;
    altAllele = altAllele.substr(found1 + 1, found2 - found1 - 1);
    std::size_t sepa = altAllele.find_first_of(":");
    if (sepa == std::string::npos) return false;
    chr2Name = altAllele.substr(0, sepa);
    endval = boost::lexical_cast<int32_t>(altAllele.substr(sepa+1));
    std::string chr1Name = bcf_hdr_id2name(hdr, rec->rid);
    if (chr1Name == chr2Name) {
      // Fix SV type
      int32_t pos = rec->pos + 1;
      if (pos > endval) return false;
      std::string svclass = "NA";
      if (_parse_bcf_string(hdr, rec, "SVCLASS", svclass)) {
	if (svclass == "DEL") {
	  svtval = "DEL";
	  ctval = "3to5";
	} else if (svclass == "DUP") {
	  svtval = "DUP";
	  ctval = "5to3";
	} else if (svclass == "h2hINV") {
	  svtval = "INV";
	  ctval = "3to3";
	} else if (svclass == "t2tINV") {
	  svtval = "INV";
	  ctval = "5to5";
	} else if (svclass == "INS") {
	  svtval = "INS";
	  ctval = "NtoN";
	  endval = rec->pos + 2;
	}
      }
      if (svtval == "BND") {
	// Still not fixed, try CT
	if (ctval == "3to5") svtval = "DEL";
	else if (ctval == "5to3") svtval = "DUP";
	else if (ctval == "3to3") svtval = "INV";
	else if (ctval == "5to5") svtval = "INV";
	else if (ctval == "NtoN") svtval = "INS";
	else return false;
      }	
      //std::cerr << std::string(rec->d.allele[1]) << ',' << chr1Name << ',' << pos << ',' << chr2Name << ',' << endval << ',' << svtval << ',' << ctval << ',' << svclass << std::endl;
    }
    return true;
  }

  inline int32_t
  deriveSvLength(bcf1_t* rec, std::string const& svtval, int32_t const endval, int32_t const svlenval) {
    if (svlenval != 0) {
      if (svlenval > 0) return svlenval;
      else return -svlenval;
    }
    else {
      if ((svtval == "DEL") || (svtval == "DUP") || (svtval == "INV")) {
	if ((endval != -1) && (endval > rec->pos)) {
	  return endval - (rec->pos + 1);
	}
      }
    }
    return svlenval;
  }

  inline double
  _sharedCarriers(std::vector<int32_t> const& gt1, std::vector<int32_t> const& gt2) {
    // Percentage of shared carriers
    uint32_t carnum = 0;
    uint32_t carshared = 0;
    for(uint32_t k = 0; k < gt1.size(); ++k) {
      if ((gt1[k] != 0) || (gt2[k] != 0)) {
	++carnum;
	if ((gt1[k] != 0) && (gt2[k] != 0)) ++carshared;
      }
    }
    return (double) (carshared) / (double) (carnum);
  }

  inline double
  nonrefGtConc(std::vector<int32_t> const& gt1, std::vector<int32_t> const& gt2) {
    uint32_t totgt = 0;
    uint32_t matchgt = 0;
    for(uint32_t k = 0; k < gt1.size(); ++k) {
      if ((gt1[k] != 0) || (gt2[k] != 0)) {
	++totgt;
	if (gt1[k] == gt2[k]) ++matchgt;
      }
    }
    return (double) (matchgt) / (double) (totgt);
  }
  
  inline double
  gtConc(std::vector<int32_t> const& gt1, std::vector<int32_t> const& gt2) {
    uint32_t totgt = 0;
    uint32_t matchgt = 0;
    for(uint32_t k = 0; k < gt1.size(); ++k) {
      ++totgt;
      if (gt1[k] == gt2[k]) ++matchgt;
    }
    return (double) (matchgt) / (double) (totgt);
  }
  
  inline void
  _makeCanonical(SV& sv) {
    if (sv.chr != sv.chr2) {
      if (sv.chr < sv.chr2) {
	int32_t tmpChr = sv.chr;
	sv.chr = sv.chr2;
	sv.chr2 = tmpChr;
	int32_t tmpPos = sv.svStart;
	sv.svStart = sv.svEnd;
	sv.svEnd = tmpPos;
	if (sv.svt == DELLY_SVT_TRANS + 2) sv.svt = DELLY_SVT_TRANS + 3;
	else if (sv.svt == DELLY_SVT_TRANS + 3) sv.svt = DELLY_SVT_TRANS + 2;
      }
    }
  }

  
  inline bool
  is_gff3(boost::filesystem::path const& f) {
    std::ifstream in(f.string().c_str());
    if (!in) return false;
    in.close();

    // Check file ending
    if (f.string().size()>=6) {
      int32_t strn = f.string().size();
      if (f.string().substr(strn - 6, 6) == "gff.gz") return true;
    }
    if (f.string().size()>=7) {
      int32_t strn = f.string().size();
      if (f.string().substr(strn - 7, 7) == "gff3.gz") return true;
    }

    // Check file content
    std::ifstream file(f.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::getline(instream, gline);
    bool gff = false;
    if ((gline.size()>=5) && (gline.substr(0,5) == "##gff")) gff = true;
    file.close();
    return gff;
  }

  inline bool
  is_gtf(boost::filesystem::path const& f) {
    std::ifstream in(f.string().c_str());
    if (!in) return false;
    in.close();

    // Check file ending
    if (f.string().size()>=6) {
      int32_t strn = f.string().size();
      if (f.string().substr(strn - 6, 6) == "gtf.gz") return true;
    }

    // Check file content
    std::ifstream file(f.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::getline(instream, gline);
    bool gtf = false;
    if ((gline.size()>=2) && (gline.substr(0,2) == "#!")) gtf = true;
    file.close();
    return gtf;
  }

  inline bool
  is_gz(boost::filesystem::path const& f) {
    std::ifstream bfile(f.string().c_str(), std::ios_base::binary | std::ios::ate);
    bfile.seekg(0, std::ios::beg);
    char byte1;
    bfile.read(&byte1, 1);
    char byte2;
    bfile.read(&byte2, 1);
    bfile.close();
    if ((byte1 == '\x1F') && (byte2 == '\x8B')) return true;
    else return false;
  }

  inline std::string
  _translateSvType(int32_t const svt) {
    if (svt == 0) return "INV";
    else if (svt == 1) return "INV";
    else if (svt == 2) return "DEL";
    else if (svt == 3) return "DUP";
    else if (svt == 4) return "INS";
    else if (svt == 9) return "CNV";
    else if (svt == 10) return "CPX";
    else if (svt == 11) return "CTX";
    else if (svt == 12) return "SCNA";
    else if (svt == 13) return "MCNV";
    else return "BND";
  }

  inline std::string
  _translateCt(int32_t const svt) {
    if (svt == 0) return "3to3";
    else if (svt == 1) return "5to5";
    else if (svt == 2) return "3to5";
    else if (svt == 3) return "5to3";
    else if (svt == 4) return "NtoN";
    else if (svt == 5) return "3to3";
    else if (svt == 6) return "5to5";
    else if (svt == 7) return "3to5";
    else if (svt == 8) return "5to3";
    else return "NA";
  }    
  
  // Decode Orientation
  inline int32_t
  _decodeOrientation(std::string const& value, std::string const& svt) {
    if ((svt == "BND") || (svt == "TRA")) {
      if (value=="NA") return DELLY_SVT_TRANS + 0;
      else if (value=="3to3") return DELLY_SVT_TRANS + 0;
      else if (value=="5to5") return DELLY_SVT_TRANS + 1;
      else if (value=="3to5") return DELLY_SVT_TRANS + 2;
      else if (value=="5to3") return DELLY_SVT_TRANS + 3;
      else return -1;
    } else {
      if (value=="NA") {
	if (svt == "DEL") return 2;
	else if (svt == "INV") return 0;
	else if (svt == "DUP") return 3;
	else if (svt == "INS") return 4;
	else if (svt == "CNV") return 9;  
	else if (svt == "CPX") return 10;
	else if (svt == "CTX") return 11;
	else if (svt == "SCNA") return 12;
	else if (svt == "MCNV") return 13;
	else return -1;
      }
      else if (value=="3to3") return 0;
      else if (value=="5to5") return 1;
      else if (value=="3to5") return 2;
      else if (value=="5to3") return 3;
      else return 4;
    }
  }

  template<typename TChrMap>
  inline int32_t
  chrMapSize(TChrMap const& chrMap) {
    int32_t maxRID = 0;
    for(typename TChrMap::const_iterator itcm = chrMap.begin(); itcm != chrMap.end(); ++itcm) {
      if (itcm->second >= maxRID) maxRID = itcm->second + 1;
    }
    return maxRID;
  }
    
  
  template<typename TConfig>
  inline int32_t
  fixChrNames(TConfig& c) {
    typedef typename TConfig::TChrMap TChrMap;
    
    // Take care of 1 vs. chr1, X vs chrX, ...
    int32_t maxRID = chrMapSize(c.nchr);
    std::vector<std::string> altName(maxRID);
    for(typename TChrMap::const_iterator itcm = c.nchr.begin(); itcm != c.nchr.end(); ++itcm) {
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
	else if (cn == "MT") altName[itcm->second] = "chrMT";
      } else if (cn.size() == 4) {
	if (cn == "chr1") altName[itcm->second] = "1";
	else if (cn == "chr2") altName[itcm->second] = "2";
	else if (cn == "chr3") altName[itcm->second] = "3";
	else if (cn == "chr4") altName[itcm->second] = "4";
	else if (cn == "chr5") altName[itcm->second] = "5";
	else if (cn == "chr6") altName[itcm->second] = "6";
	else if (cn == "chr7") altName[itcm->second] = "7";
	else if (cn == "chr8") altName[itcm->second] = "8";
	else if (cn == "chr9") altName[itcm->second] = "9";
	else if (cn == "chrX") altName[itcm->second] = "X";
	else if (cn == "chrY") altName[itcm->second] = "Y";
	else if (cn == "chrM") altName[itcm->second] = "M";
      } else if (cn.size() == 5) {
	if (cn == "chr10") altName[itcm->second] = "10";
	else if (cn == "chr11") altName[itcm->second] = "11";
	else if (cn == "chr12") altName[itcm->second] = "12";
	else if (cn == "chr13") altName[itcm->second] = "13";
	else if (cn == "chr14") altName[itcm->second] = "14";
	else if (cn == "chr15") altName[itcm->second] = "15";
	else if (cn == "chr16") altName[itcm->second] = "16";
	else if (cn == "chr17") altName[itcm->second] = "17";
	else if (cn == "chr18") altName[itcm->second] = "18";
	else if (cn == "chr19") altName[itcm->second] = "19";
	else if (cn == "chr20") altName[itcm->second] = "20";
	else if (cn == "chr21") altName[itcm->second] = "21";
	else if (cn == "chr22") altName[itcm->second] = "22";
	else if (cn == "chrMT") altName[itcm->second] = "MT";
      }
    }

    // Insert alternative names
    for(uint32_t i = 0; i < altName.size(); ++i) {
      if (!altName[i].empty()) c.nchr.insert(std::make_pair(altName[i], i));
    }
    return maxRID;
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
