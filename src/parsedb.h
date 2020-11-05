#ifndef PARSEDB_H
#define PARSEDB_H

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
  

  template<typename TConfig, typename TMap>
  inline bool
  parseDB(TConfig const& c, TMap& chrMap) {

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Parse annotation database" << std::endl;
    
    // Load bcf file
    htsFile* ifile = bcf_open(c.db.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to load " << c.db.string() << std::endl;
      return false;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    
    // Read SV information
    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t nsvlen = 0;
    int32_t* svlen = NULL;
    int32_t npos2 = 0;
    int32_t* pos2 = NULL;
    int32_t nsvt = 0;
    char* svt = NULL;
    int32_t nchr2 = 0;
    char* chr2 = NULL;
    int32_t nct = 0;
    char* ct = NULL;
    
    // Get the valid columns
    uint32_t fieldIndex = 0;
    typedef std::map<std::string, uint32_t> TColumnMap;
    TColumnMap cMap;
    cMap["chr"] = fieldIndex++;
    cMap["start"] = fieldIndex++;
    cMap["chr2"] = fieldIndex++;
    cMap["end"] = fieldIndex++;
    cMap["id"] = fieldIndex++;
    cMap["qual"] = fieldIndex++;
    cMap["pass"] = fieldIndex++;
    if (_isKeyPresent(hdr, "SVTYPE")) cMap["svtype"] = fieldIndex++;
    if (_isKeyPresent(hdr, "CT")) cMap["ct"] = fieldIndex++;
    if (_isKeyPresent(hdr, "SVLEN")) cMap["svlen"] = fieldIndex++;

    // Dump file
    boost::iostreams::filtering_ostream dumpOut;
    if (c.hasDumpFile) {
      dumpOut.push(boost::iostreams::gzip_compressor());
      dumpOut.push(boost::iostreams::file_sink(c.dumpfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
      dumpOut << "chr\tstart\tchr2\tend\tid\tqual\tpass\tsvtype\tct\tsvlen" << std::endl;
    }

    // Parse VCF records
    bcf1_t* rec = bcf_init();
    uint32_t siteCount = 0;
    int32_t lastRID = -1;
    while (bcf_read(ifile, hdr, rec) == 0) {
      // Augment chromosome map
      if (rec->rid != lastRID) {
	lastRID = rec->rid;
	std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
	if (chrMap.find(chrName) == chrMap.end()) chrMap.insert(std::make_pair(chrName, rec->rid));
      }

      // Only bi-allelic
      if (rec->n_allele != 2) continue;

      // Defaults
      bool endPresent = false;
      bool pos2Present = false;
      int32_t endsv = -1;
      int32_t svlength = -1;
      
      bcf_unpack(rec, BCF_UN_INFO);
      if (_isKeyPresent(hdr, "SVTYPE")) {
	if (bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt) > 0) {
	  // Nop
	} else continue;
      }

      ++siteCount;
      std::string ctval("NA");
      if (_isKeyPresent(hdr, "CT")) {
	if (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0) ctval = std::string(ct);
      }
      bool passSite = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);
      std::string chr2Name(bcf_hdr_id2name(hdr, rec->rid));
      if (_isKeyPresent(hdr, "CHR2")) {
	if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) chr2Name = std::string(chr2);
      }
      if (_isKeyPresent(hdr, "POS2")) {
	if (bcf_get_info_int32(hdr, rec, "POS2", &pos2, &npos2) > 0) {
	  pos2Present = true;
	}
      }
      if (_isKeyPresent(hdr, "SVLEN")) {
	if (bcf_get_info_int32(hdr, rec, "SVLEN", &svlen, &nsvlen) > 0) {
	  svlength = *svlen;
	}
      }
      if (_isKeyPresent(hdr, "END")) {
	if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) {
	  endPresent = true;
	}
      }
      if ((pos2Present) && (endPresent)) {
	if (std::string(svt) == "BND") {
	  endsv = *pos2;
	} else {
	  endsv = *svend;
	}
      }
      else if (pos2Present) endsv = *pos2;
      else if (endPresent) endsv = *svend;
      else {
	std::string refAllele = rec->d.allele[0];
	std::string altAllele = rec->d.allele[1];
	int32_t diff = refAllele.size() - altAllele.size();
	endsv = rec->pos + diff + 2;
      }

      if (c.hasDumpFile) {
	dumpOut << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << chr2Name << "\t" << endsv << "\t" << rec->d.id << "\t" << rec->qual << "\t" << passSite << "\t" << svt << "\t" << ctval << "\t" << svlength << std::endl;
      }
    }
  
    // Clean-up
    if (svend != NULL) free(svend);
    if (svlen != NULL) free(svlen);
    if (pos2 != NULL) free(pos2);
    if (svt != NULL) free(svt);
    if (chr2 != NULL) free(chr2);
    if (ct != NULL) free(ct);
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    bcf_destroy(rec);
    
    return true;
  }

}

#endif
