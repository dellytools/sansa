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


  template<typename TConfig, typename TSV, typename TMap>
  inline bool
  parseDB(TConfig& c, TSV& svs, TMap& chrMap) {

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Parse annotation database" << std::endl;
    
    // Load bcf file
    htsFile* ifile = bcf_open(c.db.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to load " << c.db.string() << std::endl;
      return false;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    
    // Open output VCF file
    htsFile *ofile = hts_open(c.annofile.string().c_str(), "wb");
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "ANNOID");
    bcf_hdr_append(hdr_out, "##INFO=<ID=ANNOID,Number=1,Type=String,Description=\"Annotation ID that links query SVs to database SVs.\">");
    if (bcf_hdr_write(ofile, hdr_out) != 0) {
      std::cerr << "Error: Failed to write BCF header!" << std::endl;
      return false;
    }

    // Temporary chr2 map until we have seen all chromosomes
    typedef std::map<std::string, int32_t> TChr2Map;
    TChr2Map chr2Map;
    
    // Parse VCF records
    bcf1_t* rec = bcf_init();
    int32_t svid = 0;
    int32_t sitecount = 0;
    int32_t lastRID = -1;
    while (bcf_read(ifile, hdr, rec) == 0) {
      int32_t startsv = rec->pos + 1;
      bool parsed = true;

      // Count records
      ++sitecount;
      
      // Augment chromosome map
      if (rec->rid != lastRID) {
	lastRID = rec->rid;
	std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
	if (chrMap.find(chrName) == chrMap.end()) chrMap.insert(std::make_pair(chrName, rec->rid));
      }

      // Only bi-allelic
      if (rec->n_allele != 2) parsed = false;

      // Unpack INFO
      bcf_unpack(rec, BCF_UN_INFO);

      // Parse INFO fields
      std::string svtval = "NA";
      if (!_parseSVTYPE(hdr, rec, svtval)) parsed = false;
      std::string ctval("NA");
      if (_parse_bcf_string(hdr, rec, "CT", ctval)) c.hasCT = true;
      std::string chr2Name(bcf_hdr_id2name(hdr, rec->rid));
      _parse_bcf_string(hdr, rec, "CHR2", chr2Name);      
      int32_t pos2val = -1;
      _parse_bcf_int32(hdr, rec, "POS2", pos2val);
      int32_t endval = -1;
      _parse_bcf_int32(hdr, rec, "END", endval);
      int32_t svlenval = -1;
      _parse_bcf_int32(hdr, rec, "SVLEN", svlenval);

      // Derive proper END and SVLEN
      int32_t endsv = deriveEndPos(rec, svtval, pos2val, endval);
      int32_t svlength = deriveSvLength(rec, svtval, endval, svlenval);
      
      // Numerical SV type
      int32_t svtint = _decodeOrientation(ctval, svtval);
      if (svtint == -1) parsed = false;
      int32_t qualval = (int32_t) (rec->qual);

      // Dump record
      //if (!parsed) std::cerr << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << chr2Name << "\t" << endsv << "\t" << rec->d.id << "\t" << qualval << "\t" << svtval << "\t" << ctval << "\t" << svtint << "\t" << svlength << std::endl;

      // Store SV
      if (parsed) {
	if (chr2Map.find(chr2Name) == chr2Map.end()) chr2Map.insert(std::make_pair(chr2Name, chr2Map.size()));
	svs.push_back(SV(rec->rid, startsv, chr2Map[chr2Name], endsv, svid, qualval, svtint, svlength));
	std::string id("id");
	std::string padNumber = boost::lexical_cast<std::string>(svid);
	padNumber.insert(padNumber.begin(), 9 - padNumber.length(), '0');
	id += padNumber;
	_remove_info_tag(hdr_out, rec, "ANNOID");
	bcf_update_info_string(hdr_out, rec, "ANNOID", id.c_str());
	bcf_write1(ofile, hdr_out, rec);
	++svid;
      }
    }

    // Remap chr names and ensure canonical order for translocations
    typedef std::map<int32_t, int32_t> TChrIdMap;
    TChrIdMap chrIdMap;
    for(typename TChr2Map::iterator itc2 = chr2Map.begin(); itc2 != chr2Map.end(); ++itc2) chrIdMap.insert(std::make_pair(itc2->second, chrMap[itc2->first]));
    for(uint32_t i = 0; i < svs.size(); ++i) {
      svs[i].chr2 = chrIdMap[svs[i].chr2];
      if (svs[i].chr != svs[i].chr2) {
	if (svs[i].chr < svs[i].chr2) {
	  int32_t tmpChr = svs[i].chr;
	  svs[i].chr = svs[i].chr2;
	  svs[i].chr2 = tmpChr;
	  int32_t tmpPos = svs[i].svStart;
	  svs[i].svStart = svs[i].svEnd;
	  svs[i].svEnd = tmpPos;
	}
      }
    }

    // Sort SVs
    sort(svs.begin(), svs.end(), SortSVs<SV>());
    
    // Statistics
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Parsed " << svid << " out of " << sitecount << " VCF/BCF records." << std::endl;
	
    // Close output VCF
    bcf_hdr_destroy(hdr_out);
    hts_close(ofile);
    
    // Close file handles
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    bcf_destroy(rec);
    
    return true;
  }

}

#endif
