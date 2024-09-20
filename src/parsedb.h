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


  template<typename TConfig, typename TSV>
  inline bool
  parseDB(TConfig& c, TSV& svs) {

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Parse SV annotation database" << std::endl;
    
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

    // Parse VCF records
    bcf1_t* rec = bcf_init();
    int32_t svid = 0;
    int32_t sitecount = 0;
    int32_t lastRID = -1;
    int32_t refIndex = -1;
    while (bcf_read(ifile, hdr, rec) == 0) {
      int32_t startsv = rec->pos + 1;
      bool parsed = true;

      // Count records
      ++sitecount;
      
      // New chromosome?
      if (rec->rid != lastRID) {
	lastRID = rec->rid;
	std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
	refIndex = c.nchr[chrName];
      }

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
      int32_t svlenval = 0;
      _parse_bcf_int32(hdr, rec, "SVLEN", svlenval);

      // Derive proper END and SVLEN
      int32_t endsv = deriveEndPos(rec, svtval, pos2val, endval);
      bool parseALTBND = parseAltBnd(hdr, rec, svtval, ctval, chr2Name, endsv);
      if (!parseALTBND) parsed = false;
      int32_t svlength = deriveSvLength(rec, svtval, endsv, svlenval);
      
      // Numerical SV type
      int32_t svtint = _decodeOrientation(ctval, svtval);
      if (svtint == -1) parsed = false;
      int32_t qualval = 0;
      if (rec->qual > 0) qualval = (int32_t) (rec->qual);

      // Dump record
      //std::cerr << parsed << "\t" << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << chr2Name << "\t" << endsv << "\t" << rec->d.id << "\t" << qualval << "\t" << svtval << "\t" << ctval << "\t" << svtint << "\t" << svlength << std::endl;

      // Store SV
      if (parsed) {
	SV dbsv = SV(refIndex, startsv, c.nchr[chr2Name], endsv, svid, qualval, svtint, svlength);
	_makeCanonical(dbsv);
	svs.push_back(dbsv);
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
    bcf_destroy(rec);
    
    // Sort SVs
    sort(svs.begin(), svs.end());
    
    // Statistics
    now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Parsed " << svid << " out of " << sitecount << " VCF/BCF records." << std::endl;
	
    // Close output VCF
    bcf_hdr_destroy(hdr_out);
    bcf_hdr_destroy(hdr);
    hts_close(ofile);
    bcf_close(ifile);

    // Build BCF index
    bcf_index_build(c.annofile.string().c_str(), 14);
    
    return true;
  }

}

#endif
