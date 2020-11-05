#ifndef QUERY_H
#define QUERY_H

#include <boost/dynamic_bitset.hpp>
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
  query(TConfig& c, TSV& svs, TMap& chrMap) {

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Query input SVs" << std::endl;
    
    // Load bcf file
    htsFile* ifile = bcf_open(c.infile.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to load " << c.infile.string() << std::endl;
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

    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.matchfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dataOut << "anno.id\tquery.chr\tquery.start\tquery.chr2\tquery.end\tquery.id\tquery.qual\tquery.svtype\tquery.svlen" << std::endl;
    
    // Parse VCF records
    bcf1_t* rec = bcf_init();
    int32_t lastRID = -1;
    int32_t sitecount = 0;
    int32_t parsedSV = 0;
    int32_t refIndex = -1;
    int32_t refIndex2 = -1;
    while (bcf_read(ifile, hdr, rec) == 0) {
      ++sitecount;
      
      // Defaults
      std::string svtval = "NA";
      bool endPresent = false;
      bool pos2Present = false;
      int32_t endsv = -1;
      int32_t startsv = rec->pos + 1;
      int32_t svlength = -1;

      // New chromosome?
      if (rec->rid != lastRID) {
	lastRID = rec->rid;
	std::string chrName = bcf_hdr_id2name(hdr, rec->rid);	
	if (chrMap.find(chrName) == chrMap.end()) refIndex = -1;
	else refIndex = chrMap[chrName];
      }
      if (refIndex == -1) continue;
      
      // Only bi-allelic
      if (rec->n_allele != 2) continue;

      // Unpack INFO
      bcf_unpack(rec, BCF_UN_INFO);

      // SVTYPE
      if (_isKeyPresent(hdr, "SVTYPE")) {
	if (bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt) > 0) {
	  svtval = std::string(svt);
	} else continue;
      }

      // CT
      std::string ctval("NA");
      if (_isKeyPresent(hdr, "CT")) {
	if (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0) {
	  ctval = std::string(ct);
	}
      }

      // CHR2
      std::string chr2Name(bcf_hdr_id2name(hdr, rec->rid));
      refIndex2 = refIndex;
      if (_isKeyPresent(hdr, "CHR2")) {
	if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) {
	  chr2Name = std::string(chr2);
	  if (chrMap.find(chr2Name) == chrMap.end()) continue;
	  else refIndex2 = chrMap[chr2Name];
	}
      }

      // POS2
      if (_isKeyPresent(hdr, "POS2")) {
	if (bcf_get_info_int32(hdr, rec, "POS2", &pos2, &npos2) > 0) {
	  pos2Present = true;
	}
      }

      // SVLEN
      if (_isKeyPresent(hdr, "SVLEN")) {
	if (bcf_get_info_int32(hdr, rec, "SVLEN", &svlen, &nsvlen) > 0) {
	  svlength = *svlen;
	}
      }

      // END
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

      // Numerical SV type
      int32_t svtint = _decodeOrientation(ctval, svtval);
      if (svtint == -1) continue;
      int32_t qualval = (int32_t) (rec->qual);

      // Any breakpoint hit?
      typename TSV::iterator itSV = std::lower_bound(svs.begin(), svs.end(), SV(refIndex, std::max(0, startsv - c.bpwindow), refIndex2, endsv), SortSVs<SV>());
      for(;itSV != svs.end(); ++itSV) {
	if (std::abs(itSV->svStart - startsv) > c.bpwindow) break;
	if (itSV->chr2 != refIndex2) continue;
	if ((c.matchSvType) && (itSV->svt != svtint)) continue;
	if (std::abs(itSV->svEnd - endsv) > c.bpwindow) continue;
	if ((itSV->svlen > 0) && (svlength > 0)) {
	  float rat = itSV->svlen / svlength;
	  if (svlength < itSV->svlen) rat = svlength / itSV->svlen;
	  if (rat < c.sizediff) continue;
	}
	std::string id("id");
        std::string padNumber = boost::lexical_cast<std::string>(itSV->id);
        padNumber.insert(padNumber.begin(), 9 - padNumber.length(), '0');
        id += padNumber;
	dataOut << id << "\t" << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << chr2Name << "\t" <<  endsv << "\t" << rec->d.id << "\t" << qualval << "\t" << svtval << "\t" << svlength << std::endl;
      }
      
      // Successful parse
      ++parsedSV;
    }

    // Statistics
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Parsed " << parsedSV << " out of " << sitecount << " VCF/BCF records." << std::endl;
	
    // Clean-up
    if (svend != NULL) free(svend);
    if (svlen != NULL) free(svlen);
    if (pos2 != NULL) free(pos2);
    if (svt != NULL) free(svt);
    if (chr2 != NULL) free(chr2);
    if (ct != NULL) free(ct);

    // Close file handles
    dataOut.pop();
    dataOut.pop();
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    bcf_destroy(rec);
    
    return true;
  }

}

#endif
