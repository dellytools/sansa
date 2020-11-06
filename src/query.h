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

    // Output file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.matchfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dataOut << "[1]ANNOID\tquery.chr\tquery.start\tquery.chr2\tquery.end\tquery.id\tquery.qual\tquery.svtype\tquery.svlen" << std::endl;
    
    // Parse VCF records
    bcf1_t* rec = bcf_init();
    int32_t lastRID = -1;
    int32_t sitecount = 0;
    int32_t parsedSV = 0;
    int32_t refIndex = -1;
    int32_t refIndex2 = -1;
    while (bcf_read(ifile, hdr, rec) == 0) {
      int32_t startsv = rec->pos + 1;
      ++sitecount;

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

      // Parse INFO fields
      std::string svtval = "NA";
      if (!_parseSVTYPE(hdr, rec, svtval)) continue;
      std::string ctval("NA");
      _parse_bcf_string(hdr, rec, "CT", ctval);
      std::string chr2Name(bcf_hdr_id2name(hdr, rec->rid));
      refIndex2 = refIndex;
      if (_parse_bcf_string(hdr, rec, "CHR2", chr2Name)) {
	if (chrMap.find(chr2Name) == chrMap.end()) continue;
	else refIndex2 = chrMap[chr2Name];
      }
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
      if (svtint == -1) continue;
      int32_t qualval = (int32_t) (rec->qual);

      // Any breakpoint hit?
      typename TSV::iterator itSV = std::lower_bound(svs.begin(), svs.end(), SV(refIndex, std::max(0, startsv - c.bpwindow), refIndex2, endsv), SortSVs<SV>());
      int32_t bestID = -1;
      float bestScore = -1;
      bool noMatch = true;
      for(; itSV != svs.end(); ++itSV) {
	int32_t startDiff = std::abs(itSV->svStart - startsv);
	if (startDiff > c.bpwindow) break;
	if (itSV->chr2 != refIndex2) continue;
	if ((c.matchSvType) && (itSV->svt != svtint)) continue;
	int32_t endDiff = std::abs(itSV->svEnd - endsv);
	if (endDiff > c.bpwindow) continue;
	if (itSV->id == -1) continue;
	float score = 0;
	if ((itSV->svlen > 0) && (svlength > 0)) {
	  float rat = itSV->svlen / svlength;
	  if (svlength < itSV->svlen) rat = svlength / itSV->svlen;
	  if (rat < c.sizediff) continue;
	  score += rat;
	}
	
	// Found match
	noMatch = false;
	if (c.bestMatch) {
	  if (c.bpwindow > 0) {
	    if (startDiff > endDiff) score += (1 - float(startDiff) / (float(c.bpwindow)));
	    else score += (1 - float(endDiff) / (float(c.bpwindow)));
	  } else score += 1;
	  if (score > bestScore) {
	    bestScore = score;
	    bestID = itSV->id;
	  }
	} else {
	  std::string id("id");
	  std::string padNumber = boost::lexical_cast<std::string>(itSV->id);
	  padNumber.insert(padNumber.begin(), 9 - padNumber.length(), '0');
	  id += padNumber;
	  dataOut << id << '\t' << bcf_hdr_id2name(hdr, rec->rid) << '\t' << (rec->pos + 1) << '\t' << chr2Name << '\t' <<  endsv << '\t' << rec->d.id << '\t' << qualval << '\t' << svtval << '\t' << svlength << std::endl;
	}
      }
      if (((c.bestMatch) && (bestID != -1)) || ((c.reportNoMatch) && (noMatch))) {
	std::string id("id");
	if (noMatch) {
	  id = "None";
	} else {
	  std::string padNumber = boost::lexical_cast<std::string>(bestID);
	  padNumber.insert(padNumber.begin(), 9 - padNumber.length(), '0');
	  id += padNumber;
	}
	dataOut << id << '\t' << bcf_hdr_id2name(hdr, rec->rid) << '\t' << (rec->pos + 1) << '\t' << chr2Name << '\t' <<  endsv << '\t' << rec->d.id << '\t' << qualval << '\t' << svtval << '\t' << svlength << std::endl;
      }
      
      // Successful parse
      ++parsedSV;
    }

    // Statistics
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Parsed " << parsedSV << " out of " << sitecount << " VCF/BCF records." << std::endl;
	
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
