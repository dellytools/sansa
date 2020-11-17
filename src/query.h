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



  template<typename TConfig, typename TGenomicRegions, typename TGeneIds>
  inline void
  geneAnnotation(TConfig const& c, TGenomicRegions const& gRegions, TGeneIds const& geneIds, int32_t const refIndex, int32_t const svStart, int32_t const refIndex2, int32_t const svEnd, std::string& featureBp1, std::string& featureBp2) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;

    // Search nearby genes
    for(uint32_t bp = 0; bp < 2; ++bp) {
      // Distance vector
      typedef std::pair<int32_t, int32_t> TDistOffset;
      typedef std::vector<TDistOffset> TDistVector;
      TDistVector dist;

      // Fetch overlapping features
      int32_t rid = refIndex;
      int32_t bpoint = svStart;
      if (bp) {
	rid = refIndex2;
	bpoint = svEnd;
      }
      for(typename TChromosomeRegions::const_iterator itg = gRegions[rid].begin(); itg != gRegions[rid].end(); ++itg) {
	if (itg->start - bpoint > c.maxDistance) break;
	if (bpoint - itg->end > c.maxDistance) continue;
	int32_t featureDist = 0;
	if (bpoint > itg->end) featureDist = itg->end - bpoint;
	if (bpoint < itg->start) featureDist = itg->start - bpoint;
	dist.push_back(std::make_pair(featureDist, itg - gRegions[rid].begin()));
      }

      // Sort by distance
      std::sort(dist.begin(), dist.end());

      // Assign gene IDs
      bool firstFeature = true;
      for(uint32_t i = 0; i < dist.size(); ++i) {
	if (!firstFeature) {
	  if (!bp) featureBp1 += ",";
	  else featureBp2 += ",";
	} else firstFeature = false;
	if (!bp) {
	  featureBp1 += geneIds[gRegions[rid][dist[i].second].lid] + '(' + boost::lexical_cast<std::string>(dist[i].first) + ';' + gRegions[rid][dist[i].second].strand + ')' ;
	} else {
	  featureBp2 += geneIds[gRegions[rid][dist[i].second].lid] + '(' + boost::lexical_cast<std::string>(dist[i].first) + ';' + gRegions[rid][dist[i].second].strand + ')' ;
	}
      }
    }
  }




  
  template<typename TConfig, typename TSV, typename TGenomicRegions, typename TGeneIds>
  inline bool
  query(TConfig& c, TSV& svs, TGenomicRegions& gRegions, TGeneIds& geneIds) {

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
    dataOut << "[1]ANNOID\tquery.chr\tquery.start\tquery.chr2\tquery.end\tquery.id\tquery.qual\tquery.svtype\tquery.ct\tquery.svlen\tquery.startfeature\tquery.endfeature" << std::endl;
    
    // Parse VCF records
    bcf1_t* rec = bcf_init();
    int32_t lastRID = -1;
    int32_t sitecount = 0;
    int32_t parsedSV = 0;
    int32_t refIndex = -1;
    while (bcf_read(ifile, hdr, rec) == 0) {
      ++sitecount;

      // New chromosome?
      if (rec->rid != lastRID) {
	lastRID = rec->rid;
	std::string chrName = bcf_hdr_id2name(hdr, rec->rid);	
	if (c.nchr.find(chrName) == c.nchr.end()) refIndex = -1;
	else refIndex = c.nchr[chrName];
      }
      if (refIndex == -1) continue;
      
      // Unpack INFO
      bcf_unpack(rec, BCF_UN_INFO);

      // Parse INFO fields
      std::string svtval = "NA";
      if (!_parseSVTYPE(hdr, rec, svtval)) continue;
      std::string ctval("NA");
      _parse_bcf_string(hdr, rec, "CT", ctval);
      std::string chr2Name(bcf_hdr_id2name(hdr, rec->rid));
      if (_parse_bcf_string(hdr, rec, "CHR2", chr2Name)) {
	if (c.nchr.find(chr2Name) == c.nchr.end()) continue;
      }
      int32_t pos2val = -1;
      _parse_bcf_int32(hdr, rec, "POS2", pos2val);
      int32_t endval = -1;
      _parse_bcf_int32(hdr, rec, "END", endval);
      int32_t svlenval = 0;
      _parse_bcf_int32(hdr, rec, "SVLEN", svlenval);
      
      // Derive proper END and SVLEN
      int32_t endsv = deriveEndPos(rec, svtval, pos2val, endval);
      bool parseALTBND = parseAltBnd(hdr, rec, svtval, ctval, chr2Name, endsv);
      if (!parseALTBND) continue;
      int32_t refIndex2 = c.nchr[chr2Name];
      int32_t svlength = deriveSvLength(rec, svtval, endsv, svlenval);

      // Numerical SV type
      int32_t svtint = _decodeOrientation(ctval, svtval);
      if (svtint == -1) continue;
      int32_t qualval = 0;
      if (rec->qual > 0) qualval = (int32_t) (rec->qual);
      //std::cerr << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << chr2Name << "\t" << endsv << "\t" << rec->d.id << "\t" << qualval << "\t" << svtval << "\t" << ctval << "\t" << svtint << "\t" << svlength << std::endl;

      // Generate query SV
      SV qsv = SV(refIndex, rec->pos + 1, refIndex2, endsv, 0, qualval, svtint, svlength);
      _makeCanonical(qsv);
            
      // Annotate genes
      std::string featureBp1 = "";
      std::string featureBp2 = "";
      if (c.gtfFileFormat != -1) geneAnnotation(c, gRegions, geneIds, qsv.chr, qsv.svStart, qsv.chr2, qsv.svEnd, featureBp1, featureBp2);
      if (featureBp1.empty()) featureBp1 = "NA";
      if (featureBp2.empty()) featureBp2 = "NA";

      // Any breakpoint hit?
      typename TSV::iterator itSV = std::lower_bound(svs.begin(), svs.end(), SV(qsv.chr, std::max(0, qsv.svStart - c.bpwindow), qsv.chr2, qsv.svEnd), SortSVs<SV>());
      int32_t bestID = -1;
      float bestScore = -1;
      bool noMatch = true;
      for(; itSV != svs.end(); ++itSV) {	
	int32_t startDiff = std::abs(itSV->svStart - qsv.svStart);
	if (startDiff > c.bpwindow) break;
	if (itSV->chr2 != qsv.chr2) continue;
	if ((c.matchSvType) && (itSV->svt != qsv.svt)) continue;
	int32_t endDiff = std::abs(itSV->svEnd - qsv.svEnd);
	if (endDiff > c.bpwindow) continue;
	if (itSV->id == -1) continue;
	float score = 0;
	if ((itSV->svlen > 0) && (qsv.svlen > 0)) {
	  float rat = itSV->svlen / qsv.svlen;
	  if (qsv.svlen < itSV->svlen) rat = qsv.svlen / itSV->svlen;
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
	  dataOut << id << '\t' << bcf_hdr_id2name(hdr, qsv.chr) << '\t' << qsv.svStart << '\t' << bcf_hdr_id2name(hdr, qsv.chr2) << '\t' <<  qsv.svEnd << '\t' << rec->d.id << '\t' << qualval << '\t' << _translateSvType(qsv.svt) << '\t' << _translateCt(qsv.svt) << '\t' << svlength << '\t' << featureBp1 << '\t' << featureBp2 << std::endl;
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
	dataOut << id << '\t' << bcf_hdr_id2name(hdr, qsv.chr) << '\t' << qsv.svStart << '\t' << bcf_hdr_id2name(hdr, qsv.chr2) << '\t' <<  qsv.svEnd << '\t' << rec->d.id << '\t' << qualval << '\t' << _translateSvType(qsv.svt) << '\t' << _translateCt(qsv.svt) << '\t' << svlength << '\t' << featureBp1 << '\t' << featureBp2 << std::endl;
      }

      // Successful parse
      ++parsedSV;
    }
    bcf_destroy(rec);

    // Statistics
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Parsed " << parsedSV << " out of " << sitecount << " VCF/BCF records." << std::endl;
	
    // Close file handles
    dataOut.pop();
    dataOut.pop();

    bcf_hdr_destroy(hdr);
    bcf_close(ifile);

    return true;
  }

}

#endif
