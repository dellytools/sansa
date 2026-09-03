#ifndef SOMATIC_H
#define SOMATIC_H

#include <iostream>
#include <set>
#include <vector>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>

#include <htslib/vcf.h>

#include "version.h"
#include "util.h"
#include "compvcf.h"

namespace sansa
{

  inline int
  somaticRun(CompvcfConfig& c) {
    // Load the VCF/BCFs
    std::vector<CompSVRecord> ponsv;
    if (!_loadCompSVs(c, c.base.string(), ponsv)) return -1;
    std::vector<CompSVRecord> tumorsv;
    if (!_loadCompSVs(c, c.vcffile.string(), tumorsv)) return -1;

    // Sort SVs
    sort(ponsv.begin(), ponsv.end());
    sort(tumorsv.begin(), tumorsv.end());

    // Compare tumor SVs against the PoN
    compareSVs(c, ponsv, tumorsv);

    // Somatic SVs are absent in PoN
    std::set<std::string> somaticIds;
    for(uint32_t j = 0; j < tumorsv.size(); ++j) {
      if (!tumorsv[j].match) somaticIds.insert(tumorsv[j].id);
    }
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << somaticIds.size() << " of " << tumorsv.size() << " tumor SVs are absent in the PoN (somatic)" << std::endl;

    // Write somatic SVs
    htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to load " << c.vcffile.string() << std::endl;
      return -1;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    if (hdr == NULL) {
      std::cerr << "Fail to read header of " << c.vcffile.string() << std::endl;
      bcf_close(ifile);
      return -1;
    }
    htsFile* ofile = hts_open(c.outfile.string().c_str(), "wb");
    if (ofile == NULL) {
      std::cerr << "Fail to open output file " << c.outfile.string() << std::endl;
      bcf_hdr_destroy(hdr);
      bcf_close(ifile);
      return -1;
    }
    if (bcf_hdr_write(ofile, hdr) != 0) {
      std::cerr << "Error: Failed to write BCF header!" << std::endl;
      hts_close(ofile);
      bcf_hdr_destroy(hdr);
      bcf_close(ifile);
      return -1;
    }
    bcf1_t* rec = bcf_init();
    while (bcf_read(ifile, hdr, rec) == 0) {
      bcf_unpack(rec, BCF_UN_STR);
      if (somaticIds.find(std::string(rec->d.id)) != somaticIds.end()) bcf_write1(ofile, hdr, rec);
    }
    bcf_destroy(rec);
    hts_close(ofile);
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);

    // Build BCF index
    bcf_index_build(c.outfile.string().c_str(), 14);
    
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  int somatic(int argc, char **argv) {
    CompvcfConfig c;

    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("pon,a", boost::program_options::value<boost::filesystem::path>(&c.base), "panel of normals (PoN) VCF/BCF file")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("somatic.bcf"), "output somatic VCF/BCF file")
      ("quality,y", boost::program_options::value<int32_t>(&c.qualthres)->default_value(200), "min. SV site quality")
      ("minsize,m", boost::program_options::value<int32_t>(&c.minsize)->default_value(0), "min. SV size")
      ("maxsize,n", boost::program_options::value<int32_t>(&c.maxsize)->default_value(250000000), "max. SV size")
      ("bpdiff,b", boost::program_options::value<int32_t>(&c.bpdiff)->default_value(1000), "max. SV breakpoint offset")
      ("sizeratio,s", boost::program_options::value<float>(&c.sizeratio)->default_value(0.5), "min. SV size ratio")
      ("divergence,d", boost::program_options::value<float>(&c.divergence)->default_value(0.3), "max. SV allele divergence")
      ("troffset", boost::program_options::value<int32_t>(&c.trOffset)->default_value(200), "tandem repeat max. breakpoint offset")
      ("trfrac", boost::program_options::value<float>(&c.trFrac)->default_value(0.25), "tandem repeat fraction")
      ("trseqid", boost::program_options::value<float>(&c.trSeqId)->default_value(0.7), "min. tandem repeat seq. identity")
      ("nosvt,t", "ignore the SV type")
      ("pass,p", "filter sites for PASS")
      ("ct,c", "require matching CT value in addition to SV type")
      ;

    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "tumor VCF/BCF file")
      ;
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("pon"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: sansa " << argv[0] << " [OPTIONS] -a <pon.bcf> <tumor.bcf>" << std::endl;
      std::cerr << visible_options << "\n";
      return 0;
    }

    // Filter for PASS
    if (vm.count("pass")) c.filterForPass = true;
    else c.filterForPass = false;

    // Ignore SV type
    if (vm.count("nosvt")) c.checkSVT = false;
    else c.checkSVT = true;

    // Check CT
    if (vm.count("ct")) c.checkCT = true;
    else c.checkCT = false;

    // Ignore AC for somatic filtering
    c.checkID = true;
    c.minac = 0;
    c.maxac = 100000000;

    // Check PoN file
    if (!(boost::filesystem::exists(c.base) && boost::filesystem::is_regular_file(c.base) && boost::filesystem::file_size(c.base))) {
      std::cerr << "PoN VCF/BCF file is missing: " << c.base.string() << std::endl;
      return 1;
    }
    // Check tumor file
    if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
      std::cerr << "Tumor VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
      return 1;
    }

    // Show cmd
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] ";
    std::cerr << "sansa ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;

    // Run somatic filtering
    return somaticRun(c);
  }

}

#endif
