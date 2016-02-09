/**
 * Aaron Fienberg
 * fienberg@uw.edu
 *
 * code for pulse analysis from UW lab daq
 */

// std includes
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <cassert>

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// project includes
#include "fitterStructs.hh"
#include "json11.hpp"

/**
 * @brief check if file exists
 *
 * @param name path to file
 * @return true if file exists, false otherwise
 */
bool exists(const std::string& name) {
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}

/**
 * @brief parse config file, build collection of fitConfigurations/
 * also construct TApplication if any drawing is to happen
 */
json11::Json parseConfig(const std::string& confFileName,
                         std::vector<digitizer>& digs);

/**
 * @brief Display a root plot of a pulse fit
 */
void displayFit(TemplateFitter& tf, TemplateFitter::Output out,
                std::vector<UShort_t> sampleTimes, std::vector<UShort_t> trace,
                detector& det);

int main(int argc, char const* argv[]) {
  std::string configfile;
  if (argc < 3) {
    std::cout << "Usage: ./pulseAnalyzer <infile> <outfile> [configfile]"
              << std::endl;
    exit(EXIT_FAILURE);
  } else if (argc == 3) {
    configfile =
        "/home/newg2/Workspace/L1Tests/fitting/config/"
        "defaultFitConfig.json";
  } else {
    configfile = argv[3];
  }

  for (auto& file : std::vector<std::string>({argv[1], configfile})) {
    if (!exists(file)) {
      std::cerr << "Error: " << file << " doesn't exist" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  std::vector<digitizer> digs;
  auto conf = parseConfig(configfile, digs);

  // setup input and output files and trees
  TFile inFile(argv[1]);
  std::unique_ptr<TTree> inTree((TTree*)inFile.Get("t"));
  inTree->SetBranchStatus("*", 0);

  TFile outf(argv[2], "recreate");
  TTree outTree("t", "t");
  for (auto& dig : digs) {
    if (dig.type == "caen1742") {
      inTree->SetBranchStatus(dig.branchName.c_str(), 1);

      inTree->SetBranchAddress(dig.branchName.c_str(),
                               &dig.daqData.system_clock);

      for (auto& det : dig.detectors) {
        outTree.Branch(
            det.name.c_str(), &det.pSum.energy,
            "energy/D:baseline/D:threeSampleAmpl/D:time/D:threeSampleTime/"
            "D:chi2/D:fitConverged/O");
      }
    }
  }

  // do the fitting
  double chi2cutoff = conf.object_items().at("chi2Cutoff").number_value();

  for (int i = conf["startEntry"].int_value(); i < inTree->GetEntries(); ++i) {
    inTree->GetEntry(i);

    for (auto& dig : digs) {
      if (dig.type == "caen1742") {
        for (auto& det : dig.detectors) {
          UShort_t* trace = dig.daqData.trace[det.conf.channel];

          std::vector<UShort_t> fitSamples(det.conf.fitLength);
          UShort_t* peakptr;
          if (det.conf.negPolarity) {
            peakptr = std::min_element(trace, trace + CAEN_1742_LN);
          } else {
            peakptr = std::max_element(trace, trace + CAEN_1742_LN);
          }

          assert(peakptr - det.conf.peakIndex >= trace);
          assert(peakptr - det.conf.peakIndex + det.conf.fitLength <=
                 trace + CAEN_1742_LN);
          std::copy(peakptr - det.conf.peakIndex,
                    peakptr - det.conf.peakIndex + det.conf.fitLength,
                    fitSamples.begin());

          auto out = det.fitter.fit(fitSamples, det.conf.peakIndex);
	  //try again if bad fit
	  if ((out.chi2 > chi2cutoff) || 
	      (det.conf.negPolarity ? (out.scales[0] > 0) : (out.scales[0] < 0))) {
	    out = det.fitter.fit(fitSamples, det.conf.peakIndex - 1);
	  }
	  //try one more time
	  if ((out.chi2 > chi2cutoff) || 
	      (det.conf.negPolarity ? (out.scales[0] > 0) : (out.scales[0] < 0))) {
	    out = det.fitter.fit(fitSamples, det.conf.peakIndex + 1);
	  }
	  //this fit has failed
	  if ((out.chi2 > chi2cutoff) || 
	      (det.conf.negPolarity ? (out.scales[0] > 0) : (out.scales[0] < 0))) {		
	    out.converged = false;
	  }
       
          double tsa =
              peakptr[0] +
              (peakptr[1] - peakptr[-1]) * (peakptr[1] - peakptr[-1]) /
                  (16.0 * peakptr[0] - 8.0 * (peakptr[1] + peakptr[-1]));
          double tst =
              peakptr - trace +
              (peakptr[1] - peakptr[-1]) /
                  (4.0 * peakptr[0] - 2.0 * (peakptr[1] + peakptr[-1]));

          det.pSum = {out.scales[0], out.pedestal, tsa - out.pedestal,
                      out.times[0] + (peakptr - trace), tst, out.chi2,
                      out.converged};

          if (det.conf.negPolarity) {
            det.pSum.energy *= -1;
            det.pSum.threeSampleAmpl *= -1;
          }

          if (det.conf.draw) {
            std::vector<UShort_t> times(fitSamples.size());
            std::iota(times.begin(), times.end(),
                      peakptr - det.conf.peakIndex - trace);
            displayFit(det.fitter, out, times, fitSamples, det);
          }
        }
      }
    }
    outTree.Fill();
  }

  outTree.Write();
  outf.Write();

  return 0;
}
