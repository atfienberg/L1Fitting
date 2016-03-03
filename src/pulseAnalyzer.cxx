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
                         std::vector<std::unique_ptr<digitizer>>& digs);

/**
 * @brief Display a root plot of a pulse fit
 */
void displayFit(TemplateFitter& tf, const TemplateFitter::Output& out,
                const std::vector<UShort_t>& sampleTimes,
                const std::vector<UShort_t>& trace, const detector& det);

void processTrace(UShort_t* trace, detector& det, std::size_t len);

int main(int argc, char const* argv[]) {
  std::string configfile;
  if (argc < 3) {
    std::cout << "Usage: ./pulseAnalyzer <infile> <outfile> [configfile]"
              << std::endl;
    exit(EXIT_FAILURE);
  } else if (argc == 3) {
    configfile =
        "/Users/atfienberg/g-2/frascati/L1Fitting/config/"
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

  std::vector<std::unique_ptr<digitizer>> digs;
  auto conf = parseConfig(configfile, digs);

  // setup input and output files and trees
  TFile inFile(argv[1]);
  std::unique_ptr<TTree> inTree((TTree*)inFile.Get("t"));
  inTree->SetBranchStatus("*", 0);

  TFile outf(argv[2], "recreate");
  TTree outTree("t", "t");
  for (auto& dig : digs) {
    if ((dig->type == "caen5730") || (dig->type == "caen5742")) {
      inTree->SetBranchStatus(dig->branchName.c_str(), 1);

      inTree->SetBranchAddress(dig->branchName.c_str(),
                               dig->getStructAddress());

      for (auto& det : dig->detectors) {
        outTree.Branch(
            det.name.c_str(), &det.pSum.energy,
            "energy/D:baseline/D:threeSampleAmpl/D:time/D:threeSampleTime/"
            "D:chi2/D:fitConverged/O:successfulFit/O");
      }
    }
  }

  for (int i = conf["startEntry"].int_value(); i < inTree->GetEntries(); ++i) {
    inTree->GetEntry(i);

    if (i % 10000 == 0) {
      std::cout << "event " << i << std::endl;
    }

    for (auto& dig : digs) {
      if ((dig->type == "caen5730") || (dig->type == "caen5742")) {
        for (auto& det : dig->detectors) {
          processTrace(dig->getTrace(det.conf.channel), det,
                       dig->getTraceLength());
        }
      }
    }
    outTree.Fill();
  }

  outTree.Write();
  outf.Write();

  return 0;
}

void processTrace(UShort_t* trace, detector& det, std::size_t len) {
  std::vector<UShort_t> fitSamples(det.conf.fitLength);
  UShort_t* peakptr;
  if (det.conf.negPolarity) {
    peakptr = std::min_element(trace, trace + len);
  } else {
    peakptr = std::max_element(trace, trace + len);
  }

  if ((*peakptr < 3000) || (*peakptr > 3450)) {
    return;
  }

  assert(peakptr - det.conf.peakIndex >= trace);
  assert(peakptr - det.conf.peakIndex + det.conf.fitLength <= trace + len);
  std::copy(peakptr - det.conf.peakIndex,
            peakptr - det.conf.peakIndex + det.conf.fitLength,
            fitSamples.begin());

  // try fit at 3 different starting points before giving up
  std::vector<int> timeOffsets = {0, 1, -1};
  TemplateFitter::Output out;
  bool successfulFit = false;
  for (std::size_t i = 0; (!successfulFit) && (i < timeOffsets.size()); ++i) {
    // for now noise is set to one here, doesn't matter as long as it's flat
    out = det.fitter.fit(fitSamples, det.conf.peakIndex + timeOffsets[i]);
    if ((std::abs(out.times[0] - det.conf.peakIndex) < det.conf.wiggleRoom) &&
        (out.converged) &&
        (det.conf.negPolarity ? (out.scales[0] < 0) : (out.scales[0] > 0))) {
      successfulFit = true;
    }
  }

  if (out.chi2 > 150) {
    // try double pulse fit
    auto correctedSamples = fitSamples;
    for (std::size_t i = 0; i < correctedSamples.size(); ++i) {
      correctedSamples[i] -=
          (out.scales[0] * det.templateSpline->Eval(i - out.times[0]) +
           out.pedestal);
    }
    double secondPeak =
        std::min_element(correctedSamples.begin(), correctedSamples.end()) -
        correctedSamples.begin();
    std::cout << "Second peak raw: " << secondPeak << std::endl;
    std::cout << "Second peak: "
              << secondPeak + (peakptr - trace - det.conf.peakIndex)
              << std::endl;
    out = det.fitter.fit(fitSamples,
                         {static_cast<double>(det.conf.peakIndex), secondPeak});

    std::vector<std::pair<double, double>> timeTries = {
        {-2, 1}, {-1, 1}, {-2, 2}, {-5, 0}, {0, 5}, {-10,0}, {0, 10}};
    for (std::size_t i = 0; (out.chi2 > 200) && (i < timeTries.size()); ++i) {
      std::cout << "second double fit try " << i << std::endl;
      out = det.fitter.fit(fitSamples,
                           {det.conf.peakIndex + timeTries[i].first,
                            det.conf.peakIndex + timeTries[i].second});
    }
  }

  double tsa = peakptr[0] +
               (peakptr[1] - peakptr[-1]) * (peakptr[1] - peakptr[-1]) /
                   (16.0 * peakptr[0] - 8.0 * (peakptr[1] + peakptr[-1]));
  double tst = peakptr - trace +
               (peakptr[1] - peakptr[-1]) /
                   (4.0 * peakptr[0] - 2.0 * (peakptr[1] + peakptr[-1]));

  det.pSum = {out.scales[0], out.pedestal, tsa - out.pedestal,
              out.times[0] + (peakptr - trace - det.conf.peakIndex), tst,
              out.chi2, out.converged, successfulFit};

  if (det.conf.negPolarity) {
    det.pSum.energy *= -1;
    det.pSum.threeSampleAmpl *= -1;
  }

  if (det.conf.draw) {
    std::vector<UShort_t> times(fitSamples.size());
    std::iota(times.begin(), times.end(), peakptr - det.conf.peakIndex - trace);
    displayFit(det.fitter, out, times, fitSamples, det);
  }
}
