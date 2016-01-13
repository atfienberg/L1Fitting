#pragma once

#include "Rtypes.h"
#include "TSpline.h"
#include "TemplateFitter.hh"

#include <memory>
#include <map>
#include <string>

#include "daqStructs.hh"

/**
 * structs used in pulse analysis program
 */

struct pulseSummary {
  Double_t energy;
  Double_t baseline;
  Double_t threeSampleAmpl;
  Double_t time;
  Double_t threeSampleTime;
  Double_t chi2;
  Bool_t fitConverged;
};

struct fitConfiguration {
  UInt_t channel;
  Double_t templateBuffer;
  Double_t templateLength;
  UInt_t fitLength;
  UInt_t peakIndex;
  Bool_t negPolarity;
  Bool_t draw;
};

struct detector {
  std::string name;
  fitConfiguration conf;
  std::unique_ptr<TSpline3> templateSpline;
  TemplateFitter fitter;
  pulseSummary pSum;
};

struct digitizer {
  caen_1742
      daqData;  // will have to be variant if more than one type of digitizer
  std::string type;
  std::string branchName;
  std::vector<detector> detectors;
};
