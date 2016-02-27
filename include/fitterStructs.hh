#pragma once

#include "Rtypes.h"
#include "TSpline.h"
#include "TemplateFitter.hh"

#include <memory>
#include <map>
#include <string>

#include "common.hh"

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
  UInt_t wiggleRoom;
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

class digitizer {
public:
  virtual UShort_t* getTrace(int i) = 0;
  virtual ULong64_t* getStructAddress() = 0;
  virtual std::size_t getTraceLength() const = 0;
  std::string type;
  std::string branchName;
  std::vector<detector> detectors;
};

class digitizerCaen5730 : public digitizer {
public:
  std::size_t getTraceLength() const { return CAEN_5730_LN; }
  UShort_t* getTrace(int i) { return data.trace[i]; }
  ULong64_t* getStructAddress() { return &data.event_index; }
private:
  daq::caen_5730 data; 
};

class digitizerCaen5742 : public digitizer {
public:
  std::size_t getTraceLength() const { return CAEN_6742_LN; }
  UShort_t* getTrace(int i) { return data.trace[i]; }
  ULong64_t* getStructAddress() { return &data.system_clock; }
private:
  daq::caen_6742 data; 
};
