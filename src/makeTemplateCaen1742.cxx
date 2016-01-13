/*
  Aaron Fienberg
  fienberg@uw.edu
  code for generating "fuzzy templates" based on digitized datasets
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "TSystem.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TSpline.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TString.h"
#include "time.h"

#include "json11.hpp"

#include "daqStructs.hh"
#include "fitterStructs.hh"

using namespace std;

namespace {
int templateLength;
int nBinsPseudoTime;
int nTimeBins;
int traceLength = CAEN_1742_LN;
int baselineFitLength;
int bufferZone;
int minPeak;
int channel;
bool negPolarity;
};

typedef struct traceSummary {
  double pseudoTime;
  int peakIndex;
  double baseline;
  double integral;
  double normalizedAmpl;
  bool bad;
} traceSummary;

traceSummary processTrace(unsigned short* trace);
vector<double> correctTrace(unsigned short* trace, traceSummary summary);
void readConfigs(const char* fitConf, const char* detectorName);
json11::Json valueFromDetectorOrDefault(const std::string& key,
                                        const json11::Json::object& detector,
                                        const json11::Json::object& def);

int main(int argc, char* argv[]) {
  // clock_t t1, t2;
  // t1 = clock();

  if (argc < 4) {
    cout << "usage: ./makeTemplate <inputfile> <outputfile> <detectorName> "
            "[fitter config]" << endl;
    return -1;
  }

  if (argc == 4) {
    readConfigs(
        "/home/newg2/Workspace/L1Tests/fitting/config/defaultFitConfig.json",
        argv[3]);
  } else {
    readConfigs(argv[4], argv[3]);
  }

  // read input file
  gSystem->Load("libTree");
  TFile infile(argv[1]);
  TTree* t = (TTree*)infile.Get("t");
  caen_1742 c;
  t->SetBranchAddress("caen_0", &c.system_clock);

  // process traces
  // cout << "Processing traces... " << endl;
  vector<traceSummary> summaries(t->GetEntries());
  TH1D pseudoTimesHist("ptimes", "ptimes", nBinsPseudoTime, 0, 1);
  TH1D normalizedMaxes("maxes", "maxes", 100, 0.0, 0.0);
  TH1D integralHist("integrals", "integrals", 100, 0.0, 0.0);
  for (int i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    summaries[i] = processTrace(c.trace[channel]);
    pseudoTimesHist.Fill(summaries[i].pseudoTime);
    normalizedMaxes.Fill(summaries[i].normalizedAmpl);
    integralHist.Fill(summaries[i].integral);
    if (i % 1000 == 0) {
      // cout << "Trace " << i << " processed." << endl;
    }
  }
  pseudoTimesHist.Scale(1.0 / pseudoTimesHist.Integral());

  // find max for fuzzy template bin range
  normalizedMaxes.Fit("gaus", "q0");
  double binRangeMax = normalizedMaxes.GetFunction("gaus")->GetParameter(1) +
                       5 * normalizedMaxes.GetFunction("gaus")->GetParameter(2);

  // create map to real time
  TGraph realTimes(0);
  realTimes.SetName("realTimeGraph");
  realTimes.SetPoint(0, 0, 0);
  for (int i = 0; i < nBinsPseudoTime; ++i) {
    realTimes.SetPoint(i, pseudoTimesHist.GetBinLowEdge(i + 2),
                       pseudoTimesHist.Integral(1, i + 1));
  }
  TSpline3 rtSpline = TSpline3("realTimeSpline", &realTimes);
  rtSpline.SetName("realTimeSpline");

  // fill the timeslices and make the master fuzzy template
  TH2D masterFuzzyTemplate =
      TH2D("masterFuzzy", "Fuzzy Template", templateLength * nTimeBins,
           -.5 - bufferZone, templateLength - .5 - bufferZone, 1000,
           -.2 * binRangeMax, binRangeMax);

  // cout << "Populating timeslices... " << endl;
  for (int i = 3; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (summaries[i].bad) {
      continue;
    }
    double realTime = rtSpline.Eval(summaries[i].pseudoTime);
    int thisSlice = static_cast<int>(realTime * nTimeBins);
    if (thisSlice == nTimeBins) --thisSlice;
    auto ctrace = correctTrace(c.trace[channel], summaries[i]);
    for (int j = 0; j < templateLength; ++j) {
      masterFuzzyTemplate.Fill(j - realTime + 0.5 - bufferZone, ctrace[j]);
    }
    if (i % 1000 == 0) {
      // cout << "Trace " << i << " placed." << endl;
    }
  }

  // step through fuzzy template to get errors and means
  // cout << "Calculating errors and means... " << endl;
  TGraphErrors masterGraph(0);
  masterGraph.SetName("masterGraph");
  TGraph errorGraph(0);
  errorGraph.SetName("errorGraph");
  TGraph errorVsMean(0);
  errorVsMean.SetName("errorVsMean");
  for (int i = 0; i < templateLength * nTimeBins; ++i) {
    TH1D* xBinHist = masterFuzzyTemplate.ProjectionY("binhist", i + 1, i + 1);
    xBinHist->Fit("gaus", "q0", "",
                  xBinHist->GetMean() - xBinHist->GetRMS() * 3,
                  xBinHist->GetMean() + xBinHist->GetRMS() * 3);
    double mean = xBinHist->GetFunction("gaus")->GetParameter(1);
    double sig = xBinHist->GetFunction("gaus")->GetParameter(2);
    errorGraph.SetPoint(i, static_cast<double>(i) / nTimeBins - bufferZone - .5,
                        sig);
    masterGraph.SetPoint(
        i, static_cast<double>(i) / nTimeBins - bufferZone - .5, mean);
    masterGraph.SetPointError(i, 0, sig);
    errorVsMean.SetPoint(i, mean, sig);
    delete xBinHist;
  }
  // cout << "Errors and Means Calculated" << endl;

  TSpline3 masterSpline("masterSpline", &masterGraph);
  masterSpline.SetName("masterSpline");
  masterSpline.SetNpx(10000);
  TSpline3 errorSpline("errorSpline", &errorGraph);
  errorSpline.SetName("errorSpline");
  errorSpline.SetNpx(10000);

  // save data
  TFile outf(argv[2], "recreate");
  rtSpline.Write();
  pseudoTimesHist.Write();
  masterFuzzyTemplate.Write();
  errorGraph.Write();
  masterGraph.Write();
  masterSpline.Write();
  errorSpline.Write();
  errorVsMean.Write();
  outf.Write();
  outf.Close();

  // finish up
  delete t;
  // t2 = clock();
  // double diff((double)t2 - (double)t1);
  // cout << "Time elapsed: " << diff / CLOCKS_PER_SEC << " seconds." << endl;
  return 0;
}

traceSummary processTrace(unsigned short* trace) {
  traceSummary results;
  results.bad = false;

  // find maximum
  int maxdex = 0;
  for (int i = 0; i < traceLength; ++i) {
    if (negPolarity) {
      maxdex = trace[i] < trace[maxdex] ? i : maxdex;
    } else {
      maxdex = trace[i] > trace[maxdex] ? i : maxdex;
    }
  }
  results.peakIndex = maxdex;

  // calculate pseudotime
  if (trace[maxdex] == trace[maxdex + 1])
    results.pseudoTime = 1;
  else {
    results.pseudoTime =
        2.0 / M_PI *
        atan(static_cast<double>(trace[maxdex - 1] - trace[maxdex]) /
             (trace[maxdex + 1] - trace[maxdex]));
  }

  if (negPolarity) {
    if (trace[maxdex] > minPeak) {
      results.bad = true;
      return results;
    }
  } else {
    if (trace[maxdex] < minPeak) {
      results.bad = true;
      return results;
    }
  }

  // get the baseline
  if (maxdex - baselineFitLength - bufferZone < 0) {
    cout << "Baseline fit walked off the end of the trace!" << endl;
    results.bad = true;
    return results;
  }
  double runningBaseline = 0;
  for (int i = 0; i < baselineFitLength; ++i) {
    runningBaseline =
        runningBaseline + trace[maxdex - bufferZone - baselineFitLength + i];
  }
  results.baseline = runningBaseline / baselineFitLength;

  // get the normalization
  if (maxdex - bufferZone + templateLength > traceLength) {
    results.bad = true;
    return results;
  }
  double runningIntegral = 0;
  for (int i = 0; i < templateLength; ++i) {
    runningIntegral =
        runningIntegral + trace[maxdex - bufferZone + i] - results.baseline;
  }
  results.integral = runningIntegral;

  results.normalizedAmpl =
      (trace[maxdex] - results.baseline) / results.integral;

  return results;
}

vector<double> correctTrace(unsigned short* trace, traceSummary summary) {
  vector<double> correctedTrace(templateLength);
  if (summary.bad) {
    for (int i = 0; i < templateLength; ++i) correctedTrace[i] = 0;
    return correctedTrace;
  }
  for (int i = 0; i < templateLength; ++i) {
    correctedTrace[i] =
        (trace[summary.peakIndex - bufferZone + i] - summary.baseline) * 1.0 /
        (summary.integral);
  }
  return correctedTrace;
}

void readConfigs(const char* fitConf, const char* detectorName) {
  // first read the templateConf
  const char* tempConfName =
      "/home/newg2/Workspace/L1Tests/fitting/config/makeTemplateConf.json";
  std::stringstream ss;
  std::ifstream configfile(tempConfName);
  ss << configfile.rdbuf();
  configfile.close();
  std::string err;
  auto confJson = json11::Json::parse(ss.str(), err);
  if (err.size() != 0) {
    std::cerr << "Parsing error for " << tempConfName << " : " << err
              << std::endl;
    exit(EXIT_FAILURE);
  }
  auto confMap = confJson.object_items();

  nBinsPseudoTime = confMap.at("nBinsPseudoTime").int_value();
  nTimeBins = confMap.at("nTimeBins").int_value();
  baselineFitLength = confMap.at("baselineFitLength").int_value();
  minPeak = confMap.at("minPeak").int_value();

  // now other info from detector conf
  ss.str("");
  configfile.open(fitConf);
  ss << configfile.rdbuf();
  configfile.close();
  confJson = json11::Json::parse(ss.str(), err);
  if (err.size() != 0) {
    std::cerr << "Parsing error for " << fitConf << " : " << err << std::endl;
    exit(EXIT_FAILURE);
  }
  confMap = confJson.object_items();

  // find detector configuration
  auto defaults = confMap.at("defaultDetector").object_items();
  json11::Json::object detector;
  bool found = false;
  for (auto dig : confMap.at("digitizers").array_items()) {
    auto digMap = dig.object_items();
    if (digMap.at("type") == "caen1742") {
      for (auto det : digMap.at("detectors").array_items()) {
        if (det.object_items().at("name") == detectorName) {
          found = true;
          detector = det.object_items();
          break;
        }
      }
    }
    if (found) {
      break;
    }
  }

  if (!found) {
    cerr << detectorName << " not in config file " << fitConf << endl;
    exit(EXIT_FAILURE);
  }

  templateLength = valueFromDetectorOrDefault("templateLength", detector,
                                              defaults).int_value();
  bufferZone = valueFromDetectorOrDefault("templateBuffer", detector, defaults)
                   .int_value();
  negPolarity = valueFromDetectorOrDefault("negPolarity", detector, defaults)
                    .bool_value();
  channel = detector.at("channel").int_value();
}