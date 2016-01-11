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
#include "TApplication.h"
#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TPaveText.h"
#include <TSystem.h>
#include <TF1.h>

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
            "D:chi2/D");
      }
    }
  }

  // do the fitting
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

          assert(peakptr - trace >= det.conf.peakIndex);
          std::copy(peakptr - det.conf.peakIndex,
                    peakptr - det.conf.peakIndex + det.conf.fitLength,
                    fitSamples.begin());

          auto out = det.fitter.fit(fitSamples, det.conf.peakIndex);

          double tsa =
              peakptr[0] +
              (peakptr[1] - peakptr[-1]) * (peakptr[1] - peakptr[-1]) /
                  (16.0 * peakptr[0] - 8.0 * (peakptr[1] + peakptr[-1]));
          double tst =
              peakptr - trace +
              (peakptr[1] - peakptr[-1]) /
                  (4.0 * peakptr[0] - 2.0 * (peakptr[1] + peakptr[-1]));

          det.pSum = {out.scales[0], out.pedestal, tsa - out.pedestal,
                      out.times[0] + (peakptr - trace), tst, out.chi2};

          if (det.conf.negPolarity){
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

/*
  Helper functions
 */

json11::Json valueFromDetectorOrDefault(const std::string& key,
                                        const json11::Json::object& detector,
                                        const json11::Json::object& def) {
  if (detector.find(key) != detector.end()) {
    return detector.at(key);
  } else {
    return def.at(key);
  }
}
json11::Json parseConfig(const std::string& confFileName,
                         std::vector<digitizer>& digs) {
  std::stringstream ss;
  std::ifstream configfile(confFileName);
  ss << configfile.rdbuf();
  configfile.close();
  std::string err;
  auto confJson = json11::Json::parse(ss.str(), err);
  if (err.size() != 0) {
    std::cerr << "Parsing error for " << confFileName << " : " << err
              << std::endl;
    exit(EXIT_FAILURE);
  }
  auto confMap = confJson.object_items();

  digs.resize(0);
  bool drawingAny = false;
  auto defaults = confJson["defaultDetector"].object_items();
  for (const auto& digEntry : confJson["digitizers"].array_items()) {
    digs.push_back(digitizer());
    auto digMap = digEntry.object_items();
    digs.back().type = digMap.at("type").string_value();
    digs.back().branchName = digMap.at("branchName").string_value();

    for (const auto detEntry : digEntry["detectors"].array_items()) {
      auto detectorMap = detEntry.object_items();

      digs.back().detectors.push_back(detector());
      struct detector& thisDetector = digs.back().detectors.back();

      thisDetector.name = detectorMap.at("name").string_value();

      TFile templateFile(
          (confMap.at("templateBaseDir").string_value() + "/" +
           detectorMap.at("templateFile").string_value()).c_str());
      thisDetector.templateSpline.reset(
          (TSpline3*)templateFile.Get("masterSpline"));

      thisDetector.conf.channel = detectorMap.at("channel").int_value();

      thisDetector.conf.templateBuffer =
          valueFromDetectorOrDefault("templateBuffer", detectorMap, defaults)
              .number_value();
      thisDetector.conf.templateLength =
          valueFromDetectorOrDefault("templateLength", detectorMap, defaults)
              .number_value();
      thisDetector.conf.fitLength =
          valueFromDetectorOrDefault("fitLength", detectorMap, defaults)
              .int_value();
      thisDetector.conf.peakIndex =
          valueFromDetectorOrDefault("peakIndex", detectorMap, defaults)
              .int_value();
      thisDetector.conf.negPolarity =
          valueFromDetectorOrDefault("negPolarity", detectorMap, defaults)
              .bool_value();
      thisDetector.conf.draw = valueFromDetectorOrDefault(
                                   "draw", detectorMap, defaults).bool_value();

      thisDetector.fitter.setTemplate(
          thisDetector.templateSpline.get(),
          -1 * thisDetector.conf.templateBuffer,
          thisDetector.conf.templateLength - thisDetector.conf.templateBuffer,
          10000);

      if ((!drawingAny)) {
        drawingAny = thisDetector.conf.draw;
      }
    }
  }

  if (drawingAny) {
    new TApplication("app", 0, nullptr);
  }

  return confJson;
}

void displayFit(TemplateFitter& tf, TemplateFitter::Output out,
                std::vector<UShort_t> sampleTimes, std::vector<UShort_t> trace,
                detector& det) {
  const int nPulses = out.times.size();

  // print to terminal
  std::cout << det.name << std::endl;

  for (int i = 0; i < nPulses; ++i) {
    std::cout << "t" << i + 1 << ": " << out.times[0] + sampleTimes[0]
              << " +/- " << sqrt(tf.getCovariance(i, i)) << std::endl;
    std::cout << "scale" << i + 1 << ": " << out.scales[0] << " +/- "
              << sqrt(tf.getCovariance(i + nPulses, i + nPulses)) << std::endl;
  }
  std::cout << "pedestal: " << out.pedestal << " +/- "
            << sqrt(tf.getCovariance(2 * nPulses, 2 * nPulses)) << std::endl;
  std::cout << "chi2: " << out.chi2 << std::endl;
  std::cout << std::endl;
  std::cout << "covariance matrix" << std::endl;
  for (int i = 0; i < 2 * nPulses + 1; ++i) {
    for (int j = 0; j < 2 * nPulses + 1; ++j) {
      std::cout << std::setw(12) << tf.getCovariance(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // make plot
  std::unique_ptr<TCanvas> c(new TCanvas((det.name + "_canvas").c_str(),
                                         (det.name + "_canvas").c_str()));

  std::unique_ptr<TGraph> g(new TGraph(0));
  g->SetTitle(det.name.c_str());
  for (size_t i = 0; i < trace.size(); ++i) {
    g->SetPoint(g->GetN(), sampleTimes[i], trace[i]);
  }

  const TSpline* tSpline = det.templateSpline.get();
  // room for up to three pulses
  auto templateFunction = [&](double* x, double* p) {
    double returnValue = p[6];
    for (int i = 0; (i < nPulses) && (i < 3); ++i) {
      if ((x[0] - p[2 * i] > -1 * det.conf.templateBuffer) &&
          (x[0] - p[2 * i] <
           det.conf.templateLength - det.conf.templateBuffer)) {
        returnValue += p[1 + 2 * i] * tSpline->Eval(x[0] - p[2 * i]);
      }
    }
    return returnValue;
  };

  std::unique_ptr<TF1> func(new TF1("fitFunc", templateFunction, sampleTimes[0],
                                    sampleTimes.back(), 7));
  func->SetLineColor(kBlack);
  func->SetParameters(std::vector<double>(7, 0).data());

  g->SetMarkerStyle(20);
  g->Draw("ap");
  g->GetXaxis()->SetRangeUser(sampleTimes[0], sampleTimes.back());
  g->GetXaxis()->SetTitle("sample number");
  g->GetYaxis()->SetTitle("ADC counts");
  g->GetYaxis()->SetTitleOffset(1.5);

  double yMin = g->GetYaxis()->GetXmin();
  double yMax = g->GetYaxis()->GetXmax();
  std::unique_ptr<TPaveText> txtbox(
      new TPaveText(15 + sampleTimes[0], yMin + (yMax - yMin) * 0.5,
                    29 + sampleTimes[0], yMin + (yMax - yMin) * 0.1));
  txtbox->SetFillColor(kWhite);
  func->SetParameter(6, out.pedestal);
  for (int i = 0; i < nPulses; ++i) {
    func->SetParameter(2 * i, out.times[i] + sampleTimes[0]);
    txtbox->AddText(Form("t_{%i}: %.3f #pm %.3f", i + 1,
                         out.times[i] + sampleTimes[0],
                         sqrt(tf.getCovariance(i, i))));
    func->SetParameter(2 * i + 1, out.scales[i]);
    txtbox->AddText(Form("E_{%i}: %.0f #pm %.0f", i + 1, out.scales[i],
                         sqrt(tf.getCovariance(nPulses + i, nPulses + i))));
  }
  txtbox->AddText(Form("pedestal: %.0f #pm %.1f", out.pedestal,
                       sqrt(tf.getCovariance(2 * nPulses, 2 * nPulses))));
  txtbox->AddText(Form("#chi^{2} / NDF : %.2f", out.chi2));

  std::vector<std::unique_ptr<TF1>> components;
  if (nPulses > 1) {
    int colors[3] = {kRed, kBlue, kMagenta + 2};
    for (int i = 0; i < nPulses; ++i) {
      components.emplace_back(new TF1("fitFunc", templateFunction, 0, 30, 7));
      components.back()->SetParameters(std::vector<double>(7, 0).data());
      components.back()->SetParameter(6, out.pedestal);
      components.back()->SetParameter(2 * i, out.times[i] + sampleTimes[0]);
      components.back()->SetParameter(2 * i + 1, out.scales[i]);
      components.back()->SetLineColor(colors[i]);
      components.back()->SetNpx(1000);
      components.back()->Draw("same");
    }
  }

  func->SetNpx(1000);
  func->Draw("same");
  func->SetLineColor(kRed);
  txtbox->Draw("same");
  c->Print((det.name + ".pdf").c_str());
  c->Write();

  c->Modified();
  c->Update();
  c->Draw();
  gSystem->ProcessEvents();

  std::cout << det.name << " displayed. Any key to move on" << std::endl;
  std::cin.ignore();
}
