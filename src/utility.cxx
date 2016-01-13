/**
 * Some utility helper functions used throughout
 */

#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TGraph.h"

#include "fitterStructs.hh"

#include "json11.hpp"

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
      new TPaveText(18 + sampleTimes[0], yMin + (yMax - yMin) * 0.5,
                    28 + sampleTimes[0], yMin + (yMax - yMin) * 0.1));
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
