#include "FunzioneDistribuzione2.hpp"  
#include <vector>                       
#include <TGraphErrors.h>               
#include <TCanvas.h>                    
#include <TStyle.h>                     
#include <TColor.h>                     
#include <TF1.h>
#include <cmath>

TGraphErrors* FunzioneDistribuzione2::ErGraficoConErrori(const std::vector<double>& mean,
                                                           const std::vector<double>& stddev,
                                                           double xmin, double xmax)
{
    int bins = mean.size();
    double binWidth = (xmax - xmin) / bins;

    auto* gr = new TGraphErrors(bins);
    for (int i = 0; i < bins; ++i)
    {
        double x = xmin + (i + 0.5) * binWidth;  // centro del bin
        gr->SetPoint(i, x, mean[i]);
        gr->SetPointError(i, 0.0, stddev[i]);
    }

    gr->SetTitle("Distribuzione con errori; x; f(x)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);

    auto* c = new TCanvas("c", "Distribuzione con errori", 800, 600);
    gr->Draw("APE");
    c->SaveAs("DistribuzioneConErrori.pdf");

    return gr;
}


double FunzioneDistribuzione2::Coseno(double* x, double* par)
{
    return pow(cos(par[0] * x[0] + par [1]), 2) + par[2];
}

// Fit del grafico con parametri B e C fissati
void FunzioneDistribuzione2::FitCoseno(TGraphErrors* gr, double b_fixed)
{
    double xmin = gr->GetX()[0];
    double xmax = gr->GetX()[gr->GetN() - 1];

    TF1* f_coseno = new TF1("f_coseno", FunzioneDistribuzione2::Coseno, xmin, xmax, 3);
    f_coseno->SetParameter(0, 5.2); // k libero
    f_coseno->SetParameter(1, 1.8); // phi libero
    f_coseno->FixParameter(2, b_fixed);     // b fissato

    gr->Fit(f_coseno, "R");

    // Disegno del fit sopra il grafico
    auto* cFit = new TCanvas("cFit", "Fit Coseno", 800, 600);
    gr->Draw("APE");
    f_coseno->SetLineColor(kRed);
    f_coseno->Draw("same");
    cFit->Update(); 
    cFit->SaveAs("FitCosenoConParametriFissi.pdf");
}

