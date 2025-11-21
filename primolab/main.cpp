#include "FunzioneDistribuzione.hpp"
#include <iostream>

int main()
{
    FunzioneDistribuzione f{5.2, 1.8, 0.2};

    auto ftodraw = f.GetFunction();
    TCanvas c1{};
    ftodraw->Draw();
    c1.SaveAs("funzione.png");

    auto htodraw = f.ErHistoGeneratore(10000, 100);
    TCanvas c2{};
    htodraw->Draw("HIST");
    c2.SaveAs("histo.png");

    auto f_normy = f.GetFunctionNormalized();
    double rmse = f.CarcolaRMSE(htodraw, f_normy, 10000);
    std::cout << "RMSE = " << rmse << std::endl;

    auto hist_normy = f.ErHistoGeneratoreNormalized(10000, 100);
    TCanvas c3{};
    hist_normy->Draw("HIST");
    f_normy->Draw("Same");
    c3.SaveAs("fhistonormalized.png");

    TCanvas c4{};
    htodraw->Draw("HIST");
    c4.SaveAs("histonormalized.png");

    int N = 10000;
    int bins = 50;
    int M = 100;
    auto [mean1, stddev1] = f.StimaIncertezza(N, bins, M);

    // Ad esempio, stampare le incertezze dei primi 5 bin
    for (int i = 0; i < 5; ++i)
    {
        std::cout << "Bin " << i + 1 << ": mean = " << mean1[i] << ", stddev = " << stddev1[i] << std::endl;
    }

    TH1D *h_smear = f.BinSmearing(bins, 0.0, 4.5, 0.05);
    TCanvas *c5 = new TCanvas("c5", "Bin Smearing Example", 800, 600);
    h_smear->Draw("HIST");
    c5->SaveAs("bin_smearing.png");

    delete h_smear;
    delete c5;

    auto [mean2, stddev2] = f.StimaIncertezzaBinSmearing(bins, M, 0.0, 4.5, 0.05);

    // 3) Stampiamo i risultati (media e deviazione standard per i primi 5 bin, per esempio)
    std::cout << "Bin\tMedia\tIncertezza (stddev)" << std::endl;
    for (int i = 0; i < 5; ++i)
    {
        std::cout << i + 1 << "\t" << mean2[i] << "\t" << stddev2[i] << std::endl;
    }

    double sigma_k = 0.05;
    double sigma_phi = 0.1;
    double sigma_b = 0.02;
    auto [mean_rig, std_rig, mean_sme, std_sme] =
        f.PropagaIncertezzeParametri(sigma_k, sigma_phi, sigma_b,
                                     5000, 50,
                                     100,
                                     50,
                                     0.0, 4.5, 0.05);

    // Stampa un riepilogo per alcuni bin
    for (int i = 0; i < 5; ++i)
    {
        std::cout << "Bin " << i + 1
                  << "  ⟨rigenerazione⟩=" << mean_rig[i]
                  << " ± " << std_rig[i]
                  << "  ⟨smearing⟩=" << mean_sme[i]
                  << " ± " << std_sme[i]
                  << std::endl;
    }

    return 0;
}
