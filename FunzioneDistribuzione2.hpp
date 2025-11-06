#include <fstream>
#include <iostream>
#include <cmath>

#include "TRandom3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include <TF1.h>

#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include <cmath>
#include <vector>

class FunzioneDistribuzione
{
private:
    double k;
    double phi;
    double b;
    TF1 *functionf; // Puntatore, così possiamo inizializzarlo dinamicamente

public:
    FunzioneDistribuzione(double k_val, double phi_val, double b_val)
        : k(k_val), phi(phi_val), b(b_val)
    {
        // Lambda che accede ai membri della classe
        functionf = new TF1("f", [this](double *x, double *)
                            { return std::pow(std::cos(k * x[0] + phi), 2) + b; }, 0, 4.5, 0); // dominio da 0 a 4.5

        functionf->SetTitle("f(x) = cos^{2}(k x + #phi) + b; x; f(x)");
        functionf->SetLineColor(kBlue);
        functionf->SetLineWidth(2);
    }

    TF1 *GetFunction() const
    {
        return functionf;
    }

    TH1D *ErHistoGeneratore(int N, int bins, double xmin = 0.0, double xmax = 4.5)
    {
        TH1D *hist = new TH1D("hist", "Generazione dati; x; Conteggi", bins, xmin, xmax);
        for (int i = 0; i < N; ++i)
        {
            double x = functionf->GetRandom(xmin, xmax);
            hist->Fill(x);
        }

        return hist;
    }

    TH1D *ErHistoGeneratoreNormalized(int N, int bins, double xmin = 0.0, double xmax = 4.5)
    {
        TH1D *hist_normy = new TH1D("hist_normy", "Generazione dati; x; Conteggi", bins, xmin, xmax);
        for (int i = 0; i < N; ++i)
        {
            double x = functionf->GetRandom(xmin, xmax);
            hist_normy->Fill(x);
        }

        hist_normy->Scale(1.0 / hist_normy->Integral("width"));
        return hist_normy;
    }

    TF1 *GetFunctionNormalized() const
    {
        double integral = functionf->Integral(0, 4.5);
        return new TF1("f_normy", [this, integral](double *x, double *)
                       { return (std::pow(std::cos(k * x[0] + phi), 2) + b) / integral; }, 0, 4.5, 0);
    }

    double CarcolaRMSE(TH1D *hist, TF1 *functionf_normy, int N)
    {
        int B = hist->GetNbinsX();
        double rmse_sum = 0.0;

        for (int i = 1; i <= B; ++i)
        {
            double x_center = hist->GetBinCenter(i);
            double expected_density = functionf_normy->Eval(x_center);

            double observed_density = hist->GetBinContent(i);
            double diff = observed_density - expected_density;
            rmse_sum += diff * diff;
        }

        return std::sqrt(rmse_sum / B);
    }

    // stima incertezza su molti esperimenti
    std::pair<std::vector<double>, std::vector<double>> StimaIncertezza(int N, int bins, int M, double xmin = 0.0, double xmax = 4.5)
    {
        // matrice con bins righe e M colonne. le righe contengono le occorrenze di ciascun bin (da 1 a bins) mentre le colonne indicano la ripetizione (da 1 a M)
        std::vector<std::vector<double>> bin_values(bins);

        // Ripetere la generazione M volte per riempire la matrice
        for (int j = 0; j < M; ++j)
        {
            TH1D *hist_tmp = ErHistoGeneratore(N, bins, xmin, xmax);

            for (int i = 1; i <= bins; ++i)
            {
                bin_values[i - 1].push_back(hist_tmp->GetBinContent(i));
            }

            delete hist_tmp; // pulizia memoria
        }

        // vengono preparati due vettori di dimensione bins che verranno usati per salvare media e dev std di ogni bin
        std::vector<double> mean(bins, 0.0);
        std::vector<double> stddev(bins, 0.0);

        for (int i = 0; i < bins; ++i)
        {
            double sum = 0.0;
            double sum_sq = 0.0;
            int count = bin_values[i].size(); // dovrebbe essere M

            // Sommi tutti i valori di quel bin sulle M generazioni e anche i quadrati di quei valori per la varianza.
            for (double val : bin_values[i])
            {
                sum += val;
                sum_sq += val * val;
            }

            mean[i] = sum / count;
            double variance = (sum_sq / count) - (mean[i] * mean[i]);
            stddev[i] = (variance > 0) ? std::sqrt(variance) : 0.0;
        }

        return std::make_pair(mean, stddev);
    }

    TH1D *BinSmearing(int bins, double xmin = 0.0, double xmax = 4.5, double rel_sigma = 0.05)
    {
        // Creiamo l'istogramma teorico dalla funzione normalizzata
        TF1 *f_norm = GetFunctionNormalized();
        TH1D *h_theoretical = new TH1D("h_theoretical", "Funzione teorica fluttuata; x; f(x)", bins, xmin, xmax);

        TRandom3 rand(0); // inizializza generatore casuale

        for (int i = 1; i <= bins; ++i)
        {
            double x_center = h_theoretical->GetBinCenter(i);
            double expected_value = f_norm->Eval(x_center);

            // deviazione standard relativa (es. 5% del valore)
            double sigma = rel_sigma * expected_value;

            // fluttuazione gaussiana
            double fluctuated = rand.Gaus(expected_value, sigma);

            // Assicuriamoci che non vada sotto zero
            if (fluctuated < 0)
                fluctuated = 0;

            h_theoretical->SetBinContent(i, fluctuated);
        }

        // Normalizziamo l’istogramma
        h_theoretical->Scale(1.0 / h_theoretical->Integral("width"));
        delete f_norm;
        return h_theoretical;
    }

    std::pair<std::vector<double>, std::vector<double>> StimaIncertezzaBinSmearing(
        int bins, int M, double xmin = 0.0, double xmax = 4.5, double rel_sigma = 0.05)
    {
        std::vector<std::vector<double>> bin_values(bins);

        for (int j = 0; j < M; ++j)
        {
            TH1D *h_smear = BinSmearing(bins, xmin, xmax, rel_sigma);
            for (int i = 1; i <= bins; ++i)
            {
                bin_values[i - 1].push_back(h_smear->GetBinContent(i));
            }
            delete h_smear;
        }

        std::vector<double> mean(bins, 0.0);
        std::vector<double> stddev(bins, 0.0);

        for (int i = 0; i < bins; ++i)
        {
            double sum = 0.0, sum_sq = 0.0;
            int count = bin_values[i].size();
            for (double val : bin_values[i])
            {
                sum += val;
                sum_sq += val * val;
            }
            mean[i] = sum / count;
            double variance = (sum_sq / count) - (mean[i] * mean[i]);
            stddev[i] = (variance > 0) ? std::sqrt(variance) : 0.0;
        }

        return std::make_pair(mean, stddev);
    }

    std::tuple<
        std::vector<double>, std::vector<double>, // mean e stddev da rigenerazione
        std::vector<double>, std::vector<double>  // mean e stddev da bin-smearing
        >
    PropagaIncertezzeParametri(
        double k_sigma, double phi_sigma, double b_sigma, // incertezze sui parametri
        int N, int bins, int M_param,                     // N eventi, numero di bin, numero di campioni di parametri
        int M_stat = 100,                                 // ripetizioni statistiche per 3.2 e 3.3
        double xmin = 0.0, double xmax = 4.5, double rel_sigma = 0.05)
    {
        TRandom3 rand(0);

        // vettori per accumulare i risultati delle medie di ogni bin
        std::vector<std::vector<double>> mean_rig(bins), mean_sme(bins);

        // ciclo sui campioni di parametri
        for (int m = 0; m < M_param; ++m)
        {
            // genera parametri fluttuati
            double k_fluc = rand.Gaus(k, k_sigma);
            double phi_fluc = rand.Gaus(phi, phi_sigma);
            double b_fluc = rand.Gaus(b, b_sigma);

            // nuova funzione con parametri fluttuati
            FunzioneDistribuzione f_flutt(k_fluc, phi_fluc, b_fluc);

            // 3.2: incertezza da rigenerazione
            auto [mean_rig_tmp, std_rig_tmp] = f_flutt.StimaIncertezza(N, bins, M_stat, xmin, xmax);

            // 3.3: incertezza da bin-smearing
            auto [mean_sme_tmp, std_sme_tmp] = f_flutt.StimaIncertezzaBinSmearing(bins, M_stat, xmin, xmax, rel_sigma);

            // accumula i valori medi (per poi stimare varianza dovuta ai parametri)
            for (int i = 0; i < bins; ++i)
            {
                mean_rig[i].push_back(mean_rig_tmp[i]);
                mean_sme[i].push_back(mean_sme_tmp[i]);
            }
        }

        // calcolo media e stddev sulle M_param variazioni dei parametri
        std::vector<double> mean_rig_tot(bins, 0.0), std_rig_tot(bins, 0.0);
        std::vector<double> mean_sme_tot(bins, 0.0), std_sme_tot(bins, 0.0);

        for (int i = 0; i < bins; ++i)
        {
            double sum_r = 0.0, sum_r2 = 0.0;
            double sum_s = 0.0, sum_s2 = 0.0;
            int count = mean_rig[i].size();

            for (int j = 0; j < count; ++j)
            {
                sum_r += mean_rig[i][j];
                sum_r2 += mean_rig[i][j] * mean_rig[i][j];
                sum_s += mean_sme[i][j];
                sum_s2 += mean_sme[i][j] * mean_sme[i][j];
            }

            mean_rig_tot[i] = sum_r / count;
            mean_sme_tot[i] = sum_s / count;

            double var_r = (sum_r2 / count) - std::pow(mean_rig_tot[i], 2);
            double var_s = (sum_s2 / count) - std::pow(mean_sme_tot[i], 2);

            std_rig_tot[i] = (var_r > 0) ? std::sqrt(var_r) : 0.0;
            std_sme_tot[i] = (var_s > 0) ? std::sqrt(var_s) : 0.0;
        }

        return std::make_tuple(mean_rig_tot, std_rig_tot, mean_sme_tot, std_sme_tot);
    }
};

#ifndef FUNZIONEDISTRIBUZIONE2_HPP
#define FUNZIONEDISTRIBUZIONE2_HPP

class FunzioneDistribuzione2 {
public:
    FunzioneDistribuzione2() = default;

    TGraphErrors* ErGraficoConErrori(const std::vector<double>& mean,
                                       const std::vector<double>& stddev,
                                       double xmin, double xmax);

   
    void FitCoseno(TGraphErrors* gr, double b_fixed);

    // Funzione seno per TF1
    static double Coseno(double* x, double* par);
};

#endif
