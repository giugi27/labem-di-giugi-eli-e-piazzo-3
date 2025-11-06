#include "FunzioneDistribuzione2.hpp"
#include "TApplication.h"
#include <fstream>
#include <iostream>


int main(int argc, char **argv)
{
    TApplication app("app", &argc, argv); // per usare ROOT grafico

    //Creiamo la funzione
    double k = 5.2, phi = 1.8, b = 0.2;
    FunzioneDistribuzione f(k, phi, b);

    //Stimiamo le incertezze da rigenerazione
    int N = 5000, bins = 50, M_stat = 100;
    auto [mean, stddev] = f.StimaIncertezza(N, bins, M_stat, 0.0, 4.5);

    //Stampiamo le incertezze
    std::ofstream file("incertezze_rigenerazione.txt");
    file << "# Bin  Mean  StdDev\n";
    for (int i = 0; i < bins; ++i)
        file << i+1 << "  " << mean[i] << "  " << stddev[i] << "\n";
    file.close();

    std::cout << "Incertezze salvate in 'incertezze_rigenerazione.txt'\n";

    //creiamo grafico con errori usando la classe
    FunzioneDistribuzione2 f_distrib;
    TGraphErrors* gr = f_distrib.ErGraficoConErrori(mean, stddev, 0.0, 4.5);

    double b_fixed = 0.2;  
    f_distrib.FitCoseno(gr, b_fixed);

    app.Run();
    return 0;
}
