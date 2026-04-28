#include <iostream>
#include <cmath>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>

double S_over_N(int N, int E) {
    if (E == N || E == -N) {
        return 0; // log(0)이 발생하지 않도록 처리
    }

    double term1 = N * std::log(N);
    double term2 = ((N + E) / 2.0) * std::log((N + E) / 2.0);
    double term3 = ((N - E) / 2.0) * std::log((N - E) / 2.0);

    return (term1 - term2 - term3) / N;
}

void draw_S_over_N_vs_E_over_N(int N=100) {
    const int points = 2 * N + 1; // E의 범위: -N에서 N까지
    double x[points], y[points];

    for (int i = 0, E = -N; E <= N; ++E, ++i) {
        x[i] = static_cast<double>(E) / N;  // E/N
        y[i] = S_over_N(N, E);              // S(N, E) / N
    }

    TCanvas *c = new TCanvas("c", "S(N,E)/N vs E/N", 800, 600);
    TGraph *graph = new TGraph(points, x, y);
    
    graph->SetTitle(Form("S(N,E)/N vs E/N for N=%d", N));
    graph->GetXaxis()->SetTitle("E/N");
    graph->GetYaxis()->SetTitle("S(N,E)/N");
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);
    
    graph->Draw("AL");
    c->Update();
    c->Draw();
}
