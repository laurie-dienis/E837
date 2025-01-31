#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TStyle.h>

void set_plot_style() {
  const Int_t NRGBs = 5;
  const Int_t NCont = 455;

  Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  Double_t red[NRGBs] = {0.00, 0.00, 0.87, 1.00, 0.51};
  Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
  Double_t blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void si_distance() {
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
   Double_t x[100], y[100];
  //Double_t x_data[10] = {5, 6, 7, 7.5, 8, 8.5, 8.75, 8.9, 9, 10};
  //Double_t y_data[10] = {0.1485, 0.1485, 0.1645,  0.1992, 0.0556,
  //               0.0331, 0.01014, -0.006, -0.01, -0.08};
  
  Double_t x_data[7] = {7.5, 8, 8.5, 8.75, 8.9, 9, 10};
  Double_t y_data[7] = {0.1992, 0.0556,
                 0.0331, 0.01014, -0.006, -0.01, -0.08};


   Int_t n = 7;
   for (Int_t i=0;i<n;i++) {
     x[i] = x_data[i];
     y[i] = y_data[i];
   }
    
   set_plot_style(); 
   auto gr = new TGraph(n,x,y);
    gr->GetXaxis()->SetTitle("Si distance [cm]");
    gr->GetYaxis()->SetTitle("Ex Mean [MeV]");
   gr->SetMarkerStyle(8); // Style de marqueur
   gr->SetMarkerSize(2); // Taille du marqueur
   gr->Draw("AP");


  TLine *line1 = new TLine(7.3, 0, 10.5, 0);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(3);
  line1->Draw("same");

}
