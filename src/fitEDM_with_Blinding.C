#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "Blinders.hh"

using std::cout;
using std::endl;
using namespace blinding;

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
//Blinding string - change this to something private if working on data. We blind run 1,2,3 all separately. 
blinding::Blinders getBlinded( ftype, "this is a test blinding string" ); 

double blinded_wiggle(double *x, double *p){
//g-2 blinded wiggle function for fitting
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double time = x[0];
    double omega = getBlinded.paramToFreq(R);
    return norm * exp(-time/life) * (1 - asym*cos(omega*time + phi));
}

double blinded_EDM(){
//EDM blinding calculator based on the same string as the g-2 blinding and a central reference value
    double du = 1.9e-19;
    double ppmshift = 4.81; //it's i^i^i^i :D
    double omega_diff = ((getBlinded.paramToFreq(ppmshift) / getBlinded.referenceValue()) - 1) / 1e-6;
    return omega_diff*du;
}

double GetDelta(double dMu) {
//convert between the du value and the tilt angle
    double aMu = 11659208.9e-10;
    double mMu = 105.6583715; // u
    double mMuKg = mMu * 1.79e-30; // kg
    double c = 299792458.; // m/s
    double cm2m = 100.0; // cm -> m
    double hbar = 1.05457e-34;
    double gmagic = std::sqrt( 1.+1./aMu );
    double beta   = std::sqrt( 1.-1./(gmagic*gmagic) );
    double alpha = 0.10; //0.13; // asymmetry factor
    double eta = ((4 * mMuKg * c * dMu)/ (hbar * cm2m) );
  
    double tan_delta = (eta * beta) / (2 * aMu);
    double delta = atan(tan_delta);
    return delta;
}

double EDMFunc( double *x, double *p )  {
//EDM function for injecting a large blinded oscillation
    double time = x[0];
    return (-p[0] * cos(p[1]* time + p[2]));
}

TGraphErrors *BlindedModulo(TProfile* gr_thetaY_mod, TF1 *blindEDMFunc) {
//combination function to add the blinding to the data
//Assumes you've used a TProfile for the data - different functions for TGraphErrors
    int n = gr_thetaY_mod->GetNbinsX();
    std::vector<double> x,y,ex,ey;

    for (int i(0); i<n; i++) {
      double time = gr_thetaY_mod->GetBinCenter(i+1);
      double theta_y = gr_thetaY_mod->GetBinContent(i+1);
      double theta_y_shift = blindEDMFunc->Eval(time);
      
      x.push_back(time);
      ex.push_back(0);
      y.push_back(theta_y + theta_y_shift);
      ey.push_back(gr_thetaY_mod->GetBinError(i+1));
   }
    return new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
}

template<class VAR>
void Enable(TTree* tree, TString name, VAR &var) { 
    std::cout<<"Activating branch "<<name<<std::endl;
    var=0;
    tree->SetBranchStatus(name,1);
    tree->SetBranchAddress(name, &var);
}



//void fitEDM_with_Blinding (){
int main(){
        TCanvas* c1 = new TCanvas("c1", "Blinded wiggle", 0., 0., 1200, 800);
        auto file = TFile::Open("/scratch/gm2/Samples/simTree.dMu.5.4e-18.root");
     
        //get tracker tree with tracker wiggle
        TTree* tree = (TNtuple*)file->Get("trackerNTup/tracker");
   
 	TH1F* t = new TH1F("time", "wiggleplot", 1000, 30, 300);

        float tracktime, mom, posY, momx,momy,momz;
        bool QualityTrack, QualityVertex;
        Int_t n = tree->GetEntries();


	tree->SetBranchStatus("*",0); 
        Enable( tree, "decayTime", tracktime);
        Enable( tree, "trackMomentum", mom);
        Enable( tree, "decayVertexMomX", momx);
        Enable( tree, "decayVertexMomY", momy);
        Enable( tree, "decayVertexMomZ", momz);

        Enable( tree, "decayVertexPosY", posY);
        Enable( tree, "passTrackQuality", QualityTrack);
        Enable( tree, "passVertexQuality", QualityVertex);

        for (int i = 0; i < n; ++i) {
	  tree->GetEntry(i);
	  //high p cut for wiggle plot
	  if (mom>1800 and tracktime >= 30e3)  t->Fill(tracktime*1e-3);
	}
	
        //define fit for blinded omega_a
        TF1 *f1 = new TF1("f1", blinded_wiggle, 30,300,5);
        f1->SetParNames("N","#tau","A","R","#phi");
        f1->SetParameters(0.8e4,6.4e1,0.4,1e-6,0);
        f1->SetLineColor(2);
        f1->SetNpx(1000);

        t->Fit(f1); 
        t->Draw();
      
        t->SetXTitle("Time [#mus]");
        t->GetYaxis()->SetTitleOffset(1.5);
        t->SetYTitle("Counts/bin [/0.27 #mus]");

        gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
        gStyle->SetOptFit(1);

	double omega = f1->GetParameter(3);
        double omega_err = f1->GetParError(3);
        double phi = f1->GetParameter(4);
        double phi_err = f1->GetParError(4);

    	double aMu = 11659208.9e-10;
    	double gmagic = std::sqrt( 1.+1./aMu );
        double alpha = 0.1;      
       
        //get blinding tilt
        double delta_blind = GetDelta(blinded_EDM());
        //boost effect
        double tan_A_edm = tan(delta_blind) / gmagic;
        //asymmetry reduction
        double A_edm = alpha*atan(tan_A_edm) * 1e3; 

        double omega_a = getBlinded.referenceValue();
        double_t Tg2 = (2*TMath::Pi())/omega_a;

        //define blinded EDM fit function using BNL omega_a
        TF1 *blindEDMFunc = new TF1("blindEDMFunc",EDMFunc,0,Tg2,3);
        blindEDMFunc->SetParameters(A_edm,omega_a,phi+(TMath::Pi()/2));

        //make data EDM modulo plot
        TProfile* edm = new TProfile("vangle","",100,0,Tg2,-100,100);

        for (int i = 0; i < n; ++i) {
                tree->GetEntry(i);
                double xzmom = sqrt(pow(momx,2)+pow(momz,2));
                double v_angle = 1e3*atan(momy/xzmom);
                if (tracktime > 30e3 and QualityVertex == 1 and tracktime < 300e3) {
		double modTime = fmod((tracktime)*1e-3,Tg2);
                double shiftedModTime = 0;
                if (modTime >= Tg2/2) shiftedModTime = modTime - Tg2/2;
                if (modTime < Tg2/2) shiftedModTime = modTime + Tg2/2;
		edm->Fill(shiftedModTime,(v_angle));
		}
	}

        //blind the data modulo plot
        TGraphErrors *blindedEDMfit = BlindedModulo(edm, blindEDMFunc);
     
      
        TF1 *f2 = new TF1("f2","[0]*sin([1]*x+[2])+[3]",0,5);
        f2->SetParameter(0,0.1);
        f2->SetParameter(3,-0.1);
	f2->FixParameter(1,omega_a);
        f2->FixParameter(2,phi-(TMath::Pi())); 
   
        blindedEDMfit->Fit("f2");

	TCanvas* c2 = new TCanvas("c2", "Volumes Hit", 0., 0., 1200, 800);
	blindedEDMfit->Draw("P");
        blindedEDMfit->SetTitle(";Time % g-2 [#mus];Average vertical angle [mrad];");
        blindedEDMfit->GetXaxis()->SetLimits(0,Tg2);
        blindedEDMfit->SetMarkerSize(3);
        blindedEDMfit->Draw();

        
        c2->Update();

	TFile * output_file = new TFile("output.root", "RECREATE");
	output_file->cd();

	c1->Write();
	c2->Write();

	cout<<"Written file "<<output_file->GetName()<<endl;
	output_file->Close();
	
	return 1;
	

}
