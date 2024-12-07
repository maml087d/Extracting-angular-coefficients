// -*- C++ -*-
#include <iostream> //debug output
#include <fstream> //writing file
#include <cmath> // Include cmath for sqrt
#include <cstdlib> // std::exit() and system()
#include <array> // arrays
#include <vector> // vectors
#include <numeric> // std::accumulate

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
  

  /// VBFZ in pp at 13 TeV
  class ATLAS_2020_I1803608 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2020_I1803608);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // for fileoutput
      // file.open("/data/horse/ws/maml087d-workspace/data/data.json", std::ios::trunc);
      
      // check if file is open
      // if (!(file.is_open())){
        // std::cout << "couldnt open file -> terminating" << std::endl;
        // std::exit(EXIT_FAILURE);
      // }

      // file << "[\n";
      // file << "[\"theta\", \"phi\"]";

      FinalState fs(Cuts::abseta < 5.0);

      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON);
      PromptFinalState muons(Cuts::abspid == PID::MUON);

      bool nocuts = getOption("CUTS") == "NO";
      if (nocuts) {
        std::cout << "no cuts" << std::endl;
        DressedLeptons dressed_electrons(photons, electrons, 0.1);
        declare(dressed_electrons, "DressedElectrons");

        DressedLeptons dressed_muons(photons, muons, 0.1);
        declare(dressed_muons, "DressedMuons");
      }
      else {
        std::cout << "cuts" << std::endl;
        Cut cuts_el = (Cuts::pT > 25*GeV) && ( Cuts::abseta < 1.37 || (Cuts::abseta > 1.52 && Cuts::abseta < 2.47) );
        Cut cuts_mu = (Cuts::pT > 25*GeV) && (Cuts::abseta < 2.4);

        DressedLeptons dressed_electrons(photons, electrons, 0.1, cuts_el);
        declare(dressed_electrons, "DressedElectrons");

        DressedLeptons dressed_muons(photons, muons, 0.1, cuts_mu);
        declare(dressed_muons, "DressedMuons");

      }

      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "Jets");

      _doControl = bool(getOption("TYPE") != "EW_ONLY");
      if (_doControl) {
        initialisePlots(SRplots, "SR");
        initialisePlots(CRAplots, "CRA");
        initialisePlots(CRBplots, "CRB");
        initialisePlots(CRCplots, "CRC");
      } else {
        initialisePlots(SRplots, "EW");
      }
      int bins = 8;
      ptbins = arange(20.0, 200.0, 5.0);
      vector<double> addptbin = arange(225.0, 400.0, 25.0);
      ptbins.insert(ptbins.end(), addptbin.begin(), addptbin.end());
      // ptbins = {20,30,45,70,100,140,200,275,400,550,1050};
      init_data(bins);
      save_pt(ptbins, "data/ptbinned/ptbins.txt");
      // save_pt(ptbins, "debugdata/ptbins.txt");
      // std::cout << "dataobject size"<< dat.size() << std::endl;
      book(_hist1, "costheta", bins, -1, 1);
      book(_hist2, "phi", bins, 0, TWOPI);
      book(_hist3, "costhetacheck", bins, -1, 1);
      book(_hist4, "xsec_bin", 1, -0.5, 0.5);
      book(_hist5, "data_2d", bins, -1, 1, bins, 0, TWOPI);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Access fiducial electrons and muons
      const Particle *l1 = nullptr, *l2 = nullptr;
      Particles muons = apply<DressedLeptons>(event, "DressedMuons").particles();
      Particles elecs = apply<DressedLeptons>(event, "DressedElectrons").particles();

      // Dilepton selection 1: =2 leptons of the same kind
      if (muons.size()+elecs.size() != 2) vetoEvent;
      if      (muons.size()==2) { l1=&muons[0]; l2=&muons[1]; }
      else if (elecs.size()==2) { l1=&elecs[0]; l2=&elecs[1]; }
      else vetoEvent;

      // Dilepton selection 2: oppostie-charge and in mass range
      if ( !oppCharge(*l1, *l2) )  vetoEvent;
      if ( !inRange((l1->mom()+l2->mom()).mass()/GeV, 81.0, 101.0) ) vetoEvent;

      // Electron-jet overlap removal (note: muons are not included in jet finding)
      // make sure jets do not overlap with an electron within DR<0.2
      Jets jets;
      for (const Jet& j : apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::absrap < 4.4)) {
        if (elecs.size() == 2 && (deltaR(j, *l1, RAPIDITY) < 0.2 || deltaR(j, *l2, RAPIDITY) < 0.2 )) {
          continue;
        }
        jets += j;
      }

      // Require 2 jets with pT > 85 and 80 GeV
      if (jets.size() < 2) vetoEvent;

      // Calculate the observables
      Variables vars(jets, l1, l2);

      // make sure neither lepton overlaps with a jet within 0.4
      for (const Jet& j : jets) {
        if (deltaR(j, *l1, RAPIDITY) < 0.4 || deltaR(j, *l2, RAPIDITY) < 0.4)  vetoEvent;
      }

      if (jets[0].pt() < 85*GeV || jets[1].pt() < 80*GeV ) vetoEvent;

      bool pass_VBFZtopo = (vars.mjj > 250*GeV && vars.Dyjj > 2.0 && vars.pTll > 20*GeV && vars.Zcent < 1.0 && vars.pTbal < 0.15);

      if (pass_VBFZtopo) {
        if      (_doControl && vars.Ngj  > 0 && vars.Zcent <  0.5) fillPlots(vars, CRAplots);
        else if (_doControl && vars.Ngj  > 0 && vars.Zcent >= 0.5) fillPlots(vars, CRBplots);
        else if (_doControl && vars.Ngj == 0 && vars.Zcent >= 0.5) fillPlots(vars, CRCplots);

        if ( vars.Ngj == 0 && vars.Zcent < 0.5 ) {
          fillPlots(vars, SRplots);
          if (vars.mjj > 1000*GeV){
		  _hist1->fill(vars.cos_theta);
		  _hist2->fill(vars.phi);
		  _hist3->fill(vars.cos_thetacheck);
		  _hist4->fill(0);
		  _hist5->fill(vars.cos_theta, vars.phi);
          fill_data(vars);
          }
		  // file << ",\n[" << vars.theta << ", " << vars.phi << "]";
        }
      }
    }


    void finalize() {
      // file << "]";
      // file.close();
      const double xsec = crossSectionPerEvent()/femtobarn;
      scalePlots(SRplots, xsec);
      scalePlots(CRAplots, xsec);
      scalePlots(CRBplots, xsec);
      scalePlots(CRCplots, xsec);
      scale(_hist1, xsec);
      scale(_hist2, xsec);
      scale(_hist3, xsec);
      scale(_hist4, xsec);
      scale(_hist5, xsec);
      scaledata(xsec);
    }

    /// @}


    /// @name Analysis helpers
    /// @{

    struct Variables {
      Variables(const vector<Jet>& jets, const Particle* l1, const Particle* l2) {
        // get the jets
        assert(jets.size()>=2);
        FourMomentum j1 = jets[0].mom(), j2 = jets[1].mom();
        pTj1 = j1.pT()/GeV; pTj2 = j2.pT()/GeV;
        assert(pTj1 >= pTj2);
        
        // build dilepton system
        FourMomentum ll = (l1->mom() + l2->mom());
        pTll = ll.pT(); mll = ll.mass();
        
        Nj = jets.size();
        Dyjj = std::abs(j1.rap() - j2.rap());
        mjj = (j1 + j2).mass();
        Dphijj = ( j1.rap() > j2.rap() ) ? mapAngleMPiToPi(j1.phi() - j2.phi()) : mapAngleMPiToPi(j2.phi() - j1.phi());
        
        Jets gjets = getGapJets(jets);
        Ngj = gjets.size();
        pTgj = Ngj? gjets[0].pT()/GeV : 0;
        
        FourMomentum vecSum = (j1 + j2 + l1->mom() + l2->mom());
        double HT = j1.pT() + j2.pT() + l1->pT() + l2->pT();
        if (Ngj) { 
          vecSum += gjets[0].mom(); 
          HT += pTgj;
        }
        pTbal = vecSum.pT() / HT;

        Zcent = std::abs(ll.rap() - (j1.rap() + j2.rap())/2) / Dyjj;

	// CS-variables:
	
	// Select the negatively charged lepton
	FourMomentum p_l = (l1->charge() < 0) ? l1->mom() : l2->mom();

	// checking with analytical formula
	FourMomentum p_lp = (l1->charge() > 0) ? l1->mom() : l2->mom();
	cos_thetacheck = costhetacheck(ll, p_l, p_lp);

	// rotation or not
	const bool z_sgn_pos = (ll.eta() > 0) ? true : false;
//	std::cout << "pt_ll = " << ll.pT() << std::endl;
//	std::cout << "pz_ll = " << ll.pz() << std::endl;
	const bool pt_null = isZero(ll.pT());

	// Lorentz transformation to CS-frame with 2 pure boosts
        LorentzTransform CS_boost1 = LorentzTransform::mkFrameTransformFromBeta(Vector3(0, 0, ll.pz() / ll.E()));
        ll = CS_boost1.transform(ll);
        LorentzTransform CS_boost2 = LorentzTransform::mkFrameTransformFromBeta(ll.pTvec() / ll.E());
	ll = CS_boost2.transform(ll);
//	std::cout << "pt_ll_boost= " << ll.pT() << std::endl;
//	std::cout << "pz_ll_boost = " << ll.pz() << std::endl;
	
	// boosting lepton to CS-frame
	FourMomentum l_boost = CS_boost2.transform(CS_boost1.transform(p_l));
	


	// rotating around y axis by pi in case of negative eta
	Matrix<4> rotmat; 
	if (!z_sgn_pos) {
		Vector<4> dia;
		dia.set(0, 1); dia.set(1, -1); dia.set(2, 1); dia.set(3, -1);
		rotmat = Matrix<4>::mkDiag(dia);
		l_boost = multiply(rotmat, l_boost);
	}

	// rotating around z axis -> y axis normalvector of the beam plane in rest frame
	if (!pt_null) {
		FourMomentum p1 = FourMomentum::mkXYZE(0,0, 6500, 6500);
		FourMomentum p2 = FourMomentum::mkXYZE(0,0,-6500, 6500);
		
		p1 = CS_boost2.transform(CS_boost1.transform(p1));
		p2 = CS_boost2.transform(CS_boost1.transform(p2));

//		std::cout << "p1 boost = " << p1 << std::endl;
//		std::cout << "p2 boost = " << p2 << std::endl;

		if (!z_sgn_pos) {
			p1 = multiply(rotmat, p1);
			p2 = multiply(rotmat, p2);
		}

		Vector3 p1p = p1.vector3();
		Vector3 p2p = p2.vector3();

		Matrix3 rot2; 
		Vector3 cross = (p1p.cross(p2p)).unit();
		rot2 = Matrix3(cross, Vector3::mkY());
		Vector3 newp1p = multiply(rot2, p1p.unit());
		Vector3 newp2p = multiply(rot2, p2p.unit());
		Vector3 yax = newp1p.cross(newp2p);
//		std::cout << "cross" << cross << std::endl;
		if (!isZero(yax.x()) || !isZero(yax.z())) std::cout << "y-ax-err" << std::endl;
		const Matrix<4> rot24 = _mkMatrix4(rot2);
		l_boost = multiply(rot24, l_boost);
	}
	else std::cout << "pT == 0" << std::endl;

//	std::cout << "-----------------------------------------------" << std::endl << std::endl;
	
	// extracting CS-frame variables
	theta = l_boost.theta();
	cos_theta = cos(theta);
	phi = l_boost.phi();
      }

      double Zcent, pTj1, pTj2, pTgj, pTll, mll, Dyjj, mjj, Dphijj, pTbal, cos_theta, phi, cos_thetacheck, theta;
      size_t Nj, Ngj;

      Jets getGapJets(const Jets& jets) {
        Jets gjets;
        if (jets.size() <= 2)  return gjets;
        FourMomentum j1 = jets[0].mom(), j2 = jets[1].mom();
        double yFwd = j1.rap(), yBwd = j2.rap();
        if (yFwd < yBwd) std::swap(yFwd,yBwd);
        for (size_t i = 2; i < jets.size(); ++i)
         if (inRange(jets[i].rap(), yBwd, yFwd)) gjets += jets[i];
        return gjets;
      }

    }; // struct variables


    struct Plots {
      string label;
      Histo1DPtr m_jj, Dphi_jj, Dy_jj, pT_ll;
    };


    void initialisePlots(Plots& plots, const string& phase_space) {
      plots.label   = phase_space;
      size_t region = 0;
      if (phase_space == "SR")   region = 4;
      if (phase_space == "CRA")  region = 8;
      if (phase_space == "CRB")  region = 12;
      if (phase_space == "CRC")  region = 16;
      book(plots.m_jj,    1 + region, 1, 1);
      book(plots.Dy_jj,   2 + region, 1, 1);
      book(plots.pT_ll,   3 + region, 1, 1);
      book(plots.Dphi_jj, 4 + region, 1, 1);
    }


    void fillPlots(const Variables& vars, Plots& plots) {
      // The mjj plot extends down to 250 GeV
      plots.m_jj->fill(vars.mjj/GeV);
      if (vars.mjj > 1000*GeV) {
        plots.Dy_jj->fill(vars.Dyjj);
        plots.Dphi_jj->fill(vars.Dphijj);
        plots.pT_ll->fill(vars.pTll/GeV);
      }
    }


    void scalePlots(Plots& plots, const double xsec) {
      scale(plots.m_jj, xsec);
      scale(plots.Dy_jj, xsec);
      scale(plots.Dphi_jj, xsec);
      scale(plots.pT_ll, xsec);
    }

    struct data{
      string label;
      std::array<double, 2> xedges;
      Histo1DPtr cos, phi, xsec;
        Histo2DPtr cph;
    };

    //void init_data(std::vector<double>& ptbins, int bins){
    void init_data(const int bins){
      // initialize data histos
      // static std::array<data, ptbins.size()-1> dat;
      dat.resize(ptbins.size()-1);
      for (size_t i = 0; i < dat.size(); ++i) {
      // for (data &dat)
          dat[i].label = "datbin" + std::to_string(i);
            dat[i].xedges = {ptbins[i], ptbins[i+1]};
            book(dat[i].cos, dat[i].label + "_cos", bins, -1, 1);
            book(dat[i].phi, dat[i].label + "_phi", bins, 0, TWOPI);
            book(dat[i].xsec, dat[i].label + "_xsec", 1, -0.5, 0.5);
            book(dat[i].cph, dat[i].label + "_2D" , bins, -1, 1, bins, 0, TWOPI);
      }

      // initialize TProfile for angcs
      for (int i = 0; i < 8; ++i) {
          book(_tpfs[i], "A_" + to_string(i), ptbins);
          book(_tpf_of[i], "overflow_ang" + to_string(i), 1, 0, 1);
      }

      // initialize overflow
      book(_hist_of1D_cos, "overflow_cos", bins, -1, 1);
      book(_hist_of1D_phi, "overflow_phi", bins, 0, TWOPI);
      book(_hist_of1D_xsec, "overflow_xsec", 1, -0.5, 0.5);
      book(_hist_of_2D, "overflow_2D", bins, -1, 1, bins, 0, TWOPI);
    }

    void fill_data(Variables& vars){
      // std::cout << "datsize: " << dat.size();
      std::array<double, 8> cs = ang_coeffs(vars.theta, vars.phi);
      // for (size_t i=0; i<cs.size(); ++i){
        // std::cout <<  "A_" << i << " = " << cs[i] << "; ";
      // }
      // std::cout << std::endl << "theta = " << vars.theta << "; phi = " << vars.phi << std::endl;
      // A0.push_back(cs[0]);
      // std::cout << "pT = " << vars.pTll << std::endl;
      if (vars.pTll > ptbins.back()) {
        // filling overflow
        // std::cout << "filling overflow .... ";
        _hist_of1D_cos->fill(vars.cos_theta);
        _hist_of1D_phi->fill(vars.phi);
        _hist_of_2D->fill(vars.cos_theta, vars.phi);
        _hist_of1D_xsec->fill(0);
        for (size_t i=0; i<8; ++i) {
          _tpf_of[i]->fill(0.5, cs[i]);
        }
        // std::cout << "filled overlow" << std::endl;
      }

      else {
          // fill binned dat object
          size_t i = 0;
          // std::cout << "filling data object: ";
          // std::cout << "ptbin[" << i << "] = " << ptbins[i];
          while (vars.pTll > ptbins[i+1]) {
            i += 1;
            // std::cout << "ptbin[" << i << "] = " << ptbins[i];
          }  // find right ptbin
          // std::cout << "i = " << i;
          if (i >= dat.size()){
            std::cout << std::endl <<"!!!!!!!!!i out of bounds; i = " << i << std::endl;
            std::exit(EXIT_FAILURE);
          }
          dat[i].cos->fill(vars.cos_theta);
          dat[i].phi->fill(vars.phi);
          dat[i].cph->fill(vars.cos_theta, vars.phi);
          dat[i].xsec->fill(0);
          // std::cout << " filled dataobject"<< std::endl;

          // fill TProfiles for each angular coeff
          for (int j = 0; j < 8; ++j){
              _tpfs[j]->fill(vars.pTll, cs[j]);
          }
      }
    }

    void scaledata(const double xsec){
      // scaling data
      for (size_t i=0; i<dat.size(); ++i) {
        scale(dat[i].cos, xsec);
        scale(dat[i].phi, xsec);
        scale(dat[i].cph, xsec);
        scale(dat[i].xsec, xsec);
      }
      // scaling overflow
      scale(_hist_of1D_cos, xsec);
      scale(_hist_of1D_phi, xsec);
      scale(_hist_of_2D, xsec);
      scale(_hist_of1D_xsec, xsec);
    }

    static double costhetacheck(const FourMomentum& ll, const FourMomentum& l1, const FourMomentum& l2) {

      double p1p, p2p, p1m, p2m, denom, cos;

      p1p = (l1.E() + l1.pz())/SQRT2;
      p2p = (l2.E() + l2.pz())/SQRT2;
      p1m = (l1.E() - l1.pz())/SQRT2;
      p2m = (l2.E() - l2.pz())/SQRT2;

      denom = sqrt(ll.mass2()+ ll.pt2())*ll.mass();

      cos = 2 * sign(ll.pz()) * (p1p * p2m - p1m * p2p) / denom;

      return cos;
    }

    static Matrix4 _mkMatrix4(const Matrix3& m3) {
      Matrix4 m4 = Matrix4::mkIdentity();
      for (size_t i = 0; i < 3; ++i) {
              for (size_t j = 0; j < 3; ++j) {
                  m4.set(i+1, j+1, m3.get(i, j));
              }
          }
          return m4;
    }

    static std::array<double, 8> ang_coeffs(const double theta, const double phi){
      double c = cos(theta);
      double a0 = ((1.0 - 3.0 * std::pow(c, 2)) / 2.0) * 20.0/3.0 + 2.0 / 3.0;
      double a1 = sin(2.0 * theta) * cos(phi) * 5.0;
      double a2 = std::pow(sin(theta), 2) * cos(2.0 * phi) * 10.0;
      double a3 = sin(theta) * cos(phi) * 4.0;
      double a4 = 4.0 * c;
      double a5 = std::pow(sin(theta), 2) * sin(2.0*phi) * 5.0;
      double a6 = 5.0 * sin(2 * theta) * sin(phi);
      double a7 = 4.0 * sin(theta) * sin(phi);

      return {a0, a1, a2, a3, a4, a5, a6, a7};
    }

    static std::vector<double> arange(const double start, const double stop, double step) {
      std::vector<double> result;
      for (double value = start; value < stop+0.00001; value += step) {
        result.push_back(value);
      }
      return result;
    }

    static void save_pt(std::vector<double>& ptbins, const string& filename){
      std::ofstream file; file.open(filename);
      if (!(file.is_open())){
        std::cout << "couldnt open" << filename <<  "-> terminating" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      for (size_t i=0; i<ptbins.size(); ++i){
        file << ptbins[i] << "\n";
      }
    }
 
    /// @}


    private:

    Plots SRplots, CRAplots, CRBplots, CRCplots;
    Histo1DPtr _hist1, _hist2, _hist3, _hist4, _hist_of1D_cos, _hist_of1D_phi, _hist_of1D_xsec;
    Histo2DPtr _hist5, _hist_of_2D;
    std::array<Profile1DPtr, 8> _tpf_of;
    std::array<Profile1DPtr, 8> _tpfs;
    std::vector<double> ptbins;
    std::vector<data> dat;
    // std::vector<double> A0;

    bool _doControl;
    // std::ofstream /*file*/;

  };


  
  RIVET_DECLARE_PLUGIN(ATLAS_2020_I1803608);

}
