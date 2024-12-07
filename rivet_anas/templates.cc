// -*- C++ -*-
#include <iostream> //debug output
#include <fstream> //writing file
#include <cmath> // Include cmath for sqrt
#include <cstdlib> // std::exit()
#include <array>
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

	  ptbins = get_ptbins("/data/horse/ws/maml087d-workspace/rivet_anas/data/ptbinned/ptbins.txt");
	  coeffs = get_coeffs("/data/horse/ws/maml087d-workspace/rivet_anas/data/ptbinned/histo_nocutscoeffs.txt", static_cast<int>(ptbins.size()));
	  coeffs_of = get_coeffs_of("/data/horse/ws/maml087d-workspace/rivet_anas/data/ptbinned/histo_nocutscoeffs_of.txt");

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
	  book(_hist1, "costheta", 40, -1, 1);
	  book(_hist2, "phi", 40, 0, TWOPI);
	  book(_hist3, "costhetacheck", 40, -1, 1);
	  book(_hist4, "xsec_bin", 1, -0.5, 0.5);
	  //initplots();
	  init_ptbintemps(static_cast<int>(ptbins.size()), 8);
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
		// for (i == 0; i < event)
          fillPlots(vars, SRplots);
		  // double weight = 1 / dxsec(coeffs, vars.theta, vars.cos_theta, vars.phi);
		  if (vars.mjj > 1000 * GeV){
		  _hist1->fill(vars.cos_theta);
		  _hist2->fill(vars.phi);
		  _hist3->fill(vars.cos_thetacheck);
		  _hist4->fill(0);
		  // filltemps("temps", vars.theta, vars.cos_theta, vars.phi);
		  fill_ptbintemps(vars);
		  }
//		  std::cout << "fills_histo" << std::endl;

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
      //scaletemps("temps", xsec);
	  scale_ptbintemps(xsec);
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

		if (!isZero(yax.x()) || !isZero(yax.z())) std::cout << "y-ax-err" << std::endl;
		const Matrix<4> rot24 = _mkMatrix4(rot2);
		l_boost = multiply(rot24, l_boost);
	}
	else std::cout << "pT == 0" << std::endl;

	
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

    // templates in struct
    struct template_histos {
	    string label;
	    Histo1DPtr t0, t1, t2, t3, t4, t5, t6, t7, t8;
    };

    struct template_histos2D {
	    string label;
	    Histo2DPtr t0, t1, t2, t3, t4, t5, t6, t7, t8;
    };

    struct templates {
	    string label;
	    template_histos cos, phi;
	    template_histos2D cph;

    };
	void init_ptbintemps (const int nptbins, const int n) {
		for (int i=0; i<nptbins-1; ++i) {
			inittemps("ptbin" + to_string(i), n);
		}
		inittemps("ptoverflow", n);
	}

	void fill_ptbintemps (Variables& vars) {
		// fill overflow
		if (vars.pTll > ptbins.back()) {
			filltemps("ptoverflow", vars.theta, vars.cos_theta, vars.phi, coeffs_of);
		}
		// fill binned tempaltes
		else{
			size_t i = 0;
			while (vars.pTll > ptbins[i+1]) i += 1;
			if (i >= coeffs.size()){
				std::cout << std::endl <<"!!!!!!!!!i out of bounds; i = " << i << std::endl;
				std::exit(EXIT_FAILURE);
			}
			filltemps("ptbin" + to_string(i), vars.theta, vars.cos_theta, vars.phi, coeffs[i]);
		}
	}

	void scale_ptbintemps (const double xsec) {
		for (size_t i=0; i<coeffs.size(); ++i){
			scaletemps("ptbin" + to_string(i), xsec);
		}
		scaletemps("ptoverflow", xsec);
	}

    void inittemps (const string& name, const int& n) {
	    _t[name].label = name;
	    _t[name].cos.label = name + "_cos";
	    _t[name].phi.label = name + "_phi";
	    _t[name].cph.label = name + "_2D";
	    
	    double startcos = -1;
	    double endcos = 1;
	    double startphi = 0;
	    double endphi = TWOPI;

	    book(_t[name].cos.t0, _t[name].cos.label + "0", n, startcos, endcos);
	    book(_t[name].cos.t1, _t[name].cos.label + "1", n, startcos, endcos);
	    book(_t[name].cos.t2, _t[name].cos.label + "2", n, startcos, endcos);
	    book(_t[name].cos.t3, _t[name].cos.label + "3", n, startcos, endcos);
	    book(_t[name].cos.t4, _t[name].cos.label + "4", n, startcos, endcos);
	    book(_t[name].cos.t5, _t[name].cos.label + "5", n, startcos, endcos);
	    book(_t[name].cos.t6, _t[name].cos.label + "6", n, startcos, endcos);
	    book(_t[name].cos.t7, _t[name].cos.label + "7", n, startcos, endcos);
	    book(_t[name].cos.t8, _t[name].cos.label + "8", n, startcos, endcos);

	    book(_t[name].phi.t0, _t[name].phi.label + "0", n, startphi, endphi);
	    book(_t[name].phi.t1, _t[name].phi.label + "1", n, startphi, endphi);
	    book(_t[name].phi.t2, _t[name].phi.label + "2", n, startphi, endphi);
	    book(_t[name].phi.t3, _t[name].phi.label + "3", n, startphi, endphi);
	    book(_t[name].phi.t4, _t[name].phi.label + "4", n, startphi, endphi);
	    book(_t[name].phi.t5, _t[name].phi.label + "5", n, startphi, endphi);
	    book(_t[name].phi.t6, _t[name].phi.label + "6", n, startphi, endphi);
	    book(_t[name].phi.t7, _t[name].phi.label + "7", n, startphi, endphi);
	    book(_t[name].phi.t8, _t[name].phi.label + "8", n, startphi, endphi);

	    book(_t[name].cph.t0, _t[name].cph.label + "0", n, startcos, endcos, n, startphi, endphi);
	    book(_t[name].cph.t1, _t[name].cph.label + "1", n, startcos, endcos, n, startphi, endphi);
	    book(_t[name].cph.t2, _t[name].cph.label + "2", n, startcos, endcos, n, startphi, endphi);
	    book(_t[name].cph.t3, _t[name].cph.label + "3", n, startcos, endcos, n, startphi, endphi);
	    book(_t[name].cph.t4, _t[name].cph.label + "4", n, startcos, endcos, n, startphi, endphi);
	    book(_t[name].cph.t5, _t[name].cph.label + "5", n, startcos, endcos, n, startphi, endphi);
	    book(_t[name].cph.t6, _t[name].cph.label + "6", n, startcos, endcos, n, startphi, endphi);
	    book(_t[name].cph.t7, _t[name].cph.label + "7", n, startcos, endcos, n, startphi, endphi);
	    book(_t[name].cph.t8, _t[name].cph.label + "8", n, startcos, endcos, n, startphi, endphi);
    }
 /*
    void filltemps(const string& name, const double theta, const double cos, const double phi) {
	    double normweight;
	    array<double, 9> pols;
	    array<double, 9> weight;

	    normweight =  1 / dxsec(coeffs, theta, cos, phi);
	    pols = polis(coeffs, theta, cos, phi);
//	    std::cout << "pols: [";
	    for (int i=0; i<9; i++){
//		    std::cout << pols[i] << ", ";
		    weight[i] = normweight * pols[i];
	    }
//	    std::cout << "]" << std::endl;
//	    std::cout << "weights" << weight << std::endl;

	    _t[name].cos.t0->fill(cos, weight[0]);
	    _t[name].cos.t1->fill(cos, weight[1]);
	    _t[name].cos.t2->fill(cos, weight[2]);
	    _t[name].cos.t3->fill(cos, weight[3]);
	    _t[name].cos.t4->fill(cos, weight[4]);
	    _t[name].cos.t5->fill(cos, weight[5]);
	    _t[name].cos.t6->fill(cos, weight[6]);
	    _t[name].cos.t7->fill(cos, weight[7]);
	    _t[name].cos.t8->fill(cos, weight[8]);
	    
	    _t[name].phi.t0->fill(phi, weight[0]);
	    _t[name].phi.t1->fill(phi, weight[1]);
	    _t[name].phi.t2->fill(phi, weight[2]);
	    _t[name].phi.t3->fill(phi, weight[3]);
	    _t[name].phi.t4->fill(phi, weight[4]);
	    _t[name].phi.t5->fill(phi, weight[5]);
	    _t[name].phi.t6->fill(phi, weight[6]);
	    _t[name].phi.t7->fill(phi, weight[7]);
	    _t[name].phi.t8->fill(phi, weight[8]);

	    _t[name].cph.t0->fill(cos, phi, weight[0]);
	    _t[name].cph.t1->fill(cos, phi, weight[1]);
	    _t[name].cph.t2->fill(cos, phi, weight[2]);
	    _t[name].cph.t3->fill(cos, phi, weight[3]);
	    _t[name].cph.t4->fill(cos, phi, weight[4]);
	    _t[name].cph.t5->fill(cos, phi, weight[5]);
	    _t[name].cph.t6->fill(cos, phi, weight[6]);
	    _t[name].cph.t7->fill(cos, phi, weight[7]);
	    _t[name].cph.t8->fill(cos, phi, double theta, double cos, double phiweight[8]);
    }*/

    void filltemps(const string& name, const double theta, const double cos, const double phi, const array<double, 9>& coeffs) {
		double normweight;
		array<double, 9> pols;
		array<double, 9> weight;

		normweight =  1 / dxsec(coeffs, theta, cos, phi);
		pols = polis(coeffs, theta, cos, phi);
		//	    std::cout << "pols: [";
		for (int i=0; i<9; i++){
			//		    std::cout << pols[i] << ", ";
			weight[i] = normweight * pols[i];
		}
		//	    std::cout << "]" << std::endl;
		//	    std::cout << "weights" << weight << std::endl;

		_t[name].cos.t0->fill(cos, weight[0]);
		_t[name].cos.t1->fill(cos, weight[1]);
		_t[name].cos.t2->fill(cos, weight[2]);
		_t[name].cos.t3->fill(cos, weight[3]);
		_t[name].cos.t4->fill(cos, weight[4]);
		_t[name].cos.t5->fill(cos, weight[5]);
		_t[name].cos.t6->fill(cos, weight[6]);
		_t[name].cos.t7->fill(cos, weight[7]);
		_t[name].cos.t8->fill(cos, weight[8]);

		_t[name].phi.t0->fill(phi, weight[0]);
		_t[name].phi.t1->fill(phi, weight[1]);
		_t[name].phi.t2->fill(phi, weight[2]);
		_t[name].phi.t3->fill(phi, weight[3]);
		_t[name].phi.t4->fill(phi, weight[4]);
		_t[name].phi.t5->fill(phi, weight[5]);
		_t[name].phi.t6->fill(phi, weight[6]);
		_t[name].phi.t7->fill(phi, weight[7]);
		_t[name].phi.t8->fill(phi, weight[8]);

		_t[name].cph.t0->fill(cos, phi, weight[0]);
		_t[name].cph.t1->fill(cos, phi, weight[1]);
		_t[name].cph.t2->fill(cos, phi, weight[2]);
		_t[name].cph.t3->fill(cos, phi, weight[3]);
		_t[name].cph.t4->fill(cos, phi, weight[4]);
		_t[name].cph.t5->fill(cos, phi, weight[5]);
		_t[name].cph.t6->fill(cos, phi, weight[6]);
		_t[name].cph.t7->fill(cos, phi, weight[7]);
		_t[name].cph.t8->fill(cos, phi, weight[8]);
	}

    void scaletemps (const string& name, const double xsec){

	    scale(_t[name].cos.t0, xsec);
	    scale(_t[name].cos.t1, xsec);
	    scale(_t[name].cos.t2, xsec);
	    scale(_t[name].cos.t3, xsec);
	    scale(_t[name].cos.t4, xsec);
	    scale(_t[name].cos.t5, xsec);
	    scale(_t[name].cos.t6, xsec);
	    scale(_t[name].cos.t7, xsec);
	    scale(_t[name].cos.t8, xsec);

	    scale(_t[name].phi.t0, xsec);
	    scale(_t[name].phi.t1, xsec);
	    scale(_t[name].phi.t2, xsec);
	    scale(_t[name].phi.t3, xsec);
	    scale(_t[name].phi.t4, xsec);
	    scale(_t[name].phi.t5, xsec);
	    scale(_t[name].phi.t6, xsec);
	    scale(_t[name].phi.t7, xsec);
	    scale(_t[name].phi.t8, xsec);

	    scale(_t[name].cph.t0, xsec);
	    scale(_t[name].cph.t1, xsec);
	    scale(_t[name].cph.t2, xsec);
	    scale(_t[name].cph.t3, xsec);
	    scale(_t[name].cph.t4, xsec);
	    scale(_t[name].cph.t5, xsec);
	    scale(_t[name].cph.t6, xsec);
	    scale(_t[name].cph.t7, xsec);
	    scale(_t[name].cph.t8, xsec);
    }

//    void initplots(){
//	    for (int i = 0; i<134; i++) {
//		    string namec = "cos";
//		    string namep = "phi";
//		    namec += std::to_string(i);
//		    namep += std::to_string(i);
//		    book(_h[namec], namec, 40, -1, 1);
//		    book(_h[namep], namep, 40, 0, TWOPI);
//		    histcos[i] = namec;
//		    histphi[i] = namep;
//	//	    std::cout << "" << << std::endl;
//	    }
//    }

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

    static double dxsec(const array<double,9>& coeffs, double theta, double cos_theta, double phi) {
	    double P8 = 1 + std::pow(cos_theta,2);
	    double P0 = (coeffs[0] * (1 - 3 * std::pow(cos_theta, 2))) / 2;
	    double P1 = coeffs[1] * sin(2 * theta) * cos(phi);
	    double P2 = coeffs[2] * std::pow(sin(theta),2) * cos(2 * phi) / 2;
	    double P3 = coeffs[3] * sin(theta) * cos(phi);
	    double P4 = coeffs[4] * cos_theta;
	    double P5 = coeffs[5] * std::pow(sin(theta), 2) * sin(2 * phi);
	    double P6 = coeffs[6] * sin(2 * theta) * sin(phi);
	    double P7 = coeffs[7] * sin(theta) * sin(phi);
	    double norm = coeffs[8];
	    
	    double dsig = 3/(16 * PI) * norm * (P0 + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8);
	    return dsig;
    }

    static array<double, 9> polis(const array<double,9>& coeffs, const double theta, const double cos_theta, const double phi) {
	    double P8 = 1 + std::pow(cos_theta,2);
	    double P0 = (coeffs[0] * (1 - 3 * std::pow(cos_theta, 2))) / 2;
	    double P1 = coeffs[1] * sin(2 * theta) * cos(phi);
	    double P2 = coeffs[2] * std::pow(sin(theta),2) * cos(2 * phi) / 2;
	    double P3 = coeffs[3] * sin(theta) * cos(phi);
	    double P4 = coeffs[4] * cos_theta;
	    double P5 = coeffs[5] * std::pow(sin(theta), 2) * sin(2 * phi);
	    double P6 = coeffs[6] * sin(2 * theta) * sin(phi);
	    double P7 = coeffs[7] * sin(theta) * sin(phi);
	    double norm = coeffs[8] * 3/(16 * PI);
	    array<double, 9> pols{{P0 * norm, P1 * norm, P2 * norm, P3 * norm, P4 * norm, P5 * norm, P6 * norm, P7 * norm, P8 * norm}};

	    return pols;
    }
    static int countlines(const string& filename){
		std::ifstream file; file.open(filename);
		int i = 0;
		string line;
		while (std::getline(file, line)){
			i += 1;
		}
		file.close();
		return i;
	}

	static std::vector<double> get_ptbins(const string& filename) {
		vector<double> ptbins;
		std::ifstream f; f.open(filename);
		if (f.is_open()) {
			string line;
			while (std::getline(f, line)){
				ptbins.push_back(std::stod(line));
			}
		}
		else {
			std::cout << "couldnt open " << filename << " -> terminating" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		f.close();
		return ptbins;
	}
	static vector<array<double, 9>> get_coeffs(const string& filename, const int nptbins) {
		std::ifstream f; f.open(filename);
		string line; vector<array<double, 9>> coeffs;

		if (f.is_open()) {
			coeffs.resize(nptbins-1);
			for (int i = 0; i<9; i++){
				for (int j = 0; j < nptbins-1; ++j){
					std::getline(f, line);
					coeffs[j][i] = std::stod(line);
				}
			}
			f.close();
			return coeffs;
		}
		else {
			std::cout << "couldnt open " << filename << " -> terminating" << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	static array<double, 9> get_coeffs_of(const string& filename) {
		std::ifstream f; f.open(filename);
		string line; array<double, 9> coeffs_of;
		if (f.is_open()) {
			for (int i = 0; i<9; i++){
				std::getline(f, line);
				coeffs_of[i] = stod(line);
			}
			f.close();
			return coeffs_of;
		}
		else {
			std::cout << "couldnt open " << filename << " -> terminating" << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}


    /// @}


    private:

      // plot declarations
      Plots SRplots, CRAplots, CRBplots, CRCplots;
      Histo1DPtr _hist1, _hist2, _hist3, _hist4;
      map<string, Histo1DPtr> _h;
      map<string, templates> _t;
            
      bool _doControl;

      // file handling + angular_coefficients
      // std::ofstream file;
      vector<array<double, 9>> coeffs;
	  array<double, 9> coeffs_of;
	  vector<double> ptbins;

  };


  
  RIVET_DECLARE_PLUGIN(ATLAS_2020_I1803608);

}
