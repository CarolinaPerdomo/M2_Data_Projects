#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TDatime.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "TString.h"
#include "TPaveText.h" // Add it by me
#include "TPaveLabel.h" // Add it by me
#include "TFile.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TLine.h"
#include "TPave.h" // Add it by me
#include "TMultiGraph.h" // <--
#include "TGraph.h" // <--
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TVirtualPad.h"
#include "TAxis.h"
#include "TGaxis.h"

#include <iostream> // <--
#include <fstream> // <--

// pT bins

const Int_t nPtBins = 16 ;
Float_t pTbins[nPtBins] = { 25 , 30 , 40 , 45 , 50 , 55 , 65 , 75 , 85 , 105 , 125 , 145 , 205 , 300 , 500 , 999 } ; //pTbins	

// eta bins

const Int_t nEtaBins = 6 ;
Float_t etaBins[nEtaBins] = { 0 , 0.6 , 1.37 , 1.52 , 1.81 , 2.37 } ; //absoetabins

void studySinglePhotons ( 
       
        Float_t _whichFraction = 1 
){

	cout << "Welcome to studySinglePhotons!" << endl ;

	cout << "Will use the following arguments : " << endl ;
        cout << "\t\t\t_whichFraction = " << _whichFraction << endl ; 

	cout << "Will  declare the variables from the SinglePhoton tree that will be analysed. " << endl ;
        Float_t y_pt ;
        Float_t y_eta ;
        Float_t y_phi ;
        Int_t y_convType ;
        Float_t y_topoetcone40 ;
        Bool_t y_IsTight ;
        Bool_t y_IsLoosePrime4 ;

	// Which conv
	
	const Int_t convorunconv = 3;
        const char *Conv[convorunconv] ;
        Conv[0]="conv" ;
        Conv[1]="uncov" ;
        Conv[2]="all" ;
	
	// Data or MC

	const Int_t dataormc = 3 ;
        const char *Type[dataormc] ;
        Type[0]="data" ;
        Type[1]="sherpa" ;
        Type[2]="pythia" ;

	cout << "Create histograms ... " << endl ;
 
        Float_t minIso = -10 ; //<--
        Float_t maxIso =  55 ; //<--

        TH1F * histoIso_loose4[nPtBins][nEtaBins][convorunconv][dataormc] ;
        TH1F * histoIso_tight[nPtBins][nEtaBins][convorunconv][dataormc] ;
	TH1F * histoIso_loose4norm[nPtBins][nEtaBins][convorunconv][dataormc] ;
	TH1F * histoIso_backgroundsub[nPtBins][nEtaBins][convorunconv][dataormc] ;	
	TH1F * histoIso_sherpanorm[nPtBins][nEtaBins][convorunconv][dataormc] ;
	TH1F * histoIso_pythianorm[nPtBins][nEtaBins][convorunconv][dataormc] ;
	TH1F * histoIso_XiSH[nPtBins][nEtaBins][convorunconv] ;
	TH1F * histoIso_XiPy[nPtBins][nEtaBins][convorunconv] ;
	TH1F * histoIso_SHmenosPy[nPtBins][nEtaBins][convorunconv] ;
	TH1F * histoalphaminSH[nEtaBins][convorunconv] ;
	TH1F * histoalphaminpy[nEtaBins][convorunconv] ;
	TH1F * _etaHisto[convorunconv][dataormc] ;
	TH1F * _ptHisto[dataormc] ;
	TH1F * histotightSH[nEtaBins][convorunconv] ;
	TH1F * histotightpy[nEtaBins][convorunconv] ;
	
	Int_t nBins = 0 ;	
	
	for ( Int_t iConv = 0 ; iConv < 3 ; iConv ++ ){
		TString _whichConv = Conv[iConv] ;
	for ( Int_t iType = 0 ; iType < 3 ; iType ++ ){
		TString _whichType = Type[iType] ;
        for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
                if ( iEta != 2 ) {
	for ( Int_t iPt = 0 ; iPt < nPtBins - 1 ; iPt ++ ){
		if ( ( iPt == 14 ) && (iEta == 3 || iEta == 4 ) ){
                nBins = 90 ;
                }
                if ( ( iPt == 14 ) && (iEta == 0 || iEta == 1 ) ){
                nBins = 90 ;
                }
                if ( ( iPt == 12 || iPt == 13 ) && (iEta == 3 || iEta == 4 ) && iConv == 0 ){
                nBins = 200 ;
                }
                if ( ( iPt == 12 || iPt == 13 ) && (iEta == 3 || iEta == 4 ) && iConv == 1 ){
                nBins = 180 ;
                }
                if ( ( iPt == 12 || iPt == 13 ) && (iEta == 0 || iEta == 1 ) ){
                nBins = 160 ;
                }
                if ( iPt == 8 || iPt == 9 || iPt == 10 || iPt == 11 ){
                nBins = 1500 ;
                }
                if ( iPt == 3 || iPt == 4 || iPt == 5 || iPt == 6 || iPt == 7 ){
                nBins = 1200 ;
                }
                if ( iPt == 0 || iPt == 1 || iPt == 2 ){
                nBins = 700 ;
                }
	
			// Create loose histo
			TString histoName = "histoIso_loose4_" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			TString histoTitle = " Loose'4 histogram for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons in " ;
			histoTitle += Type[iType] ;
                        histoIso_loose4[iPt][iEta][iConv][iType] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        TString histoXTitle = "(Loose'4 , pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_loose4[iPt][iEta][iConv][iType]->SetXTitle( histoXTitle ) ;
                        histoIso_loose4[iPt][iEta][iConv][iType]->SetYTitle( "Entries/0.5 GeV" ) ;
                       	// Create tight histos
			histoName = "histoIso_tight_" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
                        histoName += "_" ;
                        histoName += iType ;
			histoTitle = " Tight histogram (in red ) for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons in " ;
			histoTitle += Type[iType] ;
                        histoIso_tight[iPt][iEta][iConv][iType] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "(Tight , pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_tight[iPt][iEta][iConv][iType]->SetXTitle( histoXTitle ) ;
                        histoIso_tight[iPt][iEta][iConv][iType]->SetYTitle( "Entries/0.5 GeV" ) ;
			// Create loose norm histo
			histoName = "histoIso_loose4norm_" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " Loose'4 normalice histogram for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons in " ;
			histoTitle += iType ;
                        histoIso_loose4norm[iPt][iEta][iConv][iType] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "(Loose'4 norm, pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_loose4norm[iPt][iEta][iConv][iType]->SetXTitle( histoXTitle ) ;
                        histoIso_loose4norm[iPt][iEta][iConv][iType]->SetYTitle( "Entries/0.5 GeV" ) ;
			// Create background subtraction histo
			histoName = "histoIso_backgroundsub_" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " Background subtraction histogram (black) for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons in " ;
			histoTitle += " data. " ;
                        histoIso_backgroundsub[iPt][iEta][iConv][iType] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "(Background sub, pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_backgroundsub[iPt][iEta][iConv][iType]->SetXTitle( histoXTitle ) ;
                        histoIso_backgroundsub[iPt][iEta][iConv][iType]->SetYTitle( "Entries/0.5 GeV" ) ;
			// Create norm sherpa histo
			histoName = "histoIso_sherpanorm_" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " Normalize sherpa histogram for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons in " ;
			histoTitle += iType ;
                        histoIso_sherpanorm[iPt][iEta][iConv][iType] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "(Norm Sherpa histo, pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_sherpanorm[iPt][iEta][iConv][iType]->SetXTitle( histoXTitle ) ;
                        histoIso_sherpanorm[iPt][iEta][iConv][iType]->SetYTitle( "Entries/0.5 GeV" ) ;
			// Create norm pythia histo
			histoName = "histoIso_pythianorm" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " Normalize sherpa histogram for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons in " ;
			histoTitle += iType ;
                        histoIso_pythianorm[iPt][iEta][iConv][iType] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "(Pithya norm histo, pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_pythianorm[iPt][iEta][iConv][iType]->SetXTitle( histoXTitle ) ;
                        histoIso_pythianorm[iPt][iEta][iConv][iType]->SetYTitle( "Entries/0.5 GeV" ) ;	
			//HistotightSHERPA
			histoName = "histoIso_tightSHERPA" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " IsoTight fraction_Sherpa vs pT for iEta: [" ;
			histoTitle += etaBins[iEta] ;
			histoTitle += " : " ;
			histoTitle += etaBins[iEta+1] ;
			histoTitle += " ] and " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons. " ;
                        histotightSH[iEta][iConv] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = " pT " ; 
                        histotightSH[iEta][iConv]->SetXTitle( histoXTitle ) ;
                        histotightSH[iEta][iConv]->SetYTitle( "Fraction of Events < 5GeV for (tight) photons" ) ;
			//Histotightpythia
			histoName = "histoIso_tightpythia" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " IsoTight fraction_pythia vs pT for iEta: [" ;
			histoTitle += etaBins[iEta] ;
			histoTitle += " : " ;
			histoTitle += etaBins[iEta+1] ;
			histoTitle += " ] and " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons. " ;
                        histotightpy[iEta][iConv] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = " pT " ; 
                        histotightpy[iEta][iConv]->SetXTitle( histoXTitle ) ;
                        histotightpy[iEta][iConv]->SetYTitle( "Fraction of Events < 5GeV for (tight) photons" ) ;
			//Eta histogram for the photon
			histoName = "histoeta_" ;	
			histoName += iPt ;
			histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType ;
			histoTitle = " Eta histogram for " ;
			histoTitle += Conv[iConv] ;
			histoTitle += " photons in "
			histoTitle += _whichType ;	
			_etaHisto[iConv][iType] = new TH1F( histoName , histoTitle , 500 , -2.5 , 2.5 ) ;
			histoxTitle = " eta " ;
			_etaHisto[iConv][iType]->SetXTitle( histoXTitle ) ;	
			//pt histogram for the photon
			histoName = "histoIso_tightpT" ;	
			histoName += iPt ;
			histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType ;
			histoTitle = " (tight) pT histogram in " ;
			histoTitle += _whichType ;
			_ptHisto[iType] = new TH1F( histoName , histoTitle , 500 , 0 , 600 ) ;
			histoxTitle = " pT " ;
		}
        }
        }
	}
	}

	cout << "Will create some first histograms" << endl ;

	for ( Int_t iType = 0 ; iType < 3 ; iType ++ ){
		TString _whichType = Type[iType] ;
		cout << " In " << _whichType << endl ;
		TString fileName = "NONE" ;
        	if ( _whichType == Type[0] ) fileName = "data_AllYear_v03" ;
        	if ( _whichType == Type[1] ) fileName = "ShPt15_inf_mc15b_v03" ;
		if ( _whichType == Type[2] ) fileName = "PyPt17_inf_mc15b_v03" ;
     		if ( fileName =="NONE" ) {
                	cout << "ERROR : file name does not match \"data\" , \"sherpa\" or \"pythia\" ! " << endl ;
                	return ;
        	}
		
		char filePath[200] ;
        	sprintf( filePath , "SinglePhoton_v3/%s.root" , fileName.Data() ) ;
        	cout << "Will Open a rootfile in this path : " << filePath << endl ;
	
        	TFile * openFile = TFile::Open( filePath , "READ" ) ;
        	openFile->ls() ; //Open

        	TTree * myTree = ( TTree * ) openFile->Get( "SinglePhoton" ) ;

        	cout << "Opened a Tree that contains " << myTree->GetEntries() << " entries" << endl ;

       		myTree->SetBranchAddress( "y_pt" , & y_pt ) ;
        	myTree->SetBranchAddress( "y_eta" , & y_eta ) ;
        	myTree->SetBranchAddress( "y_phi" , & y_phi ) ;
        	myTree->SetBranchAddress( "y_convType" , & y_convType ) ;
        	myTree->SetBranchAddress( "y_topoetcone40" , & y_topoetcone40 ) ;
        	myTree->SetBranchAddress( "y_IsTight" , & y_IsTight ) ;
        	myTree->SetBranchAddress( "y_IsLoosePrime4" , & y_IsLoosePrime4 ) ;	

		TCanvas * etaconvdata = new TCanvas ( "etaconvdata" , "etaconvdata" ,  800 , 800 ) ;
	
		TString histoNamedataetaconv = "histoeta_" ;
                histoNamedataetaconv += iType ;
                histoNamedataetaconv += "_for_convphotons" ;
		TString histoTitledataetaconv = " Eta histogram for " ;
		histoTitledataetaconv += Type[0] ;
		histoTitledataetaconv += " with convphotons." ;	
		TH1F * _etaHistoconv = new TH1F ( histoNamedataetaconv , histoTitledataetaconv , 500 , -2.5 , 2.5 ) ;			
	
		myTree->Draw("y_eta>>_etaHistoconv","y_convType==0") ;
		etaconvdata->Update() ;

		TString etaconvdpdf = "etaconv_" ;
		etaconvdpdf = iType ;
		etaconvdpdf = ".pdf" ;	
		etaconvdata->Print( etaconvdpdf ) ;
	
		TCanvas * etaunconvdata = new TCanvas ( "etaunconvdata" , "etaunconvdata" ,  800 , 800 ) ;
	
		TString histoNamedataetaunconv = "histoeta_" ;
                histoNamedataetaunconv += iType ;
                histoNamedataetaunconv += "_for_unconvphotons" ;
		TString histoTitledataetaunconv = " Eta histogram for " ;
		histoTitledataetaunconv += Type[0] ;
		histoTitledataetaunconv += " with unconvphotons." ;	
		TH1F * _etaHistounconv = new TH1F ( histoNamedataetaunconv , histoTitledataetaunconv , 500 , -2.5 , 2.5 ) ;

		myTree->Draw("y_eta>>_etaHistounconv","y_convType==1") ;
		etaunconvdata->Update() ;

		TString etaunconvdpdf = "etaunconv_" ;
		etaunconvdpdf = iType ;
		etaunconvdpdf = ".pdf" ;	
		etaunconvdata->Print( etaunconvdpdf ) ;
	
		TCanvas * ptdata = new TCanvas ( "ptdata" , "ptdata" ,  800 , 800 ) ;
	
		TString histoNamedatapt = "histopt_" ;
                histoNamedatapt += iType ; 
		TString histoTitledatapt = " pt histogram for " ;
		histoTitledatapt += Type[0] ;	
		TH1F * _ptHisto = new TH1F ( histoNamedatapt , histoTitledatapt , 500 , 0 , 600 ) ;	
		
		myTree->Draw("y_pt>>_ptHisto") ;	
		ptdata->Update() ;
		
		TString ptdpdf = "pt_" ;
		ptdpdf = iType ;
		ptdpdf = ".pdf" ;	
		ptdata->Print( ptdpdf ) ;

		TCanvas * etaconvSH = new TCanvas ( "etaconvSH" , "etaconvSH" ,  800 , 800 ) ;
	
		TString histoNameSherpaetaconv = "histoeta_" ;
                histoNameSherpaetaconv += iType ;
                histoNameSherpaetaconv += "_for_convphotons" ;
		TString histoTitleSherpaetaconv = " Eta histogram for " ;
		histoTitleSherpaetaconv += Type[1] ;
		histoTitleSherpaetaconv += " with conv photons." ;
		TH1F * _etaHistoconvSH = new TH1F ( histoNameSherpaetaconv , histoTitleSherpaetaconv , 500 , -2.5 , 2.5 ) ;
		
		myTree->Draw("y_eta>>_etaHistoconvSH","y_convType==0") ;
		etaconvSH->Update() ;

		TString etaconvdpdf = "etaconv_" ;
		etaconvdpdf = iType ;
		etaconvdpdf = ".pdf" ;	
		etaconvSH->Print( etaconvdpdf ) ;

		TCanvas * etaunconvSH = new TCanvas ( "etaunconvSH" , "etaunconvSH" ,  800 , 800 ) ;
	
		TString histoNameSherpaetaunconv = "histoeta_" ;
                histoNameSherpaetaunconv += iType ;
                histoNameSherpaetaunconv += "_for_unconvphotons" ;
		TString histoTitleSherpaetaunconv = " Eta histogram for " ;
		histoTitleSherpaetaunconv += Type[1] ;
		histoTitleSherpaetaunconv += " with unconv photons." ;
		TH1F * _etaHistounconvSH = new TH1F ( histoNameSherpaetaconv , histoTitleSherpaetaconv , 500 , -2.5 , 2.5 ) ;
		
		myTree->Draw("y_eta>>_etaHistounconvSH","y_convType==1") ;			
		etaunconvSH->Update() ;

		TString etaunconvdpdf = "etaunconv_" ;
		etaunconvdpdf = iType ;
		etaunconvdpdf = ".pdf" ;	
		etaunconvSH->Print( etaunconvdpdf ) ;

		TCanvas * ptSH = new TCanvas ( "ptSH" , "ptSH" ,  800 , 800 ) ;

		TString histoNameSherpapt = "histopt_" ;
                histoNameSherpapt += iType ; 
		TString histoTitleSherpapt = " pt histogram for " ;
		histoTitleSherpapt += Type[1] ;	
		TH1F * _ptHistoSH = new TH1F ( histoNameSherpapt , histoTitleSherpapt , 500 , 0 , 600 ) ;
		
		myTree->Draw("y_pt>>_ptHistoSH") ;	
		ptSH->Update() ;
		
		TString ptdpdf = "pt_" ;
		ptdpdf = iType ;
		ptdpdf = ".pdf" ;	
		ptSH->Print( ptdpdf ) ;

		TCanvas * etaconvpy = new TCanvas ( "etaconvpy" , "etaconvpy" ,  800 , 800 ) ;

		TString histoNamepythiaetaconv = "histoeta_" ;
                histoNamepythiaetaconv += iType ;
                histoNamepythiaetaconv += "_for_convphotons" ;
		TString histoTitlepythiaetaconv = " Eta histogram for " ;
		histoTitlepythiaetaconv += Type[2] ;
		histoTitlepythiaetaconv += " with convphotons." ;	
		TH1F * _etaHistoconvpy = new TH1F ( histoNamepythiaetaconv , histoTitlepythiaetaconv , 500 , -2.5 , 2.5 ) ;	
		
		myTree->Draw("y_eta>>_etaHistoconvpy","y_convType==0") ;
		etaconvpy->Update() ;

		TString etaconvdpdf = "etaconv_" ;
		etaconvdpdf = iType ;
		etaconvdpdf = ".pdf" ;	
		etaconvpy->Print( etaconvdpdf ) ;

		TCanvas * etaunconvpy = new TCanvas ( "etaunconvpy" , "etaunconvpy" ,  800 , 800 ) ;

		TString histoNamepythiaetaunconv = "histoeta_" ;
                histoNamepythiaetaunconv += iType ;
                histoNamepythiaetaunconv += "_for_unconvphotons" ;
		TString histoTitlepythiaetaunconv = " Eta histogram for " ;
		histoTitlepythiaetaunconv += Type[2] ;
		histoTitlepythiaetaunconv += " with unconvphotons." ;	
		TH1F * _etaHistounconvpy = new TH1F ( histoNamepythiaetaunconv , histoTitlepythiaetaunconv , 500 , -2.5 , 2.5 ) ;	
		
		myTree->Draw("y_eta>>_etaHistounconvpy","y_convType==1") ;
		etaunconvpy->Update() ;

		TString etaunconvdpdf = "etaunconv_" ;
		etaunconvdpdf = iType ;
		etaunconvdpdf = ".pdf" ;	
		etaunconvpy->Print( etaunconvdpdf ) ;

		TCanvas * ptpy = new TCanvas ( "ptpy" , "ptpy" ,  800 , 800 ) ;

		TString histoNamepythiapt = "histopt_" ;
                histoNamepythiapt += iType ;
		TString histoTitlepythiapt = " pt histogram for " ;
		histoTitlepythiapt += Type[2] ;
		TH1F * _ptHistopy = new TH1F ( histoNamepythiapt , histoTitlepythiapt , 500 , 0 , 600 ) ;
		
		myTree->Draw("y_pt>>_ptHistopy") ;
		ptpy->Update() ;
		
		TString ptdpdf = "pt_" ;
		ptdpdf = iType ;
		ptdpdf = ".pdf" ;	
		ptpy->Print( ptdpdf ) ;

		cout << "Will loop over events in the tree" << endl ;	
	
		for ( Int_t iEntry = 0 ; iEntry < _whichFraction * myTree->GetEntries() ; iEntry ++ ) {
	               	myTree->GetEntry ( iEntry ) ;	
			
			// Print out just a few events	

			Bool_t verbose = kFALSE ;
       			if ( iEntry < 10 ) verbose = kTRUE ;
            		if ( iEntry%2000000 == 0 ) verbose = kTRUE ;
               		if ( gRandom->Uniform(1000000) < 2 ) verbose = kTRUE ;
               									
			if ( verbose ) {
              			cout << "\t Read event " << iEntry << endl ;
                    		cout << "\t\t y_pt = " << y_pt << " ; y_eta = " << y_eta << " ; y_phi = " << y_phi << " ; y_convType = " << y_convType << endl ;
                		}
		
			for ( Int_t iConv = 0 ; iConv < 3 ; iConv ++ ){
               			Bool_t selectEvent = kTRUE ;
		
				Int_t iPt = - 999 ;
               			for ( Int_t jPt = 0 ; jPt < nPtBins - 1 ; jPt ++ ){
                     			if ( pTbins[jPt] < y_pt && pTbins[jPt+1] > y_pt ) iPt = jPt ;
               			}
              			if ( iPt == - 999 ) selectEvent = kFALSE ;

               			Int_t iEta = - 999 ;
              			Float_t absEta = TMath::Abs( y_eta ) ;
               			for ( Int_t jEta = 0 ; jEta < nEtaBins - 1 ; jEta ++ ){
                       			if ( etaBins[jEta] < absEta && etaBins[jEta+1] > absEta ) iEta = jEta ;
               			}	
               			if ( iEta == 2 || iEta == - 999 ) selectEvent = kFALSE ;
			
				TString _whichConv = Conv[iConv] ;		
		
				if ( _whichConv == Conv[0] && y_convType == 0 ) selectEvent = kFALSE ; 	
				if ( _whichConv == Conv[1] && y_convType != 0 ) selectEvent = kFALSE ; 	
	
				if ( verbose ) cout << "\t\t selectEvent = " << selectEvent << endl;
				if ( selectEvent ) {
                       			if ( verbose ) cout << "\t\t\t pT,|eta| = " << y_pt << " , " << absEta << " corresponds to bins = " << iPt << " , " << iEta << " for " << Conv[iConv] << " photons in " << Type[iType] << " . " << endl ;
                     			if ( y_IsTight ) histoIso_tight[iPt][iEta][iConv][iType]->Fill( y_topoetcone40 ) ;
                    			if ( ! y_IsTight && y_IsLoosePrime4 ) histoIso_loose4[iPt][iEta][iConv][iType]->Fill( y_topoetcone40 ) ;
               			}
       			} //For the iConv

		} //For the iEntry
		
		cout << "Done with the loop !" << endl ;
	
		openFile->Close() ;

		cout << " Close file for " << _whichType << endl ;

	} //For the iType	

	Int_t Internalsuffering = 0 ;

	//For tight MC
	Float_t vector_tightSH[nPtBins][nEtaBins][convorunconv - 1] ; // For the TGraph MC tight vs pt
	Float_t vector_tightpy[nPtBins][nEtaBins][convorunconv - 1] ; //this is for later
	Float_t vector_tightSHerror[nPtBins][nEtaBins][convorunconv - 1] ;
	Float_t vector_tightpyerror[nPtBins][nEtaBins][convorunconv - 1] ;

        for ( Int_t iConv = 0 ; iConv < 2 ; iConv ++ ){
	for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
                if ( iEta != 2 ) {	
	for ( Int_t iPt = 0 ; iPt < nPtBins - 1 ; iPt ++ ){		
		if ( ( iPt == 14 ) && (iEta == 3 || iEta == 4 ) ){
                nBins = 90 ;
                }
                if ( ( iPt == 14 ) && (iEta == 0 || iEta == 1 ) ){
                nBins = 90 ;
                }
                if ( ( iPt == 12 || iPt == 13 ) && (iEta == 3 || iEta == 4 ) && iConv == 0 ){
                nBins = 200 ;
                }
                if ( ( iPt == 12 || iPt == 13 ) && (iEta == 3 || iEta == 4 ) && iConv == 1 ){
                nBins = 180 ;
                }
                if ( ( iPt == 12 || iPt == 13 ) && (iEta == 0 || iEta == 1 ) ){
                nBins = 160 ;
                }
                if ( iPt == 8 || iPt == 9 || iPt == 10 || iPt == 11 ){
                nBins = 1500 ;
                }
                if ( iPt == 3 || iPt == 4 || iPt == 5 || iPt == 6 || iPt == 7 ){
                nBins = 1200 ;
                }
                if ( iPt == 0 || iPt == 1 || iPt == 2 ){
                nBins = 700 ;
                }
	cout << " Now, as an example for iPt, iEta : " << iPt << " , " << iEta << " . It will show some of the histograms. " << endl ;
	
		Float_t hightisolationregion = 20 ;// for the region above 20 GeV
		Float_t fondo = 5 ;
	
		// For normalice looseprime4           
		if ( histoIso_loose4[iPt][iEta][iConv][0]->Integral() > 0 && histoIso_tight[iPt][iEta][iConv][0]->Integral() > 0) {
               		Float_t tightIntegral = 0 ;
       			Float_t loose4Integral = 0 ;
			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
                       		Float_t xBin = histoIso_tight[iPt][iEta][iConv][0]->GetBinCenter( iBin ) ; //<---
                       		if ( xBin > hightisolationregion ) {
					tightIntegral += histoIso_tight[iPt][iEta][iConv][0]->GetBinContent( iBin ) ; //<---		
					loose4Integral += histoIso_loose4[iPt][iEta][iConv][0]->GetBinContent( iBin ) ; //<--	
                       		}
               		}
			
			cout << " For " << Conv[iConv] << " photons with iPt, iEta " << iPt << " , " << iEta <<endl ;
			cout << "\t *loose'4 histo in " << Type[0] << " : the events above " << hightisolationregion << " GeV are " << loose4Integral << endl ;
			cout << "\t *tight histo in " << Type[0] << " : the events above " << hightisolationregion << " Gev are " << tightIntegral << endl ;
	
			Float_t frac = tightIntegral/loose4Integral ; //<-- OJO
			
			// Let's fill loose'4 norm
			
			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
                       		
				Float_t yBinloose4 = histoIso_loose4[iPt][iEta][iConv][0]->GetBinContent( iBin ) ;	
				histoIso_loose4norm[iPt][iEta][iConv][0]->SetBinContent( iBin , frac*yBinloose4 ) ;
	
				Float_t sBinloose4 = histoIso_loose4[iPt][iEta][iConv][0]->GetBinError( iBin ) ;
				histoIso_loose4norm[iPt][iEta][iConv][0]->SetBinError( iBin , frac*sBinloose4 ) ;
				
			}
	
      			//histoIso_loose4norm[iPt][iEta][iConv][0]->Draw("SAME") ;	
			
			Float_t loose4normover20Integral = 0 ;
			Float_t backgroundsubover20Integral = 0 ;	
			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
                       		Float_t xBin = histoIso_tight[iPt][iEta][iConv][0]->GetBinCenter( iBin ) ;
                       		if ( xBin > hightisolationregion ) { 
					loose4normover20Integral += histoIso_loose4norm[iPt][iEta][iConv][0]->GetBinContent( iBin ) ;	
                       		}
				Float_t yBintight = histoIso_tight[iPt][iEta][iConv][0]->GetBinContent( iBin ) ;
				Float_t yBinloosenorm = histoIso_loose4norm[iPt][iEta][iConv][0]->GetBinContent( iBin ) ;
				histoIso_backgroundsub[iPt][iEta][iConv][0]->SetBinContent( iBin , yBintight - yBinloosenorm ) ;
				
				Float_t sBintight = histoIso_tight[iPt][iEta][iConv][0]->GetBinError( iBin );
				Float_t sBinloosenorm = histoIso_loose4norm[iPt][iEta][iConv][0]->GetBinError ( iBin ) ;
				histoIso_backgroundsub[iPt][iEta][iConv][0]->SetBinError( iBin , TMath::Sqrt((sBintight*sBintight) + (sBinloosenorm*sBinloosenorm))) ;
				if ( xBin > hightisolationregion ){
					backgroundsubover20Integral += histoIso_backgroundsub[iPt][iEta][iConv][0]->GetBinContent( iBin );
				}
			}
			
			// Draw the normalize

			TCanvas * normalize1 = new TCanvas( "normalize1" , "normalize1" , 1200 , 1200 ) ;

			Internalsuffering ++ ;
			normalize1->cd( Internalsuffering ) ;
			histoIso_backgroundsub[iPt][iEta][iConv][0]->SetLineColor(kBlack) ;
			histoIso_backgroundsub[iPt][iEta][iConv][0]->Draw() ;
			TLine * line = new TLine(-10 , 0 , 45 , 0 ) ;
  			line->SetLineColor(kOrange);
  			line->Draw("SAME");
			normalize1->Update() ;
	
			cout << "\t *Area for normalize loose'4 over 20 GeV: " << loose4normover20Integral <<  endl ;

			cout << "\t *Area for data tight background subtracted over 20 Gev: " << backgroundsubover20Integral << endl ;

			cout << "\t\t **The fraction of events above 20 GeV for the tight histogram and the loose'4 is " << frac << endl; 
			
			//For normalize the MC

			Float_t bgsubIntegral = 0 ;
       			Float_t tightSHIntegral = 0 ;
			Float_t tightPyIntegral = 0 ;

			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
				tightSHIntegral += histoIso_tight[iPt][iEta][iConv][1]->GetBinContent( iBin ) ;	
       					
				tightPyIntegral += histoIso_tight[iPt][iEta][iConv][2]->GetBinContent( iBin ) ;
				
				bgsubIntegral += histoIso_backgroundsub[iPt][iEta][iConv][0]->GetBinContent( iBin ) ;
			}
	
			Float_t tight5SHIntegral = 0 ;
			Float_t tight5pyIntegral = 0 ;	
			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
                       		Float_t xSHBin = histoIso_tight[iPt][iEta][iConv][1]->GetBinCenter( iBin ) ;
                       		Float_t xpyBin = histoIso_tight[iPt][iEta][iConv][2]->GetBinCenter( iBin ) ;
				if ( xSHBin < fondo ) { 
					tight5SHIntegral += histoIso_tight[iPt][iEta][iConv][1]->GetBinContent( iBin ) ;	
                       		}
				if ( xpyBin < fondo ) { 
					tight5pyIntegral += histoIso_tight[iPt][iEta][iConv][2]->GetBinContent( iBin ) ;	
                       		}	
			}

			cout << "\t *The background substracted in " << Type[0] << " have " << bgsubIntegral << " events. <-- " << endl ;
			
			cout << "\t *Tight Monte Carlo in " << Type[1] << " have " << tightSHIntegral << endl ;

			Float_t FractionSH = tight5SHIntegral / histoIso_tight[iPt][iEta][iConv][1]->Integral() ; //Filling the arrays
			vector_tightSH[iPt][iEta][iConv] = FractionSH ;
			Float_t ErrorSH = TMath::Sqrt( FractionSH * ( 1.0 - FractionSH ) / histoIso_tight[iPt][iEta][iConv][1]->Integral() ) ;
			vector_tightSHerror[iPt][iEta][iConv] = ErrorSH ;

			cout << "\t *Tight Monte Carlo in " << Type[2] << " have " << tightPyIntegral << endl ;
		
			Float_t Fractionpy = tight5pyIntegral / histoIso_tight[iPt][iEta][iConv][1]->Integral() ;
			vector_tightpy[iPt][iEta][iConv] = Fractionpy ;
			Float_t Errorpy = TMath::Sqrt( Fractionpy * ( 1.0 - Fractionpy ) / histoIso_tight[iPt][iEta][iConv][2]->Integral() ) ;
			vector_tightpyerror[iPt][iEta][iConv] = Errorpy ;

			histotightSH[iEta][iConv]->SetBinContent( iPt + 1 , vector_tightSH[iPt][iEta][iConv] ) ;
			histotightSH[iEta][iConv]->SetBinError( iPt +1 , vector_tightSHerror[iPt][iEta][iConv] ) ;
			histotightSH[iEta][iConv]->SetXTitle( "pT" ) ;
			histotightSH[iEta][iConv]->SetYTitle( " Fraction of Events < 5GeV for (tight) photons[GeV]" ) ;

			histotightpy[iEta][iConv]->SetBinContent( iPt + 1 , vector_tightpy[iPt][iEta][iConv] ) ;
			histotightpy[iEta][iConv]->SetBinError( iPt +1 , vector_tightpyerror[iPt][iEta][iConv] ) ;
			histotightpy[iEta][iConv]->SetXTitle( "pT" ) ;
			histotightpy[iEta][iConv]->SetYTitle( " Fraction of Events < 5GeV for (tight) photons[GeV]" ) ;
		
			//Continues				
			Float_t fracSH = bgsubIntegral/tightSHIntegral ; //<-- OJO
			Float_t fracPy = bgsubIntegral/tightPyIntegral ; //<-- OJO
				
			cout << "\t\t **The fraction of events for the background substracted histogram and tight " << Type[1] << " is " << fracSH << endl; 
			cout << "\t\t **The fraction of events for the background substracted histogram and tight " << Type[2] << " is " << fracPy << endl;
			
			// Let's fill sherpanorm
				
			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
                       		
				Float_t yBinsherpa = histoIso_tight[iPt][iEta][iConv][1]->GetBinContent( iBin ) ;	
				histoIso_sherpanorm[iPt][iEta][iConv][1]->SetBinContent( iBin , fracSH*yBinsherpa ) ;	

				Float_t sBinsherpa = histoIso_tight[iPt][iEta][iConv][1]->GetBinError( iBin ) ;
				histoIso_sherpanorm[iPt][iEta][iConv][1]->SetBinError( iBin , fracSH*sBinsherpa ) ;			

			}
			
			// Let's fill pythianorm
			
			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
                       		
				Float_t yBinpythia = histoIso_tight[iPt][iEta][iConv][2]->GetBinContent( iBin ) ;	
				histoIso_pythianorm[iPt][iEta][iConv][2]->SetBinContent( iBin , fracPy*yBinpythia ) ;	

				Float_t sBinpythia = histoIso_tight[iPt][iEta][iConv][2]->GetBinError( iBin ) ;
				histoIso_pythianorm[iPt][iEta][iConv][2]->SetBinError( iBin , fracPy*sBinpythia ) ;
				
			}

			//Check normalization	

			Float_t sherpanormIntegral = 0 ;	
			Float_t pythianormIntegral = 0 ;
			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
				sherpanormIntegral += histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinContent( iBin ) ;
				pythianormIntegral += histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinContent( iBin ) ;
			}
		
			cout << " \t *Tight normalize, Monte Carlo in " << Type[1] << " have " << sherpanormIntegral << " events. <-- " << endl ;
			cout << " \t *Tight normalize, Monte Carlo in " << Type[2] << " have " << pythianormIntegral << " events. <-- " << endl ;

			//Some Values

			Float_t meanSHnorm = histoIso_sherpanorm[iPt][iEta][iConv][1]->GetMean() ;
			Float_t meanPynorm = histoIso_pythianorm[iPt][iEta][iConv][2]->GetMean() ;
			Float_t rmsSHnorm = histoIso_sherpanorm[iPt][iEta][iConv][1]->GetRMS() ;
			Float_t rmsPynorm =  histoIso_pythianorm[iPt][iEta][iConv][2]->GetRMS() ;

			cout << "\t * Some values: " << endl ;

			cout << "\t\t **Tight normalize, Monte Carlo in " << Type[1] << " have a mean value of " << meanSHnorm << " and the RMS is " << rmsSHnorm << endl ;
			cout << "\t\t **Tight normalize, Monte Carlo in " << Type[2] << " have a mean value of " << meanPynorm << " and the RMS is " << rmsPynorm << endl ;

			cout << "\t\t\t ** BIN WIDTH in SH " << histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinWidth(1) << " , in py " << histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinWidth(1) << endl ;

			histoIso_sherpanorm[iPt][iEta][iConv][1]->SetLineColor(kRed) ;
			histoIso_sherpanorm[iPt][iEta][iConv][1]->Draw("SAME") ;	
			normalize1->Update() ;				
			histoIso_pythianorm[iPt][iEta][iConv][2]->SetLineColor(kBlue) ;
			histoIso_pythianorm[iPt][iEta][iConv][2]->Draw("SAME") ;	
   			normalize1->Update() ;
		
			TString normalizationpdf = "nBins_" ;
			normalizationpdf += nBins ;
			normalizationpdf += "_normalization_data_MC_in_iPt,iEta:";
			normalizationpdf += iPt ;
			normalizationpdf += "," ;
			normalizationpdf += iEta ;
			normalizationpdf += "For_" ;
			normalizationpdf += Conv[iConv] ;
			normalizationpdf += "photons" ;
			normalizationpdf += "_binwidth_" ;
			normalizationpdf += histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinWidth(1) ;
			normalizationpdf += ".pdf" ;
			normalize1->Print( normalizationpdf ) ;

		}//thisbigif

	}
		}
	}
	}

  /*      for ( Int_t iConv = 0 ; iConv < 2 ; iConv ++ ){
                for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
                        if ( iEta != 2 ) {
 
                                TCanvas * tightSHCanvas = new TCanvas( "tightSHCanvas" , "tightSHCanvas" , 1700 , 1700 ) ;                                   
                                //->SetMaximum(1.6) ;
                                //->SetMinimum(-1) ;
                                histotightSH[iEta][iConv]->Draw() ;
                                tightSHCanvas->Update() ;

                                TString tightSHprint = "TightSHfraction_" ;
                                tightSHprint += Conv[iConv] ;
                                tightSHprint += "photons_Etabin" ;
                                tightSHprint += iEta ; 
                                tightSHprint += ".pdf" ;
                                tightSHCanvas->Print(tightSHprint) ; //<--

                                TCanvas * tightpyCanvas = new TCanvas( "tightpyCanvas" , "tightpyCanvas" , 1700 , 1700 ) ;

                                //->SetMaximum(1.6) ;
                                //->SetMinimum(-1) ;
                                histotightpy[iEta][iConv]->Draw() ;
                                tightpyCanvas->Update() ;

                                TString tightpyprint = "Tightpyfraction_" ;
                                tightpyprint += Conv[iConv] ;
                                tightpyprint += "photons_Etabin" ;
                                tightpyprint += iEta ; 
                                tightpyprint += ".pdf" ;
                                tightpyCanvas->Print(tightpyprint) ; //<--
                        }
                }
        }
*/	 
};

	
