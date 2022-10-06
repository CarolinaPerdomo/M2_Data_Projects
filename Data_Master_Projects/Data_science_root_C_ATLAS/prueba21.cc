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

using namespace std;

Float_t computeChi2( TH1F* histo0 , TH1F * histo1 , Float_t xMin = -99999.0f , Float_t xMax = -99999.0f ) {

        /* Input arguments : 
		* 2 histograms histo0 as reference and histo1 as test.
		* xMin and xMax are the intervals to be used to compute Chi2. 
        // Return argument : Chi2 value. */

	// Check that the two histos are compatible (same binning...)
	
	Int_t nBins0 = histo0->GetNbinsX() ;
	Float_t min0 = histo0->GetBinLowEdge( 1 ) ;
	Float_t max0 = histo0->GetBinLowEdge( nBins0 ) + histo0->GetBinWidth( nBins0 ) ;
	Int_t nBins1 = histo1->GetNbinsX() ;
	Float_t min1 = histo1->GetBinLowEdge( 1 ) ;
	Float_t max1 = histo1->GetBinLowEdge( nBins1 ) + histo1->GetBinWidth( nBins1 ) ;
	Float_t integral0 = histo0->Integral() ;
	Float_t integral1 = histo1->Integral() ;

	Bool_t canComputeChi2 = kTRUE ;
	if ( nBins0 != nBins1 || min0 != min1 || max0 != max1 ) canComputeChi2 = kFALSE ;
	if ( ! canComputeChi2 ) {
		cout << "ERROR : histograms have different binnings, cannot computeChi2 chi2 !!! " << endl ;
	return -99999 ;
        }

        // check the histos are not empty

	if ( integral0 <= 0 || integral1 <= 0 ) canComputeChi2 = kFALSE ;
	if ( ! canComputeChi2 ) {
		cout << "ERROR : one of the histograms is empty ! " << endl ;
		return -99999 ;
	}

        // Check whether the complete histogram range will be used
        
	if ( xMin == -99999.0f && xMax == -99999.0f ) {
		xMin = min0 ;
		xMax = max0 ;
        }

	// Normalize histo1 to histo0
	Float_t normalizationFactor = integral0 / integral1 ;

	Float_t Chi2 = 0.0f ;
	for ( Int_t iBin = 1 ; iBin <= nBins1 ; iBin ++ ) {
		Float_t xBin = histo0->GetBinLowEdge( iBin ) ;
                Float_t dBin = histo0->GetBinWidth( iBin ) ;
		if ( xBin > xMin && xBin + dBin < xMax ) {// Only compute Chi2 in the xMin:xMax interval
			Float_t yBin0 = histo0->GetBinContent( iBin ) ;
			Float_t sBin0 = histo0->GetBinError( iBin ) ;
                        Float_t yBin1 = histo1->GetBinContent( iBin ) ;
                        Float_t sBin1 = histo1->GetBinError( iBin ) ;
                        Float_t Delta2 = TMath::Power( yBin0 - normalizationFactor * yBin1 , 2 ) ;
                        Float_t Error2 =  TMath::Power( sBin0 , 2 ) + TMath::Power(  normalizationFactor * sBin1 , 2 ) ;
                        if ( Error2 > 0.0f ) Chi2 += Delta2 / Error2 ;
                }
        }

        return Chi2 ;
}

TH1F * shiftHistogram( TH1F * histoInput , Int_t nBinsShift = 0 ) {
    
	/* Input arguments :
		*the histogram to be shifted
		*the shift to be applied (in number of bins)
      
	Return argument
		*a shifted histogram */

    //Check the binning of input histogram
    Int_t nBins = histoInput->GetNbinsX() ;
    Float_t xMin = histoInput->GetBinCenter( 1 )  - 0.5f * histoInput->GetBinWidth( 1 ) ;
    Float_t xMax = histoInput->GetBinCenter( nBins ) + 0.5f * histoInput->GetBinWidth( nBins ) ;

    TString histoName = histoInput->GetName() ;
    TString histoTitle = histoInput->GetTitle() ;

    //Define new names for the new histogram to create, add "_Shift_" and nBinsShift
    histoName += "_shift_" ;
    if ( nBinsShift >= 0 ) {
        histoName += "p" ;
        histoName += nBinsShift ;
    } else {
        histoName += "m" ;
        histoName += -nBinsShift ;

    }
    histoTitle += "_shift_" ;
    if ( nBinsShift >= 0 ) {
        histoTitle += "p" ;
        histoTitle += nBinsShift ;
    } else {
        histoTitle += "m" ;
        histoTitle += -nBinsShift ;

    }


    // Create a new histo with the same binning as the input one
    TH1F * histoOutput = new TH1F( histoName , histoTitle , nBins , xMin , xMax ) ;

    // Loop over bins
    for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) { 
	Int_t jBin = iBin + nBinsShift ;
        if ( jBin >=1 && jBin <= nBins ) {
            Float_t yBin = histoInput->GetBinContent( iBin ) ;// Check bin content and error of input histo
            Float_t sBin = histoInput->GetBinError( iBin ) ;
            histoOutput->SetBinContent( jBin , yBin ) ;// Fill the new histo nBinsShift away of the input bin
            histoOutput->SetBinError( jBin , sBin ) ;
        }
    }

    return histoOutput ; 
}

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

        Int_t nBins = 680 ; //<--
        Float_t minIso = -10 ; //<--
        Float_t maxIso =  45 ; //<--

        TH1F * histoIso_loose4[nPtBins][nEtaBins][convorunconv][dataormc] ;
        TH1F * histoIso_tight[nPtBins][nEtaBins][convorunconv][dataormc] ;
	TH1F * histoIso_loose4norm[nPtBins][nEtaBins][convorunconv][dataormc] ;
	TH1F * histoIso_backgroundsub[nPtBins][nEtaBins][convorunconv][dataormc] ;	
	TH1F * histoIso_sherpanorm[nPtBins][nEtaBins][convorunconv][dataormc] ;
	TH1F * histoIso_pythianorm[nPtBins][nEtaBins][convorunconv][dataormc] ;
	TH1F * histoIso_XiSH[nPtBins][nEtaBins][convorunconv] ;
	TH1F * histoIso_XiPy[nPtBins][nEtaBins][convorunconv] ;
	TH1F * histoIso_SHmenosPy[nPtBins][nEtaBins][convorunconv] ;
	TH1F * histoIso_test[nPtBins][nEtaBins][convorunconv] ;

	for ( Int_t iConv = 0 ; iConv < 3 ; iConv ++ ){
		TString _whichConv = Conv[iConv] ;
	for ( Int_t iType = 0 ; iType < 3 ; iType ++ ){
		TString _whichType = Type[iType] ;
	for ( Int_t iPt = 0 ; iPt < nPtBins - 1 ; iPt ++ ){
        for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
                if ( iEta != 2 ) {
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
			// Create Sherpa Xi histo
			histoName = "histoIso_XiSH" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " Distribution for Data-Sherpa (red), Data-Pythia (blue) and Sherpa-Pythia (black ) histogram for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
                        histoIso_XiSH[iPt][iEta][iConv] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "Data-Sherpa histo, pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_XiSH[iPt][iEta][iConv]->SetXTitle( histoXTitle ) ;
                        histoIso_XiSH[iPt][iEta][iConv]->SetYTitle( "Entries subtraction" ) ;	
			// Create Pythia Xi histo
			histoName = "histoIso_XiPy" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " Distribution for Data-Sherpa (red), Data-Pythia (blue) and Sherpa-Pythia (black ) histogram for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
                        histoIso_XiPy[iPt][iEta][iConv] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "(Data-Pythia histo, pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_XiPy[iPt][iEta][iConv]->SetXTitle( histoXTitle ) ;
                        histoIso_XiPy[iPt][iEta][iConv]->SetYTitle( "Entries subtraction" ) ;
			// Create MC resta histo
			histoName = "histoIso_MCsub" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType;
			histoTitle = " Distribution for Sherpa-Pythia histograms for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
                        histoIso_SHmenosPy[iPt][iEta][iConv] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "(Sherpa-Pythia distributio histo, pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_SHmenosPy[iPt][iEta][iConv]->SetXTitle( histoXTitle ) ;
                        histoIso_SHmenosPy[iPt][iEta][iConv]->SetYTitle( "Sherpa normalized- Pythia normalized" ) ;	
			// Create test
			histoName = "histoIso_thetesthisto" ;
                        histoName += iPt ;
                        histoName += "_" ;
                        histoName += iEta ;
			histoName += "_" ;
			histoName += iConv ;
			histoName += "_" ;
			histoName += iType ;
			histoTitle = " Xi function for the iPt, iEta bin " ;
			histoTitle += iPt ;
			histoTitle += " , " ;
			histoTitle += iEta ;
			histoTitle += " and for " ;
			histoTitle += Conv[iConv] ;
                        histoIso_test[iPt][iEta][iConv] = new TH1F( histoName , histoTitle , nBins , minIso , maxIso ) ;
                        histoXTitle = "(test histogram, pT,eta = " ;
                        histoXTitle += iPt ;
                        histoXTitle += " , " ;
                        histoXTitle += iEta ;
                        histoXTitle += ") topoETcone40 [GeV]" ;
                        histoIso_test[iPt][iEta][iConv]->SetXTitle( histoXTitle ) ;
                        histoIso_test[iPt][iEta][iConv]->SetYTitle( "Bins" ) ;		
		}
        }
        }
	}
	} 

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
	
       	Int_t iEta = 3 ;		

	/*TCanvas * ptCanvas1 = new TCanvas( "ptCanvas1" , "ptCanvas1" , 6000 , 1600 ) ;
	ptCanvas1->Divide( 5 , 3 ) ; //RowxLine
	
	TCanvas * ptCanvas2 = new TCanvas( "ptCanvas2" , "ptCanvas2" , 6000 , 1600 ) ;
	ptCanvas2->Divide( 5 , 3 );

	TCanvas * ptCanvas3 = new TCanvas( "ptCanvas3" , "ptCanvas3" , 6000 , 1600 ) ;
	ptCanvas3->Divide( 5 , 3 ) ; //RowxLine
	
	TCanvas * ptCanvas4 = new TCanvas( "ptCanvas4" , "ptCanvas4" , 6000 , 1600 ) ;
	ptCanvas4->Divide( 5 , 3 ) ;
	
	TCanvas * ptCanvas5 = new TCanvas( "ptCanvas5" , "ptCanvas5" , 6000 , 1600 ) ;
	ptCanvas5->Divide( 5 , 3 ) ;

	TCanvas * ptCanvas6 = new TCanvas( "ptCanvas6" , "ptCanvas6" , 6000 , 1600 ) ;
	ptCanvas6->Divide( 5 , 3 ) ;

	TCanvas * ptCanvas7 = new TCanvas( "ptCanvas7" , "ptCanvas7" , 6000 , 1600 ) ;
	ptCanvas7->Divide( 5 , 3 ) ;

	TCanvas * ptCanvas8 = new TCanvas( "ptCanvas8" , "ptCanvas8" , 6000 , 1600 ) ;
	ptCanvas8->Divide( 5 , 3 ) ;
*/
	TCanvas * normalize1 = new TCanvas( "normalize1" , "normalize1" , 1200 , 480 ) ;
	normalize1->Divide( 1 , 2 ); /*	

       	Int_t iPad1 = 0 ;

	Int_t iPad2 = 0 ;
	Int_t iPad3 = 0 ;
	Int_t iPad4 = 0 ;
	Int_t iPad5 = 0 ;
	Int_t iPad6 = 0 ;
	Int_t iPad7 = 0 ;
	Int_t iPad8 = 0 ;

	cout << " As example, it is going to show for differet iPt the tight and loose'4 histograms for data " << endl ;
	
	for ( Int_t iPt = 0 ; iPt < nPtBins - 1 ; iPt ++ ){
 	
			// For tight conv
			iPad1 ++ ;
			ptCanvas1->cd( iPad1 ) ;
			histoIso_tight[iPt][iEta][0][0]->SetLineColor(kBlack) ;
			//gPad->SetLogy(1) ;
                	histoIso_tight[iPt][iEta][0][0]->Draw("e") ;
			ptCanvas1->Update() ;

			// For looseprime4 conv
			iPad2 ++ ;
			ptCanvas2->cd( iPad2 ) ;
			histoIso_loose4[iPt][iEta][0][0]->SetLineColor(kBlack) ;
			histoIso_loose4[iPt][iEta][0][0]->Draw("e") ;
			ptCanvas2->Update() ;
	
			// For tight unconv
			iPad3 ++ ;
			ptCanvas3->cd( iPad3 ) ;
			histoIso_tight[iPt][iEta][1][0]->SetLineColor(kBlack) ;
			//gPad->SetLogy(1) ;
                	histoIso_tight[iPt][iEta][1][0]->Draw("e") ;
			ptCanvas3->Update() ;

			// For looseprime4 unconv
			iPad4 ++ ;
			ptCanvas4->cd( iPad4 ) ;
			histoIso_loose4[iPt][iEta][1][0]->SetLineColor(kBlack) ;
			histoIso_loose4[iPt][iEta][1][0]->Draw("e") ;
			ptCanvas4->Update() ;
	
			// For tight conv SH
			iPad5 ++ ;
			ptCanvas5->cd( iPad5 ) ;
			histoIso_tight[iPt][iEta][0][1]->SetLineColor(kRed) ;
                	histoIso_tight[iPt][iEta][0][1]->Draw("e") ;
			ptCanvas5->Update() ;
	
			// For tight unconv SH
			iPad6 ++ ;
			ptCanvas6->cd( iPad6 ) ;
			histoIso_tight[iPt][iEta][1][1]->SetLineColor(kRed) ;
                	histoIso_tight[iPt][iEta][1][1]->Draw("e") ;
			ptCanvas6->Update() ;

			// For tight conv Py
			iPad7 ++ ;
			ptCanvas7->cd( iPad7 ) ;
			histoIso_tight[iPt][iEta][0][2]->SetLineColor(kBlue) ;
                	histoIso_tight[iPt][iEta][0][2]->Draw("e") ;
			ptCanvas7->Update() ;
	
			// For tight unconv Py
			iPad8 ++ ;
			ptCanvas8->cd( iPad8 ) ;
			histoIso_tight[iPt][iEta][1][2]->SetLineColor(kBlue) ;	
                	histoIso_tight[iPt][iEta][1][2]->Draw("e") ;
			ptCanvas8->Update() ;
	}*/
	
	Int_t Internalsuffering = 0 ;
	
	Int_t iPt = 10 ;

	cout << " Now, as an example for iPt, iEta : " << iPt << " , " << iEta << " . It will show some of the histograms. " << endl ;

	TLegend * leg = new TLegend ( 0.1 , 0.7 , 0.48 , 0.9 ) ;
	//TLegend * leg2 = new TLegend ( 0.1 , 0.7 , 0.48 , 0.9 ) ;

	for ( Int_t iConv = 0 ; iConv < 2 ; iConv ++ ){	

		Float_t hightisolationregion = 20 ;// for the region above 20 GeV
	
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
			
			cout << " For " << Conv[iConv] << " photons: " << endl ;
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

			Internalsuffering ++ ;
			normalize1->cd( Internalsuffering ) ;
			histoIso_backgroundsub[iPt][iEta][iConv][0]->SetLineColor(kBlack) ;
			histoIso_backgroundsub[iPt][iEta][iConv][0]->Draw() ;
			TLine * line = new TLine(-10 , 0 , 45 , 0 ) ;
  			line->SetLineColor(kOrange);
  			line->Draw("SAME");
			//Float_t ymax = histoIso_backgroundsub[iPt][iEta][iConv][0]->GetMaximum() ;
			//TLine * line2 = new TLine( 20 , 0 , 20 , ymax ) ;
			//line2->SetLineColor(kGreen) ;
			//line2->Draw("SAME") ;
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

			cout << "\t *The background substracted in " << Type[0] << " have " << bgsubIntegral << " events. <-- " << endl ;
			
			cout << "\t *Tight Monte Carlo in " << Type[1] << " have " << tightSHIntegral << endl ;

			cout << "\t *Tight Monte Carlo in " << Type[2] << " have " << tightPyIntegral << endl ;
			
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

			histoIso_sherpanorm[iPt][iEta][iConv][1]->SetLineColor(kRed) ;
			histoIso_sherpanorm[iPt][iEta][iConv][1]->Draw("SAME") ;	
			leg->SetHeader("Legend") ;
   			leg->AddEntry( histoIso_sherpanorm[iPt][iEta][1][1] ," Normalize Sherpa Histogram ","l") ;	
   			leg->Draw("SAME") ;
			normalize1->Update() ;	
			
			histoIso_pythianorm[iPt][iEta][iConv][2]->SetLineColor(kBlue) ;
			histoIso_pythianorm[iPt][iEta][iConv][2]->Draw("SAME") ;	
   			leg->AddEntry( histoIso_pythianorm[iPt][iEta][1][2] ," Normalize Pythia Histogram ","l") ;	
   			leg->Draw("SAME") ;
			normalize1->Update() ;

		}//thisbigif
	}//For the iConv

	// Data-MC study
	
	cout << " Now, we are going to plot the distributions for Data-MC for these histograms in the previous given iPt, iEta bins as an example. " << endl ;

	TCanvas * sub1 = new TCanvas( "sub1" , "sub1" , 2400 , 1000 ) ;
	sub1->Divide( 1 , 2 );

	TCanvas * sub2 = new TCanvas( "sub2" , "sub2" , 2400 , 1000 ) ;
	sub2->Divide( 1 , 2 ) ;

	Int_t iPad9 = 0 ;

	Int_t iPad10 = 0 ;

	for ( Int_t iConv = 0 ; iConv < 2 ; iConv ++ ){
		if ( histoIso_sherpanorm[iPt][iEta][iConv][1]->Integral() > 0 && histoIso_pythianorm[iPt][iEta][iConv][2]->Integral() > 0 && histoIso_backgroundsub[iPt][iEta][iConv][0]->Integral() > 0 ) {
		
			for ( Int_t iBin = 1 ; iBin <= nBins ; iBin ++ ) {
				Float_t yBinbgsub = histoIso_backgroundsub[iPt][iEta][iConv][0]->GetBinContent( iBin ) ;
				Float_t yBinSHnorm = histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinContent( iBin ) ;
				Float_t yBinPynorm = histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinContent( iBin ) ;	
				
				histoIso_XiSH[iPt][iEta][iConv]->SetBinContent( iBin , (yBinbgsub - yBinSHnorm) ) ;
				histoIso_XiPy[iPt][iEta][iConv]->SetBinContent( iBin , (yBinbgsub - yBinPynorm) ) ; 			
				histoIso_SHmenosPy[iPt][iEta][iConv]->SetBinContent( iBin , (yBinSHnorm - yBinPynorm) ) ;	
					
				Float_t sBinbgsub = histoIso_backgroundsub[iPt][iEta][iConv][0]->GetBinError( iBin ) ;
				Float_t sBinSHnorm = histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinError( iBin ) ;
				Float_t sBinPynorm = histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinError( iBin ) ;

				histoIso_XiSH[iPt][iEta][iConv]->SetBinError( iBin , TMath::Sqrt((sBinbgsub*sBinbgsub) + (sBinSHnorm*sBinSHnorm)) ) ;
				histoIso_XiPy[iPt][iEta][iConv]->SetBinError( iBin , TMath::Sqrt((sBinbgsub*sBinbgsub) + (sBinPynorm*sBinPynorm)) ) ;
			
				histoIso_SHmenosPy[iPt][iEta][iConv]->SetBinError( iBin ,TMath::Sqrt((sBinSHnorm*sBinSHnorm) + (sBinPynorm*sBinPynorm)) ) ;
			}
		
		iPad9 ++ ;
		sub1->cd( iPad9 ) ;
		histoIso_XiSH[iPt][iEta][iConv]->SetLineColor(kRed) ;
		histoIso_XiSH[iPt][iEta][iConv]->Draw() ;
		TLine * line = new TLine(-10 , 0 , 45 , 0 ) ;
  		line->SetLineColor(kOrange);
  		line->Draw("SAME");	
		sub1->Update() ;
		histoIso_XiPy[iPt][iEta][iConv]->SetLineColor(kBlue) ;
		histoIso_XiPy[iPt][iEta][iConv]->Draw("SAME") ;
		histoIso_SHmenosPy[iPt][iEta][iConv]->SetLineColor(kBlack) ;
		histoIso_SHmenosPy[iPt][iEta][iConv]->Draw("SAME") ;
		sub1->Update() ;
		
		iPad10 ++ ;
		sub2->cd( iPad10 ) ;
		histoIso_SHmenosPy[iPt][iEta][iConv]->SetLineColor(kBlack) ;
		histoIso_SHmenosPy[iPt][iEta][iConv]->Draw() ;
		TLine * line2 = new TLine(-10 , 0 , 45 , 0 ) ;
  		line2->SetLineColor(kOrange);
  		line2->Draw("SAME");		
		sub2->Update() ;
		}
	}

	cout << " Then, let's plot the function that depends of the shift. " << endl ;

	// Test for calculing Xi with the Data and the MC	
		
	Int_t iPad11 = 0 ;
	Int_t iPad12 = 0 ;	

	// Check that both histograms have the same binning	
		
	TCanvas * xi = new TCanvas( "xi" , "xi" , 1600 , 800 );
	xi->Divide( 1 , 2 ) ;	
	
	TCanvas * xi2 = new TCanvas( "xi2" , "xi2" , 1600 , 800 );
	xi2->Divide( 1 , 2 ) ;	

	for ( Int_t iConv = 0 ; iConv < 2 ; iConv ++ ){
		
		Float_t Chi2shtest = computeChi2( histoIso_backgroundsub[iPt][iEta][iConv][0] , histoIso_sherpanorm[iPt][iEta][iConv][1] ) ;
		Float_t Chi2pytest = computeChi2( histoIso_backgroundsub[iPt][iEta][iConv][0] , histoIso_pythianorm[iPt][iEta][iConv][2] ) ;
		
		cout << " For " << Conv[iConv] << " photons: " << endl;
		cout << "\t *Data histogram and normalize Sherpa histogram have as a test Chi2 = " << Chi2shtest << endl ;
		cout << "\t *Data histogram and the normalize Pythia histogram have as a test Chi2 = " << Chi2pytest << endl ;
	
		/* In the lines 83, 84 and 85 is define the histoIsomin, histoIsomax and the nBins,
		 just for remember they are: histoMin= -10 , histoMax= 45 and nBins= 170 */

		// Now scan alpha values

		cout << " The content of the bin with more events in Sherpa is: " << histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinContent(histoIso_sherpanorm[iPt][iEta][iConv][1]->GetMaximumBin()) << endl ;
		cout << " The content of the bin with more events in pythia is: " << histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinContent(histoIso_pythianorm[iPt][iEta][iConv][2]->GetMaximumBin()) << endl ;

		Float_t alpha_min = -20 * histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinWidth(1) ; // alpha is two times the bin width
	        Float_t alpha_max = 20 * histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinWidth(1);		

		cout << " The bin width in Sherpa is " << histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinWidth(1) << " GeV " << endl ;
		cout << " The bin width in Pythia is " << histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinWidth(1) << " GeV " << endl ;

		const Int_t nSteps = 41 ;		
       		
		cout << " Will scan " << nSteps << " alpha values from " << alpha_min << " to " << alpha_max << endl ;

		Float_t vector_alphaSH[nSteps] ;// These arrays are needed to construct a TGraph
		Float_t vector_alphapy[nSteps] ;
		Float_t vector_Chi2sh[nSteps] ;
		Float_t vector_Chi2py[nSteps] ;		
		
		Float_t chi2minSH_guess = 9.99e99 ;// Dummy variables to store the "best" scanned values (ad to be used later for the fit)
		Float_t chi2minpy_guess = 9.99e99 ;
		Float_t alpha_guessSH = 0.0f ;	
		Float_t alpha_guesspy = 0.0f ;

		Int_t howmanybins = -20 ;		

		for ( Int_t iAlpha = 0 ; iAlpha < nSteps ; iAlpha ++ ) {// loop over alpha values 
                
			vector_alphaSH[iAlpha] = howmanybins * histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinWidth(1) ;
			vector_alphapy[iAlpha] = howmanybins * histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinWidth(1) ;

      			// Shift histos
			
			TH1F * shiftHSH = shiftHistogram( histoIso_sherpanorm[iPt][iEta][iConv][1] , howmanybins ) ;
			TH1F * shiftHpy = shiftHistogram( histoIso_pythianorm[iPt][iEta][iConv][2] , howmanybins ) ;	
 			
			// Compute Chi2 for the alpha histogram
 			
			Float_t Chi2sh = computeChi2( histoIso_backgroundsub[iPt][iEta][iConv][0] , shiftHSH ) ;	
			Float_t Chi2py = computeChi2( histoIso_backgroundsub[iPt][iEta][iConv][0] , shiftHpy ) ;              		
	
			cout << "\t\t\t ( Alpha GeV , Chi2sh ) = " << howmanybins * shiftHSH->GetBinWidth(1) << " , " << Chi2sh << endl ;	
			cout << "\t\t\t ( Alpha GeV , Chi2py ) = " << howmanybins * shiftHpy->GetBinWidth(1) << " , " << Chi2py << endl ;
			
			vector_Chi2sh[iAlpha] = Chi2sh ;// Store Chi2 value in array
			vector_Chi2py[iAlpha] = Chi2py ;
 
                	if ( Chi2sh < chi2minSH_guess ) {
                      		chi2minSH_guess = Chi2sh ;
                      		alpha_guessSH = howmanybins * histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinWidth(1) ;
               		}

			if ( Chi2py < chi2minpy_guess ) {
                      		chi2minpy_guess = Chi2sh ;
                      		alpha_guesspy = howmanybins * histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinWidth(1) ;
               		}
 		
		howmanybins += 1 ;
                 	
		} //alphafor
		
		// TGraph Chi2 vs alpha

		TGraph * gr_Chi2alphaSH = new TGraph( nSteps , vector_alphaSH , vector_Chi2sh ) ;
		TGraph * gr_Chi2alphapy = new TGraph( nSteps , vector_alphapy , vector_Chi2py ) ;		
		
		// Parabolic Function

		TF1 * fitChi2SH = new TF1( "fitChi2SH" , "[0]+TMath::Power((x-[1])/[2],2)" , alpha_min , alpha_max ) ;
		fitChi2SH->SetParName( 0 , "Chi2minSH" ) ;
		fitChi2SH->SetParName( 1 , "FitValueSH" ) ;
		fitChi2SH->SetParName( 2 , "FitErrorSH" ) ;

		TF1 * fitChi2py = new TF1( "fitChi2py" , "[0]+TMath::Power((x-[1])/[2],2)" , alpha_min , alpha_max ) ;
		fitChi2py->SetParName( 0 , "Chi2minpy" ) ;
		fitChi2py->SetParName( 1 , "FitValuepy" ) ;
		fitChi2py->SetParName( 2 , "FitErrorpy" ) ;

		// Estimation for the parameter CHECK CHECK CHECK CHECK COME BACK HERE AGAIN

		Float_t deltaChi2SH = vector_Chi2sh[0] - chi2minSH_guess ;
		if ( deltaChi2SH < 0 ) deltaChi2SH = - deltaChi2SH ;
		Float_t error_guessSH = ( alpha_guessSH - alpha_min ) / TMath::Sqrt( deltaChi2SH ) ;
		fitChi2SH->SetParameter( 0 , chi2minSH_guess ) ;
		fitChi2SH->SetParameter( 1 , alpha_guessSH ) ;
		fitChi2SH->SetParameter( 2 , error_guessSH ) ;

		cout << "Fit test in Sherpa using chi2min = " << chi2minSH_guess << " , alpha = " << alpha_guessSH << " +/- " << error_guessSH <<  " for " << Conv[iConv] << " photons. " << endl ;
		gr_Chi2alphaSH->Fit( "fitChi2SH" ) ; 

		Float_t deltaChi2py = vector_Chi2py[0] - chi2minpy_guess ;
		if ( deltaChi2py < 0 ) deltaChi2py = - deltaChi2py ;
		Float_t error_guesspy = ( alpha_guesspy - alpha_min ) / TMath::Sqrt( deltaChi2py ) ;
		fitChi2py->SetParameter( 0 , chi2minpy_guess ) ;
		fitChi2py->SetParameter( 1 , alpha_guesspy ) ;
		fitChi2py->SetParameter( 2 , error_guesspy ) ;

		cout << "Fit test in Pythia using chi2min = " << chi2minpy_guess << " , alpha = " << alpha_guesspy << " +/- " << error_guesspy << " for " << Conv[iConv] << " photons. " << endl ;
		gr_Chi2alphapy->Fit( "fitChi2py" ) ;

		// Rerun the fit over a restricted range (only if it is necessary)

		Float_t alpha_fitSH = fitChi2SH->GetParameter( 1 ) ;
		Float_t error_fitSH = fitChi2SH->GetParameter( 2 ) ;
		
		Bool_t redoFitSH = kFALSE ;

		if ( alpha_fitSH - 20 * error_fitSH > alpha_min ){
			alpha_min = alpha_fitSH - 20 * error_fitSH ;
			redoFitSH = kTRUE ;
		}

		if ( alpha_fitSH + 20 * error_fitSH < alpha_max ){
			alpha_max = alpha_fitSH + 20 * error_fitSH ;
			redoFitSH = kTRUE ;
		}
		if ( redoFitSH) gr_Chi2alphaSH->Fit( "fitChi2SH" , "" , "" , alpha_min , alpha_max ) ;

		Float_t alpha_fitpy = fitChi2py->GetParameter( 1 ) ;
		Float_t error_fitpy = fitChi2py->GetParameter( 2 ) ;
		
		Bool_t redoFitpy = kFALSE ;

		if ( alpha_fitpy - 20 * error_fitpy > alpha_min ){
			alpha_min = alpha_fitpy - 20 * error_fitpy ;
			redoFitpy = kTRUE ;
		}

		if ( alpha_fitpy + 20 * error_fitpy < alpha_max ){
			alpha_max = alpha_fitpy + 20 * error_fitpy ;
			redoFitpy = kTRUE ;
		}
		if ( redoFitpy) gr_Chi2alphapy->Fit( "fitChi2py" , "" , "" , alpha_min , alpha_max ) ;

		// Draw	
		
		iPad11 ++ ;
		xi->cd ( iPad11 ) ;

		TString histoTitle1 = "  Chi2(alpha) vs alpha for " ;	
		histoTitle1 += Conv[iConv] ;
		histoTitle1 += " photons between data and Sherpa. " ;
		gr_Chi2alphaSH->SetTitle( histoTitle1 ) ;	
		gr_Chi2alphaSH->SetMarkerColor( kBlue ) ;
		gr_Chi2alphaSH->SetMarkerStyle( 7 ) ;
		gr_Chi2alphaSH->SetMarkerColor( 4 ) ;
		gr_Chi2alphaSH->GetXaxis()->SetTitle( "" ) ;
		gr_Chi2alphaSH->GetYaxis()->SetTitle( "# chi^2 value2 " ) ;
		gr_Chi2alphaSH->Draw( "AP" ) ;
		xi->Update() ;

		iPad12 ++ ;
		xi2->cd ( iPad12 ) ;

		TString histoTitle2 = "  Chi2(alpha) vs alpha for " ;	
		histoTitle2 += Conv[iConv] ;
		histoTitle2 += " photons between data and Pythia. " ;
		gr_Chi2alphapy->SetTitle( histoTitle2 ) ;	
		gr_Chi2alphapy->SetMarkerColor( kBlue ) ;
		gr_Chi2alphapy->SetMarkerStyle( 7 ) ;
		gr_Chi2alphapy->SetMarkerColor( 4 ) ;
		gr_Chi2alphapy->GetXaxis()->SetTitle( "" ) ;
		gr_Chi2alphapy->GetYaxis()->SetTitle( "# chi^2 value2 " ) ;
		gr_Chi2alphapy->Draw( "AP" ) ;
		xi2->Update() ;

		// José's print out
		
		cout << "The Chi2 vs. alpha trend was fitted with a parabolic function with parameters : " << endl ;
	        cout << "\tChi2Min2 = " << fitChi2SH->GetParameter( 0 ) << endl ;
       		cout << "\tAlpha = " << fitChi2SH->GetParameter( 1 ) << endl ;
       	 	cout << "\tAlpha_error = " << fitChi2SH->GetParameter( 2 ) << endl ;

		cout << "The Chi2 vs. alpha trend was fitted with a parabolic function with parameters : " << endl ;
	        cout << "\tChi2Min2 = " << fitChi2py->GetParameter( 0 ) << endl ;
       		cout << "\tAlpha = " << fitChi2py->GetParameter( 1 ) << endl ;
       	 	cout << "\tAlpha_error = " << fitChi2py->GetParameter( 2 ) << endl ;	
		
		cout << " This should appeard 2 times " << endl ;

	}//for conv

	// Shifted histo for check

	Int_t iPad13 = 0 ;
	Int_t iPad14 = 0 ;
	Int_t iPad15 = 0 ;
	Int_t iPad16 = 0 ;

	Int_t howmanybins2 = -20 ;
	Int_t howmanybins3 = -20 ;
	const Int_t nSteps = 51 ;

	TCanvas * shiftedhistosSHconv = new TCanvas( "shiftedhistosSHconv" , "shiftedhistosSHconv" , 2000 , 2000 ) ;
	shiftedhistosSHconv->Divide( 5 , 5 ) ;

	TCanvas * shiftedhistospyconv = new TCanvas( "shiftedhistospyconv" , "shiftedhistospyconv" , 2000 , 2000 ) ;
	shiftedhistospyconv->Divide( 5 , 5 ) ;

	TCanvas * shiftedhistosSHunconv = new TCanvas( "shiftedhistosSHunconv" , "shiftedhistosSHunconv" , 2000 , 2000 ) ;
	shiftedhistosSHunconv->Divide( 5 , 5 ) ;

	TCanvas * shiftedhistospyunconv = new TCanvas( "shiftedhistospyunconv" , "shiftedhistospyunconv" , 2000 , 2000 ) ;
	shiftedhistospyunconv->Divide( 5 , 5 ) ;

	for ( Int_t iAlpha = 0 ; iAlpha < nSteps ; iAlpha ++ ) {// loop over alpha values 
		
		iPad13 ++ ;

		TH1F * shiftHSH = shiftHistogram( histoIso_sherpanorm[iPt][iEta][0][1] , howmanybins2 ) ;	
 				
		shiftedhistosSHconv->cd( iPad13 ) ;

		TString histoTitle3 = "Sherpa shift(red) respect data(black)_" ;
                histoTitle3 += howmanybins2 ;
		histoTitle3 += "alpha GeV in " ;	
                histoTitle3 += Conv[0] ;
                histoTitle3 += " photons." ;
                histoIso_backgroundsub[iPt][iEta][0][0]->SetTitle( histoTitle3 ) ;
                histoIso_backgroundsub[iPt][iEta][0][0]->DrawNormalized() ;
		shiftedhistosSHconv->Update() ;		
	
		shiftHSH->SetLineColor( kRed ) ;
                shiftHSH->DrawNormalized( "same" ) ;
                shiftedhistosSHconv->Update() ;
                 
		howmanybins2 += 1 ;	
	} //iAlphafor 	

	for ( Int_t iAlpha = 0 ; iAlpha < nSteps ; iAlpha ++ ) {// loop over alpha values 
		
		iPad14 ++ ;

		TH1F * shiftHpy = shiftHistogram( histoIso_pythianorm[iPt][iEta][0][2] , howmanybins3 ) ;	
 				
		shiftedhistospyconv->cd( iPad14 ) ;

		TString histoTitle4 = "Pythia shift(red) respect data(black)_" ;
               	histoTitle4 += howmanybins3 ;
		histoTitle4 += "alpha GeV in " ;	
                histoTitle4 += Conv[0] ;
                histoTitle4 += " photons." ;
                histoIso_backgroundsub[iPt][iEta][0][0]->SetTitle( histoTitle4 ) ;
                histoIso_backgroundsub[iPt][iEta][0][0]->DrawNormalized() ;
		shiftedhistospyconv->Update() ;		
	
		shiftHpy->SetLineColor( kRed ) ;
                shiftHpy->DrawNormalized( "same" ) ;
                shiftedhistospyconv->Update() ;
                 
		howmanybins3 += 1 ;	
	} //iAlphafor

	Int_t howmanybins4 = -20 ;
	Int_t howmanybins5 = -20 ;

	for ( Int_t iAlpha = 0 ; iAlpha < nSteps ; iAlpha ++ ) {// loop over alpha values 
		
		iPad15 ++ ;

		TH1F * shiftHSH = shiftHistogram( histoIso_sherpanorm[iPt][iEta][1][1] , howmanybins4 ) ;	
 				
		shiftedhistosSHunconv->cd( iPad15 ) ;

		TString histoTitle5 = "Sherpa shift(red) respect data(black)_" ;
                histoTitle5 += howmanybins4 ;
		histoTitle5 += "alpha GeV in " ;	
                histoTitle5 += Conv[1] ;
                histoTitle5 += " photons." ;
                histoIso_backgroundsub[iPt][iEta][1][0]->SetTitle( histoTitle5 ) ;
                histoIso_backgroundsub[iPt][iEta][1][0]->DrawNormalized() ;
		shiftedhistosSHunconv->Update() ;		
	
		shiftHSH->SetLineColor( kRed ) ;
                shiftHSH->DrawNormalized( "same" ) ;
                shiftedhistosSHunconv->Update() ;
                 
		howmanybins4 += 1 ;	
	} //iAlphafor 	

	for ( Int_t iAlpha = 0 ; iAlpha < nSteps ; iAlpha ++ ) {// loop over alpha values 
		
		iPad16 ++ ;

		TH1F * shiftHpy = shiftHistogram( histoIso_pythianorm[iPt][iEta][1][2] , howmanybins5 ) ;	
 				
		shiftedhistospyunconv->cd( iPad16 ) ;

		TString histoTitle6 = "Pythia shift(red) respect data(black)_" ;
                histoTitle6 += howmanybins5 ;	
                histoTitle6 += "alpha GeV in " ;
                histoTitle6 += Conv[1] ;
                histoTitle6 += " photons." ;
                histoIso_backgroundsub[iPt][iEta][1][0]->SetTitle( histoTitle6 ) ;
                histoIso_backgroundsub[iPt][iEta][1][0]->DrawNormalized() ;
		shiftedhistospyunconv->Update() ;		
	
		shiftHpy->SetLineColor( kRed ) ;
                shiftHpy->DrawNormalized( "same" ) ;
                shiftedhistospyunconv->Update() ;
                 
		howmanybins5 += 1 ;	
	} //iAlphafor

	/*// alpha min vs pt in function of eta

	Float_t vector_ptbins[nPtBins - 1] ;
	Float_t vector_alphamin[nPtBins - 1] ;

	const Int_t Steps = nPtBins - 1 ;

	Float_t vector_alphaminSH = 0.0f ;	
	Float_t vector_alphaminpy = 0.0f ;
	Float_t chi2minSH_guess = 9.99e99 ;
	Float_t chi2minpy_guess = 9.99e99 ;

	Int_t binsshift = -4 ;

	for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
                if ( iEta != 2 ) {
		for ( Int_t iPt = 0 ; iPt < nPtBins - 1 ; iPt ++ ){
			vector_ptbins[iPt] = iPt ;
			for ( Int_t iConv = 0 ; iConv < 2 ; iConv ++ ){		

			for ( Int_t iAlpha = 0 ; iAlpha < Steps ; iAlpha ++ ) {// loop over alpha values 
                
      			// Shift histos
			
			TH1F * shiftHSH = shiftHistogram( histoIso_sherpanorm[iPt][iEta][iConv][1] , binsshift ) ;
			TH1F * shiftHpy = shiftHistogram( histoIso_pythianorm[iPt][iEta][iConv][2] , binsshift ) ;	
 			
			// Compute Chi2 for the alpha histogram
 			
			Float_t Chi2sh = computeChi2( histoIso_backgroundsub[iPt][iEta][iConv][0] , shiftHSH ) ;	
			Float_t Chi2py = computeChi2( histoIso_backgroundsub[iPt][iEta][iConv][0] , shiftHpy ) ;              		
	 
                	if ( Chi2sh < chi2minSH_guess ) {
                      		chi2minSH_guess = Chi2sh ;
                      		vector_alphaminSH[iAlpha] = binsshift * histoIso_sherpanorm[iPt][iEta][iConv][1]->GetBinWidth(1) ;
               		}

			if ( Chi2py < chi2minpy_guess ) {
                      		chi2minpy_guess = Chi2sh ;
                      		vector_alphaminpy[iAlpha] = binsshift * histoIso_pythianorm[iPt][iEta][iConv][2]->GetBinWidth(1) ;
               		}
 		
		binsshift += 1 ;
                 	
			} //alphafor		

		// TGraph Alphamin vs pt

		TGraph * gr_alphaminptSH = new TGraph( Steps , vector_alphaminSH , vector_ptbins ) ;
		TGraph * gr_alphaminptpy = new TGraph( Steps , vector_alphaminpy , vector_ptbins ) ;

		cout << " Where the alphamins are " << endl ;
		cout << "\t * For iPt " << iPt << endl ;
		cout << "\t\t ** For iEta " << iEta << endl ;
		cout << "\t\t *** For iConv " << iConv << endl ;
		cout << " \t\t\t First " << vector_alphaminSH[0] << " in Sherpa. " << endl ;
		cout << " \t\t\t First " << vector_alphaminpy[0] << " in pythia. " << endl;
		cout << " \t\t\t Last " << vector_alphaminSH[nPtbins - 1] << " in Sherpa. " << endl ;
		cout << " \t\t\t Last " << vector_alphaminpy[nPtbins - 1] << " in pythia. " << endl ;
		
			}
		}
		}
	}*/ 

};

	
