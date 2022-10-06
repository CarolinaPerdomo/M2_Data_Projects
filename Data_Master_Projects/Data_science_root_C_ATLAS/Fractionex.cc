{

	gROOT->LoadMacro( "studySinglePhotons.cc" ) ; //Call1stProgram

	const Int_t nEtaBins = 6 ;

	Float_t etaBins[nEtaBins] = { 0 , 0.6 , 1.37 , 1.52 , 1.81 , 2.37 } ; //absoetabins

	const Int_t nPtBins = 16 ;

	Float_t pTbins[nPtBins] = { 25 , 30 , 40 , 45 , 50 , 55 , 65 , 75 , 85 , 105 , 125 , 145 , 205 , 300 , 500 , 999 } ; //pTbins

	Float_t myFraction[nPtBins][nEtaBins] ;
        Float_t myError[nPtBins][nEtaBins] ; //Errors
        Float_t myFr[1] ;
        Float_t myErr[1] ; //And for eta?? heuheuheu

       	for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
		if ( iEta != 2 ) {
			for ( Int_t iPt = 0 ; iPt < nPtBins - 1 ; iPt ++ ){
				myFr[0] = -9999 ;
				myErr[0] = -9999 ; 
				studySinglePhotons( myFr , myErr , pTbins[iPt] , pTbins[iPt+1] , etaBins[iEta] , etaBins[iEta+1] ) ;
				cout << "After calling studySinglePhotons for iPt,iEta = " << iPt << " , " << iEta << " now myFraction = " << myFr[0] << " +/- " << myErr[0] << endl ;
				myFraction[iPt][iEta] = myFr[0] ;
				myError[iPt][iEta] = myErr[0] ;
			}
		}
	}
	
	TH1F * histoFrac[nEtaBins] ;
       	for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
		if ( iEta != 2 ) {
			TString histoName = "histoFrac_" ;
			histoName += iEta ;
			TString histoTitle = "Fraction of events vs pT for eta bin " ;
			histoTitle += iEta ;
			histoFrac[iEta] = new TH1F( histoName , histoTitle , nPtBins - 1 , pTbins ) ;
		}
	}

//Something wrong in this loop 
        for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
		if ( iEta != 2 ) {
			cout << "For eta [ " << etaBins[iEta] << " : " << etaBins[iEta+1] << "]" << endl ;
			for ( Int_t iPt = 0 ; iPt < nPtBins - 1 ; iPt ++ ) {
				cout << "*)pT bin [ " <<  pTbins[iPt] << " : " << pTbins[iPt+1] << " ] found fraction = " << myFraction[iPt][iEta] << " +/- " << myError[iPt][iEta] <<  endl ;
				histoFrac[iEta]->SetBinContent( iPt + 1 , myFraction[iPt][iEta] ) ;
				histoFrac[iEta]->SetBinError( iPt + 1 , myError[iPt][iEta] ) ;
				histoFrac[iEta]->SetXTitle( "pT" ) ;
				histoFrac[iEta]->SetYTitle( "Fraction of Events < 5GeV for (tight) photons" ) ;
			}
		}
	 }

        TCanvas * aCanvas = new TCanvas() ;
	Int_t iPad = 0 ;
	aCanvas->Divide( 2 , 2 ) ;
        for ( Int_t iEta = 0 ; iEta < nEtaBins - 1 ; iEta ++ ){
		if ( iEta != 2 ) {
			iPad ++ ;
			aCanvas->cd( iPad ) ;
			histoFrac[iEta]->Draw() ;
		}
	}
	

        aCanvas->Print( "FractionvspT.png" ) ;
}
