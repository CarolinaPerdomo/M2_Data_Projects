{

	gROOT->LoadMacro( "studySinglePhotons.cc" ) ; //Call1stProgram

	Float_t eta[7] = { 0 , 0.6 , 1.37 , 1.52 , 1.81 , 2.37 , 99999 } ; //absoetabins

	Float_t pTbins[16] = { 25 , 30 , 40 , 45 , 50 , 55 , 65 , 75 , 85 , 105 , 125 , 145 , 205 , 300 , 500 , 99999 } ; //pTbins

	Float_t myFraction[14] ;
        Float_t myError[14] ; //Errors
        Float_t myFr[1] ;
        Float_t myErr[1] ; //And for eta?? "< heuheuheu

	//for ( Int_t iPt = 0 ; iPt < 5 ; iPt ++ ) {
       	for ( Int_t jeta = 0 ; jeta < 5 ; jeta ++ ){
		for ( Int_t iPt = 0 ; iPt < 14 ; iPt ++ ){
			myFr[0] = -9999 ;
                	myErr[0] = -9999 ; 
                		studySinglePhotons( myFr , myErr , pTbins[iPt] , pTbins[iPt+1] , eta[jeta] , eta[jeta+1] ) ;
                cout << "After calling studySinglePhotons, now myFraction = " << myFr[0] << " +/- " << myErr[0] << endl ;
                myFraction[iPt] = myFr[0] ;
                myError[iPt] = myErr[0] ;
		}
	}
	TH1F * histoFrac = new TH1F( "histoFrac" , "histoFrac" , 100 , 0.5 , 100.5 ) ;

        for ( Int_t jeta = 0 ; jeta < 5 ; jeta ++ ){
		cout << "For eta [ " << eta[jeta] << " : " << eta[jeta+1] << "]" << endl ;
			for ( Int_t iPt = 0 ; iPt < 14 ; iPt ++ ) {
				cout << "pT bin [ " <<  pTbins[iPt] << " : " << pTbins[iPt+1] << " ] found fraction = " << myFraction[iPt] << " +/- " << myError[iPt] <<  endl ;

			histoFrac->SetBinContent( pTbins[iPt] , myFraction[iPt] ) ;
			histoFrac->SetBinError( pTbins[iPt] , myError[iPt] ) ;
       			}
	 }

        TCanvas * aCanvas = new TCanvas() ;
        histoFrac->Draw() ;

        aCanvas->Print( "FractionvspT.png" ) ;
}
