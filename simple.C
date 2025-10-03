// Author: Ron Belmont
// Date: 2025-09-30

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>

#include <sys/time.h>

#include "test.h"
#include "flow_functions.h"

//#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;



const float pi = 3.141592653589793;



Long64_t nevents = 0;
Long64_t ntracks = 0;




// Main part of program
int main(int argc, char *argv[])
{

  //Char_t inFile[100];
  //Char_t outFile[100];
  char inFile[100];
  char outFile[100];

  if(argc==1)
    {
      cout<<"Now beginning program"<<endl;
      cout<<"Program name is "<<argv[0]<<endl;
      cout<<"Please enter input file name"<<endl;
      cin>>inFile;
      cout<<"Input file is "<<inFile<<endl;
      cout<<"Please enter output file name"<<endl;
      cin>>outFile;
      cout<<"Output file is "<<outFile<<endl;
    }
  else if(argc==3)
    {
      strcpy(inFile,argv[1]);
      strcpy(outFile,argv[2]);
      cout<<"Now beginning program"<<endl;
      cout<<"Program name is "<<argv[0]<<endl;
      cout<<"Input file is "<<inFile<<endl;
      cout<<"Output file is "<<outFile<<endl;
    }
  else
    {
      cout<<"Wrong number of input arguments"<<endl;
      cout<<"This program takes 0 or 2 arguments"<<endl;
      cout<<"With 0 arguments it prompts the user for the file list and output root file"<<endl;
      cout<<"With 2 arguments the first is the file list and the second is the output root file"<<endl;
      return 1;
    }

  // ----------------------------

  struct timeval Time;

  gettimeofday(&Time,0);
  int begintime = Time.tv_sec;
  cout<<"begintime is "<<begintime<<endl;

  // ----------------------------

  TFile *mData = new TFile(outFile,"recreate"); // declare output file



  // -------------------------- //
  // --- Declare histograms --- //
  // -------------------------- //

  TH1D *hhcent = new TH1D("hhcent","cent",100,-0.5,99.5);
  TH1D *hhmbdz = new TH1D("hhmbdz","mbdz",100,-50,50);

  TH1D* th1d_sepd_all_e = new TH1D("th1d_sepd_all_e","",100,-1,9);
  TH1D* th1d_sepd_all_r = new TH1D("th1d_sepd_all_r","",100,0,100);
  TH1D* th1d_sepd_all_phi = new TH1D("th1d_sepd_all_phi","",630,-6.3,6.3);
  TH1D* th1d_sepd_all_arm = new TH1D("th1d_sepd_all_arm","",5,-2.5,2.5);
  TH1D* th1d_sepd_all_good = new TH1D("th1d_sepd_all_good","",5,-2.5,2.5);


  // ---------------------------- //
  // --- Done with Histograms --- //
  // ---------------------------- //



  // --- Now read in the pDSTs listed in the input files

  int ifile=0;
  char filename[100];
  ifstream fin(inFile);
  if(!fin)
    {
      cout<<"list input error: file does not exist "<<endl;
      return 1;
    }
  else cout << "Successfully opened " << inFile << endl;
  while(fin.getline(filename,100))
    {

      ifile++;
      cout<<ifile<<" "<<filename<<endl;

      TFile *f = TFile::Open(filename);
      if(!f)
	{
	  cout<<"pDST input error: file does not exist "<<endl;
	  continue;
	}
      else cout << "Successfully opened " << filename << endl;

      //nevents += (Long64_t)((TH1F *)f->Get("hcent"))->GetEntries();

      //TTree *t=(TTree *)f->Get("hadrontree");
      TTree *t=(TTree *)f->Get("T");
      if(!t)
	{
	  cout<<"pDST input error: cannot find tree "<<endl;
	  continue;
	}
      else cout << "Successfully got the tree" << endl;


      int nevt = (int)t->GetEntries(); // number of events in tree
      test *tree = new test(t); // pointer to tree
      for ( int ievt = 0; ievt < nevt; ++ievt ) // loop over events
	{

          // very stupid event counter
          ++nevents;

	  tree->GetEntry(ievt);

	  float mbdz = tree->zvrtx;
	  float cent = tree->centrality;

	  hhcent->Fill(cent);
	  hhmbdz->Fill(mbdz);

          //if ( ievt < 100 ) cout << cent << endl;
          //if ( cent > 0 ) cout << cent << endl;

          std::vector<double> all_phi;
          std::vector<double> south_phi;
          std::vector<double> north_phi;
          std::vector<std::pair<double,double>> all_phi_e;
          std::vector<std::pair<double,double>> south_phi_e;
          std::vector<std::pair<double,double>> north_phi_e;
          for ( int iepd = 0; iepd < 744; ++iepd )
            {
              float e = tree->sepd_energy[iepd];
              float r = tree->sepd_radius[iepd];
              float phi = tree->sepd_phi[iepd];
              int arm = tree->sepd_arm[iepd];
              int good = tree->sepd_isgood[iepd];

              th1d_sepd_all_e->Fill(e);
              th1d_sepd_all_r->Fill(r);
              th1d_sepd_all_phi->Fill(phi);
              th1d_sepd_all_arm->Fill(arm);
              th1d_sepd_all_good->Fill(good);

              all_phi.push_back(phi);
              all_phi_e.push_back(std::make_pair(phi,e));
              if ( arm == 0 )
                {
                  south_phi.push_back(phi);
                  south_phi_e.push_back(std::make_pair(phi,e));
                }
              if ( arm == 1 )
                {
                  north_phi.push_back(phi);
                  north_phi_e.push_back(std::make_pair(phi,e));
                }

            }

          TComplex all_Q2 = flow_functions::get_flow_vector(all_phi,2);
          TComplex south_Q2 = flow_functions::get_flow_vector(south_phi,2);
          TComplex north_Q2 = flow_functions::get_flow_vector(north_phi,2);
          TComplex all_weighted_Q2 = flow_functions::get_weighted_flow_vector(all_phi_e,2);
          TComplex south_weighted_Q2 = flow_functions::get_weighted_flow_vector(north_phi_e,2);
          TComplex north_weighted_Q2 = flow_functions::get_weighted_flow_vector(north_phi_e,2);

          if ( ievt < 10 )
            {
              cout << "all Q2 " << all_Q2 << endl;
              cout << "south Q2 " << south_Q2 << endl;
              cout << "north Q2 " << north_Q2 << endl;
              cout << "all weighted Q2 " << all_weighted_Q2 << endl;
              cout << "south weighted Q2 " << south_weighted_Q2 << endl;
              cout << "north weighted Q2 " << north_weighted_Q2 << endl;
            }

	} // End of event loop

      t->Delete();
      delete tree;
      f->Close();
      delete f;

    } // End of pDST loop

  // write the histograms to file and close
  mData->Write();
  mData->Close();

  cout<<"Number of events: "<<nevents<<endl;

  gettimeofday(&Time,0);
  int endtime = Time.tv_sec;
  //cout<<"endtime is "<<endtime<<endl;

  int tdiff = endtime-begintime;

  cout<<"End of program."<<endl;
  cout<<"Execution time: "<<tdiff<<" seconds"<<endl;

  exit(0);

} // end of main

