// Author: Ron Belmont
// Date: 2009-07-28

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>

#include <sys/time.h>

#include "test.h"

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

  TH1F *hhcent = new TH1F("hhcent","cent",100,0,100); // doesn't work, event bias
  TH1F *hhbbcz = new TH1F("hhbbcz","bbcz",100,-50,50);



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


      int n=(int)t->GetEntries(); // number of events in tree
      //hadrontree *ktree = new hadrontree(t); // pointer to tree
      test *ktree = new test(t); // pointer to tree
      for(int i=0;i<n;i++) // loop over events
	{

	  ktree->GetEntry(i);

	  float bbcz = ktree->zvrtx;
	  float cent = ktree->centrality;

	  hhcent->Fill(cent);
	  hhbbcz->Fill(bbcz);

	} // End of event loop

      t->Delete();
      delete ktree;
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

