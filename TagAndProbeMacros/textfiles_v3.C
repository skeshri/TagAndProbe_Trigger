#include <iostream>

void textfiles_v3(){

  string file;
  std::ofstream outfile;

  outfile.open("Ele35_pt_eta_efficiency_withSys_Run2017B.txt");
  outfile<<std::fixed;

  ifstream nominal;
  nominal.open("Ele35_pt_eta_EGM_2017B_nominal_efficiency.txt");

  ifstream zmassup;
  zmassup.open("Ele35_pt_eta_EGM_2017B_Zmass_up_efficiency.txt");

  ifstream zmassdown;
  zmassdown.open("Ele35_pt_eta_EGM_2017B_Zmass_down_efficiency.txt");

  ifstream tagup;
  tagup.open("Ele35_pt_eta_EGM_2017B_TagPt_up_efficiency.txt");

  ifstream tagdown;
  tagdown.open("Ele35_pt_eta_EGM_2017B_TagPt_up_efficiency.txt");


//  const int nptbins = 5;
//  const int netabins = 5;
//  double ptmin[nptbins+1] = {10.,20.,30.,40.,50.,2000.};
//  double etamin[netabins+1] = {-2.5,-2.0,-1.566,-1.4442,-0.8,0,0.8,1.4442,1.566,2.0,2.5};
// double etamin[netabins+1] = {0,0.8,1.4442,1.566,2.0,2.5};  

  double absdiffLO,absdiff_fitR,absdiffExp,absdiffCB,absdiffTagS;
  while(!nominal.eof()){

    double detamin, detamax, dptmin, dptmax, eff, errlow, errhigh;
    double zupdetamin, zupdetamax, zupdptmin, zupdptmax, zupeff, zuperrlow, zuperrhigh;
    double zdowndetamin, zdowndetamax, zdowndptmin, zdowndptmax, zdowneff, zdownerrlow, zdownerrhigh;
    double tagdowndetamin, tagdowndetamax, tagdowndptmin, tagdowndptmax, tagdowneff, tagdownerrlow, tagdownerrhigh;
    double tagupdetamin, tagupdetamax, tagupdptmin, tagupdptmax, tagupeff, taguperrlow, taguperrhigh;

    double defftagup, defftagdown, deffzup, deffzdown;
    double toterrup, toterrdown;

    nominal>>detamin>>detamax>>dptmin>>dptmax>>eff>>errlow>>errhigh;
    zmassup>>zupdetamin>>zupdetamax>>zupdptmin>>zupdptmax>>zupeff>>zuperrlow>>zuperrhigh;
    zmassdown>>zdowndetamin>>zdowndetamax>>zdowndptmin>>zdowndptmax>>zdowneff>>zdownerrlow>>zdownerrhigh;
    tagup>>tagupdetamin>>tagupdetamax>>tagupdptmin>>tagupdptmax>>tagupeff>>taguperrlow>>taguperrhigh;
    tagdown>>tagdowndetamin>>tagdowndetamax>>tagdowndptmin>>tagdowndptmax>>tagdowneff>>tagdownerrlow>>tagdownerrhigh;    



     defftagup = fabs(eff-tagupeff);
     defftagdown = fabs(eff-tagdowneff);
     deffzup = fabs(eff-zupeff);                    
     deffzdown = fabs(eff-zdowneff);
 
    toterrup = sqrt(errhigh*errhigh + defftagup*defftagup + deffzup*deffzup);
    toterrdown = sqrt(errlow*errlow + defftagdown*defftagdown + deffzdown*deffzdown);

outfile<<setprecision(3)<<detamin<<"\t"<<detamax<<"\t"<<dptmin<<"\t"<<dptmax<<"\t"<<eff<<"\t"<<toterrdown<<"\t"<<toterrup<<endl;


/*
//    mcfile>>mcptmin>>mcptmax>>mcetamin>>mcetamax>>mceff>>mcerr;
     mcfile>>mceff>>mcerr;    
    //cout<<mcptmin<<"  "<<mcptmax<<"  "<<mcetamin<<"  "<<mcetamax<<"  "<<mceff<<"  "<<mcerr<<endl;
 //   cout << mceff << " " << mcerr << endl;    

//    datafileLO>>dptminLO>>dptmaxLO>>detaminLO>>detamaxLO>>deffLO>>derrLO;
    //cout<<dptmin<<"  "<<dptmax<<"  "<<detamin<<"  "<<detamax<<"  "<<deff<<"  "<<derr<<endl;
    
    mcfileLO>>mceffLO>>mcerrLO;    
    //cout<<mcptmin<<"  "<<mcptmax<<"  "<<mcetamin<<"  "<<mcetamax<<"  "<<mceff<<"  "<<mcerr<<endl;
    
    datafilefitR>>dptminfitR>>dptmaxfitR>>detaminfitR>>detamaxfitR>>defffitR>>derrfitR;
    //cout<<dptmin<<"  "<<dptmax<<"  "<<detamin<<"  "<<detamax<<"  "<<deff<<"  "<<derr<<endl; 

    datafileExp>>dptminExp>>dptmaxExp>>detaminExp>>detamaxExp>>deffExp>>derrExp;
    //cout<<dptmin<<"  "<<dptmax<<"  "<<detamin<<"  "<<detamax<<"  "<<deff<<"  "<<derr<<endl;
    
    datafileCB>>dptminCB>>dptmaxCB>>detaminCB>>detamaxCB>>deffCB>>derrCB;
    

    datafileTagS>>dptminTagS>>dptmaxTagS>>detaminTagS>>detamaxTagS>>deffTagS>>derrTagS;
//    mcfileTagS>>mcptminTagS>>mcptmaxTagS>>mcetaminTagS>>mcetamaxTagS>>mceffTagS>>mcerrTagS;
i*/
/*
    absdiffLO = fabs((deff/mceff) - (deffLO/mceffLO));
    absdiff_fitR = fabs((deff/mceff) - (defffitR/mcefffitR));
    absdiffExp = fabs((deff/mceff) - (deffExp/mceff));
    absdiffCB = fabs((deff/mceff) - (deffCB/mceff));
    absdiffTagS = fabs((deff/mceff) - (deffTagS/mceff));
*/
/*
//minEta   maxEta   minPt   maxPt   effData    statError effMC  statError   systBkgShape    systSigShape   systFitRange systNLOvsLO(or whatever 2 MCs you have)    systPU    systTagSelection

//    outfile<<setprecision(3)<<dptmin<<"\t"<<dptmax<<"\t"<<detamin<<"\t"<<detamax<<"\t"<<deff<<"\t"<<derr<<"\t"<<mceff<<"\t"<<mcerr<<"\t"<<absdiffCB<<"\t"<<absdiffExp<<"\t"<<absdiff_fitR<<"\t"<<absdiffLO<<"\t"<<absdiffTagS<<"\t"<<"-1"<<endl;

 //    outfile<<setprecision(3)<<detamin<<"\t"<<detamax<<"\t"<<dptmin<<"\t"<<dptmax<<"\t"<<deff<<"\t"<<derr<<"\t"<<mceff<<"\t"<<mcerr<<"\t"<<absdiffExp << "\t"<<absdiffCB<<"\t"<<absdiffLO << "\t" << absdiffTagS << "\t" << "-1" << "\t" << "-1" << "\t" << absdiff_fitR<<endl;

//    outfile<<setprecision(3)<<detamin<<"\t"<<detamax<<"\t"<<dptmin<<"\t"<<dptmax<<"\t"<<deff<<"\t"<<derr<<"\t"<<mceff<<"\t"<<mcerr<<"\t"<<absdiffExp << "\t"<<absdiffCB<<"\t"<< "-1" << "\t" << absdiffTagS << "\t" << "-1" << "\t" << "-1" << "\t" << absdiff_fitR<<endl;

//outfile<<setprecision(3)<<detamin<<"\t"<<detamax<<"\t"<<dptmin<<"\t"<<dptmax<<"\t"<<deff<<"\t"<<derr<<"\t"<<mceff<<"\t"<<mcerr<<"\t"<< deffExp << "\t"<<deffCB<<"\t"<< mceffLO << "\t" << deffTagS << "\t" << "-1" << "\t" << "-1" << "\t" << defffitR <<endl;

outfile<<setprecision(3)<<detamin<<"\t"<<detamax<<"\t"<<dptmin<<"\t"<<dptmax<<"\t"<<deff<<"\t"<<derr<<"\t"<<mceff<<"\t"<<mcerr << "\t" << deffExp << "\t"<<deffCB<<"\t"<< mceffLO << "\t" << deffTagS << "\t" << "-1" << "\t" << "-1" << "\t" << defffitR <<endl;



//cout << dptmin << "  " << dptmax << endl;

//   outfile<<setprecision(3)<<sfetamin<<"\t"<<sfetamax<<"\t"<<sfptmin<<"\t"<<sfptmax<<"\t"<<sf<<"\t"<<sferr<<endl;
 */

  }

}
