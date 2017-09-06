#include <Rcpp.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>

#include "msdef.h"
#include "msstringconversion.h"
#include "mshrmassintensmaputil.h"
#include "mslrmassintensmap.h"
#include "msconvert.h"
#include "mselementmultmaputil.h"
#include "mselementmultmapfilter.h"
#include "mselementintervalmaputil.h"
#include "msexception.h"
#include "msaddion.h"
#include "msgenformbymassalgo.h"
#include "msisotope.h"
#include "mscomparelr.h"
#include "mshrmatchformula.h"

#include "wrap.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List GenFormMatchIsotopeMsMs_R(
    Rcpp::DataFrame ms,
    Nullable<Rcpp::DataFrame> msms= R_NilValue,
    Rcpp::Nullable<double> m = R_NilValue,
    Rcpp::Nullable<Rcpp::List> settings = R_NilValue
    )
{

  Rcpp::List output;
  Rcpp::List _settings;
  if(settings.isNotNull())
  {
    _settings = Rcpp::List(settings);
  }

  string strEl,strEIM,strMs,strIon="M+H",strMsMs,strMsMv="nsae";
  bool3 b3EvenValSum=BOOL3_PERH;
  int iMinDiffMaxVal=numeric_limits<int>::min();
  int iMinDiffAtomCount=numeric_limits<int>::min();
  int iDbeExcess=3;
  int iFragMinDiffMaxVal=0;
  int iFragMinDiffAtomCount=0;
  bool bHCFilter=false,bAnalyze=false,bLoss=false;
  bool bAllowRadicalIons=false,bTex=false;
  bool bHighMass=true,bStripCalc=false;
  bool bPercent=false,bWriteCalcMass=false,bWriteDBE=false;
  bool bObligatoryHeteroAtom=false,bMultipleValencies=false;
  bool bOrganicCompound=false;
  bool bIntens=false, bShowRef=true;
  double dMass=0.0,dPPM=5.0;
  double dMinMsMv=-numeric_limits<double>::max();
  double dMinMsMsMv=-numeric_limits<double>::max();
  double dMinCombMv=-numeric_limits<double>::max();
  int iPpmAccept=2,iPpmReject=4,iExponent=5;
  string strWeightMass="false",strWeightIntens="false";
  //ostream *pOst=NULL,*pOutMs=NULL,*pOutMsMs=NULL,*pOutCleanMsMs=NULL;

  enum SortMethod { SortPpm, SortMsMv, SortMsMsMv, SortCombMv, SortUndefined };
  const char* pSortValue[SortUndefined]={"ppm","msmv","msmsmv","combmv"};
  unsigned int iSortMethod=SortUndefined;


  // // TODO: convert output
  // if((it=mapArgValue.find("out"))!=mapArgValue.end())
  // { pOst=it->second.empty()?&cout:new ofstream(it->second.c_str()); mapArgValue.erase(it); }
  //
  // if((it=mapArgValue.find("oms"))!=mapArgValue.end())
  // { pOutMs=it->second.empty()?&cout:new ofstream(it->second.c_str()); mapArgValue.erase(it); }
  //
  // if((it=mapArgValue.find("omsms"))!=mapArgValue.end())
  // { pOutMsMs=it->second.empty()?&cout:new ofstream(it->second.c_str()); mapArgValue.erase(it); }
  //
  // if((it=mapArgValue.find("oclean"))!=mapArgValue.end())
  // { pOutCleanMsMs=it->second.empty()?&cout:new ofstream(it->second.c_str()); mapArgValue.erase(it); }

  if(m.isNotNull())
    dMass = NumericVector(m)[0];

  if(_settings.containsElementNamed("ppm"))
    dPPM = _settings["ppm"];

  if(_settings.containsElementNamed("acc"))
    iPpmAccept = _settings["acc"];

  if(_settings.containsElementNamed("rej"))
    iPpmReject = _settings["rej"];

  if(_settings.containsElementNamed("thms"))
    dMinMsMv = _settings["thms"];
  if(_settings.containsElementNamed("thmsms"))
    dMinMsMsMv = _settings["thmsms"];
  if(_settings.containsElementNamed("thcomb"))
    dMinCombMv = _settings["thcomb"];

  if(_settings.containsElementNamed("el"))
    strEl = as<string>(_settings["el"]);

  // this is slightly different from original: OC can now be True of False, not just Set
  if(_settings.containsElementNamed("oc"))
    bOrganicCompound = _settings["oc"];

  if(_settings.containsElementNamed("ff"))
    strEIM = as<string>(_settings["ff"]);

  if(_settings.containsElementNamed("het"))
    bObligatoryHeteroAtom = _settings["het"];

  // vsp is some kind of tristate bool that is PERH if not set.
  if(_settings.containsElementNamed("vsp"))
  {
    if(LogicalVector::is_na(_settings["vsp"]))
    {
      b3EvenValSum = BOOL3_PERH;
    }
    else
    {
      bool vsp = _settings["vsp"];
      b3EvenValSum = vsp?BOOL3_TRUE:BOOL3_FALSE;
    }
  }

  if(_settings.containsElementNamed("vsm2mv"))
    iMinDiffMaxVal = _settings["vsm2mv"];

  if(_settings.containsElementNamed("vsm2ap2"))
    iMinDiffAtomCount = _settings["vsm2ap2"];

  if(_settings.containsElementNamed("exist"))
  {
    b3EvenValSum = BOOL3_TRUE;
    iMinDiffMaxVal=iMinDiffAtomCount=0;
    string exist = _settings["exist"];
    bMultipleValencies = (exist == "mv");
  }

  if(_settings.containsElementNamed("hcf"))
    bHCFilter = _settings["hcf"];
  if(_settings.containsElementNamed("oei"))
    bAllowRadicalIons = _settings["oei"];

  if(_settings.containsElementNamed("analyze"))
    bAnalyze = _settings["analyze"];
  if(_settings.containsElementNamed("loss"))
    bLoss = _settings["loss"];
  if(_settings.containsElementNamed("sc"))
    bStripCalc = _settings["sc"];
  if(_settings.containsElementNamed("wi"))
    strWeightIntens = as<string>(_settings["wi"]);
  if(_settings.containsElementNamed("wm"))
    strWeightMass = as<string>(_settings["wm"]);

  if(_settings.containsElementNamed("ion"))
    strIon = as<string>(_settings["ion"]);
  if(_settings.containsElementNamed("exp"))
    iExponent = _settings["exp"];

  if(_settings.containsElementNamed("dbeexc"))
    iDbeExcess = _settings["dbeexc"];

  if(_settings.containsElementNamed("ivsm2mv"))
    iFragMinDiffMaxVal = _settings["ivsm2mv"];
  if(_settings.containsElementNamed("ivsm2ap2"))
    iFragMinDiffAtomCount = _settings["ivsm2ap2"];

  if(_settings.containsElementNamed("msmv"))
    strMsMv = as<string>(_settings["msmv"]);


  if(_settings.containsElementNamed("tex"))
    bTex = _settings["tex"];

  if(_settings.containsElementNamed("pc"))
    bPercent = _settings["pc"];

  if(_settings.containsElementNamed("dbe"))
    bWriteDBE = _settings["dbe"];


  if(_settings.containsElementNamed("cm"))
    bWriteCalcMass = _settings["cm"];

  if(_settings.containsElementNamed("intens"))
    bIntens = _settings["intens"];

  // this is not needed here, it's for displaying the help file
  // // "ref" is originally "noref" andwould shut OFF showing ref.
  // if(_settings.containsNamedElement("ref"))
  //   bShowRef = as(_settings["ref"]);



  if(_settings.containsElementNamed("msmv"))
    strMsMv = as<string>(_settings["msmv"]);
  // if((it=mapArgValue.find("sort"))!=mapArgValue.end())
  // {
  if(_settings.containsElementNamed("sort"))
  {
    string sort = as<string>(_settings["sort"]);
    if(sort == "")
      iSortMethod = SortCombMv;
    else
      for(iSortMethod=0;iSortMethod!=SortUndefined;iSortMethod++)
        if(sort==pSortValue[iSortMethod]) break;

        if(iSortMethod==SortUndefined)
        {
          stop("invalid value for key 'sort'");
        }
  }


  // TODO: this is not a bad idea!
  // Remove assigned arguments from the list.
  // if(!mapArgValue.empty())
  // {
  //   cerr << "Error:\tunknown key(s)\n\t ";
  //   for(it=mapArgValue.begin();it!=mapArgValue.end();it++)
  //     cerr << it->first;
  //   cerr << "\n"; GenFormMatchIsotopeMsMsUsage(strProgName); exit(1);
  //
  // }


  LrMassIntensMap lrmCalc,lrmExp;
  HrMassIntensMap hrmMs,hrmMsMs,hrmMsMsNorm,hrmMsMsClean;

  // ifstream in(strMs.c_str());
  // if(!in)
  // { cerr << "Error opening " << strMs << "\n"; exit(1); }
  // while(in)
  // {
  //   double dMass,dIntens;
  //   in >> dMass >> dIntens;
  //
  //   if(in)
  //     hrmMs[dMass]=dIntens;
  // }

  hrmMs = toMap(ms);

  Init(lrmExp,hrmMs);

  //Normalize to intensity sum = 1;
  //this is important for the ShimadzuSimilarityIndex
  lrmExp.Normalize(1.0,false);
  hrmMs.Normalize(1.0,false);


  if(dMass==0.0)
  {
    if(hrmMs.GetBasePeak()!=hrmMs.end())
      dMass=hrmMs.GetBasePeak()->first;
    else
    { stop("Error searching ms base peak"); }
  }

  if(msms.isNotNull())
  {
    Rcpp::DataFrame msmsTemp = Rcpp::DataFrame(msms);
    hrmMsMs = toMap(msmsTemp);
  }

  hrmMsMs.Normalize();
  hrmMsMsNorm=hrmMsMs;
  WeightSpectrum(hrmMsMs,strWeightMass,strWeightIntens,iExponent);
  hrmMsMs.Normalize();
  //if(pOutMsMs) *pOutMsMs << hrmMsMs;

  AddIon AI;
  try{AI=AddIon(strIon);}
  catch(MsError& ME)
  {
    stop(ME.what());
  }
  dMass=AI.CalcMolMass(dMass);

  double dMinMass=dMass*(1.0-dPPM/1000000.0),dMaxMass=dMass*(1.0+dPPM/1000000.0);
  unsigned long iTotalCount=0,iValidCount=0,iHCFCount=0,iIsotopeCount=0,iMsMsCount=0,iCombCount=0;
  ElementIntervalMap EIM;

  if(!strEIM.empty())
  {
    try{Init(EIM,strEIM);}
    catch(MsError& ME)
    {
      stop(ME.what());
    }
  }
  else
  {
    Interval I(0,numeric_limits<ElementMult>::max());
    if(strEl.empty()) strEl="CHBrClFINOPSSi";

    try{Init(EIM,strEl);}
    catch(MsError& ME)
    {
      stop(ME.what());
    }

    EIM.SetInterval(I);

    if(bOrganicCompound)
      EIM.SetInterval(e_C,Interval(1,numeric_limits<ElementMult>::max()));
  }

  // raise lower limits
  if(!EIM.RaiseLowerLimits(AI.GetFormulaMinus()))
  {
    stop("Invalid  element bounds when considering ionization type"); //" " << EIM);
  }

  // in order to write explained peaks and explaining formula to output
  map<double,list<ElementMultMap> > mapMassEMM;

  // in order to write candidate formulas sorted by MV
  multimap<double,Rcpp::List> mapSort;

  //	CpuTime	time; time.Start();

  GenFormByMassAlgo<double> BFBMA(EIM,dMinMass,dMaxMass);

  // MST: here the output variables are defined
  StringVector fo_name;
  NumericVector fo_DBE;
  NumericVector fo_CalcMass;
  NumericVector fo_ppm;
  NumericVector fo_msmv;
  NumericVector fo_msmsmv;
  NumericVector fo_cmv;
  // They should later go in a data frame
  List fo_msms;
  // this list is a list of dataframes, where all the analyzed msms go

  for(BFBMA.begin();!BFBMA.end();BFBMA.operator++())
  {
    ElementMultMap emmIon,emmMol;

    if(bObligatoryHeteroAtom)
    {
      BFBMA.GetElementMultMap(emmMol);
      if(!HasHeteroAtom(emmMol))
        continue;
    }

    iTotalCount++;

    if(emmMol.empty()&&bMultipleValencies)
      BFBMA.GetElementMultMap(emmMol);

    bool bValid=bMultipleValencies?
    IsGraphicalMultVal(emmMol):
      BFBMA.IsValid(b3EvenValSum,iMinDiffMaxVal,iMinDiffAtomCount);

    if(bValid)
    {
      iValidCount++;

      if(emmMol.empty())
        BFBMA.GetElementMultMap(emmMol);

      if(bHCFilter)
      { if(HeuerdingClercFilter(emmMol)==false) continue; }
      iHCFCount++;

      AI.CalcEMM(emmMol,emmIon);
      Init(lrmCalc,emmIon);

      if(bHighMass)
        Init(lrmExp,hrmMs,GetNominalMass<int>(emmIon)-GetNominalMass<double>(emmIon)+0.5);

      //cout << setprecision(12) << GetNominalMass<double>(emmIon)-(c_dMassElectron*iCharge) << "\n";

      if(bStripCalc)
      {
        //remove Peaks in calculated spectrum that are not present in Exp spectrum
        lrmCalc.erase(lrmCalc.begin(),lrmCalc.lower_bound(lrmExp.begin()->first));
        lrmCalc.erase(lrmCalc.upper_bound(lrmExp.rbegin()->first),lrmCalc.end());
      }

      //cout << lrmCalc << '\n' << lrmExp;

      double dMs=1.0,dMsMs=0.0;

      if(strMsMv=="ndp")
        dMs-=CompareNormalizedDotProduct(lrmExp,lrmCalc);
      else
        if(strMsMv=="nsse")
          dMs-=CompareNormalizedSumSquares(lrmExp,lrmCalc);
        else
        {
          //LrMassIntensMap hrmMsLr;

          //Init(hrmMsLr,ShiftMass(hrmMs,HrMassIntensMap(),GetNominalMass<int>(emmIon)-GetNominalMass<double>(emmIon)));
          dMs-=CompareShimadzuSimilarityIndex(lrmExp,lrmCalc);

          //	dMs-=CompareShimadzuSimilarityIndexSpecial(lrmExp,lrmCalc);
        }

        if(dMs<dMinMsMv) continue;
        iIsotopeCount++;

        if(hrmMsMs.size())
          dMsMs=ComputeHrMatchvalue(hrmMsMs,emmIon,iPpmAccept,iPpmReject,
                                    iDbeExcess,iFragMinDiffMaxVal,iFragMinDiffAtomCount,
                                    bAllowRadicalIons,AI.GetCharge()>0,NULL,(bAnalyze)?&mapMassEMM:NULL);

        if(dMsMs<dMinMsMsMv) continue;
        iMsMsCount++;

        double dMvComb=dMsMs*dMs;

        if(dMsMs*dMs<dMinCombMv) continue;
        iCombCount++;

        if(bPercent)
        {
          dMs*=100.0;dMsMs*=100.0;dMvComb*=100.0;
        }

        if(true)
        {
          //	Use this code, if all three MS matchvalues shall be printed
          //				pOst->setf(ios::fixed,ios::floatfield);
          //				double dMassCalc=GetNominalMass<double>(emmMol);
          //				*pOst << setw(15) << setiosflags(ios_base::left) << (bTex?TexName(emmMol):ChemName(emmMol)) << (bTex?"&\t":"\t")
          //					  << setiosflags(ios_base::right) << setw(6) << setprecision(1) << CalcDiffPPM(dMassCalc,dMass)
          //					  << (bTex?"&\t":"\t") << setprecision(6) << 1.0-CompareNormalizedDotProduct(lrmExp,lrmCalc)
          //					  << (bTex?"&\t":"\t") << setprecision(6) << 1.0-CompareShimadzuSimilarityIndex(lrmExp,lrmCalc)
          //					  << (bTex?"&\t":"\t") << setprecision(6) << 1.0-CompareNormalizedSumSquares(lrmExp,lrmCalc) << "\n" << resetiosflags(ios::adjustfield);

          int iWritePrec=bPercent?3:5;



          double dMassCalc=GetNominalMass<double>(emmMol);
          double dPpmCalc=CalcDiffPPM(dMassCalc,dMass);

          fo_name.push_back(bTex?TexName(emmMol):ChemName(emmMol));
          //sStream << setw(15) << setiosflags(ios_base::left) << (bTex?TexName(emmMol):ChemName(emmMol)) << (bTex?"&\t":"\t");
          // if(bWriteDBE)
            fo_DBE.push_back(GetDBE(emmMol));
            //sStream << setiosflags(ios_base::right) << setw(5) << setprecision(1) << GetDBE(emmMol) << (bTex?"&\t":"\t");
          // if(bWriteCalcMass)
            fo_CalcMass.push_back(AI.CalcIonMass(dMassCalc));
          //  sStream << setiosflags(ios_base::right) << setw(10) << setprecision(5) << AI.CalcIonMass(dMassCalc) << (bTex?"&\t":"\t");
          fo_ppm.push_back(dPpmCalc);
          fo_msmv.push_back(dMs);
          fo_msmsmv.push_back(dMsMs);
          fo_cmv.push_back(dMvComb);

          // sStream << setiosflags(ios_base::right) << setw(6) << setprecision(1) << dPpmCalc
          //         << (bTex?"&\t":"\t") << setiosflags(ios_base::right) << setw(7) << setprecision(iWritePrec) << dMs;
          // if(hrmMsMs.size())
          //   sStream << (bTex?"&\t":"\t") << setiosflags(ios_base::right) << setw(7) << setprecision(iWritePrec) << dMsMs
          //           << (bTex?"&\t":"\t") << setiosflags(ios_base::right) << setw(7) << setprecision(iWritePrec) << dMvComb;
          //   sStream << (bTex?"\\\\ \\hline\n":"\n") << resetiosflags(ios::adjustfield);
          //   *pOst << sStream.str();

          //TODO:: handle sorting, probably on the R side
            // if(iSortMethod!=SortUndefined)
            //   switch(iSortMethod)
            //   {
            //   case SortPpm   : mapSort.insert(make_pair(-fabs(dPpmCalc),formulaOut)); break;
            //   case SortMsMv  : mapSort.insert(make_pair(dMs,formulaOut)); break;
            //   case SortMsMsMv: mapSort.insert(make_pair(dMsMs,formulaOut)); break;
            //   case SortCombMv: mapSort.insert(make_pair(dMvComb,formulaOut)); break;
            //   }

            if(bAnalyze)
            {

              NumericVector mz;
              NumericVector intens;
              NumericVector dbe;
              StringVector loss;
              StringVector formula;
              NumericVector ppm;
              NumericVector calcMass;

              // This iterates through the masses in the spectrum
              for(map<double,list<ElementMultMap> >::const_iterator itMass=mapMassEMM.begin();itMass!=mapMassEMM.end();itMass++)
              {
                // ie this is run once per mass

                const list<ElementMultMap>& listEMM=itMass->second;
                double _mz = itMass->first;
                double _intens = hrmMsMsNorm.GetIntens(itMass->first);

               // These two lines kinda returned the "title" for the spectrum, i.e. mz / intensity.

                //*pOst << setw(10) << setiosflags(ios::fixed) << setprecision(5) << itMass->first << '\t';
                //if(bIntens)
                //  *pOst << setw(10) << setiosflags(ios::fixed) << setprecision(8) << hrmMsMsNorm.GetIntens(itMass->first) << '\t';

                // for now just:
                // TODO: CONTINUE HERE

                for(list<ElementMultMap>::const_iterator itEMM=listEMM.begin();itEMM!=listEMM.end();itEMM++)
                {
                  // this iterates through the formulas matching to the mass
                  double dIonMassCalc=GetNominalMass<double>(*itEMM)-(c_dMassElectron*AI.GetCharge());
                  ElementMultMap emmLoss;
                  string strLoss;
                  if(bLoss)
                  {
                    emmLoss=emmIon-*itEMM;
                    strLoss=emmLoss.empty()?"":"-"+ChemName(emmLoss);
                  }
                  mz.push_back(_mz);
                  intens.push_back(_intens);
                  loss.push_back(strLoss);
                  formula.push_back(ChemName(*itEMM));
                  dbe.push_back(CalcDBE(*itEMM));
                  calcMass.push_back(dIonMassCalc);
                  ppm.push_back(CalcDiffPPM(dIonMassCalc, _mz));

//
//                   if(itEMM!=listEMM.begin())
//                     *pOst << "          \t" << (bIntens?"          \t":"");
//                   *pOst << setw(15) << setiosflags(ios_base::left) << (bLoss?strLoss:ChemName(*itEMM))  << '\t';
//                   if(bWriteDBE)
//                     *pOst << setiosflags(ios_base::right) << setw(5) << setprecision(1) << CalcDBE(*itEMM) << (bTex?"&\t":"\t");
                  // if(bWriteCalcMass)
                  //   *pOst << setw(10) << setiosflags(ios::fixed) << setprecision(5) << dIonMassCalc << '\t';
                  // *pOst << setw(6) << setprecision(1) << CalcDiffPPM(dIonMassCalc,itMass->first) << '\n' << resetiosflags(ios::adjustfield);
                }
              }
              DataFrame df = DataFrame::create(
                Named("mz") =  mz,
                Named("int") = intens,
                Named("dbe") = dbe,
                Named("loss") = loss,
                Named("formula") = formula,
                Named("dppm") = ppm,
                Named("calcMass") = calcMass);
              fo_msms.push_back(df);
            }

            // TODO: see if I want to support this and what this is good for
              // if(pOutCleanMsMs)
              //   for(map<double,list<ElementMultMap> >::const_iterator itMass=mapMassEMM.begin();itMass!=mapMassEMM.end();itMass++)
              //     hrmMsMsClean.AddPeak(itMass->first,hrmMsMsNorm.GetIntens(itMass->first));
        }
    }
  }

  // TODO: see above
  // if(pOutCleanMsMs) *pOutCleanMsMs << hrmMsMsClean;

  // TODO: is this needed? (@Emma)
  // if(bHCFilter)
  //   cout << iCombCount << "/" << iMsMsCount << "/" << iIsotopeCount << "/" << iHCFCount << "/" << iValidCount << "/" << iTotalCount
  //        << " (final/MSMS/MS-/HC-filter/valid/total) formula(s)\n";// in " << setprecision(1) << time.GetSecs() << "s\n";
  //   else
  //     cout << iCombCount << "/" << iMsMsCount << "/" << iIsotopeCount << "/" <<  iValidCount << "/" << iTotalCount
  //          << " (final/MSMS/MS-filter/valid/total) formula(s)\n";// in " << setiosflags(ios::fixed) << setprecision(1) << time.GetSecs() << "s\n";

      // if(iSortMethod!=SortUndefined)
      //   for(multimap<double,string>::reverse_iterator itSort=mapSort.rbegin();itSort!=mapSort.rend();itSort++)
      //     *pOst << itSort->second;

      // if(pOst!=NULL&&pOst!=&cout) delete pOst;
      // if(pOutMs!=NULL&&pOutMs!=&cout) delete pOutMs;
      // if(pOutMsMs!=NULL&&pOutMsMs!=&cout) delete pOutMsMs;
      // if(pOutCleanMsMs!=NULL&&pOutCleanMsMs!=&cout) delete pOutCleanMsMs;

      DataFrame formulas = DataFrame::create(
        Named("formula") = fo_name,
        Named("dbe")= fo_DBE,
        Named("calcMass") = fo_CalcMass,
        Named("ppm") =  fo_ppm,
        Named("msmv")=  fo_msmv,
        Named("msmsmv") = fo_msmsmv,
        Named("cmv") =  fo_cmv
      );

      Rcpp::List results = List::create(
          Named("formulas") = formulas,
          Named("msms") = fo_msms
      );

      return(results);

}
