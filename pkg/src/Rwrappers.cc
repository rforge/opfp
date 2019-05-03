/* Dear emacs, please set -*- compile-command: "R CMD INSTALL .." -*- */

/* February 2014 Guillem Rigaill <rigaill@evry.inra.fr> 

   This file is part of the R package fpop

   fpop is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License (LGPL) as published by
   the Free Software Foundation; either version 2.1 of the License, or
   (at your option) any later version.
   
   opfp is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public License
   along with opfp; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <list>
#include "BinSeg_MultiDim.h"
#include "colibri.h"
#include <R_ext/Rdynload.h>

void colibri_op_R_c (double *profil, int *nbi, double *lambda_, double *mini, double *maxi, int *origine,
		     double *cout_n){
  colibri_op_c (profil, nbi, lambda_, mini, maxi, origine, cout_n);
}

void colibri_op_R_c_analysis (double *profil, int *nbi, double *lambda_, double *mini, double *maxi, int *origine,
			      double *cout_n, int *nbcandidate){
  colibri_op_c_analysis (profil, nbi, lambda_, mini, maxi, origine, cout_n, nbcandidate);
}

void BinSeg_interface
(double * dataVec_, int * Kmax_, int * n_, int * P_,
 int * Ruptures_, double * RupturesCost_){
  int n = *n_;
  int P = *P_;
  int Kmax = *Kmax_;
  double ** data = new double*[n];
  for(int i=0; i< n; i++)
    data[i] = new double[P];
  for(int i =0; i< n; i++)
    for(int j=0; j < P; j++)
      data[i][j] = dataVec_[i+ n*j];
  BinSeg_MultiDim BinSegRun(data, n, P, 2*Kmax + 10);
  BinSegRun.Initialize(Kmax);
  int i=0;
  for
    (std::list<int>::iterator iter = BinSegRun.Ruptures.begin();
     iter != BinSegRun.Ruptures.end(); iter++){
    Ruptures_[i] = *iter;
    i++;
  }
  i=0;
  for
    (std::list<double>::iterator iter = BinSegRun.RupturesCost.begin();
     iter != BinSegRun.RupturesCost.end(); iter++){
    RupturesCost_[i] = *iter;
    i++;
  }
  // delete
  for(int i = 0; i<n; i++)
    delete[] data[i];
  delete[] data;	
}


R_CMethodDef cMethods[] = {
  {"BinSeg_interface",(DL_FUNC) &BinSeg_interface, 6},
  {"colibri_op_R_c",(DL_FUNC) &colibri_op_R_c, 7},
  {"colibri_op_R_c_analysis",(DL_FUNC) &colibri_op_R_c_analysis, 8},
  {NULL, NULL, 0}
};

extern "C" {
  void R_init_fpop(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    //R_useDynamicSymbols call says the DLL is not to be searched for
    //entry points specified by character strings so .C etc calls will
    //only find registered symbols.
    R_useDynamicSymbols(info, FALSE);
  }
}
