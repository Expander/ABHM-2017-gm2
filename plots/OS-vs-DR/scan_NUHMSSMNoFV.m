Get["models/NUHMSSMNoFV/NUHMSSMNoFV_librarylink.m"];
Install["/home/avoigt/packages/FeynHiggs-2.13.0/build/MFeynHiggs"];

invalid;

LinearRange[start_, stop_, steps_] :=
    Table[start + i/steps (stop - start), {i, 0, steps}];

LogRange[start_, stop_, steps_] :=
    Module[{i, result = {}},
           For[i = 0, i <= steps, i++,
               result = AppendTo[result, Exp[Log[start] + (Log[stop] - Log[start]) i / steps]];
              ];
           result
          ];

RunNUHMSSMNoFV[TB_, MS_, MR_] :=
    Module[{handle, spectrum, obs},
           handle = FSNUHMSSMNoFVOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
                   maxIterations -> 0,                (* FlexibleSUSY[1] *)
                   solver -> 0,                       (* FlexibleSUSY[2] *)
                   calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 2,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 2,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> 0,                (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 2,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> 0          (* MODSEL[12] *)
               },
               fsSMParameters -> {
                   alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
                   GF -> 1.16637*^-5,      (* SMINPUTS[2] *)
                   alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
                   MZ -> 91.1876,          (* SMINPUTS[4] *)
                   mbmb -> 4.18,           (* SMINPUTS[5] *)
                   Mt -> 173.34,           (* SMINPUTS[6] *)
                   Mtau -> 1.777,          (* SMINPUTS[7] *)
                   Mv3 -> 0,               (* SMINPUTS[8] *)
                   MW -> 80.385,           (* SMINPUTS[9] *)
                   Me -> 0.000510998902,   (* SMINPUTS[11] *)
                   Mv1 -> 0,               (* SMINPUTS[12] *)
                   Mm -> 0.1056583715,     (* SMINPUTS[13] *)
                   Mv2 -> 0,               (* SMINPUTS[14] *)
                   md2GeV -> 0.00475,      (* SMINPUTS[21] *)
                   mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
                   ms2GeV -> 0.104,        (* SMINPUTS[23] *)
                   mcmc -> 1.27,           (* SMINPUTS[24] *)
                   CKMTheta12 -> 0,
                   CKMTheta13 -> 0,
                   CKMTheta23 -> 0,
                   CKMDelta -> 0,
                   PMNSTheta12 -> 0,
                   PMNSTheta13 -> 0,
                   PMNSTheta23 -> 0,
                   PMNSDelta -> 0,
                   PMNSAlpha1 -> 0,
                   PMNSAlpha2 -> 0,
                   alphaEm0 -> 1/137.035999074,
                   Mh -> 125.09
               },
               fsModelParameters -> {
                   TanBeta -> TB,
                   Qin -> 454.7,
                   M1 -> 140,
                   M2 -> MS,
                   M3 -> 300,
                   AtIN -> 0,
                   AbIN -> 0,
                   AtauIN -> 0,
                   AcIN -> 0,
                   AsIN -> 0,
                   AmuonIN -> 0,
                   AuIN -> 0,
                   AdIN -> 0,
                   AeIN -> 0,
                   MuIN -> -160,
                   mA2IN -> MS^2,
                   ml11IN -> MS,
                   ml22IN -> MS,
                   ml33IN -> MS,
                   me11IN -> MR,
                   me22IN -> MR,
                   me33IN -> MR,
                   mq11IN -> MS,
                   mq22IN -> MS,
                   mq33IN -> MS,
                   mu11IN -> MS,
                   mu22IN -> MS,
                   mu33IN -> MS,
                   md11IN -> MS,
                   md22IN -> MS,
                   md33IN -> MS
               }
           ];
           spectrum = FSNUHMSSMNoFVCalculateSpectrum[handle];
           observables = FSNUHMSSMNoFVCalculateObservables[handle];
           (* Print[spectrum]; *)
           FSNUHMSSMNoFVCloseHandle[handle];
           If[spectrum === $Failed || observables === $Failed,
              $Failed,
              Join[spectrum, observables, {FHgm2 -> CalcFHgm2[spectrum]}]
             ]
          ];

GetPar[spec_, par__] :=
    GetPar[spec, #]& /@ {par};

GetPar[spec_, par_] :=
    If[spec =!= $Failed, (par /. spec), invalid];

GetPar[spec_, par_[n__?IntegerQ]] :=
    If[spec =!= $Failed, (par /. spec)[[n]], invalid];

RunNUHMSSMNoFVMh[args__] :=
    Module[{spec = RunNUHMSSMNoFV[args]},
           If[spec === $Failed,
              invalid,
              GetPar[NUHMSSMNoFV /. spec, Pole[M[hh]][1]]
             ]
          ];

RunNUHMSSMNoFVAMU[args__] :=
    Module[{spec = RunNUHMSSMNoFV[args]},
           If[spec === $Failed,
              Array[invalid&, 5],
              {
                  GetPar[spec, FlexibleSUSYObservable`aMuon],
                  GetPar[spec, FlexibleSUSYObservable`aMuonUncertainty],
                  GetPar[spec, FlexibleSUSYObservable`aMuonGM2Calc],
                  GetPar[spec, FlexibleSUSYObservable`aMuonGM2CalcUncertainty],
                  GetPar[spec, FHgm2]
              }
             ]
          ];

CalcFHgm2[spec_] :=
    Module[{},
           FHSetFlags[4,0,0,2,0,2,3,1,1,0];
           FHSetSMPara[137.035999074, 127.916, 0.1184, 1.16637*^-5,
                       0.000510998902, 0.0024, 0.00475, 0.1056583715, 1.27, 0.104, 1.777, 4.18,
                       80.385, 91.1876, 0, 0,
                       0, 0, 0, 0];
           FHSetPara[1,
                     173.34,
                     GetPar[NUHMSSMNoFV /. spec, vu] / GetPar[NUHMSSMNoFV /. spec, vd],
                     GetPar[NUHMSSMNoFV /. spec, M[Ah][2]],
                     0 MHp,
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, ml2[3,3]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, me2[3,3]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, mq2[3,3]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, mu2[3,3]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, md2[3,3]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, ml2[2,2]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, me2[2,2]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, mq2[2,2]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, mu2[2,2]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, md2[2,2]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, ml2[1,1]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, me2[1,1]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, mq2[1,1]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, mu2[1,1]],
                     Sqrt @ GetPar[NUHMSSMNoFV /. spec, md2[1,1]],
                     GetPar[NUHMSSMNoFV /. spec, \[Mu]],
                     GetPar[NUHMSSMNoFV /. spec, T[Ye][3,3]] / GetPar[NUHMSSMNoFV /. spec, Ye[3,3]],
                     GetPar[NUHMSSMNoFV /. spec, T[Yu][3,3]] / GetPar[NUHMSSMNoFV /. spec, Yu[3,3]],
                     GetPar[NUHMSSMNoFV /. spec, T[Yd][3,3]] / GetPar[NUHMSSMNoFV /. spec, Yd[3,3]],
                     GetPar[NUHMSSMNoFV /. spec, T[Ye][2,2]] / GetPar[NUHMSSMNoFV /. spec, Ye[2,2]],
                     GetPar[NUHMSSMNoFV /. spec, T[Yu][2,2]] / GetPar[NUHMSSMNoFV /. spec, Yu[2,2]],
                     GetPar[NUHMSSMNoFV /. spec, T[Yd][2,2]] / GetPar[NUHMSSMNoFV /. spec, Yd[2,2]],
                     GetPar[NUHMSSMNoFV /. spec, T[Ye][1,1]] / GetPar[NUHMSSMNoFV /. spec, Ye[1,1]],
                     GetPar[NUHMSSMNoFV /. spec, T[Yu][1,1]] / GetPar[NUHMSSMNoFV /. spec, Yu[1,1]],
                     GetPar[NUHMSSMNoFV /. spec, T[Yd][1,1]] / GetPar[NUHMSSMNoFV /. spec, Yd[1,1]],
                     GetPar[NUHMSSMNoFV /. spec, MassB],
                     GetPar[NUHMSSMNoFV /. spec, MassWB],
                     GetPar[NUHMSSMNoFV /. spec, MassG],
                     GetPar[NUHMSSMNoFV /. spec, SCALE],
                     GetPar[NUHMSSMNoFV /. spec, SCALE],
                     GetPar[NUHMSSMNoFV /. spec, SCALE]];
           (* Print @ FHRetrieveSMPara[]; *)
           (* Print @ FHRetrievePara[]; *)
           (* Print @ FHGetSMPara[]; *)
           (* Print @ FHGetPara[]; *)
           gm2 /. FHConstraints[]
          ];

RunFH[TB_, MS_, MR_] :=
    Module[{},
           FHSetFlags[4,0,0,2,0,2,3,1,1,0];
           FHSetSMPara[137.035999074, 127.916, 0.1184, 1.16637*^-5,
                       0.000510998902, 0.0024, 0.00475, 0.1056583715, 1.27, 0.104, 1.777, 4.18,
                       80.385, 91.1876, 0, 0,
                       0, 0, 0, 0];
           FHSetPara[1,
                     173.34, TB, MS, 0 MHp,
                     MS, MR, MS, MS, MS,
                     MS, MR, MS, MS, MS,
                     MS, MR, MS, MS, MS,
                     -160,
                     0 Atau, 0 At, 0 Ab, 0 Amu, 0 Ac, 0 As, 0 Ae, 0 Au, 0 Ad,
                     140, MS, 300,
                     454.7, 454.7, 454.7];
           (* Print @ FHRetrieveSMPara[]; *)
           (* Print @ FHRetrievePara[]; *)
           (* Print @ FHGetSMPara[]; *)
           (* Print @ FHGetPara[]; *)
           gm2 /. FHConstraints[]
          ];

data = {N[#], Sequence @@ RunNUHMSSMNoFVAMU[50, #, 200], RunFH[50, #, 200]}& /@ LogRange[800,2000,60]

Export["scan_MS_OS-vs-DR_splitting.dat", data];

data = {N[#], Sequence @@ RunNUHMSSMNoFVAMU[50, #, #], RunFH[50, #, #]}& /@ LogRange[800,2000,60]

Export["scan_MS_OS-vs-DR.dat", data];

(* Off[General::stop]; *)

(* Print[10^10 {Sequence @@ RunNUHMSSMNoFVAMU[50, 2000, 2000], RunFH[50, 2000, 2000]}]; *)
