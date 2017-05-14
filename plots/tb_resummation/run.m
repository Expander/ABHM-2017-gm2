Install["bin/gm2calc.mx"];
Install["/home/avoigt/packages/FeynHiggs-2.13.0/build/MFeynHiggs"];

LinearRange[start_, stop_, steps_] :=
    Table[start + i/steps (stop - start), {i, 0, steps}];

LogRange[start_, stop_, steps_] :=
    Module[{i, result = {}},
           For[i = 0, i <= steps, i++,
               result = AppendTo[result, Exp[Log[start] + (Log[stop] - Log[start]) i / steps]];
              ];
           result
          ];

GM2CalcSetFlags[
    loopOrder -> 2,
    tanBetaResummation -> True,
    forceOutput -> False];

GM2CalcSetSMParameters[
    alphaMZ -> 0.00775531,  (* 1L *)
    alpha0 -> 0.00729735,   (* 2L *)
    alphaS -> 0.1184,       (* 2L *)
    MW -> 80.385,           (* 1L *)
    MZ -> 91.1876,          (* 1L *)
    MT -> 173.34,           (* 2L *)
    mbmb -> 4.18,           (* 2L *)
    ML -> 1.777,            (* 2L *)
    MM -> 0.1056583715];    (* 1L *)

Calc[MS_, TanB_, Xt_] :=
    Module[{},
           {amu, Damu} /. GM2CalcAmuGM2CalcScheme[
               MAh    -> MS,     (* 2L *)
               TB     -> TanB,   (* 1L *)
               Mu     -> MS,     (* 1L *)
               MassB  -> MS,     (* 1L *)
               MassWB -> MS,     (* 1L *)
               MassG  -> MS,     (* 2L *)
               mq2    -> MS^2 IdentityMatrix[3], (* 2L *)
               ml2    -> MS^2 IdentityMatrix[3], (* 1L *)
               mu2    -> MS^2 IdentityMatrix[3], (* 2L *)
               md2    -> MS^2 IdentityMatrix[3], (* 2L *)
               me2    -> MS^2 IdentityMatrix[3], (* 2L *)
               Au      -> {{MS/TanB, 0, 0         },
                           {0, MS/TanB, 0         },
                           {0, 0, MS/TanB + Xt MS}}, (* 2L *)
               Ad     -> MS TanB IdentityMatrix[3],  (* 2L *)
               Ae     -> MS TanB IdentityMatrix[3],  (* 1L *)
               Q      -> MS]                         (* 2L *)
          ];

CalcTB[TanB_, resum_] :=
    Module[{},
           GM2CalcSetFlags[
               loopOrder -> 2,
               tanBetaResummation -> resum,
               forceOutput -> False];
           {amu, Damu} /. GM2CalcAmuGM2CalcScheme[
               MAh    -> 1500,
               TB     -> TanB,
               Mu     -> 30000,
               MassB  -> 1000,
               MassWB -> -30000,
               MassG  -> 2000,
               mq2    -> 3000^2 IdentityMatrix[3],
               ml2    -> {{1000^2,0,0}, {0,1000^2,0}, {0,0,3000^2}},
               mu2    -> 3000^2 IdentityMatrix[3],
               md2    -> 3000^2 IdentityMatrix[3],
               me2    -> {{1000^2,0,0}, {0,1000^2,0}, {0,0,3000^2}},
               Au     -> 0 IdentityMatrix[3],
               Ad     -> 0 IdentityMatrix[3],
               Ae     -> 0 IdentityMatrix[3],
               Q      -> 866.360379]
          ];

RunFH[TB_] :=
    Module[{res},
           FHSetFlags[4,0,0,2,0,2,3,1,1,0];
           FHSetSMPara[137.035999074, 127.916, 0.1184, 1.16637*^-5,
                       0.000510998902, 0.0024, 0.00475, 0.1056583715, 1.27, 0.104, 1.777, 4.18,
                       80.385, 91.1876, 0, 0,
                       0, 0, 0, 0];
           FHSetPara[1,
                     173.34, TB, 1500, 0 MHp,
                     1000, 1000, 3000, 3000, 3000,
                     1000, 1000, 3000, 3000, 3000,
                     1000, 1000, 3000, 3000, 3000,
                     30000,
                     0 Atau, 0 At, 0 Ab, 0 Amu, 0 Ac, 0 As, 0 Ae, 0 Au, 0 Ad,
                     1000, -30000, 2000,
                     866.360379, 866.360379, 866.360379];
           (* Print @ FHRetrieveSMPara[]; *)
           (* Print @ FHRetrievePara[]; *)
           (* Print @ FHGetSMPara[]; *)
           (* Print @ FHGetPara[]; *)
           res = FHConstraints[];
           If[res === Indeterminate || Head[res] === FHError, Indeterminate, gm2 /. res]
          ];

Nsteps = 100;

(* data = {#, Sequence @@ Calc[#, 20, 0]}& /@ LogRange[173.34, 10^4, Nsteps]; *)
(* Export["MS_TB-20_Xt-0.dat", data]; *)

(* data = {#, Sequence @@ Calc[#, 20, Sqrt[6]]}& /@ LogRange[173.34, 10^4, Nsteps]; *)
(* Export["MS_TB-20_Xt-sqrt6.dat", data]; *)

data = {N[#], Sequence @@ CalcTB[#,False], Sequence @@ CalcTB[#,True], RunFH[#]}& /@ LogRange[1, 10^3, Nsteps];
Export["TB.dat", data];
