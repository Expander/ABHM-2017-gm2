Install["bin/gm2calc.mx"];

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

Nsteps = 100;

(* data = {#, Sequence @@ Calc[#, 20, 0]}& /@ LogRange[173.34, 10^4, Nsteps]; *)
(* Export["MS_TB-20_Xt-0.dat", data]; *)

(* data = {#, Sequence @@ Calc[#, 20, Sqrt[6]]}& /@ LogRange[173.34, 10^4, Nsteps]; *)
(* Export["MS_TB-20_Xt-sqrt6.dat", data]; *)

data = {N[#], Sequence @@ CalcTB[#,False], Sequence @@ CalcTB[#,True]}& /@ LogRange[1, 10^3, Nsteps];
Export["TB.dat", data];
