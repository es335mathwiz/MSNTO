$Path=PrependTo[$Path,"../../ProtectedSymbols"]


BeginPackage["MSNTO`",{"ProtectedSymbols`","JLink`"}]
$gitDir=
  If[$OperatingSystem==="MacOSX",
     "~/git/",
     "/msu/scratch2/m1gsa00/git/"]
  Print["em:",$gitDir];  
pickBestMSNTO::usage="pickBestMSNTO[{maxWidth_?NumberQ,thePts_?MatrixQ,clustPts:{{_?MatrixQ,{}}..},abl_?aBestListQ,rb:{{_?aRangeSpecQ,_?aBestQ}..}}]"
MSNTOMinimizer::usage="MSNTOMinimizer[{maxWidth_?NumberQ,thePts_?MatrixQ,clustPts:{{_?MatrixQ,{}}..},abl_?aBestListQ,rb:{{_?aRangeSpecQ,_?aBestQ}..}}]"
MSNTO::usage="MSNTO[theFunc_,lowerUpper:{{_?NumberQ,_?NumberQ}..},numPts_Integer,netGen_Symbol,beta_?NumberQ,rhoValMult_?NumberQ,cntrctnFactor_?NumberQ]"
iterMSNTOFunc::usage="iterMSNTOFunc[theFunc_,numPts_Integer,netGen_Symbol,beta_?NumberQ,rhoValMult_?NumberQ,cntrctnFactor_?NumberQ]"


Nied::usage="Nied[x_, y_]"
Sobol::usage="Sobol[x_, y_]"
SuiteTore::usage="SuiteTore[x_, y_, z_]"
RandSuiteTore::usage="RandSuiteTore[x_, y_]"
Begin["`Private`"]

MSNTOMinimizer[{maxWidth_?NumberQ,thePts_?MatrixQ,clustPts:{{_?MatrixQ,{}}..},abl_?aBestListQ,rb:{{_?aRangeSpecQ,_?aBestQ}..}}]:=
With[{choices=Last/@abl},
With[{ord=Ordering[choices]},abl[[ord[[1]],1]]]]



pickBestMSNTO[{maxWidth_?NumberQ,thePts_?MatrixQ,clustPts:{{_?MatrixQ,{}}..},
abl_?aBestListQ,rb:{{_?aRangeSpecQ,_?aBestQ}..}}]:=
With[{choices=Last/@abl},Min[choices]]


MSNTO[theFunc_,listBounds:{{
{{{_?NumberQ,_?NumberQ}..},discreteVals:{{_Integer..}...}},{_?VectorQ,(_?NumberQ|Infinity)}}..},numPts_Integer,
netGen_Symbol,beta_?NumberQ,rhoValMult_?NumberQ,cntrctnFactor_?NumberQ]:=
Module[{theMeld=meldRes[Map[MSNTO[theFunc,#,numPts,netGen,beta,rhoValMult,cntrctnFactor]&,
listBounds]]},Print["MSNTO:",listBounds];
theMeld
]
sortByClusterBest[
clstrBestVals:{{_?VectorQ,(_?NumberQ|Infinity)}..},
aBndry:{{{{{_?NumberQ,_?NumberQ}..},{}},
{_?VectorQ,(_?NumberQ|Infinity)}}..}]:=
Module[{$maxClustersToKeep=Max[10,Length[Kernels[]]]},
With[{theOrds=Ordering[clstrBestVals,Min[Length[aBndry],$maxClustersToKeep],
(#1[[-1]]<#2[[-1]])&]},
aBndry[[theOrds]]]]

 
 tackOnDiscrete[nonDPts:{_?VectorQ..},{}]:=nonDPts


tackOnDiscrete[nonDPts:{_?VectorQ..},DPts:{_?VectorQ..}]:=
Flatten[tackOnOneNonD[#,DPts]&/@nonDPts,1]

tackOnOneNonD[aNonDPt_?VectorQ,DPts:{_?VectorQ..}]:=
Join[aNonDPt,#]&/@DPts                    

prepDiscrete[{}]:={}
prepDiscrete[discreteVals:{{_Integer..}..}]:=
    With[{pre=Outer[List,Sequence @@discreteVals]},Flatten[pre,Depth[pre]-3]]


                    
MSNTO[myFunc_,
{aCluster:{lowerUpper:{{_?NumberQ,_?NumberQ}..},discreteVals:{{_Integer..}...}},previousBest:{_?VectorQ,(_?NumberQ|Infinity)}},
numPts_Integer,netGen_Symbol,beta_?NumberQ,rhoValMult_?NumberQ,cntrctnFactor_?NumberQ]:=
Module[{},Print["enter MSNTO"];
SetSharedVariable[ii];ii=0;
With[{vecLen=Length[lowerUpper],maxWidth=chkBounds[lowerUpper]},
With[{thePts=tackOnDiscrete[unitToActual[lowerUpper,netGen[numPts,vecLen]],
prepDiscrete[discreteVals]]},kernel0Print["MSNTO about to parallelmap"];Print["thePts:",{numPts,vecLen,thePts}];
Print["about to test javaInitDone=",javaInitDone];
If[javaInitDone=!=True,kernel0Print["doJavaInit[]"];doJavaInit[],kernel0Print["don't,doJavaInit[]"],kernel0Print["confused don't doJavaInit[]"]];
DistributeDefinitions[myFunc,actualErrs];
theFVals=(kernel0Print["what"];ParallelMap[
(myFunc @@ #)&, thePts]);Print["theFVals=",theFVals//InputForm];
UninstallJava[];javaInitDone=False;
With[{ords=Ordering[theFVals],theDisp=dispersion[thePts]},
Print["MSNTO done parallelmap"];
With[{toRetain=Ceiling[beta*numPts]},
With[{ordPts=thePts[[ords]],ordFVals=theFVals[[ords]]},
With[{keepPts=ordPts[[Range[toRetain]]],
keepFVals=ordFVals[[Range[toRetain]]]},
With[{clusters=getAllClusters[rhoValMult*theDisp,keepPts,aCluster]},
     Print["theclustersnow:",{rhoValMult,theDisp,keepPts,aCluster,clusters}];
With[{bestVals=ParallelMap[thisClusterBest[#,keepPts,keepFVals,previousBest]&,
clusters]},
With[{newBounds=
MapThread[{
contractBounds[#1[[1]],aCluster,cntrctnFactor,#2],#}&,{bestVals,clusters}]},
{maxWidth,thePts,clusters,bestVals,
sortByClusterBest[bestVals,newBounds]}]]]]]]]]]]


thisClusterBest[aCluster:{{_?VectorQ..},_},keepPts_?MatrixQ,keepFVals_?VectorQ,
previousBest:{_?VectorQ,previousFVal:(_?NumberQ|Infinity)}]:=
With[{theLoc=Min[Flatten[Position[keepPts,#]&/@aCluster[[1]]][[1]]]},
If[previousFVal<=keepFVals[[theLoc]],previousBest,
{keepPts[[1]],keepFVals[[1]]},previousBest]]




bestResNest[{__,{_,_,_,bVals:{{_?VectorQ,_?NumberQ}..},_}}]:=
With[{theOrds=Ordering[Last/@bVals,1]},{theOrds,bVals,bVals[[theOrds,1]][[1]]}]

iterMSNTOFunc[myFunc_,numPts_Integer,
netGen_Symbol,beta_?NumberQ,rhoValMult_?NumberQ,cntrctnFactor_?NumberQ]:=
Module[{},
Function @@ {{listBounds},
MSNTO[myFunc,listBounds,numPts,netGen,beta,rhoValMult,cntrctnFactor]}
]

chkBounds[lowerUpper:{{_?NumberQ,_?NumberQ}..}]:=
Max[Abs[(#[[2]]-#[[1]])/2&/@lowerUpper]]


aRangeSpecQ[xx_]:=MatchQ[xx,{{_?VectorQ..},{}}]
aClusterQ[xx_]:=MatchQ[xx,{_?aRangeSpecQ..}]
aBestQ[xx_]:=MatchQ[xx,{_?VectorQ,(_?NumberQ|Infinity)}]
aBestListQ[xx_]:=MatchQ[xx,{_?aBestQ..}]


meldRes[theRes:{{_?NumberQ,{_?VectorQ..},_?aClusterQ,
_?aBestListQ,
{{_?aRangeSpecQ,_?aBestQ}..}}..}]:=
With[{theMax=Max[theRes[[All,1]]],
thePts=Flatten[theRes[[All,2]],1],
theClusters=Flatten[theRes[[All,3]],1],
theBestVals=Flatten[theRes[[All,4]],1],
theNewBounds=Flatten[theRes[[All,5]],1]},
{theMax,thePts,theClusters,theBestVals,theNewBounds}]


(*meldRes[xx___]:=Throw[{"meldRes bad Args",xx}]*)

unitToActual[lowerUpper:{{_?NumberQ,_?NumberQ}..},theVals_?MatrixQ]:=
Map[Function[xx,
MapThread[aUnitToActual,{lowerUpper,xx}]],theVals]

aUnitToActual[{aa_?NumberQ,bb_?NumberQ},aVal_?NumberQ]:=
Module[{theRes=aa+(bb-aa)*aVal},
theRes]

f1=Function[{xx,yy},-2*E^(-(xx^2+(yy^2)/2))]
f2=Function[{xx,yy},-(0.8*Cos[3*Pi*xx]+0.6*Cos[3*Pi*yy]-xx^2-2*yy^2+4)]
fu=Function[{xx,yy},-(1+.5*(Cos[4*Pi*(xx-.5)]+Cos[4*Pi*(yy-.5)])+
Abs[0.5*Cos[10*xx]*Cos[10*yy]])]

(*
https://mathematica.stackexchange.com/questions/83762/heuristic-method-to-create-low-discrepancy-sequences
*)



RandSuiteTore[x_, y_] :=
  Block[
        {n, dim, Stock, Stock1, step1, step2, step3, step4, Res},

        (* Initialisation  *)
         n = 1;
         dim = y;
         Stock = {};
         Stock1 = {};
        (* Boucle *)
         While[
                n <= dim,
                step1 = RandomReal[{0, 10^2}];
                step2 = NextPrime[step1, 1];
                step3 = Table[                                           
                              N@FractionalPart[i*Sqrt@step2],
                              {i, 1, x, 1}
                             ];
                Stock = Insert[Stock, step3, -1];
                Stock1 = Insert[Stock1, step2, -1];
                n++
              ];
          step4 = Partition[MapThread[## &, Stock], Length@Stock];
          Res = {step4, Stock1}
       ];

SuiteTore[x_, y_, z_] :=
  Block[
        {n, dim, Stock, Stock1, step1, step2, step3, step4, Res},
         (* Initialisation  *)
         n = 1;
         dim = y;
         Stock = {};
         Stock1 = {};
         (* Boucle *)
         While[
               n <= dim,
               step2 = z[[n]];
               step3 = Table[                                       
                             N@FractionalPart[i*Sqrt@step2],
                             {i, 1, x, 1}
                            ];
               Stock = Insert[Stock, step3, -1];
               Stock1 = Insert[Stock1, step2, -1];
               n++
              ];
        Res = Partition[MapThread[## &, Stock], Length@Stock]
       ];


Sobol[x_, y_] :=
  Block[
        {step1, step2},
         step1 = ToString@RandomReal[];
         step2 =
                BlockRandom[
                            SeedRandom[
                                       step1, 
                                       Method -> {"MKL", Method -> {"Sobol", "Dimension" -> y}}
                                      ];
                            RandomReal[{0, 1}, {x, y}]
                           ]
               ];

Nied[x_, y_] :=
  Block[
        {step1, step2},
         step1 = ToString@RandomReal[];
         step2 =
                 BlockRandom[
                             SeedRandom[
                                        step1,  
                                        Method -> {"MKL", Method -> {"Niederreiter", "Dimension" -> y}}
                                       ];

                             RandomReal[{0, 1}, {x, y}]
                            ]
               ];
(*
VanDerCorput[x_, y_,base_] := Table[
                              FromDigits[{Reverse[IntegerDigits[i, base]], 0}, y],
                              {i, x}
                             ];
myHalton[xx_,yy_]:=
N[RandHalton[xx,yy][[1]]]

RandHalton[x_, y_] :=
  Block[
        {step1, step2, step3, Res},
         step1 = Table[RandomReal[{0, 10^2}], {y}];
         step2 = NextPrime[step1, 1];
         step3 = Table[
                       VanDerCorput[x, i,Prime[Mod[i,x]]],
                       {i, step2}
                      ];
         Res = {Partition[MapThread[## &, step3], Length@step3], step2}
       ];

Halton[x_, y_, z_] :=
  Block[
        {step2, step3, Res},
         step2 = z;
         step3 = Table[
                       VanDerCorput[x, i,Prime[Mod[i,x]]],
                       {i, step2}
                      ];
         Res = Partition[MapThread[## &, step3], Length@step3]
       ];
*)
(*
https://mathematica.stackexchange.com/questions/634/error-checking-and-trapping-techniques-with-throw-and-catch
*)
General::interr = 
 "The function `1` failed due to an internal error. The failure occured in function `2`";



ClearAll[setConsistencyChecks];
Attributes[setConsistencyChecks] = {Listable};
setConsistencyChecks[function_Symbol, failTag_] :=
    function[___] := Throw[$Failed, failTag[function]];


ClearAll[catchInternalError];
Attributes[catchInternalError] = {HoldAll};
catchInternalError[code_, f_, failTag_] :=
  Catch[code, _failTag,
    Function[{value, tag},
      Message[General::interr , Style[f, Red], Style[First@tag, Red]];
      value]]; 

dispersion[theMat_?MatrixQ]:=
With[{theMinNorms=doAPoint[#,theMat]&/@Range[Length[theMat]]},
Max[theMinNorms]]
(*
dispersion[theMat_?MatrixQ]:=
With[{thptNow=First[theMat],restOfMat=Rest[theMat]},
With[{minNow=doAPoint[ptNow,restOfMat]},
Max[minNow,dispersion[restOfMat]]]]
*)
dispersion[{}]:=0


doAPoint[idx_Integer,theMat_?MatrixQ]:=
Min[Norm[theMat[[idx]]-#]&/@Drop[theMat,{idx}]]

doAPoint[thePoint_?VectorQ,{}]:=0



maxACluster[rhoVal_?NumberQ,remaining_?MatrixQ]:=
FixedPoint[growCluster,{{},rhoVal,remaining}]


maxACluster[rhoVal_?NumberQ,{}]:={}


getAllClusters[rhoVal_?NumberQ,remaining_?MatrixQ,
{lowerUpper:{{_?NumberQ,_?NumberQ}..},discreteVals:{{_Integer..}...}}]:=
	With[{ndRange=Range[Length[lowerUpper]]},
With[{
theRes=
NestWhileList[maxACluster[rhoVal,#[[-1]]]&,{remaining[[All,ndRange]]},(#[[-1]]=!={})&]},
     Print["gac1:theRes:",{theRes,ndRange,First /@ Drop[theRes,1],rhoVal,lowerUpper}];
     With[{nextLowerUppers=expandSingletons[#[[All,ndRange]],rhoVal,lowerUpper]&/@(
	First /@ Drop[theRes,1])},Print["gac2:",{rhoVal,remaining,lowerUpper,theRes,nextLowerUppers,discreteVals,{#[[ndRange]],discreteVals}&/@nextLowerUppers}];
				  {#,discreteVals}&/@nextLowerUppers]]]




End[]
EndPackage[]
