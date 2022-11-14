(* ::Package:: *)

(* ::Title:: *)
(*Spin-Weighted Spheroidal Functions*)


(* ::Section:: *)
(*Documentation *)


(* ::Section::Closed:: *)
(*Begin SWSpheroidal Package*)


BeginPackage["SWSpheroidal`"]


Unprotect[SWSphDebug];
SWSphDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[SWSphDebug,Unprotect["SWSpheroidal`*"];Unprotect["SWSpheroidal`Private`*"]];
Protect[SWSphDebug];


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*Spectral Solvers*)


SpinWeightedSpheroidal::usage=
	"SpinWeightedSpheroidal[m,s,c,N] solves the N-dimensional discrete "<>
	"approximation for the spin-weighted functions of spin-weight s, "<>
	"azimuthal index m, and oblateness parameter c.  The solution returns "<>
	"a list containing {values, vectors}, where values is a list of all "<>
	"the eigenvalues and vectors is a list of all the corresponding spectral "<>
	"coefficients.  These lists are sorted in order of the real part of the "<>
	"eigenvalues."


AngularSpectralRoot::usage=
	"AngularSpectralRoot[s,m,c,Alm,N] solves the N-dimensional discrete "<>
	"approximation for the spin-weighted functions of spin-weight s, "<>
	"azimuthal index m, and oblateness parameter c.  The solution returns "<>
	"the eigenvalue closest to Alm, and the corresponding spectral coefficients."


AngularSpectralRootIndex::usage=
	"AngularSpectralRoot[s,m,c,index,N] solves the N-dimensional discrete "<>
	"approximation for the spin-weighted functions of spin-weight s, "<>
	"azimuthal index m, and oblateness parameter c.  The solution returns "<>
	"the eigenvalue and spectral coefficients for the index-th eigenvalue."


(* ::Subsection::Closed:: *)
(*Normalization and Visualization*)


SWSFfixphase::usage=
"SWSFfixphase[m,s,L,SWdat]\n"<>
"\t m : azimuthal index\n"<>
"\t s : spin weight\n"<>
"\t L : eigenvalue index of SWdat : L = l - Max(|m|,|s|)\n"<>
"\t SWdat : eigenvector from AngularSectralRoot or AngularSpectralRootIndex\n"<>
"SWSFfixphase[m,s,L,SWdat] returns a complex phases correction that should "<>
"multiply SWdat in order to produce a properly phase-fixed spin-weighted "<>
"spheriodal function.\n\n"<>
"Options:\n"<>
"\t FixAt \[Rule] Provides preference for location used to set the phase. Defaults to Null.\n"<>
"\t\t\t If the SWSF is non-zero at one of the endpoints, the phase is set to\n"<>
"\t\t\t zero there.  If both endpoints are non-zero, FixAt must be set to\n"<>
"\t\t\t either +1 or -1.\n"<>
"\t InfoLevel \[Rule] Amount of diagnostic info to print. Defaults to 1.\n"<>
"\t ChopLevel \[Rule] Chops the evaluation of the SWSF at x=\[PlusMinus]1 to ChopLevel.\n"<>
"\t\t\t Defaults to 0."



SWSFvalues::usage=
"SWSFvalues[m,s,SWdat]\n"<>
"\t m : azimuthal index\n"<>
"\t s : spin weight\n"<>
"\t SWdat : eigenvector from AngularSectralRoot or AngularSpectralRootIndex\n"<>
"SWSFvalues[m,s,SWdat] returns a pair of lists {x,SWSF}.  "<>
"The first list x=Cos[theta] are the locations where the SWSF is evaluated.  "<>
"The second list contains the complex values of the SWSF at each value of x.\n\n"<>
"Options:\n"<>
"\t PlotPoints \[Rule] Number of vaulues of x. Defaults to 100."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[InfoLevel,FixAt];


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Spin-Weighted Spheroidal Functions : Spectral Method*)


(* ::Text:: *)
(*Evaluate the Spin-weighted Spheroidal Functions using a spectral method.*)


(* ::Subsection::Closed:: *)
(*Matrix Coefficients :*)


Flms[l_/;IntegerQ[2l],m_/;IntegerQ[2m],s_/;IntegerQ[2s]]:=Flms[l,m,s]=
									Piecewise[{{0,l<0},
												{Sqrt[(l+m+1)(l-m+1)/((2l+3)(2l+1))],s==0}},
												Sqrt[(l+m+1)(l-m+1)/((2l+3)(2l+1)) (l+s+1)(l-s+1)/((l+1)(l+1))]];


Glms[l_/;IntegerQ[2l],m_/;IntegerQ[2m],s_/;IntegerQ[2s]]:=Glms[l,m,s]=
									Piecewise[{{0,l<1}},
												Sqrt[(l+m)(l-m)/((2l+1)(2l-1)) (l+s)(l-s)/l^2]];


Hlms[l_/;IntegerQ[2l],m_/;IntegerQ[2m],s_/;IntegerQ[2s]]:=Hlms[l,m,s]=
									Piecewise[{{0,l==0||s==0}},
												-m s/(l(l+1))];


Alms[l_/;IntegerQ[2l],m_/;IntegerQ[2m],s_/;IntegerQ[2s]]:=Alms[l,m,s]=
									Flms[l+1,m,s]Flms[l,m,s];


Dlms[l_/;IntegerQ[2l],m_/;IntegerQ[2m],s_/;IntegerQ[2s]]:=Dlms[l,m,s]=
									Hlms[l+1,m,s]Flms[l,m,s]+Hlms[l,m,s]Flms[l,m,s];


Blms[l_/;IntegerQ[2l],m_/;IntegerQ[2m],s_/;IntegerQ[2s]]:=Blms[l,m,s]=
									Glms[l+1,m,s]Flms[l,m,s]+Glms[l,m,s]Flms[l-1,m,s] +Hlms[l,m,s]^2;


Elms[l_/;IntegerQ[2l],m_/;IntegerQ[2m],s_/;IntegerQ[2s]]:=Elms[l,m,s]=
									Glms[l,m,s]Hlms[l-1,m,s]+Glms[l,m,s]Hlms[l,m,s];


Clms[l_/;IntegerQ[2l],m_/;IntegerQ[2m],s_/;IntegerQ[2s]]:=Clms[l,m,s]=
									Glms[l,m,s]Glms[l-1,m,s];


(*If[!SWSphDebug,Protect[Alms,Blms,Clms,Dlms,Elms,Flms,Glms,Hlms]]; Don't Protect! *)


(* ::Subsection::Closed:: *)
(*Construct and solve the spectral matrix for Spin-weighted Oblate Spheroidal functions with azimuthal index m*)


Mat[i_Integer,j_Integer,m_/;IntegerQ[2m],s_/;IntegerQ[2s],c_?NumberQ]:=
Module[{lmin,l},
	lmin=Max[Abs[m],Abs[s]];
	l=lmin+(j-1);
	Return[Piecewise[
				{{l(l+1)-s(s+1)- c^2Blms[l,m,s]+2c s Hlms[l,m,s],i==j},
				{-c^2 Dlms[l,m,s]+2c s Flms[l,m,s],i-1==j},
				{-c^2Alms[l,m,s],i-2==j},
				{-c^2 Elms[l,m,s]+2c s Glms[l,m,s],i+1==j},
				{-c^2Clms[l,m,s],i+2==j}},
				0]];
]


SpinWeightedSpheroidal[m_/;IntegerQ[2m],s_/;IntegerQ[2s],c_?NumberQ,N_Integer] := 
Module[{unsorted,sortorder,i,j},
	unsorted=Eigensystem[Table[Mat[i,j,m,s,c],{i,1,N},{j,1,N}]];
	sortorder= Ordering[Re[unsorted[[1]]]];
	{unsorted[[1]][[sortorder]],unsorted[[2]][[sortorder]]}
]/;IntegerQ[m+s]
SpinWeightedSpheroidal[m_,s_,c_?NumberQ,N_Integer]:=Print["m and s must both be either integer or half-integer values."]


AngularSpectralRoot[s_,m_,c_?NumberQ,Alm_?NumberQ,N_Integer]:=
Module[{sol,diff,index},
	sol=SpinWeightedSpheroidal[m,s,c,N];
	diff = Abs[Alm - sol[[1]]];
	index=Position[diff,Min[diff]][[1,1]];
	{sol[[1,index]],N,sol[[2,index]]}
]


AngularSpectralRootIndex[s_,m_,c_?NumberQ,index_Integer,N_Integer]:=
Module[{sol},
	sol=SpinWeightedSpheroidal[m,s,c,N];
	{sol[[1,index]],N,sol[[2,index]]}
]


If[!QNMDebug,Protect[Mat,SpinWeightedSpheroidal,AngularSpectralRoot,AngularSpectralRootForl]];


(* ::Section::Closed:: *)
(*Normalization and Visualization*)


Options[SWSFfixphase]={InfoLevel->1,FixAt->Null[],ChopLevel->0};


SWSFfixphase[m_/;IntegerQ[2m],s_/;IntegerQ[2s],La_,SWdat_List,opts:OptionsPattern[]]:=
Module[{NC,lmin,WDplus,WDzero,WDminus,scaledcoefs,SWSFplus,SWSFzero,SWSFminus,
		il=OptionValue[InfoLevel],cl=OptionValue[ChopLevel],phase},
	NC=Length[SWdat];
	lmin=Max[Abs[m],Abs[s]];
	WDplus=ParallelTable[WignerD[{j-1+lmin,m,-s},0,0,0],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	WDzero=ParallelTable[WignerD[{j-1+lmin,m,-s},0,\[Pi]/2,0],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	WDminus=ParallelTable[WignerD[{j-1+lmin,m,-s},0,\[Pi],0],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	scaledcoefs=(-1)^m Sqrt[\[Pi]] ParallelTable[Sqrt[2(j-1+lmin)+1]SWdat[[j]],{j,1,NC}];
	SWSFplus = WDplus . scaledcoefs;
	SWSFzero = WDzero . scaledcoefs;
	SWSFminus = WDminus . scaledcoefs;
	If[il>0,
		Print["Max coeff : ",MaximalBy[SWdat,Abs]," at ",Position[SWdat,MaximalBy[SWdat,Abs][[1]]]];
		Print["SWSF[+1] = ",SWSFplus];
		Print["SWSF[0] = ",SWSFzero];
		Print["SWSF[-1] = ",SWSFminus];
	];
	If[Chop[SWSFplus,cl]==0,
		If[Chop[SWSFminus,cl]==0,
			(*phase=Exp[\[ImaginaryI](-Arg[SWSFzero])];If[il>0,Print["Phase set at 0"]]*)
			phase=1;If[il>0,Print["Phase IGNORED at 0"]],
			phase=Exp[I(-Arg[(-1)^(La+Max[Abs[m],Abs[s]])SWSFminus])];If[il>0,Print["Phase set at -1"]],
			Print["Logic error 1"];Abort[]
		],
		If[Chop[SWSFminus,cl]==0,
			phase=Exp[I(-Arg[(-1)^m SWSFplus])];If[il>0,Print["Phase set at +1"]],
			Switch[OptionValue[FixAt],
				-1,phase=Exp[I(-Arg[(-1)^(La+Max[Abs[m],Abs[s]])SWSFminus])],
				1,phase=Exp[I(-Arg[(-1)^m SWSFplus])],
				_,Print["SWSF non-zero at both ends: set FixAt=\[PlusMinus]1"];Abort[]
			],
			Print["Logic error 2"];Abort[]
		],
		Print["Logic error 3"];Abort[]
	];
	phase
]/;IntegerQ[m+s]
SWSFfixphase[m_,s_,La_,SWdat_List,opts:OptionsPattern[]]:=Print["m and s must both be either integer or half-integer values."]


Options[SWSFvalues]={PlotPoints->100};


SWSFvalues[m_/;IntegerQ[2m],s_/;IntegerQ[2s],SWdat_List,opts:OptionsPattern[]]:=
Module[{NC,x,theta,Ntheta,lmin,Matdlx,SWSF,choplev,npoints=OptionValue[PlotPoints]},
	NC=Length[SWdat];
	x=Table[x,{x,-1,1,1/(npoints-1)}];
	theta=ArcCos[#]&/@x;
	Ntheta=Length[theta];
	lmin=Max[Abs[m],Abs[s]];
	Matdlx = ParallelTable[N[WignerD[{j-1+lmin,m,-s},0,theta[[k]],0]],{k,1,Ntheta},{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	SWSF = (-1)^m Sqrt[\[Pi]] Matdlx . ParallelTable[Sqrt[2(j-1+lmin)+1]SWdat[[j]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	{x,SWSF}
]


(* ::Section::Closed:: *)
(*End of SWSpheroidal Package*)


End[] (* `Private` *)


EndPackage[]
