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
"SpinWeightedSpheroidal[m,s,c,N] "<>
"gives a list {\!\(\*StyleBox[\"values\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"vectors\",\nFontSlant->\"Italic\"]\)} from solving the "<>
"\!\(\*StyleBox[\"N\", \"TI\"]\)-dimensional approximate discrete eigensystem for "<>
"the spin-weighted spheroidal functions with spin-weight "<>
"\!\(\*StyleBox[\"s\", \"TI\"]\), azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), and oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\)."


AngularSpectralRoot::usage=
"AngularSpectralRoot[s,m,c,\!\(\*SubscriptBox[\(\[Lambda]\), \(0\)]\),N] "<>
"gives a list {\!\(\*StyleBox[\"\[Lambda]\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"N\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"vector\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"index\",\nFontSlant->\"Italic\"]\)} from solving the "<>
"\!\(\*StyleBox[\"N\", \"TI\"]\)-dimensional approximate discrete eigensystem for "<>
"the spin-weighted spheroidal functions with spin-weight "<>
"\!\(\*StyleBox[\"s\", \"TI\"]\), azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), and oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\). The "<>
"solution returns the eigenvalue closest to \!\(\*StyleBox[SubscriptBox[StyleBox[\"\[Lambda]\", \"TI\"], \"0\"], \"TI\"]\), and the "<>
"corresponding spectral coefficients.  "<>
"\!\(\*StyleBox[\"index\", \"TI\"]\) denotes the position of the eigensolution in the full list."


AngularSpectralRootIndex::usage=
"AngularSpectralRoot[s,m,c,index,N] "<>
"gives a list {\!\(\*StyleBox[\"\[Lambda]\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"N\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"vector\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"index\",\nFontSlant->\"Italic\"]\)} from solving the "<>
"\!\(\*StyleBox[\"N\", \"TI\"]\)-dimensional approximate discrete eigensystem for "<>
"the spin-weighted spheroidal functions with spin-weight "<>
"\!\(\*StyleBox[\"s\", \"TI\"]\), azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), and oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\). The "<>
"solution returns the eigenvalue and spectral coefficients "<>
"for the \!\(\*StyleBox[\"index\",\nFontSlant->\"Italic\"]\)-th eigensolution."


(* ::Subsection::Closed:: *)
(*Normalization and Visualization*)


SWSFfixphase::usage=
"SWSFfixphase[m,s,L,SWdat] "<>
"gives a complex phases correction that should multiply "<>
"\!\(\*StyleBox[\"SWdat\", \"TI\"]\) in order to produce a properly phase-fixed "<>
"spin-weighted spheroidal function with spin-weight \!\(\*StyleBox[\"s\", \"TI\"]\), "<>
"azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\)."


SLCorrectExtremumDiscontinuity::usage=
"SLCorrectExtremumDiscontinuity[m,s,L,SWdat,SLPhase] "<>
"returns a list of complex phase values which have been corrected to "<>
"to remove phase discontinuities of \[Pi] due to extrema crossing x=0."


SLCorrectInitialDiscontinuity::usage=
"SLCorrectInitialDiscontinuity[m,s,L,SWdat,SLPhase,\[Phi]dir] "<>
"returns a 2-element list {update,phase}.  If update is True, then phase "<>
"contains the complex phase correction which should replace the first "<>
"element in SLPhase."


SWSFvalues::usage=
"SWSFvalues[m,s,SWdat] "<>
"evaluates the spin-weighted spheroidal function with spin-weight "<>
"\!\(\*StyleBox[\"s\", \"TI\"]\) and azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), "<>
"represented by the spectral expansion coefficients "<>
"\!\(\*StyleBox[\"SWdat\", \"TI\"]\)\!\(\*StyleBox[\".\", \"TI\"]\)"


SWSFRealPath::usage=
"SWSFRealPath[m,s,SWdat] "<>
"returns the complex path between the poles at \!\(\*
StyleBox[\"z\", \"TI\"]\)=\[PlusMinus]1 along which the spin-weighted "<>
"spheroidal function with spin-weight \!\(\*
StyleBox[\"s\", \"TI\"]\) and azimuthal index \!\(\*
StyleBox[\"m\", \"TI\"]\), represented by the "<>
"spectral expansion coefficients \!\(\*
StyleBox[\"SWdat\", \"TI\"]\), is real."


(* ::Subsection::Closed:: *)
(*Confluent Heun Solvers*)


SWSpheroidalWronskianRoot::usage=
SpinWeightedSpheroidal::usage=
"SWSpheroidalWronskianRoot[s,m,c,Alm] "<>
"returns an eigenvalue near \!\(\*StyleBox[\"Alm\", \"TI\"]\) associated with a spin-weighted spheroidal "<>
"function with spin-weight \!\(\*StyleBox[\"s\", \"TI\"]\), "<>
"azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), and oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\)."


SWSFHeunPhases::usage=
"SWSFHeunPhases[s,m,l,c,Alm] "<>
"returns a pair of phase factors which will implement the Cook-Wang SL phase choice "<>
"for the spin-weighted spheroidal function with spin-weight \!\(\*StyleBox[\"s\", \"TI\"]\), "<>
"azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\), "<>
"and eigenvalue \!\(\*StyleBox[\"Alm\", \"TI\"]\)."


SWSFHeunNorms::usage=
"SWSFHeunNorms[s,m,c,Alm] "<>
"returns a pair of normalizatin factors which will properly normalize the "<>
"spin-weighted spheroidal function with spin-weight \!\(\*StyleBox[\"s\", \"TI\"]\), "<>
"azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\), "<>
"and eigenvalue \!\(\*StyleBox[\"Alm\", \"TI\"]\)."


SWSFHeun::usage-
"SWSFHeun[norm,s,m,c,Alm,x] "<>
"evaluates the spin-weighted spheroidal function with spin-weight \!\(\*StyleBox[\"s\", \"TI\"]\), "<>
"azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\), "<>
"and eigenvalue \!\(\*StyleBox[\"Alm\", \"TI\"]\) at location \!\(\*StyleBox[\"x\", \"TI\"]\) "
"in the range \!\(\*StyleBox[\"-1\[LessEqual]x\[LessEqual]1\", \"TI\"]\)-1\[LessEqual]x\[LessEqual]1. "<>
"\!\(\*StyleBox[\"norm\", \"TI\"]\) is a 2-element list which "<>
"provides the desired normalizatin and phase choice for the function."


SWSFHeunExpansionCoefs::usage=
"SWSFHeunExpansionCoefs[s,m,l,c,Alm] "<>
"returns a list of complex expansion coefficients necessary to represent the "<>
"spin-weighted spheroidal function with spin-weight \!\(\*StyleBox[\"s\", \"TI\"]\), "<>
"azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\), "<>
"and eigenvalue \!\(\*StyleBox[\"Alm\", \"TI\"]\) as a linear combination of sphin-weighted "<>
"spherical functions."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[PhaseChoice,SphericalLimit,Simple];


Protect[PlotFlips]


Protect[PathStart,PlotStart,PrintPoleValues,StepSize,\[Phi]guess]


Protect[ChopDelta,CoefDelta]


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
	sortorder= Ordering[Abs[unsorted[[1]]+s]]; (* to ensure symmetrical ordering with s->-s *)
	{unsorted[[1]][[sortorder]],unsorted[[2]][[sortorder]]}
]/;IntegerQ[m+s]
SpinWeightedSpheroidal[m_,s_,c_?NumberQ,N_Integer]:=Print["m and s must both be either integer or half-integer values."]


AngularSpectralRoot[s_,m_,c_?NumberQ,Alm_?NumberQ,N_Integer]:=
Module[{sol,diff,index},
	sol=SpinWeightedSpheroidal[m,s,c,N];
	diff = Abs[Alm - sol[[1]]];
	index=Position[diff,Min[diff]][[1,1]];
	{sol[[1,index]],N,sol[[2,index]],index}
]


AngularSpectralRootIndex[s_,m_,c_?NumberQ,index_Integer,N_Integer]:=
Module[{sol},
	sol=SpinWeightedSpheroidal[m,s,c,N];
	{sol[[1,index]],N,sol[[2,index]],index}
]


If[!QNMDebug,Protect[Mat,SpinWeightedSpheroidal,AngularSpectralRoot,AngularSpectralRootForl]];


(* ::Section::Closed:: *)
(*Normalization and Visualization*)


Options[SWSFfixphase]={PhaseChoice->SphericalLimit,Tolerance->10^(-10)};


SWSFfixphase[m_/;IntegerQ[2m],s_/;IntegerQ[2s],La_,SWdat_List,opts:OptionsPattern[]]:=
Module[{NC,lmin,SphericalVal,SphericalD,WDzero,WDzeroplus,WDzerominus,WDzeroD,scaledcoefs,SWSFzero,SWSFzeroD,phase,index,hoint=1,\[Delta]=OptionValue[Tolerance]},
	SWSFfixphase::badmethod="`1` is an invalid Method";
	SWSFfixphase::notpossible="Spherical value and derivative both zero.";
	If[!IntegerQ[m],hoint=-I]; (* Choice to make half-odd integer Spherical harmonics real *)
	NC=Length[SWdat];
	lmin=Max[Abs[m],Abs[s]];
	Switch[OptionValue[PhaseChoice],
	Simple,
		phase=Exp[I(-Arg[SWdat[[La+1]]])],
	SphericalLimit,
		SphericalVal=hoint (-1)^s Sqrt[La+lmin+1/2]*WignerD[{La+lmin,-m,s},0,\[Pi]/2,0];
		scaledcoefs=hoint (-1)^s ParallelTable[Sqrt[j-1+lmin+1/2]SWdat[[j]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
		WDzero=ParallelTable[WignerD[{j-1+lmin,-m,s},0,\[Pi]/2,0],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
		SWSFzero = WDzero . scaledcoefs;
		(* In general, we choose the phase so the function is real at x=0 *)
		phase=Exp[I(\[Pi]-Arg[SWSFzero])];
		If[Sign[SphericalVal]!=0,
			If[Chop[SWSFzero,\[Delta]]==0,
			(* Case where function vanishes, but spherical limit does not. *)
				WDzeroplus=ParallelTable[Sqrt[(j-1+lmin-s)*(j-1+lmin+s+1)]*If[s>=(j-1+lmin),0,WignerD[{j-1+lmin,-m,s+1},0,\[Pi]/2,0]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
				WDzerominus=ParallelTable[Sqrt[(j-1+lmin+s)*(j-1+lmin-s+1)]*If[-s>=(j-1+lmin),0,WignerD[{j-1+lmin,-m,s-1},0,\[Pi]/2,0]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
				WDzeroD=(-1/2)(WDzeroplus-WDzerominus);
				SphericalD=hoint (-1)^s Sqrt[La+lmin+1/2]WDzeroD[[La+1]];
				SWSFzeroD=scaledcoefs . WDzeroD;
				(* Since the function is zero at x=0, we can demand the derivative be real there *)
				(* and have a specified slope *)
				If[Chop[SWSFzeroD,\[Delta]]==0||Chop[SphericalD,\[Delta]]==0,
				(* This phase choice will only be used in the case where the function     *)
				(* and its derivative (or the spherical limit derivative) are effectively *)
				(* zero at x=0.  This cannot happen exactly, but can occur in the         *)
				(* asymptotic regime where expnential damping can occur at x=0.           *)
					phase=Exp[I(-Arg[SWdat[[La+1]]])],(* When all else fails, use the simple choice *)
				(* normal case, fix the derivative to be real *)
					phase=Exp[I(\[Pi]-Arg[SWSFzeroD])];
					If[Sign[Re[phase SWSFzeroD]]!= Sign[SphericalD],phase=-phase]
				],
			(* general case *)
				If[Sign[Re[phase SWSFzero]]!=Sign[SphericalVal],phase=-phase]
			]
		,(*Case when spherical value is zero.*)
			WDzeroplus=ParallelTable[Sqrt[(j-1+lmin-s)*(j-1+lmin+s+1)]*If[s>=(j-1+lmin),0,WignerD[{j-1+lmin,-m,s+1},0,\[Pi]/2,0]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
			WDzerominus=ParallelTable[Sqrt[(j-1+lmin+s)*(j-1+lmin-s+1)]*If[-s>=(j-1+lmin),0,WignerD[{j-1+lmin,-m,s-1},0,\[Pi]/2,0]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
			WDzeroD=(-1/2)(WDzeroplus-WDzerominus);
			SphericalD=hoint (-1)^s Sqrt[La+lmin+1/2]WDzeroD[[La+1]];
			SWSFzeroD=scaledcoefs . WDzeroD;
			If[Chop[SWSFzeroD,\[Delta]]==0,
			(* Case when derivative vanishes, but spherical limit does not. *)
				If[Chop[SWSFzero,\[Delta]]==0,
				(* This phase choice will only be used in the case where the function  *)
				(* and its derivative are effectively zero at x=0.  This cannot happen *)
				(* exactly, but can occur in the asymptotic regime where expnential    *)
				(* damping can occur at x=0.                                           *)
					phase=Exp[I(-Arg[SWdat[[La+1]]])],(* When all else fails, use the simple choice *)
				(* Since the derivative vanishes, the value should be non-zero and the phase *)
				(* has been chosen to make it real at x=0, so we can specify its sign        *)
					index=La+Max[Abs[m],Abs[s]]+m;
					If[OddQ[index],--index];
					If[Sign[Re[phase SWSFzero]]!= Sign[(-1)^(index/2)],phase=-phase];
				],
			(* general case *)
				phase=Exp[I(\[Pi]-Arg[SWSFzeroD])];
				If[Sign[Re[phase SWSFzeroD]]!=Sign[SphericalD],phase=-phase]
			];
			If[Sign[SphericalD]==0,Message[SWSFfixphase::notpossible];Abort[]];
		],
	_,Message[SWSFfixphase::badmethod,OptionValue[Method]];Abort[] 
	];
	phase
]/;IntegerQ[m+s]
SWSFfixphase[m_,s_,La_,SWdat_List,opts:OptionsPattern[]]:=Print["m and s must both be either integer or half-integer values."]


Options[SLCorrectExtremumDiscontinuity]={PlotFlips->False};


SLCorrectExtremumDiscontinuity[m_/;IntegerQ[2m],s_/;IntegerQ[2s],La_Integer,SWdat_List,SLPhase_List,opts:OptionsPattern[]]:=
Module[{slphase,CZPhase,PhaseDiff,flips,avg,err,flipvals,flipon,flipindex,dropped=0},
	SLCorrectExtremumDiscontinuity::SLerror="Should only need to phase flip if spherical limit vanishes at x=0.  "
		<>"Possible vansishing of fiducial CZ expansion coefficient!";
	CZPhase=Table[Exp[I(-Arg[SWdat[[ind,La+1]]])],{ind,1,Length[SWdat]}];
	PhaseDiff=Arg[Exp[I(Arg[SLPhase ]-Arg[CZPhase ])]];
	flips=Abs[Abs[Drop[PhaseDiff,1]-Drop[PhaseDiff,-1]]-\[Pi]]-\[Pi];
	avg=(Total[flips]-Max[flips])/(Length[flips]-1);
	err=Sqrt[Total[Abs[flips-avg]^2]]/Length[flips];
	flips =Abs[flips+\[Pi]];
	Print["Average diff : ",avg," : err : ",err];
	If[OptionValue[PlotFlips],Print[ListPlot[flips,PlotRange->All]]];
	flipvals=Select[flips,Abs[#]<5err&];
	If[Length[flipvals]==0,Return[SLPhase]];
	If[s==0,Return[SLPhase]]; (* Symmetry prevents problems *)
	If[WignerD[{La+Max[Abs[m],Abs[s]],m,-s},0,\[Pi]/2,0]!=0,Message[SLCorrectExtremumDiscontinuity::SLerror]];
	slphase=SLPhase;
	Do[
		flipindex=Position[flips,flipon][[1,1]];
		Print["Flip index : ",flipindex+dropped];
		slphase=Join[Take[slphase,flipindex+dropped],-Drop[slphase,flipindex+dropped]];
		flips=Drop[flips,flipindex];
		dropped +=flipindex;
	,{flipon,flipvals}];
	slphase
]/;IntegerQ[m+s]
SLCorrectExtremumDiscontinuity[m_,s_,La_,SWdat_,SLPhase_,opts:OptionsPattern[]]:=Print["m and s must both be either integer or half-integer values."]


SLCorrectInitialDiscontinuity[m_/;IntegerQ[2m],s_/;IntegerQ[2s],La_Integer,SWdat_List,SLPhase_List,\[Phi]dir_]:=
Module[{CZPhase,PhaseDiff,lmin,phase0,extrap,x,error},
	SLCorrectInitialDiscontinuity::tooshort="SWdat and SLPhase must be at least 4 elemnts long.";
	SLCorrectInitialDiscontinuity::discontinuous="The computed phase is not continuous.";
	SLCorrectInitialDiscontinuity::phasecorrected="The expansion coefficents seem to already be phase corrected.";
	If[Min[Length[SWdat],Length[SLPhase]]<5,Message[SLCorrectInitialDiscontinuity::tooshort];Abort[]];
	If[SWdat[[1,La+1]]!=1,Message[SLCorrectInitialDiscontinuity::phasecorrected];Abort[]];
	CZPhase=Table[Exp[I(-Arg[SWdat[[ind,La+1]]])],{ind,1,Min[4,Length[SWdat],Length[SLPhase]]}];
	PhaseDiff=Arg[Exp[I(Arg[Take[SLPhase,Min[4,Length[SWdat],Length[SLPhase]]]]-Arg[CZPhase ])]];
	lmin=Max[Abs[m],Abs[s]];
	If[WignerD[{La+lmin,-m,s},0,\[Pi]/2,0]!=0,Return[{False,1}],If[s==0,Return[{False,1}]]];
	phase0=Exp[-I \[Phi]dir];
	If[Abs[Arg[phase0]-PhaseDiff[[2]]]>Abs[Arg[phase0]+PhaseDiff[[2]]],phase0=-phase0];
	extrap=Fit[Transpose[{{1,2,3},Drop[PhaseDiff,1]}],{1,x,x^2},x]/.x->0;
	error=Abs[extrap-Arg[phase0]];
	If[error>10^(-12) && error>5Abs[PhaseDiff[[2]]-PhaseDiff[[4]]],Message[SLCorrectInitialDiscontinuity::discontinuous];Abort[]];
	{True,phase0}
]/;IntegerQ[m+s]
SLCorrectInitialDiscontinuity[m_,s_,La_,SWdat_,SLPhase_,\[Phi]dir_]:=Print["m and s must both be either integer or half-integer values."]


Options[SWSFvalues]={PlotPoints->100};


SWSFvalues[m_/;IntegerQ[2m],s_/;IntegerQ[2s],SWdat_List,opts:OptionsPattern[]]:=
Module[{NC,x,theta,Ntheta,lmin,Matdlx,SWSF,choplev,hoint=1,npoints=OptionValue[PlotPoints]},
	If[!IntegerQ[m],hoint=I];
	NC=Length[SWdat];
	x=Table[x,{x,-1,1,1/(npoints-1)}];
	theta=ArcCos[#]&/@x;
	Ntheta=Length[theta];
	lmin=Max[Abs[m],Abs[s]];
	Matdlx = ParallelTable[N[WignerD[{j-1+lmin,-m,s},0,theta[[k]],0]],{k,1,Ntheta},{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	SWSF = hoint (-1)^s Matdlx . ParallelTable[Sqrt[j-1+lmin+1/2]SWdat[[j]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	{x,SWSF}
]


Options[SWSFRealPath]=Union[{PathStart->Automatic,MaxSteps->5000,StepSize->1/50,PlotStart->False,\[Phi]guess->Automatic,PrintPoleValues->False}]


SWSFRealPath[m_/;IntegerQ[2m],s_/;IntegerQ[2s],SWdat_List,opts:OptionsPattern[]]:=
Module[{NC,lmin,Matdlx,SWSF,z,z0,zf,\[Delta]z,\[Phi],\[Phi]g,\[Phi]0,phase,zlist,SWSFlist,count=0,hoint=1,
		max\[Delta]r=OptionValue[StepSize],maxcount=OptionValue[MaxSteps],
		plot=OptionValue[PlotStart]},
	SWSFRealPath::PathStart="Invalid PathStart : `1`";
	If[!IntegerQ[m],hoint=I];
	NC=Length[SWdat];
	lmin=Max[Abs[m],Abs[s]];
	Matdlx = ParallelTable[WignerD[{j-1+lmin,-m,s},0,ArcCos[z],0],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	SWSF = hoint (-1)^s Matdlx . Table[Sqrt[j-1+lmin+1/2]SWdat[[j]],{j,1,NC}];
	phase=1;
	Switch[OptionValue[PathStart],
		Automatic,
			If[Abs[N[SWSF/.z->-1]]>10^(-10),phase=Exp[-I Arg[N[SWSF/.z->-1]]];z0=-1;zf=1;\[Phi]g=0,
				If[Abs[N[SWSF/.z->+1]]>10^(-10),phase=Exp[-I Arg[N[SWSF/.z->+1]]];z0=1;zf=-1;\[Phi]g=\[Pi],
				z0=-1;zf=1;\[Phi]g=0]
			],
		+1,
			phase=Exp[-I Arg[N[SWSF/.z->+1]]];z0=+1;zf=-1;\[Phi]g=\[Pi],
		-1,
			phase=Exp[-I Arg[N[SWSF/.z->-1]]];z0=-1;zf=+1;\[Phi]g=0,
		_,Message[SWSFRealPath::PathStart,OptionValue[PathStart]];Abort[]
	];
	If[OptionValue[\[Phi]guess]==Automatic,Null[],\[Phi]g=OptionValue[\[Phi]guess],\[Phi]g=OptionValue[\[Phi]guess]];
	SWSF*=phase; \[Phi]0=\[Phi]g;
	If[OptionValue[PrintPoleValues]==True,
		Print["SWSF(-1) = ",N[SWSF/.z->-1]];
		Print["SWSF(+1) = ",N[SWSF/.z->1]];
		Print["z0 = ",z0," : zf = ",zf," : \[Phi]g = ",\[Phi]g];
	];
	zlist={z0};
	SWSFlist={Chop[N[SWSF/.z->z0]]};
	While[Abs[z0-zf]>10^(-8)&& count++<maxcount,
		\[Delta]z=Min[max\[Delta]r,Abs[z0-zf]/2];
		If[plot && count==1,Print[Plot[Im[SWSF/.z->z0+\[Delta]z Exp[I \[Phi]]],{\[Phi],0,2\[Pi]}]]];
		\[Phi]g=FindRoot[Im[SWSF/.z->z0+\[Delta]z Exp[I \[Phi]]]==0,{\[Phi],\[Phi]g}][[1,2]];
		z0+=\[Delta]z Exp[I \[Phi]g];
		AppendTo[zlist,z0];AppendTo[SWSFlist,Chop[N[SWSF/.z->z0]]];
	];
	If[zf==-1,zlist=Reverse[zlist];SWSFlist=Reverse[SWSFlist]];
	{zlist,SWSFlist}
]


(* ::Section::Closed:: *)
(*Spin-Weighted Spheroidal Functions : Confluent Heun Method*)


(* ::Text:: *)
(*Evaluate the Spin-weighted Spheroidal Functions using confluent Heun functions.  Roots of the Wronskian determine values of the Eigenvalues, \!\( *)
(*\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)*)
(*\(\*SubscriptBox[\(A\), \(\[ScriptL], m\)]\)\)[c], of the angular Teukolsky which allow local Frobenius solutions to be regular simultaneously at both regular singular points.*)


(* ::Subsection::Closed:: *)
(*Spin-weighted spheroidal functions which are regular at either x=\[PlusMinus]1:*)


(* ::Text:: *)
(*Each version of the spin-weighted spheroidal function is well behaved for certain values of m+s or m-s.*)


(* ::Subsubsection::Closed:: *)
(*Regular solutions at x=+1, valid for -1<x<=1:*)


SWSFHeun[1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1-x)^((m+s)/2) (1+x)^((m-s)/2) HeunC[Alm+c^2+2 c (m+1)-m (m+1)+s(s+1),4 c (m-s+1),m+s+1,m-s+1,4 c,(1-x)/2]/;(m+s>0)&&(m-s>=0)


SWSFHeun[1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1+x)^-s HeunC[Alm+c^2-2 c (-1+s)+2 s,c (4-8 s),1,1-2 s,4 c,(1-x)/2]/;(m+s==0)&&(m-s>=0)


SWSFHeun[1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1-x)^((m+s)/2) (1+x)^((s-m)/2) HeunC[Alm+c^2+2c (m+1),4 c,m+s+1,s-m+1,4 c,(1-x)/2]/;(m+s>0)&&(m-s<=0)


SWSFHeun[1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1+x)^s HeunC[Alm+c (2+c+2 m),4 c,1,1+2 s,4 c,(1-x)/2]/;(m+s==0)&&(m-s<=0)


SWSFHeun[1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1-x)^(-((m+s)/2)) (1+x)^((s-m)/2) HeunC[Alm+c^2-2 c (m+2 s-1)-m (m-1)+s(s+1),-4 c (m+s-1),-m-s+1,s-m+1,4 c,(1-x)/2]/;(m+s<0)&&(m-s<=0)


SWSFHeun[1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1-x)^(-((m+s)/2)) (1+x)^((m-s)/2) HeunC[Alm+c^2-2c ( m+2 s-1)+2 s,-4 c (2 s-1),-m-s+1,m-s+1,4 c,(1-x)/2]/;(m+s<0)&&(m-s>=0)


(* ::Subsubsection::Closed:: *)
(*Regular solutions at x=-1, valid for -1<=x<1:*)


SWSFHeun[-1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1-x)^((m+s)/2) (1+x)^((m-s)/2) HeunC[Alm+c^2-2 c (m-2 s+1)-m (m+1)+s(s+1),-4 c (m-s+1),m-s+1,m+s+1,-4 c,(1+x)/2]/;(m+s>=0)&&(m-s>0)


SWSFHeun[-1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1-x)^s HeunC[Alm+c (-2+c+2 s),-4 c,1,1+m+s,-4 c,(1+x)/2]/;(m+s>=0)&&(m-s==0)


SWSFHeun[-1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1-x)^(-((m+s)/2)) (1+x)^((m-s)/2) HeunC[Alm+c^2-2c (m-2 s+1)+2 s,4 c (2 s-1),m-s+1,-m-s+1,-4 c,(1+x)/2]/;(m+s<=0)&&(m-s>0)


SWSFHeun[-1][s_,m_,c_,Alm_,x_]:=E^(-c x) (1-x)^-s HeunC[Alm+c^2+2 c (-1+s)+2 s,4 c (-1+2 s),1,1-2 s,-4 c,(1+x)/2]/;(m+s<=0)&&(m-s==0)








(* ::Subsection::Closed:: *)
(*Roots of the spin-weighted spheroidal function Wronskian:*)


(* ::Text:: *)
(*Each version of the Wronskian is well behaved for certain values of m+s and m-s.*)


SWSpheroidalWronskian[s_,m_,c_,Alm_]:=1/2 (HeunC[Alm+c^2-m (1+m)-2 c (1+m-2 s)+s+s^2,-4 c (1+m-s),1+m-s,1+m+s,-4 c,1/2] HeunCPrime[Alm+c^2+2 c (1+m)-m (1+m)+s+s^2,4 c (1+m-s),1+m+s,1+m-s,4 c,1/2]+HeunC[Alm+c^2+2 c (1+m)-m (1+m)+s+s^2,4 c (1+m-s),1+m+s,1+m-s,4 c,1/2] HeunCPrime[Alm+c^2-m (1+m)-2 c (1+m-2 s)+s+s^2,-4 c (1+m-s),1+m-s,1+m+s,-4 c,1/2])/;(m+s>=0)&&(m-s>=0)


SWSpheroidalWronskian[s_,m_,c_,Alm_]:=1/2 (HeunC[Alm+c (2+c+2 m),4 c,1+m+s,1-m+s,4 c,1/2] HeunCPrime[Alm+c (-2+c+2 m),-4 c,1-m+s,1+m+s,-4 c,1/2]+HeunC[Alm+c (-2+c+2 m),-4 c,1-m+s,1+m+s,-4 c,1/2] HeunCPrime[Alm+c (2+c+2 m),4 c,1+m+s,1-m+s,4 c,1/2])/;(m+s>=0)&&(m-s<=0)


SWSpheroidalWronskian[s_,m_,c_,Alm_]:=1/2 (HeunC[Alm+2 s+c (-2+c-2 m+4 s),4 c (-1+2 s),1+m-s,1-m-s,-4 c,1/2] HeunCPrime[Alm+c (2+c-2 m-4 s)+2 s,c (4-8 s),1-m-s,1+m-s,4 c,1/2]+HeunC[Alm+c (2+c-2 m-4 s)+2 s,c (4-8 s),1-m-s,1+m-s,4 c,1/2] HeunCPrime[Alm+2 s+c (-2+c-2 m+4 s),4 c (-1+2 s),1+m-s,1-m-s,-4 c,1/2])/;(m+s<=0)&&(m-s>=0)


SWSpheroidalWronskian[s_,m_,c_,Alm_]:=1/2 (HeunC[Alm+c^2+m-m^2+s+s^2-2 c (-1+m+2 s),-4 c (-1+m+s),1-m-s,1-m+s,4 c,1/2] HeunCPrime[Alm+(-2+c) c+m+2 c m-m^2+s+s^2,4 c (-1+m+s),1-m+s,1-m-s,-4 c,1/2]+HeunC[Alm+(-2+c) c+m+2 c m-m^2+s+s^2,4 c (-1+m+s),1-m+s,1-m-s,-4 c,1/2] HeunCPrime[Alm+c^2+m-m^2+s+s^2-2 c (-1+m+2 s),-4 c (-1+m+s),1-m-s,1-m+s,4 c,1/2])/;(m+s<=0)&&(m-s<=0)


Options[SWSpheroidalWronskianRoot]=Union[{Quiet->FindRoot::precw},Options[FindRoot]];
SWSpheroidalWronskianRoot[s_/;IntegerQ[2s],m_/;IntegerQ[2m],c_?NumberQ,Almguess_?NumberQ,opts:OptionsPattern[]]:=Module[{Alm,quiet=OptionValue[Quiet]},
Alm/.Quiet[FindRoot[SWSpheroidalWronskian[s,m,c,Alm],{Alm,Almguess},Evaluate[FilterRules[{opts},Options[FindRoot]]]],Evaluate[quiet]]
]/;IntegerQ[m+s]


(* ::Subsection::Closed:: *)
(*Evaluate the spin-weighted spheroidal functions:*)


Options[SWSFHeunPhases]={ChopDelta->10^-10};


SWSFHeunPhases[s_/;IntegerQ[2s],m_/;IntegerQ[2m],l_/;IntegerQ[2l],c_?NumberQ,Alm_?NumberQ,OptionsPattern[]]:=
Module[{x,sphericalval,Dsphericalval,valp,valm,Dvalp,Dvalm,phases,signs,index,hoint=1,\[Epsilon]=OptionValue[ChopDelta]},
	If[!IntegerQ[m],hoint=-I]; (* Choice to make half-odd integer Spherical harmonics real *)
	sphericalval=hoint (-1)^s Sqrt[l+1/2]WignerD[{l,-m,s},0,\[Pi]/2,0];
	valp=SWSFHeun[1][s,m,c,Alm,0];
	valm=SWSFHeun[-1][s,m,c,Alm,0];
	(* general case phase *)
	phases=Exp[I(\[Pi]-{Arg[valm],Arg[valp]})];
	If[Sign[sphericalval]!=0,
		If[Chop[valp,\[Epsilon]]==0,
		(* Case where function vanishes, but spherical limit does not. *)
			Dsphericalval=hoint ((-1)^s)/2(-Sqrt[(l-s)(l+s+1)(l+1/2)]WignerD[{l,-m,s+1},0,\[Pi]/2,0]+Sqrt[(l+s)(l-s+1)(l+1/2)]WignerD[{l,-m,s-1},0,\[Pi]/2,0]);
			Dvalp=D[SWSFHeun[1][s,m,c,Alm,x],x]/.x->0;
			If[Chop[Dvalp,\[Epsilon]]==0||Chop[Dsphericalval,\[Epsilon]]==0,
			(* This shouldn't happen, but could in the asymptotic regime where where exponential damping can occur *)
				Return[{Indeterminate,Indeterminate}];
			,(* normal case, fix the derivative to be real *)
				Dvalm=D[SWSFHeun[-1][s,m,c,Alm,x],x]/.x->0;
				phases=Exp[I(\[Pi]-{Arg[Dvalm],Arg[Dvalp]})];
				signs=Sign[phases {Dvalm,Dvalp}];
				phases*=signs Sign[Dsphericalval];
			]
			,(* general case: Check the signs for value *)
				signs=Sign[phases {valm,valp}];
				phases*=signs Sign[sphericalval];
		]
	,(*Case when spherical value is zero.*)
		Dsphericalval=hoint (-1)^s/2(-Sqrt[(l-s)(l+s+1)(l+1/2)]WignerD[{l,-m,s+1},0,\[Pi]/2,0]+Sqrt[(l+s)(l-s+1)(l+1/2)]WignerD[{l,-m,s-1},0,\[Pi]/2,0]);
		Dvalp=D[SWSFHeun[1][s,m,c,Alm,x],x]/.x->0;
		If[Chop[Dvalp,\[Epsilon]]==0,
			If[Chop[valp,\[Epsilon]]==0,
			(* This shouldn't happen, but could in the asymptotic regime where where exponential damping can occur *)
				Return[{Indeterminate,Indeterminate}];
			,(* Since the derivative vanishes, the value should be non-zero and the phase has been chosen to make it real at x=0, so we can specify its sign *)
				index=l+m;If[OddQ[index],--index];
				signs=Sign[phases {valm,valp}];
				phases*=signs Sign[(-1)^(index/2)];
			]
		,(* general case: Check the signs for derivative *)
			Dvalm=D[SWSFHeun[-1][s,m,c,Alm,x],x]/.x->0;
			phases=Exp[I(\[Pi]-{Arg[Dvalm],Arg[Dvalp]})];
			signs=Sign[phases {Dvalm,Dvalp}];
			phases*=signs Sign[Dsphericalval];
		]
	];
	phases
]/;IntegerQ[m+s]&&IntegerQ[l+s]


Options[SWSFHeunNorms]=Union[{Quiet->NIntegrate::precw},Options[NIntegrate]];


SWSFHeunNorms[s_/;IntegerQ[2s],m_/;IntegerQ[2m],c_?NumberQ,Alm_?NumberQ,opts:OptionsPattern[]]:= 
Module[{x,pintm,nintm,pinte,ninte,ratio,quiet=OptionValue[Quiet]},
	pintm=Quiet[NIntegrate[Abs[SWSFHeun[1][s,m,c,Alm,x]]^2,{x,-3/4,3/4},Evaluate[FilterRules[{opts},Options[NIntegrate]]]],Evaluate[quiet]];
	pinte=Quiet[NIntegrate[Abs[SWSFHeun[1][s,m,c,Alm,x]]^2,{x,3/4,1},Evaluate[FilterRules[{opts},Options[NIntegrate]]]],Evaluate[quiet]];
	nintm=Quiet[NIntegrate[Abs[SWSFHeun[-1][s,m,c,Alm,x]]^2,{x,-3/4,3/4},Evaluate[FilterRules[{opts},Options[NIntegrate]]]],Evaluate[quiet]];
	ninte=Quiet[NIntegrate[Abs[SWSFHeun[-1][s,m,c,Alm,x]]^2,{x,-1,-3/4},Evaluate[FilterRules[{opts},Options[NIntegrate]]]],Evaluate[quiet]];
	ratio=pintm/nintm;
	1/Sqrt[{ninte+nintm+pinte/ratio,ninte ratio+pintm+pinte}]
]/;IntegerQ[m+s]


SWSFHeun[norm_List,s_/;IntegerQ[2s],m_/;IntegerQ[2m],c_?NumberQ,Alm_?NumberQ,x_?NumberQ]:=
If[x<0,norm[[1]]SWSFHeun[-1][s,m,c,Alm,x],norm[[2]]SWSFHeun[1][s,m,c,Alm,x]]/;IntegerQ[m+s]


Options[SWSFHeunExpansionCoefs]=Union[{CoefDelta->10^-15,ChopDelta->10^-10,Quiet->NIntegrate::precw},Options[NIntegrate]];


SWSFHeunExpansionCoefs[s_/;IntegerQ[2s],m_/;IntegerQ[2m],l_/;IntegerQ[2l],c_?NumberQ,Alm_?NumberQ,opts:OptionsPattern[]]:=
Module[{norm,phase,ip,coefs={},stop=False,normsum=0,L,lmin,Npara=10,ind,keep,lastsmall=False,\[Epsilon]=OptionValue[CoefDelta],quiet=OptionValue[Quiet]},
	norm=SWSFHeunNorms[s,m,c,Alm,Evaluate[FilterRules[{opts},Options[SWSFHeunNorms]]]]; 
	phase=SWSFHeunPhases[s,m,l,c,Alm,Evaluate[FilterRules[{opts},Options[SWSFHeunPhases]]]];
	If[NumericQ[phase[[1]]],norm*=phase];(* SL phase computed *)
	L=0;lmin=Max[Abs[s],Abs[m]];
	While[Abs[normsum-1]>\[Epsilon] || !stop,
		ip=Parallelize[Table[Quiet[NIntegrate[Conjugate[SWSphericalS[s,m,lmin+Lt,x]]SWSFHeun[norm,s,m,c,Alm,x],{x,-1,1},Evaluate[FilterRules[{opts},Options[NIntegrate]]]],Evaluate[quiet]],{Lt,L,L+Max[l-lmin+3-L,Npara-1]}]];
		ind=FirstPosition[Reverse[ip],_?(Abs[#]>\[Epsilon]&)][[1]];
		keep=Min[Length[ip],If[IntegerQ[ind],Length[ip]-ind+3,ind=Npara+2;2]];
		If[lastsmall&&keep<=2,keep=1;stop=True];
		ip=Take[ip,keep];
		If[ind>=2,lastsmall=True];
		normsum+=Total[Abs[ip]^2];
		coefs = Join[coefs,ip];
		If[Npara>1&&ind>2,stop=True];
		L+=keep
	];
	If[!NumberQ[phase[[1]]],(* Indeterminant phase, impose Cook-Zalutskiy phase choice *)
		phase = Exp[-I Arg[coefs[[l-lmin+1]]]];
		coefs*=phase;
	];
	Return[coefs]
]/;IntegerQ[m+s]&&IntegerQ[l+s]


(* ::Section::Closed:: *)
(*End of SWSpheroidal Package*)


End[] (* `Private` *)


EndPackage[]
