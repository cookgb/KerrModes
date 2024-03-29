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
"the spin-weighted spheroidal functions with spin-weight "
"\!\(\*StyleBox[\"s\", \"TI\"]\), azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), and oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\)."


AngularSpectralRoot::usage=
"AngularSpectralRoot[s,m,c,\!\(\*SubscriptBox[\(\[Lambda]\), \(0\)]\),N] "<>
"gives a list {\!\(\*StyleBox[\"\[Lambda]\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"N\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"vector\",\nFontSlant->\"Italic\"]\),\!\(\*StyleBox[\"index\",\nFontSlant->\"Italic\"]\)} from solving the "<>
"\!\(\*StyleBox[\"N\", \"TI\"]\)-dimensional approximate discrete eigensystem for "<>
"the spin-weighted spheroidal functions with spin-weight "<>
"\!\(\*StyleBox[\"s\", \"TI\"]\), azimuthal index \!\(\*StyleBox[\"m\", \"TI\"]\), and oblateness parameter \!\(\*StyleBox[\"c\", \"TI\"]\). The "<>
"solution returns the eigenvalue closest to \!\(\*StyleBox[SubscriptBox[StyleBox[\"\[Lambda]\", \"TI\"], \"0\"], \"TI\"]\), and the "<>
"corresponding spectral coefficients.  "<>
"\!\(\*
StyleBox[\"index\", \"TI\"]\) denotes the position of the eigensolution in full list."


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
(*Reserved Globals*)


Protect[PhaseChoice,SphericalLimit,Simple];


Protect[PathStart,PlotStart,PrintPoleValues,StepSize,\[Phi]guess]


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


Options[SWSFfixphase]={PhaseChoice->SphericalLimit};


SWSFfixphase[m_/;IntegerQ[2m],s_/;IntegerQ[2s],La_,SWdat_List,opts:OptionsPattern[]]:=
Module[{NC,lmin,SphericalVal,SphericalD,WDzero,WDzeroplus,WDzerominus,WDzeroD,scaledcoefs,SWSFzero,SWSFzeroD,phase,hoint=1},
	SWSFfixphase::badmethod="`1` is an invalid Method";
	If[!IntegerQ[m],hoint=-I];
	NC=Length[SWdat];
	lmin=Max[Abs[m],Abs[s]];
	Switch[OptionValue[PhaseChoice],
	Simple,
		phase=Exp[I(-Arg[SWdat[[La+1]]])],
	SphericalLimit,
		SphericalVal=hoint (-1)^m*Sqrt[La+lmin+1/2]*WignerD[{La+lmin,m,-s},0,\[Pi]/2,0];
		scaledcoefs=hoint (-1)^m ParallelTable[Sqrt[j-1+lmin+1/2]SWdat[[j]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
		WDzero=ParallelTable[N[WignerD[{j-1+lmin,m,-s},0,\[Pi]/2,0]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
		SWSFzero = WDzero . scaledcoefs;
		If[Sign[SphericalVal]!=0,
			phase=Exp[I(\[Pi]-Arg[SWSFzero])];
			If[Sign[Re[phase SWSFzero]]!=Sign[SphericalVal],phase=-phase]
		,(*Case when spherical value is zero.*)
			If[Chop[SWSFzero]==0&&s==0,(* special case where symmetry causes problems *)
				phase=Exp[I(-Arg[SWdat[[La+1]]])];
			, (* general case *)
				phase = Exp[I(\[Pi]-Arg[SWSFzero])];
			];
			WDzeroplus=ParallelTable[Sqrt[(j-1+lmin-s)*(j-1+lmin+s+1)]*If[s>=(j-1+lmin),0,N[WignerD[{j-1+lmin,m,-s-1},0,\[Pi]/2,0]]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
			WDzerominus=ParallelTable[Sqrt[(j-1+lmin+s)*(j-1+lmin-s+1)]*If[-s>=(j-1+lmin),0,N[WignerD[{j-1+lmin,m,-s+1},0,\[Pi]/2,0]]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
			WDzeroD=(-1/2)*(WDzeroplus-WDzerominus);
			SphericalD=hoint (-1)^m*Sqrt[La+lmin+1/2]*WDzeroD[[La+1]];
			SWSFzeroD=phase*(scaledcoefs . WDzeroD);
			If[Sign[Re[SWSFzeroD]]!=Sign[SphericalD],phase=-phase];
			If[Sign[SphericalD]==0,Print["Value and derivative both zero."];Abort[]];
		],
	_,Message[SWSFfixphase::badmethod,OptionValue[Method]];Abort[] 
	];
	phase
]/;IntegerQ[m+s]
SWSFfixphase[m_,s_,La_,SWdat_List,opts:OptionsPattern[]]:=Print["m and s must both be either integer or half-integer values."]


Options[SWSFvalues]={PlotPoints->100};


SWSFvalues[m_/;IntegerQ[2m],s_/;IntegerQ[2s],SWdat_List,opts:OptionsPattern[]]:=
Module[{NC,x,theta,Ntheta,lmin,Matdlx,SWSF,choplev,hoint=1,npoints=OptionValue[PlotPoints]},
	If[!IntegerQ[m],hoint=-I];
	NC=Length[SWdat];
	x=Table[x,{x,-1,1,1/(npoints-1)}];
	theta=ArcCos[#]&/@x;
	Ntheta=Length[theta];
	lmin=Max[Abs[m],Abs[s]];
	Matdlx = ParallelTable[N[WignerD[{j-1+lmin,m,-s},0,theta[[k]],0]],{k,1,Ntheta},{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	SWSF = hoint (-1)^m Matdlx . ParallelTable[Sqrt[j-1+lmin+1/2]SWdat[[j]],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	{x,SWSF}
]


Options[SWSFRealPath]=Union[{PathStart->Automatic,MaxSteps->5000,StepSize->1/50,PlotStart->False,\[Phi]guess->Automatic,PrintPoleValues->False}]


SWSFRealPath[m_/;IntegerQ[2m],s_/;IntegerQ[2s],SWdat_List,opts:OptionsPattern[]]:=
Module[{NC,lmin,Matdlx,SWSF,z,z0,zf,\[Delta]z,\[Phi],\[Phi]g,\[Phi]0,phase,zlist,SWSFlist,count=0,hoint=1,
		max\[Delta]r=OptionValue[StepSize],maxcount=OptionValue[MaxSteps],
		plot=OptionValue[PlotStart]},
	SWSFRealPath::PathStart="Invalid PathStart : `1`";
	If[!IntegerQ[m],hoint=-I];
	NC=Length[SWdat];
	lmin=Max[Abs[m],Abs[s]];
	Matdlx = ParallelTable[WignerD[{j-1+lmin,m,-s},0,ArcCos[z],0],{j,1,NC},DistributedContexts->{"SWSpheroidal`Private`"}];
	SWSF = hoint (-1)^m Matdlx . Table[Sqrt[j-1+lmin+1/2]SWdat[[j]],{j,1,NC}];
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
(*End of SWSpheroidal Package*)


End[] (* `Private` *)


EndPackage[]
