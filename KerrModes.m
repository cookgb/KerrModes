(* ::Package:: *)

(* ::Title:: *)
(*Modes of Kerr*)


(* ::Section::Closed:: *)
(*Documentation *)


(* ::Item:: *)
(*Overview*)


(* ::ItemParagraph:: *)
(*The KerrModes Mathematica package provides functions to compute and plot the various modes of gravitational (|s|=2), electromagnetic (|s|=1), and scalar (s=0) perturnbations of the Kerr geometry.*)


(* ::ItemParagraph:: *)
(*Modes sequences for a specific (l,m) mode and overtone n are parameterized by the dimensionless spin parameter a \[Element][0,1).  Sequences usually start with a=0 using initial guesses stored in tables for the appropriate non-rotating (Schwarzschild) modes.*)


(* ::Item:: *)
(*Mode Type*)


(* ::ItemParagraph:: *)
(*This package is intended to be included in a "wrapper" package that supplies the definitions necessary to compute a specific type of mode: QNM, Subscript[TTM, L], Subscript[TTM, R]*)


(* ::Section::Closed:: *)
(*Begin KerrModes Package*)


BeginPackage["KerrModes`",{"SWSpheroidal`"}]


Unprotect[KerrModeDebug];
KerrModeDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[KerrModeDebug,Unprotect["KerrModes`*"];Unprotect["KerrModes`Private`*"]];
Protect[KerrModeDebug];


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[SpinWeight,ModePrecision,RadialCFDepth,RadialCFMinDepth,RadialDebug,RadialRelax,
        JacobianStep,Root\[Epsilon],SchDebug];


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Radial Equation : Modified Leaver' s Method*)


(* ::Subsection::Closed:: *)
(*Newton' s Method for finding roots of Radial Equation*)


Options[RadialLentzStep]={RadialDebug->0}


RadialLentzStep[n_Integer,s_Integer,m_Integer,
				a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ,
				\[Omega]step_Real|\[Omega]step_Rational|\[Omega]step_Integer,
				Nrcf_Integer,Nm_Integer,\[Epsilon]_Integer,OptionsPattern[]]:= 
Module[{\[Omega]r,\[Omega]i,Almr,sol,solpr,\[CapitalDelta]sol,pcount,pflag,\[CapitalDelta]\[Omega]r,\[CapitalDelta]\[Omega]i,RedRd\[Omega]r,RedRd\[Omega]i,ImdRd\[Omega]r,ImdRd\[Omega]i,\[Delta]\[Omega],WorkPrec,
		radialdebug=OptionValue[RadialDebug]},
	sol=RadialCFRem[n,s,m,a,Alm,\[Omega],Nrcf];
	\[CapitalDelta]\[Omega]r=Re[\[Omega]]\[Omega]step;\[CapitalDelta]\[Omega]i=Im[\[Omega]]\[Omega]step;
	\[Omega]r=Re[\[Omega]](1+\[Omega]step)+I Im[\[Omega]];
	\[Omega]i=Re[\[Omega]]+I Im[\[Omega]](1+\[Omega]step);
	Almr=AngularSpectralRoot[s,m,a*\[Omega]r,Alm,Nm][[1]];
	pcount=0;pflag=False;
	While[True,
		solpr=RadialCFRem[n,s,m,a,Almr,\[Omega]r,Nrcf];
		\[CapitalDelta]sol=solpr[[1]]-sol[[1]];
		Off[Precision::mnprec,Accuracy::mnprec];
		If[MyPrecision[\[CapitalDelta]sol]==Precision[\[CapitalDelta]sol],
			On[Precision::mnprec,Accuracy::mnprec];Break[],
			$MinPrecision+=4,
			$MinPrecision+=4
		];
		On[Precision::mnprec,Accuracy::mnprec];
		pflag=True;
		If[++pcount==10,Print["Excessive increase in Precision, Abort"];Abort[]];
		sol=RadialCFRem[n,s,m,a,Alm,\[Omega],Nrcf];
	];
	If[pflag,Print["Set $MinPrecision = ",$MinPrecision]];
	RedRd\[Omega]r=Re[\[CapitalDelta]sol]/\[CapitalDelta]\[Omega]r;
	ImdRd\[Omega]r=Im[\[CapitalDelta]sol]/\[CapitalDelta]\[Omega]r;
	ImdRd\[Omega]i=RedRd\[Omega]r;
	RedRd\[Omega]i=-ImdRd\[Omega]r;
	If[Abs[sol[[1]]]==0,
		Return[{Solve[{\[Delta]\[Omega]r==0,\[Delta]\[Omega]i==0},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol,{{RedRd\[Omega]r,RedRd\[Omega]i},{ImdRd\[Omega]r,ImdRd\[Omega]i}}}];
	]; (* Avoid error case *)
	If[radialdebug>4,
		Print[RedRd\[Omega]r,"\[Delta]\[Omega]r + ",RedRd\[Omega]i,"\[Delta]\[Omega]i==",-Re[sol[[1]]]];
		Print[ImdRd\[Omega]r,"\[Delta]\[Omega]r + ",ImdRd\[Omega]i,"\[Delta]\[Omega]i==",-Im[sol[[1]]]]
	];
	If[Accuracy[RedRd\[Omega]r]<=0,RedRd\[Omega]r=0];
	If[Accuracy[RedRd\[Omega]i]<=0,RedRd\[Omega]i=0];
	If[Accuracy[ImdRd\[Omega]r]<=0,ImdRd\[Omega]r=0];
	If[Accuracy[ImdRd\[Omega]i]<=0,ImdRd\[Omega]i=0];
	If[RedRd\[Omega]r==0 && RedRd\[Omega]i==0,
		Return[{Solve[{\[Delta]\[Omega]r==2 10^(\[Epsilon]+2),\[Delta]\[Omega]i==-3 10^(\[Epsilon]+2)},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol,Null}];
	]; (* Avoid error case *)
	If[ImdRd\[Omega]r==0 && ImdRd\[Omega]i==0,
		Return[{Solve[{\[Delta]\[Omega]r==-2 10^(\[Epsilon]+2),\[Delta]\[Omega]i==3 10^(\[Epsilon]+2)},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol,Null}];
	]; (* Avoid error case *)
	WorkPrec = 2Max[$MinPrecision,Quiet[Precision[sol[[1]]]]];
	\[Delta]\[Omega]=Solve[{RedRd\[Omega]r \[Delta]\[Omega]r + RedRd\[Omega]i \[Delta]\[Omega]i==-Re[sol[[1]]], ImdRd\[Omega]r \[Delta]\[Omega]r + ImdRd\[Omega]i \[Delta]\[Omega]i==-Im[sol[[1]]]},{\[Delta]\[Omega]r,\[Delta]\[Omega]i},WorkingPrecision->WorkPrec][[1]];
	If[radialdebug>3,Print["RadialLentzStep : ",\[Delta]\[Omega]," : ",sol]];
	If[radialdebug>5,Print["Precision - \[Delta]\[Omega] : ",Precision[\[Delta]\[Omega][[1]]]," " ,Precision[\[Delta]\[Omega][[2]]]]];
	{\[Delta]\[Omega],sol,{{RedRd\[Omega]r,RedRd\[Omega]i},{ImdRd\[Omega]r,ImdRd\[Omega]i}}}
]


Options[RadialLentzRoot]=Union[Options[RadialLentzStep],{RadialRelax->1,JacobianStep->-10,Root\[Epsilon]->Null[]}];


RadialLentzRoot[n_Integer,s_Integer,m_Integer,
				a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ,
				Nrcf_Integer,Nm_Integer,\[Epsilon]_Integer,
				Radius_Real|Radius_Rational|Radius_Integer,
				opts:OptionsPattern[]]:= 
Module[{sol1,\[Omega]root=\[Omega],\[Delta]\[Omega]1,\[Delta]\[Omega]2,Almc,Ninv,iteration=0,slow=False,convcount=0,
		refliptotal=0,reflipcount=0,lastrealsign=0,
		largeroot=0,largerootratio=0,largerootcount=0,largeroottotal=0,
		jacobianmatrix=Null,
		\[Omega]step=10^(OptionValue[JacobianStep]),\[Epsilon]root=OptionValue[Root\[Epsilon]],
		radialdebug=OptionValue[RadialDebug],
		radialrelax=Rationalize[OptionValue[RadialRelax]]},
	If[Not[NumberQ[\[Epsilon]root]],\[Epsilon]root=\[Epsilon]];
	lastrealsign=Sign[Re[\[Omega]]];
	Almc=AngularSpectralRoot[s,m,a*\[Omega]root,Alm,Nm][[1]];
	sol1=RadialLentzStep[n,s,m,a,Almc,\[Omega]root,\[Omega]step,Nrcf,Nm,\[Epsilon],FilterRules[{opts},Options[RadialLentzStep]]];
	\[Delta]\[Omega]1 = sol1[[1,1,2]]+I sol1[[1,2,2]];
	If[Not[NumberQ[\[Delta]\[Omega]1]],Print["RadialLentzStep failed, returning ",sol1];Abort[]];
	Ninv=n;
	If[radialdebug>0,Print["\[Delta]\[Omega]= ",\[Delta]\[Omega]1," root= ",Abs[sol1[[2,1]]]]];
	jacobianmatrix=sol1[[3]];
	While[Abs[\[Delta]\[Omega]1]>10^\[Epsilon] || Abs[sol1[[2,1]]]>10^\[Epsilon]root,
		If[++iteration > 40,
			Return[{{False,slow,jacobianmatrix},{\[Omega]root,Ninv,sol1[[2,2]],\[Epsilon],Abs[\[Delta]\[Omega]1]}}]
		];
		\[Delta]\[Omega]2=\[Delta]\[Omega]1;
		\[Delta]\[Omega]1=If[Abs[\[Delta]\[Omega]1]>Radius,Radius \[Delta]\[Omega]1/Abs[\[Delta]\[Omega]1],\[Delta]\[Omega]1];
		\[Omega]root+= radialrelax \[Delta]\[Omega]1; (* NOTE: UNDER RELAXATION *)
		If[Abs[sol1[[2,1]]]>10,
			largerootratio=largeroot;
			largeroot=Abs[sol1[[2,1]]];
			largerootratio/=largeroot;
			++largeroottotal;
			If[1/4<largerootratio<4,
				If[++largerootcount>5, (* root remains large and doesn't change much *)
					Return[{{False,True,jacobianmatrix}}];
				],
				If[largeroottotal>20, (* root is large too often *)
					Return[{{False,True,jacobianmatrix}}];
				];
				largerootcount=0;
			],
			largeroot=0;
			largerootcount=0;
		];
		If[Sign[Re[\[Omega]root]]!=lastrealsign,
			lastrealsign*=-1;
			++refliptotal;
			If[++reflipcount>5, (* real part of \[Omega] is repeatedly flipping sign *)
				Return[{{False,True,jacobianmatrix}}];
			],
			If[refliptotal>20,  (* real part of \[Omega] is flipping sign intermitently *)
				Return[{{False,True,jacobianmatrix}}];
			];
			reflipcount=0;
		];
		If[radialdebug>1,Print["\[Omega]= ",\[Omega]root]];
		Almc=AngularSpectralRoot[s,m,a*\[Omega]root,Alm,Nm][[1]];
		sol1=RadialLentzStep[Ninv,s,m,a,Almc,\[Omega]root,\[Omega]step,Nrcf,Nm,\[Epsilon],FilterRules[{opts},Options[RadialLentzStep]]];
		\[Delta]\[Omega]1 = sol1[[1,1,2]]+I sol1[[1,2,2]];
		If[Not[NumberQ[\[Delta]\[Omega]1]],Print["RadialLentzStep failed, returning ",sol1];Abort[]];
		If[radialdebug>0,Print["\[Delta]\[Omega]= ",\[Delta]\[Omega]1," root= ",Abs[sol1[[2,1]]]]];
		If[Abs[\[Delta]\[Omega]1]/Abs[\[Delta]\[Omega]2]>1/2,
		If[++convcount>5,slow=True],convcount=0,slow=False];
		If[radialdebug>2,Print["Conv.Rate: ",Abs[\[Delta]\[Omega]2]/Abs[\[Delta]\[Omega]1]," : ",slow]]; 
		If[Head[sol1[[3]]]==List,jacobianmatrix=sol1[[3]]];
	];
	\[Delta]\[Omega]1=If[Abs[\[Delta]\[Omega]1]>Radius,Radius \[Delta]\[Omega]1/Abs[\[Delta]\[Omega]1],\[Delta]\[Omega]1];
	\[Omega]root+= radialrelax \[Delta]\[Omega]1; (* NOTE: UNDER RELAXATION *)
	{{True,slow,jacobianmatrix},{\[Omega]root,Ninv,Nrcf,\[Epsilon],Abs[\[Delta]\[Omega]1]}}
]


If[!modeDebug,Protect[RadialLentzStep,RadialLentzRoot]];


(* ::Subsection::Closed:: *)
(*Evaluate nth inversion of the Radial Equation' s continued fraction equation*)


(* ::Subsubsection::Closed:: *)
(*"Bottom-up" evaluation at the Nmax element with remainder approximation*)


RadialCFRem[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:= 
Module[{Rem,i,func,t},
	If[n>Nmax,Print["inversion greater than CF depth"];Abort[]];
	Rem=If[Nmax>0,RadialCFRemainder[s,m,a,Alm,\[Omega],Nmax][[1]],-1];
(*Print["Start RadialCFRem, Remainder = ",Rem];*)
	func=Simplify[{\[Alpha]r[i,s,m,a,Alm,\[Omega]],\[Beta]r[i,s,m,a,Alm,\[Omega]],\[Gamma]r[i,s,m,a,Alm,\[Omega]]}];
	t=Table[func,{i,0,Nmax}];
	{Evaluate[t[[n+1,2]]+t[[n+1,3]]Fold[-#2[[1]]/(#2[[2]]+#2[[3]]#1)&,0,Take[t,n]]]+
		t[[n+1,1]]Fold[-#2[[3]]/(#2[[2]]+#2[[1]]#1)&,Rem,Reverse[Drop[t,n+1]]],Nmax}
]


Options[TestRadialCFConvergence]=Options[RadialLentzRoot];


TestRadialCFConvergence[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
						Alm_?NumberQ,\[Omega]_?NumberQ,Nrcf_Integer,jacobian_,\[Epsilon]_Integer,Nrcfmin_Integer,
						Nm_Integer,\[Alpha]_Real|\[Alpha]_Rational|\[Alpha]_Integer,opts:OptionsPattern[]]:= 
Module[{N1,N2,Rem,CFval,CFval1,CFval2,cfpow,newNrcf,saveNrcf,diff,diffh,diffl,sol,cfpowcut=-2},
	If[Head[jacobian]==List,Null[],Null[],
		Print["WARNING: Jacobian not set when Nradialnew needs computation!"];
		Return[{Max[Nrcfmin,Ceiling[3/2 Nrcf]],0,0,0,Null[]}]
	];
	N1=Ceiling[2*Nrcf/3];
	N2=Ceiling[3*Nrcf/2];
	Rem=RadialCFRemainder[s,m,a,Alm,\[Omega],Nrcf];
	CFval=RadialCFRem[n,s,m,a,Alm,\[Omega],Nrcf][[1]];
	CFval1=RadialCFRem[n,s,m,a,Alm,\[Omega],N1][[1]];
	CFval2=RadialCFRem[n,s,m,a,Alm,\[Omega],N2][[1]];
	diff=Abs[CFval-CFval2];
	diffh=Abs[CFval1-CFval2];
	If[diff==0,
(*Print["Debug 1: diff==0"];*)
		If[-\[Epsilon]<=IntegerPart[Accuracy[diff]-(Log10[Det[jacobian]]/2)],
			Return[{Max[Nrcfmin,Ceiling[1/2 Nrcf]],diff,CFval,Rem,Null[]}],
			Print["WARNING: \[CapitalDelta]CF=0 testing RCF Depth with Acc : 10^(-",
					Ceiling[Accuracy[diff]-(Log10[Det[jacobian]]/2)],")"];
			Return[{Nrcf,diff,CFval,Rem,Null[]}];
		];
	];
	cfpow=(Log10[diff]-Log10[diffh])/Log10[3/2];
(*Print["Debug 2: cfpow = ",cfpow];*)
	If[cfpow>-1/2,(* Untrusted Sloap *)
		Print["Debug 3: cfpow>-1/2 (diff = ",diff,", diffh = ",diffh];
		Return[{Max[Nrcfmin,Ceiling[3/2 Nrcf]],diff,CFval,Rem,Null[]}]
	];
	newNrcf=Max[Nrcfmin,Ceiling[Nrcf/10],Ceiling[Nrcf (Sqrt[Det[jacobian]]10^\[Epsilon]/diff)^(1/cfpow)]];
(*Print["Debug 4: newNrcf = ",newNrcf," jacobian : ",Sqrt[Det[jacobian]]," diff : ",diff];*)
	If[newNrcf>=2Nrcf/3 || (cfpow<cfpowcut && newNrcf>Max[Nrcfmin,Ceiling[Nrcf/10]]),
		(*Print["Debug 5"];*)
		Return[{newNrcf,diff,CFval,Rem,cfpow}]
	];
	newNrcf=Max[Nrcfmin,N1];
	sol=RadialLentzRoot[n,s,m,a,Alm,\[Omega],newNrcf,Nm,\[Epsilon],10^(-3),
							RadialRelax->\[Alpha],FilterRules[{opts},Options[RadialLentzRoot]]];
	If[sol[[1,1]] && !sol[[1,2]],
		If[Log10[Abs[\[Omega]-sol[[2,1]]]]<\[Epsilon],newNrcf=Max[Nrcfmin,N1],newNrcf=Nrcf,newNrcf=Nrcf],
		newNrcf=Nrcf,newNrcf=Nrcf;
	];
(*Print["Debug 9: sol = ",sol];*)
Print["Warning: resetting newNrcf = ",newNrcf];
	{newNrcf,diff,CFval,Rem,cfpow}
]


If[!modeDebug,Protect[RadialCF,TestRadialCFConvergence]];


(* ::Section::Closed:: *)
(*Kerr Modes methods*)


(* ::Subsection::Closed:: *)
(*Utility routines*)


MyPrecision[x_?NumberQ]:=Module[{saveprecision,returnprecision},
	saveprecision=$MinPrecision;
	$MinPrecision=0;
	returnprecision=Precision[x];
	$MinPrecision=saveprecision;
	Return[returnprecision];
]


(* ::Section:: *)
(*Initial Guesses*)


Options[SchGuess]={SpinWeight->Null[]};


SchGuess[l_Integer,n_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],guess},
	SchGuess::spinweight="SpinWeight has not been properly set!";
	SchGuess::argl="The order l is set to `1`, but must be \[GreaterEqual] |`2`|";
	SchGuess::argn="The overtone n is set to `1`, but must be \[GreaterEqual] 0";
	SchGuess::guess="Need initial guess for `1`[`2`,`3`]";
	SchGuess::count="Overtone count in `1`[`2`]=`3` is incorrect";
	If[OptionValue[SpinWeight]==Null[],Message[SchGuess::spinweight];Abort[]];
	If[l<Abs[s],Message[SchGuess::argl,l,s];Abort[]];
	If[n<0,Message[SchGuess::argn,n];Abort[]];
	If[Head[SchTable[Abs[s], 0]]==SchTable,
			Message[SchGuess::guess,SchTable,Abs[s],0];Abort[]];
	If[Head[SchTable[Abs[s], 1]]==SchTable,
			Message[SchGuess::guess,SchTable,Abs[s],1];Abort[]];
	If[Head[SchTable[Abs[s]]]==SchTable,
			Message[SchGuess::count,SchTable,Abs[s],Null];Abort[]];
	If[SchTable[Abs[s]]<2,
			Message[SchGuess::count,SchTable,Abs[s],SchTable[Abs[s]]];Abort[]];
	If[l>Abs[s],
		If[Head[SchTable[Abs[s]+1, 0]]==SchTable,
				Message[SchGuess::guess,SchTable,Abs[s]+1,0];Abort[]];
		If[Head[SchTable[Abs[s]+1, 1]]==SchTable,
				Message[SchGuess::guess,SchTable,Abs[s]+1,1];Abort[]];
		If[Head[SchTable[Abs[s]+1]]==SchTable,
				Message[SchGuess::count,SchTable,Abs[s]+1,Null];Abort[]];
		If[SchTable[Abs[s]+1]<2,
				Message[SchGuess::count,SchTable,Abs[s]+1,SchTable[Abs[s]+1]];Abort[]];
	];
	If[n<SchTable[l],
		(* Initial guess is in Table *)
		SchTable[l,n],
		If[n>SchTable[l],Message[SchGuess::count,SchTable,l,SchTable[l]];Abort[]];
		If[l>=Abs[s]+1 && n<2,
			guess=2 SchGuess[l-1,n,FilterRules[{opts},Options[SchGuess]]]-
							SchGuess[l-2,n,FilterRules[{opts},Options[SchGuess]]];
			guess[[1]]=Abs[Re[guess[[1]]]]+I Im[guess[[1]]];
			guess[[2]]=n; (* Use preferred inversion *)
			guess[[3]]=0;guess[[4]]=0;guess[[5]]=0;
			SchTable[l,n]=guess,
			guess=2 SchGuess[l,n-1,FilterRules[{opts},Options[SchGuess]]]-
							SchGuess[l,n-2,FilterRules[{opts},Options[SchGuess]]];
			guess[[1]]=Abs[Re[guess[[1]]]]+I Im[guess[[1]]];
			guess[[2]]=n; (* Use preferred inversion *)
			SchTable[l,n]=guess
		],
		Message[SchGuess::count,SchTable,l,SchTable[l]];
		Abort[];
	]
]


Options[SchwarzschildMode]=Union[{SpinWeight->Null[],SchDebug->0,
								RadialCFMinDepth->300,RadialCFDepth->1,ModePrecision->24},
								Options[TestRadialCFConvergence]];


SchwarzschildMode[l_Integer,n_Integer,opts:OptionsPattern[]] :=
Module[{s=OptionValue[SpinWeight],debug=OptionValue[SchDebug],
		rcfdepth=OptionValue[RadialCFDepth],
		Nrcf=300,testinv0,testinv1,Ninv,rcferr,rcfpow,nrpow,jacobianmatrix,Nradialnew,
		notconverged,guess,sol0,sol1,\[Epsilon]=-14,
		RCFmin=OptionValue[RadialCFMinDepth],
		precision=OptionValue[ModePrecision]},
	SchwarzschildMode::spinweight="SpinWeight has not been properly set!";
	SchwarzschildMode::argl="The order l is set to `1`, but must be \[GreaterEqual] |`2`|";
	SchwarzschildMode::argn="The overtone n is set to `1`, but must be \[GreaterEqual] 0";
	SchwarzschildMode::convergence="Solution failed to converge for (`1`,`2`)";
	If[OptionValue[SpinWeight]==Null[],Message[SchwarzschildMode::spinweight];Abort[]];
	If[l<Abs[s],Message[SchwarzschildMode::argl,l,s];Abort[]];
	If[n<0,Message[SchwarzschildMode::argn,n];Abort[]];
	$MinPrecision = precision;
	If[Head[SchTable[l]]==SchTable,SchTable[l]=0];
	If[SchTable[l]<n,SchwarzschildMode[l,n-1,FilterRules[{opts},Options[SchwarzschildMode]]]];
	Print["Computing (l=", l, ",n=", n, ")"];
	sol1=SchGuess[l, n,FilterRules[{opts},Options[SchGuess]]];
	Ninv = n;
	If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
	If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
	Nrcf=Max[Nrcf,RCFmin];
	notconverged = True;
	While[notconverged,
		sol1=RadialLentzRoot[Ninv,s,0,0,l(l+1)-s(s+1),SetPrecision[sol1[[1]],precision],Nrcf,l+2,\[Epsilon],10^(-3),FilterRules[{opts},Options[RadialLentzRoot]]];
		If[sol1[[1, 1]],
			jacobianmatrix=sol1[[1,3]];
			sol1=sol1[[2]];
			rcferr=TestRadialCFConvergence[Ninv,s,0,0,l(l+1)-s(s+1),SetPrecision[sol1[[1]],precision],Nrcf,jacobianmatrix,\[Epsilon],RCFmin,l+2,1,FilterRules[{opts},Options[TestRadialCFConvergence]]];
			Nradialnew=rcferr[[1]];
			If[debug>0,Print["RadialConverg : ",rcferr]];
			If[(Nradialnew>Nrcf),
				Nrcf=Nradialnew;
				If[debug>0,Print["Increase Nradial to ",Nrcf]];
				,
				notconverged=False;
			],
			Message[SchwarzschildMode::convergence,l,n];
			Abort[];
		];
	];
	SchTable[l,n]=sol1;
	SchTable[l]=Max[SchTable[l],n+1];
	Print[sol1];
]


(* ::Section::Closed:: *)
(*End of KerrModes Package*)


End[] (* `Private` *)


EndPackage[]
