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
(*Eigenvalue Solvers*)


RadialLentzRoot::usage=
	"RadialLentzRoot[n,s,m,a,Alm,\[Omega],N,\[Epsilon],Radius] solves the radial Teukolsky "<>
	"equation with spin-weight s, 'magnetic' index m, dimensionless angular "<>
	"momentum a, and separation constant Alm.  The solution is obtained by "<>
	"finding the root of the \!\(\*SuperscriptBox[\(n\), \(th\)]\) inversion "<>
	"of the associated continued fraction equation.  Newton's method is used "<>
	"and \[Omega] is taken as the initial guess.  The continued fraction is evalueated "<>
	"'bottom up' starting with the approximate remainder for the "<>
	"\!\(\*SuperscriptBox[\(N\), \(th\)]\) term.  Newton's method terminates "<>
	"when the corrections are smaller than \!\(\*SuperscriptBox[\(10\),\(\[Epsilon]\)]\).  "<>
	"Radius sets an upper limit to the magnitude of the Newton correction "<>
	"during each iteration.\n\n"<>
	"Options:\n"<>
	"\t RadialDebug\[Rule]0 : Integer\n"<>
	"\t\t Verbosity of debugging output.  Increasing value increases verbosity.\n"<>
	"\t RadialRelax\[Rule]1: Rational(or Integer)\n"<>
	"\t\t Under-relaxation parameter.\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t Log10 of relative step size used in numerical evaluation of the \n"<>
	"\t\t Jacobian in Newton's method."


ModeSolution::usage=
	"QNMSolution[n,s,l,m,a,\[Omega]g,Almg,\[Epsilon],relax,Nrcf,Nm,\[Omega]0,Alm0,rl,rt] finds a "<>
	"solution of the coupled radial and angular Teukolsky equations with "<>
	"spin-weight s, 'magnetic' index m, and dimensionless angular momentum a.  "<>
	"\[Omega]g and Almg are initial guesses for the frequency and separation constant, "<>
	"which are assumed to be associated with azimuthal index l and overtone "<>
	"index n.\n\n"<>
	"\tThe solution is found by iteration, calling AngularSpectralRoot and "<>
	"RadialLentzRoot until the magnitude of the change produced by either "<>
	"is less than \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  'relax' is an initial "<>
	"under-relaxation parameter used in updating the current guesses. For the "<>
	"radial equation, the continued fraction is truncated at the "<>
	"\!\(\*SuperscriptBox[\(Nrcf\),\(th\)]\) term.  For the angular equation, "<>
	"Nm sets the size of the spectral approximation matrix.\n\n"<>
	"\t\[Omega]0, Alm0, rl, and rt are used to specify solution windows centered "<>
	"around the initial guesses \[Omega]g and Almg.  \[Omega]0 and Alm0 represent the prior "<>
	"solutions in a sequence of solutions.  |\[Omega]0-\[Omega]g| and |Alm0-Almg| set a "<>
	"length scale d for each solution.  The window for each solution is a "<>
	"portion of an anulus centered on the prior solution.  The difference "<>
	"between the inner and outer rings of the anulus is 2*rl*d with the "<>
	"guess solution centered between.  The width of the wedge is 2*rt*d "<>
	"along the outer ring.  Solutions falling outside the solution window "<>
	"for either quantity are rejected.  If rl=0 or rt=0, the solution is "<>
	"always accepted.\n\n"<>
	"Options:\n"<>
	"\t SolutionDebug\[Rule]0 : Integer\n"<>
	"\t\t Verbosity of debugging output for solutions iteration.  Increasing\n"<>
	"\t\t value increases verbosity.\n"<>
	"\t RadialDebug\[Rule]0 : Integer\n"<>
	"\t\t Verbosity of debugging output.  Increasing value increases verbosity.\n"<>
	"\t RadialRelax\[Rule]1: Rational(or Integer)\n"<>
	"\t\t Under-relaxation parameter.\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t Log10 of relative step size used in numerical evaluation of the \n"<>
	"\t\t Jacobian in Newton's method."


(* ::Subsection::Closed:: *)
(*Plotting Routines*)


PlotModeFunction::usage=
"Arguments: n,s,m,a,Alm,\[Omega],Nrcf,Nm\n"<>
"\t n -> Overtone level (integer).\n"<>
"\t s -> SpinWeight (integer).\n"<>
"\t m -> Azimuthal index (integer).\n"<>
"\t a -> Dimensionless spin parameter(rational or integer).\n"<>
"\t Alm -> Angular separation constant (integer).\n"<>
"\t Nrcf -> How deep into the continued fraction (integer).\n"<>
"\t Nm -> Number of azimuthal indecies(integer).\n\n"<>
"When in ContinuedFractionMode, package will require input for all the arguments and will use SelectMode\n"<> 
"\t to replace ModeFunction with RadialCFRem.\n"<>
"When in PolynomialMode, n and Nrcf are not used in this plot function. Put dummy value in place.\n"<>
"PolynomialMode will use SelectMode to replace Modefunction with Starobinsky.\n"


PlotModeFunctionL::usage=
"Arguments: n,s,m,a,Alm,\[Omega],Nrcf,Nm\n"<>
"\t n -> Overtone level (integer).\n"<>
"\t s -> SpinWeight (integer).\n"<>
"\t m -> Azimuthal index (integer).\n"<>
"\t a -> Dimensionless spin parameter(rational or integer).\n"<>
"\t Alm -> Angular separation constant (integer).\n"<>
"\t L index -> Index into the Angular separation constants (integer).\n"<>
"\t Nrcf -> How deep into the continued fraction (integer).\n"<>
"\t Nm -> Number of azimuthal indecies(integer).\n\n"<>
"When in ContinuedFractionMode, package will require input for all the arguments and will use SelectMode\n"<> 
"\t to replace ModeFunction with RadialCFRem.\n"<>
"When in PolynomialMode, n and Nrcf are not used in this plot function. Put dummy value in place.\n"<>
"PolynomialMode will use SelectMode to replace Modefunction with Starobinsky.\n"


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[PolynomialMode,ContinuedFractionMode,RunCFConvergence];


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
	RadialLentzStep::increase="Excessive increase in Precision, Abort";
	RadialLentzStep::debug4="RadialLentzStep :`1` : `2`";
	RadialLentzStep::debug5Re="`1` \[Delta]\[Omega]r + `2` \[Delta]\[Omega]i==`3`";
	RadialLentzStep::debug5Im="`1` \[Delta]\[Omega]r + `2` \[Delta]\[Omega]i==`3`";
	RadialLentzStep::debug6="Precision - \[Delta]\[Omega] : `1`, `2`";
	RadialLentzStep::pflag="Set $MinPrecision =`1`";
	sol=ModeFunction[n,s,m,a,Alm,\[Omega],Nrcf];
	\[CapitalDelta]\[Omega]r=Re[\[Omega]]\[Omega]step;\[CapitalDelta]\[Omega]i=Im[\[Omega]]\[Omega]step;
	\[Omega]r=Re[\[Omega]](1+\[Omega]step)+I Im[\[Omega]];
	\[Omega]i=Re[\[Omega]]+I Im[\[Omega]](1+\[Omega]step);
	Almr=AngularSpectralRoot[s,m,a*\[Omega]r,Alm,Nm][[1]];
	pcount=0;pflag=False;
	While[True,
		solpr=ModeFunction[n,s,m,a,Almr,\[Omega]r,Nrcf];
		\[CapitalDelta]sol=solpr[[1]]-sol[[1]];
		Off[Precision::mnprec,Accuracy::mnprec];
		If[MyPrecision[\[CapitalDelta]sol]==Precision[\[CapitalDelta]sol],
			On[Precision::mnprec,Accuracy::mnprec];Break[],
			$MinPrecision+=4,
			$MinPrecision+=4
		];
		On[Precision::mnprec,Accuracy::mnprec];
		pflag=True;
		If[++pcount==10,Message[RadialLentzStep::increase];Abort[]];
		sol=ModeFunction[n,s,m,a,Alm,\[Omega],Nrcf];
	];
	If[pflag,Print[Style[StringForm[RadialLentzStep::pflag,$MinPrecision],{Medium,Darker[Blue]}]]];
	RedRd\[Omega]r=Re[\[CapitalDelta]sol]/\[CapitalDelta]\[Omega]r;
	ImdRd\[Omega]r=Im[\[CapitalDelta]sol]/\[CapitalDelta]\[Omega]r;
	ImdRd\[Omega]i=RedRd\[Omega]r;
	RedRd\[Omega]i=-ImdRd\[Omega]r;
	If[Abs[sol[[1]]]==0,
		Return[{Solve[{\[Delta]\[Omega]r==0,\[Delta]\[Omega]i==0},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol,{{RedRd\[Omega]r,RedRd\[Omega]i},{ImdRd\[Omega]r,ImdRd\[Omega]i}}}];
	]; (* Avoid error case *)
	If[radialdebug>4,
		Print[Style[StringForm[RadialLentzStep::debug5Re,RedRd\[Omega]r,RedRd\[Omega]i,-Re[sol[[1]]]],{Medium,Darker[Yellow]}]];
		Print[Style[StringForm[RadialLentzStep::debug5Im,ImdRd\[Omega]r,ImdRd\[Omega]i,-Im[sol[[1]]]],{Medium,Darker[Yellow]}]];
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
	If[radialdebug>3,Print[Style[StringForm[RadialLentzStep::debug4,\[Delta]\[Omega],sol],{Medium, Darker[Orange]}]]];
	If[radialdebug>5,Print[Style[StringForm[RadialLentzStep::debug6,Precision[\[Delta]\[Omega][[1]]],Precision[\[Delta]\[Omega][[2]]]],{Medium,Darker[Cyan]}]]];
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
	RadialLentzRoot::fail="RadialLentzStep failed, returning `1`";
	RadialLentzRoot::debug1="\[Delta]\[Omega]=`1`  root=`2` ";
	RadialLentzRoot::debug2="\[Omega]=`1` ";
	RadialLentzRoot::debug3="Conv.Rate: `1` : `2`";
	If[Not[NumberQ[\[Epsilon]root]],\[Epsilon]root=\[Epsilon]];
	lastrealsign=Sign[Re[\[Omega]]];
	Almc=AngularSpectralRoot[s,m,a*\[Omega]root,Alm,Nm][[1]];
	sol1=RadialLentzStep[n,s,m,a,Almc,\[Omega]root,\[Omega]step,Nrcf,Nm,\[Epsilon],FilterRules[{opts},Options[RadialLentzStep]]];
	\[Delta]\[Omega]1 = sol1[[1,1,2]]+I sol1[[1,2,2]];
	If[Not[NumberQ[\[Delta]\[Omega]1]],Message[RadialLentzStep::fail,sol1];Abort[]];
	Ninv=n;
	If[radialdebug>0,Print[Style[StringForm[RadialLentzRoot::debug1,\[Delta]\[Omega]1,Abs[sol1[[2,1]]]],{Medium, Darker[Green,0.5]}]]];
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
		(*If[radialdebug>1,Print["\[Omega]= ",\[Omega]root]];*)
		If[radialdebug>1,Print[Style[StringForm[RadialLentzRoot::debug2 ,\[Omega]root],{Medium,Darker[Blue,0.5]}]]];
		Almc=AngularSpectralRoot[s,m,a*\[Omega]root,Alm,Nm][[1]];
		sol1=RadialLentzStep[Ninv,s,m,a,Almc,\[Omega]root,\[Omega]step,Nrcf,Nm,\[Epsilon],FilterRules[{opts},Options[RadialLentzStep]]];
		\[Delta]\[Omega]1 = sol1[[1,1,2]]+I sol1[[1,2,2]];
		(*If[Not[NumberQ[\[Delta]\[Omega]1]],Print["RadialLentzStep failed, returning ",sol1];Abort[]];*)
		If[Not[NumberQ[\[Delta]\[Omega]1]],Message[RadialLentzStep::fail,sol1];Abort[]];
		(*If[radialdebug>0,Print["\[Delta]\[Omega]= ",\[Delta]\[Omega]1," root= ",Abs[sol1[[2,1]]]]];*)
		If[radialdebug>0,Print[Style[StringForm[RadialLentzRoot::debug1,\[Delta]\[Omega]1,Abs[sol1[[2,1]]]],{Medium, Darker[Green,0.5]}]]];
		If[Abs[\[Delta]\[Omega]1]/Abs[\[Delta]\[Omega]2]>1/2,
		If[++convcount>5,slow=True],convcount=0,slow=False];
		(*If[radialdebug>2,Print["Conv.Rate: ",Abs[\[Delta]\[Omega]2]/Abs[\[Delta]\[Omega]1]," : ",slow]]; *)
		If[radialdebug>2,Print[Style[StringForm[RadialLentzRoot::debug3,Abs[\[Delta]\[Omega]2]/Abs[\[Delta]\[Omega]1],slow],{Medium,Darker[Magenta,0.4]}]]]; 
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
	RadialCFRem::inversion="inversion greater than CF depth";
	If[n>Nmax,Message[RadialCFRem::inversion];Abort[]];
	Rem=If[Nmax>0,RadialCFRemainder[s,m,a,Alm,\[Omega],Nmax][[1]],-1];
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
	TestRadialCFConvergence::notset="WARNING: Jacobian not set when Nradialnew needs computation!";
	TestRadialCFConvergence::accuracyceiling="WARNING: \[CapitalDelta]CF=0 testing RCF Depth with Acc : 10^(-`1`) ";
	TestRadialCFConvergence::diff="WARNING: cfpow>-1/2 (diff =`1` , diffh =`2` )";
	If[Head[jacobian]==List,Null[],Null[],
		Message[TestRadialCFConvergence::notset];
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
		If[-\[Epsilon]<=IntegerPart[Accuracy[diff]-(Log10[Det[jacobian]]/2)],
			Return[{Max[Nrcfmin,Ceiling[1/2 Nrcf]],diff,CFval,Rem,Null[]}],
			Message[TestRadialCFConvergence::accuracyceiling,
					Ceiling[Accuracy[diff]-(Log10[Det[jacobian]]/2)],")"];
			Return[{Nrcf,diff,CFval,Rem,Null[]}];
		];
	];
	cfpow=(Log10[diff]-Log10[diffh])/Log10[3/2];
	If[cfpow>-1/2,(* Untrusted Slope *)
		Print[Style[StringForm[TestRadialCFConvergence::diff,diff,diffh],{Medium,Magenta}]];
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


(* ::Subsection::Closed:: *)
(*Starobinsky Constant*)


Starobinsky[s_Integer,m_Integer,a_Rational|a_Integer,
			Alm_?NumberQ,\[Omega]_?NumberQ]:= 
Module[{\[Lambda],starob},
		If [s <= 0,
			\[Lambda]=Alm + a^2 \[Omega]^2 -2 m a \[Omega],
			\[Lambda]=Alm + a^2 \[Omega]^2 -2 m a \[Omega]+2 s
		];
		If [Abs[s]==2,
			starob=\[Lambda]^2(\[Lambda]+2)^2+8 \[Lambda] a \[Omega](6(a \[Omega]+m)-5\[Lambda](a \[Omega]-m))+144\[Omega]^2(1+a^2(a \[Omega]-m)^2),
			If [Abs[s]==1,
				starob=\[Lambda]^2-4 a \[Omega](a \[Omega]-m),Null[],Null[]],
			Null[]
		];
	{starob,0}
]



(* ::Section::Closed:: *)
(*Kerr Modes methods*)


(* ::Subsection::Closed:: *)
(*Iterative simultaneous solution of radial & angular Teukolsky equations*)


Options[Set\[CapitalDelta]a]={Min\[CapitalDelta]alevel->1,Max\[CapitalDelta]alevel->4,\[Omega]poly->Null[],Max\[CapitalDelta]\[Phi]->0.0005,Max\[CapitalDelta]\[Omega]->0.01};


Options[ModeSolution]=Union[{SolutionDebug->0,NoNeg\[Omega]->False,RadialCFMinDepth->300,QNMPrecision->24,
							SolutionSlow->10,SolutionOscillate->10,SolutionIter->50,
							RadialCFDigits->8,RCFPower->Null[]},Options[RadialLentzRoot2]];


ModeSolution[n_Integer,s_Integer,l_Integer,m_Integer,
			a_Rational|a_Integer,\[Omega]g_?NumberQ,Almg_?NumberQ,\[Epsilon]_Integer,
			relax_Real|relax_Rational|relax_Integer,
			Nrcf_Integer,Nm_Integer,
			\[Omega]0_?NumberQ,Alm0_?NumberQ,
			rl_Real|rl_Rational|rl_Integer,
			rt_Real|rt_Rational|rt_Integer,
			opts:OptionsPattern[]]:=
Module[{c,old\[Omega],oldAlm,radialsol,angularsol,lmin,lmax,Nradial,Nmatrix,
		\[CapitalDelta]\[Omega],\[CapitalDelta]\[Omega]2,iteration=0,nrpow,rcfpow,NradFlag=False,jacobianmatrix=Null,invJacobian,
		converged=False,count,expconv,rcferr,\[Alpha],inversion,i,invcount,slowcount,
		radialfail=0,slowcount2=0,oscillate=0,\[Epsilon]2=\[Epsilon],Nradialnew,rcfpower=0,testNrcf,
		nonegfreq=OptionValue[NoNeg\[Omega]],
		slowval=OptionValue[SolutionSlow],oscval=OptionValue[SolutionOscillate],
		iterval=OptionValue[SolutionIter],jacobianstep=OptionValue[JacobianStep],
		solutiondebug=OptionValue[SolutionDebug],RCFmin=OptionValue[RadialCFMinDepth],
		precision=OptionValue[QNMPrecision],RCFdigits=OptionValue[RadialCFDigits],
		RCFpowoverride=OptionValue[RCFPower]},

	lmin = Max[Abs[m],Abs[s]];
	lmax = Max[l+Ceiling[Nm/2],lmin+Nm-1];
	Nmatrix=lmax-lmin+1;
	Nradialnew=Nradial=Nrcf;
	If[solutiondebug>1,Print["Initial Nradial : ",Nradial]];
	old\[Omega]=\[Omega]g;oldAlm=Almg;
	inversion=n;
	count=0;
	\[Alpha]=relax; (* under-relaxation parameter *)
	If[solutiondebug>0,Print["\[Omega]g = ",\[Omega]g," : Almg = ",Almg]];
	If[solutiondebug>3,Print["\[Omega]0 = ",\[Omega]0," : Alm0 = ",Alm0]];
	c=a \[Omega]g;
	angularsol=AngularSpectralRoot[s,m,c,oldAlm,Nmatrix];
	expconv=Max[Take[Abs[angularsol[[3]]],-2]];
	While[expconv >= 10^\[Epsilon]2 && ++count<5,  (* Make sure Spectral resolution is good before doing lots of work *)
		angularsol=AngularSpectralRoot[s,m,c,oldAlm,++Nmatrix];
		If[solutiondebug>4,Print["Expconv = ",expconv," : N = ",Nmatrix]];
		expconv=Max[Take[Abs[angularsol[[3]]],-2]];
	];
	count=0;
	While[Not[converged],
		If[(++iteration>iterval && \[Epsilon]2==\[Epsilon] && \[Alpha]==1) || (\[Alpha]<1 && iteration > (3/5)iterval/\[Alpha]),
			If[\[Alpha]>0.05,
				\[Alpha]*=3/5;iteration=0;
				Print["Reduce Under-relaxation parameter to \[Alpha] = ",N[\[Alpha]]],
				Print["Too Many Iterations"];
				(*Print["a=",a," \[Omega]=",old\[Omega]," Alm=",oldAlm]*);
				Return[{False}];
			];
		];
		radialsol = RadialLentzRoot[inversion,s,m,a,oldAlm,old\[Omega],Nradial,Nmatrix,\[Epsilon]2,10^(-3),
										RadialRelax->\[Alpha],FilterRules[{opts},Options[RadialLentzRoot]]];
		If[Not[radialsol[[1,1]]] && Length[radialsol]==1, (* Re[\[Omega]] flipping sign, solution fails *)
			Return[{False,\[Alpha],Nradialnew,{a,Null[],Null[]}}]
		];
		invcount=0;
		slowcount=0;
		If[Head[radialsol[[1,3]]]==List,jacobianmatrix=radialsol[[1,3]]];
		While[Not[radialsol[[1,1]]]&&Not[radialsol[[1,2]]], 
			(* Radial solution not converged, but convergence not slow *)
			If[++slowcount>=20,Break[]];
			(* Print["Looping with non-slow radial convergence"]; *)
			radialsol = RadialLentzRoot[inversion,s,m,a,oldAlm,radialsol[[2,1]],Nradial,Nmatrix,\[Epsilon]2,10^(-3),
											RadialRelax->\[Alpha],FilterRules[{opts}, Options[RadialLentzRoot]]];
			If[Not[radialsol[[1,1]]] && Length[radialsol]==1, (* Re[\[Omega]] flipping sign, solution fails *)
				Return[{False,\[Alpha],Nradialnew,{a,Null[],Null[]}}]
			];
			If[Head[radialsol[[1,3]]]==List,jacobianmatrix=radialsol[[1,3]]];
		];
		If[Not[radialsol[[1,1]]],++radialfail;(*Print["WARNING: RadialLentzRoot failed to converged"]*)];
		If[Not[NumberQ[radialsol[[2,1]]]],Print["Failure: Solution of RadialLentzRoot not a number"];Abort[]];
		c=a radialsol[[2,1]];
		angularsol=AngularSpectralRoot[s,m,c,oldAlm,Nmatrix];
		\[CapitalDelta]\[Omega]2=\[CapitalDelta]\[Omega];
		\[CapitalDelta]\[Omega]=radialsol[[2,5]];
		If[TrueQ[nonegfreq] && Re[radialsol[[2,1]]]<0,
				radialsol[[2,1]]=-Conjugate[radialsol[[2,1]]];
				angularsol[[1]]=Conjugate[angularsol[[1]]];
		];
		old\[Omega]=radialsol[[2,1]];oldAlm=angularsol[[1]];
		If[solutiondebug>2,Print["a = ",a," \[Omega] = ",old\[Omega]," alm = ",oldAlm]];
		If[solutiondebug>1||NumberQ[RCFpowoverride],Print["\[CapitalDelta]\[Omega] = ",\[CapitalDelta]\[Omega]]];
		If[\[CapitalDelta]\[Omega]>0,
			If[Abs[\[CapitalDelta]\[Omega]2]/Abs[\[CapitalDelta]\[Omega]]<1,
				If[++slowcount2>slowval,
					(* Try under-relaxation if solution is slow *)
					slowcount2=0;radialfail=0;
					\[Alpha]*=3/5;iteration=0;
					Print["Persistent slow convergence: \[Alpha]= ",N[\[Alpha]]];
					If[\[Alpha]<0.05,Return[{False}]]
				], 
				If[slowcount2 >0,slowcount2=0;
					(* Try lower precision if persisitent oscillations *)
					If[++oscillate>oscval,
						Print["Persistent oscillations: Abort"];Abort[]
						(*
						Print["Persistent oscillations: \[Epsilon]= ",\[Epsilon]2+1];
						oscillate=0;radialfail=0;
						If[\[Epsilon]2-\[Epsilon]<2,++\[Epsilon]2;iteration/=2,Return[{False}]]
						*)
					]
				]
			]
		];
		If[radialfail>10,
			(*
			radialfail=0;slowcount2=0;oscillate=0;
			Print["Persistent failure of RadialLentzRoot: \[Epsilon]= ",\[Epsilon]2+1];
			If[\[Epsilon]2-\[Epsilon]<2,++\[Epsilon]2;iteration/=2,Return[{False}]]
			*)
			Print["Persistent failure of RadialLentzRoot: Abort"];Abort[]
		];
		If[\[CapitalDelta]\[Omega]<10^\[Epsilon]2,
			If[count==0,
				rcfpower=0;
				rcferr=TestRadialCFConvergence[inversion,s,m,a,oldAlm(1+10^(\[Epsilon]+4)),old\[Omega](1-10^(\[Epsilon]+4)),
						Nradial,jacobianmatrix,\[Epsilon]2,RCFmin,Nmatrix,\[Alpha],FilterRules[{opts},Options[TestRadialCFConvergence]]];
				Nradialnew=rcferr[[1]];rcfpower=rcferr[[5]];
				If[solutiondebug>5,
					Print["Jacobian matrix : ",jacobianmatrix];
					Print["rcferr : ",rcferr];
				];
				If[solutiondebug>1,Print["Nradial : ",Nradial," ; Nradialnew : ",Nradialnew]];
				If[Nradialnew>Nradial,
					If[Nradialnew>(11/10)Nradial,
						NradFlag=True;count=0;iteration/=2;slowcount2=0;oscillate=0;
						If[solutiondebug>4,Print["Increase Nradial to ",Nradialnew]];
					];
					Nradial=Nradialnew;
					,
					If[solutiondebug>4 && Nradialnew<Nradial,Print["Decrease Nradial to ",Nradialnew]];
				];
			];

			expconv=Max[Take[Abs[angularsol[[3]]],-2]];
			(*expconv=0;*)
			If[expconv < 10^\[Epsilon]2,
				If[!NradFlag,++count];
				If[count>1&&radialsol[[1,1]],converged=True],
				If[solutiondebug>4,Print["Expconv = ",expconv," : count = ",count," : N = ",Nmatrix+1]];
				Nmatrix+=1;count=0;iteration/=2;slowcount2=0;oscillate=0
			];
			NradFlag=False,
			count=0
		];
	];
	If[rl!=0 &&rt!=0,Options[ModeSolution]=Union[{SolutionDebug->0,NoNeg\[Omega]->False,RadialCFMinDepth->300,QNMPrecision->24,
							SolutionSlow->10,SolutionOscillate->10,SolutionIter->50,
							RadialCFDigits->8,RCFPower->Null[]},Options[RadialLentzRoot2]];
		If[Not[SolutionWindow[\[Omega]0,\[Omega]g,old\[Omega],rl,rt,True]],
			Print["\[Omega] = ",old\[Omega]," outside solution window"];
			Return[{False,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix]}],angularsol}}]
		];
		If[Not[SolutionWindow[Alm0,Almg,oldAlm,rl,rt,True]],
			Print["Alm outside solution window"];
			Return[{False,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix]}],angularsol}}]
		];
	];
(*If[Nradialnew==Nrcf,Print["WARNING: Nradialnew not reset from Nrcf"]];*)
	{True,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix],$MinPrecision}],angularsol}}
]


(* ::Subsection::Closed:: *)
(*Stepsize and solution validation routines*)


SolutionWindow[F0_?NumberQ,Fg_NumberQ,Fs_NumberQ,
				rl_Real|rl_Rational|rl_Integer,
				rt_Real|rt_Rational|rt_Integer,
				plot_/;plot \[Element] Booleans]:=
Module[{Fn,Ff,Fp,Fm,Fsv,tp,tm,rad},
(*
F0: Prior solution;
Fg: Guessed solution;
Fs: Numerical solution, should be close to guess;
rl: Fraction of distance between F0 and Fg for length of error wedge (longitudinal ratio);
rt: Fraction of distance between F0 and Fg for width of error wedge (transverse ratio);
*)
	If[rt==0||rl==0,Return[True]];
	Fn=Fg-rl(Fg-F0);Ff=Fg+rl(Fg-F0);
	Fp=Ff-F0+I rt(Fg-F0);Fm=Ff-F0-I rt(Fg-F0);
	Fsv=Fs-F0;Fn=Abs[Fn-F0];Ff=Abs[Ff-F0];
	tp=Re[Fsv]Im[Fp]-Im[Fsv]Re[Fp];
	tm=Re[Fsv]Im[Fm]-Im[Fsv]Re[Fm];
	Fsv=Abs[Fsv];
	If[Fsv<Fn,(* Near fail *)Null,
	If[Fsv>Ff,(* Far fail *)Null,
	If[tp<0,(* Lateral fail + *)Null,
	If[tm>0,(* Lateral fail - *)Null,
	Return[True]]]]];
	If[plot,
		rad=Abs[Fg-F0]/50;
		Print[Show[Graphics[{Red,Circle[{Re[F0],Im[F0]},rad]}],Graphics[{Blue,Circle[{Re[Fg],Im[Fg]},rad]}],Graphics[{Black,Circle[{Re[Fp+F0],Im[Fp+F0]},rad]}],Graphics[{Green,Circle[{Re[Fm+F0],Im[Fm+F0]},rad]}],Graphics[{Black,Disk[{Re[Fs],Im[Fs]},rad]}]]]
	];
	Return[False];
]


(* ::Subsection::Closed:: *)
(*Utility routines*)


MyPrecision[x_?NumberQ]:=Module[{saveprecision,returnprecision},
	saveprecision=$MinPrecision;
	$MinPrecision=0;
	returnprecision=Precision[x];
	$MinPrecision=saveprecision;
	Return[returnprecision];
]


(* ::Section::Closed:: *)
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
	SchwarzschildMode::debug0a="RadialConverg : `1`";
	SchwarzschildMode::debug0b="Increase Nradial to `1`";
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
			If[RunCFConvergence,
				rcferr=TestRadialCFConvergence[Ninv,s,0,0,l(l+1)-s(s+1),SetPrecision[sol1[[1]],precision],Nrcf,jacobianmatrix,\[Epsilon],RCFmin,l+2,1,FilterRules[{opts},Options[TestRadialCFConvergence]]];
				Nradialnew=rcferr[[1]];
				If[debug>0,Print[Style[StringForm[SchwarzschildMode::debug0a,rcferr],{Medium,Purple}]]];
				If[(Nradialnew>Nrcf),
					Nrcf=Nradialnew;
					If[debug>0,Print[Style[StringForm[SchwarzschildMode::debug0b,Nrcf],{Medium,Darker[Purple]}]]];
					,
					notconverged=False;
				];
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
(*Graphics*)


PlotModeFunction[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			Alm_?NumberQ,\[Omega]_?NumberQ,Nrcf_Integer,Nm_Integer]:= 
Module[{Alm\[Omega]},
	Alm\[Omega]=AngularSpectralRoot[s,m,a*\[Omega],Alm,Nm][[1]];
     ModeFunction[n,s,m,a,Alm\[Omega],\[Omega],Nrcf][[1]]
]


PlotModeFunctionL[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			lindex_Integer,\[Omega]_?NumberQ,Nrcf_Integer,Nm_Integer]:= 
Module[{Alm\[Omega]},
	Alm\[Omega]=AngularSpectralRootIndex[s,m,a*\[Omega],lindex,Nm][[1]];
	ModeFunction[n,s,m,a,Alm\[Omega],\[Omega],Nrcf][[1]]
]


Options[PlotSchModes]=Union[{PlotTable->Null[]},Options[ListPlot]];


PlotSchModes[l_Integer,opts:OptionsPattern[]]:=
Module[{ptable=OptionValue[PlotTable]},
	PlotSchModes::plottable="Invalid PlotTable : `1`";
	If[ptable==Null[],Message[PlotSchModes::plottable,Null[]];Abort[]];
	plist=Table[{Re[ptable[l,n][[1]]],-Im[ptable[l,n][[1]]]},{n,0,ptable[l]-1}];
	mlist=Table[{-Re[ptable[l,n][[1]]],-Im[ptable[l,n][[1]]]},{n,0,ptable[l]-1}];
	ListLinePlot[{plist,mlist},PlotMarkers->Automatic,FilterRules[{opts},Options[ListLinePlot]]]
]


(* ::Section::Closed:: *)
(*End of KerrModes Package*)


End[] (* `Private` *)


EndPackage[]
