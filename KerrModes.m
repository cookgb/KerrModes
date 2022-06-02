(* ::Package:: *)

(* ::Title:: *)
(*Modes of Kerr*)


(* ::Section::Closed:: *)
(*Documentation*)


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
	"ModeSolution[n,s,l,m,a,\[Omega]g,Almg,\[Epsilon],relax,Nrcf,Nm,\[Omega]0,Alm0,rl,rt] finds a "<>
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
(*Utility Routines*)


KerrModeMakeMultiplet::usage=
	"KerrModeMakeMultiplet[l,m,n] converts an Integer overtone index 'n' for the "<>
	"existing Kerr QNM sequence with mode (l,m) into an overtone multiplet index.  "<>
	"By default, the sequence number is set to 0.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  An "<>
	"overtone multiplet index is a 2 element list {n,mult}, where 'n' is the Integer "<>
	"overtone number, and the sequence number 'mult' is in the range 0,1,...,(Nmult-1), "<>
	"with 'Nmult' the number of sequences with the same overtone index.\n\n"<>
	"Options:\n"<>
	"\t OTmultiple\[Rule]0\n"<>
	"\t\t Sequence number to be used in creating the overtone multiplet index."


ShortenModeSequence::usage=
	"ShortenModeSequence[l,m,n,N] removes the first N elements of the Mode sequence (l,m,n)\n"<>
	"if N>0, and removes the last N elements if N<0.\n\n"<>
	"Options:\n"<>
	"\t ShortenBy\[Rule]Drop : Drop, Take\n"<>
	"\t\t If 'Take' is chosen, the first N elemets are kept if N>0, and the last N elements\n"<>
	"\t\t are kept if N<0.\n"


(* ::Subsection::Closed:: *)
(*Plotting Routines*)


SchwarzschildOmega::usage=
	"SchwarzschildOmega[l,m,n] returns the conjucate of the Schwarzschild mode "<>
	"frequency for mode l and overtone n."


ModePlotOmega::usage=
"ModePlotOmega[l,m,n] or ModePlotOmega[l,n]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n\n"<>
"ModePlotOmega[l,m,n] plots the single specified sequence as -Im[\[Omega]] vs Re[\[Omega]].  "<>
"The overtone index can be either an integer or an overtone multiplet.\n"<>
"ModePlotOmega[l,n] plots all sequences in memory (all values of m) "<>
"as -Im[\[Omega]] vs Re[\[Omega]].  The overtone index must be an integer.  If overtone "<>
"multiplets exist, they must be specified via the OTmultiplet option.\n"<>
"Both plot the sequence and place markers at interals of 0.05 in a.\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden.\n"<>
"\t OTmultiplet \[Rule] Null[] : List of {m,\!\(\*SubscriptBox[\(n\), \(OT\)]\)} pairs, where m designates an "<>
"aximultal index having \n"<>
"\t\t\t an overtone multiplet, and \!\(\*SubscriptBox[\(n\), \(OT\)]\) is the number of multiplets.\n"<>
"\t All options available to ListPlot or ListLinePlot except for PlotMarkers."


ModePlotA::usage=
"ModePlotA[l,m,n] or ModePlotA[l,n]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n\n"<>
"ModePlotA[l,m,n] plots the single specified sequence as Im[\!\(\*SubscriptBox[\(A\), \(lm\)]\)] vs Re[\!\(\*SubscriptBox[\(A\), \(lm\)]\)].  "<>
"The overtone index can be either an integer or an overtone multiplet.\n"<>
"ModePlotOmega[l,n] plots all sequences in memory (all values of m) "<>
"as Im[\!\(\*SubscriptBox[\(A\), \(lm\)]\)] vs Re[\!\(\*SubscriptBox[\(A\), \(lm\)]\)].  The overtone index must be an integer.  If overtone "<>
"multiplets exist, they must be specified via the OTmultiplet option.\n"<>
"Both plot the sequence and place markers at interals of 0.05 in a.\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden.\n"<>
"\t OTmultiplet \[Rule] Null[] : List of {m,\!\(\*SubscriptBox[\(n\), \(OT\)]\)} pairs, where m designates an "<>
"aximultal index having \n"<>
"\t\t\t an overtone multiplet, and \!\(\*SubscriptBox[\(n\), \(OT\)]\) is the number of multiplets.\n"<>
"\t All options available to ListPlot or ListLinePlot except for PlotMarkers."


KerrOmegaListS::usage=
"KerrOmegaListS[l,m,n]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n\n"<>
"KerrOmegaListS creates a short list of {Re[\[Omega]],-Im[\[Omega]]} pairs for "<>
"modes along the specified sequence. Only modes where a is a multiple "<>
"of 0.05 are included.  (For a~1, the mode closest to a=1 is included.)\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden."


KerrOmegaList::usage=
"KerrOmegaList[l,m,n]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n\n"<>
"KerrOmegaList creates a list of {Re[\[Omega]],-Im[\[Omega]]} pairs for "<>
"modes along the specified sequence.\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden."


KerrAListS::usage=
"KerrOmegaListS[l,m,n]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n\n"<>
"KerrAListS creates a short list of {Re[\!\(\*SubscriptBox[\(A\), \(lm\)]\)],Im[\!\(\*SubscriptBox[\(A\), \(lm\)]\)]} pairs for "<>
"modes along the specified sequence. Only modes where a is a multiple "<>
"of 0.05 are included.  (For a~1, the mode closest to a=1 is included.)\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden."


KerrAList::usage=
"KerrOmegaList[l,m,n]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n\n"<>
"KerrAList creates a list of {Re[\!\(\*SubscriptBox[\(A\), \(lm\)]\)],Im[\!\(\*SubscriptBox[\(A\), \(lm\)]\)]} pairs for "<>
"modes along the specified sequence.\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden."


KerraOmegaListS::usage=
"KerraOmegaListS[l,m,n,ReIm]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n"<>
"\t ReIm : Re, Im, or Abs\n"<>
"KerraOmegaList creates a short list of either {a,Re[\[Omega]]}, {a,-Im[\[Omega]]}, or "<>
"{a,Abs[\[Omega]]} pairs for "<>
"modes along the specified sequence.  Which sequence is produced is governed "<>
"by the choice of ReIm. Only modes where a is a multiple "<>
"of 0.05 are included.  (For a~1, the mode closest to a=1 is included.)\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden."


KerraOmegaList::usage=
"KerraOmegaList[l,m,n,ReIm]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n"<>
"\t ReIm : Re, Im, or Abs\n"<>
"KerraOmegaList creates a short list of either {a,Re[\[Omega]]}, {a,-Im[\[Omega]]}, or "<>
"{a,Abs[\[Omega]]} pairs for "<>
"modes along the specified sequence.  Which sequence is produced is governed "<>
"by the choice of ReIm.\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden."


KerraAListS::usage=
"KerraAListS[l,m,n,ReIm]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n"<>
"\t ReIm : Re, Im, or Abs\n"<>
"KerraAList creates a short list of either {a,Re[\!\(\*SubscriptBox[\(A\), \(lm\)]\)]} "<>
"{a,Im[\!\(\*SubscriptBox[\(A\), \(lm\)]\)]}, or {a,Abs[\!\(\*SubscriptBox[\(A\), \(lm\)]\)]} pairs for "<>
"modes along the specified sequence.  Which sequence is produced is governed "<>
"by the choice of ReIm. Only modes where a is a multiple "<>
"of 0.05 are included.  (For a~1, the mode closest to a=1 is included.)\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden."


KerraAList::usage=
"KerraAList[l,m,n,ReIm]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n"<>
"\t ReIm : Re, Im, or Abs\n"<>
"KerraAList creates a short list of either {a,Re[\!\(\*SubscriptBox[\(A\), \(lm\)]\)]} "<>
"{a,Im[\!\(\*SubscriptBox[\(A\), \(lm\)]\)]}, or {a,Abs[\!\(\*SubscriptBox[\(A\), \(lm\)]\)]} pairs for "<>
"modes along the specified sequence.  Which sequence is produced is governed "<>
"by the choice of ReIm.\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden."


ModePlotOmegaTones::usage=
	"ModePlotOmegaTones[l,m] plots the Kerr mode frequency for all modes (l,m) "<>
	"and overtone n.  The imaginary axis is inverted.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.\n\n"<>
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t SpinWeight\[Rule] Defaults to values set by SetSpinWeight.  Can be overridden\n"<>
	"\t OTmultiple\[Rule]{}\n"<>
	"\t\t List of overtone multiplets.  An overtone multiplet is a List {n,Nmult},\n"<>
	"\t\t where 'n' is the overtone index,and 'Nmult' is the number of sequences\n"<>
	"\t\t with the same overtone index."


PlotModeFunction::usage=
"PlotModeFunction[n,s,m,a,Alm,\[Omega],Nrcf,Nm]\n"<>
"\t n : Overtone level (integer).\n"<>
"\t s : SpinWeight (integer).\n"<>
"\t m : Azimuthal index (integer).\n"<>
"\t a : Dimensionless spin parameter(rational or integer).\n"<>
"\t Alm : Approximate angular separation constant (complex number).\n"<>
"\t \[Omega] : Angular Frequency (complex number).\n"<>
"\t Nrcf : Deepth into the continued fraction (integer).\n"<>
"\t Nm : Size of the angular spectral matrix (integer).\n\n"<>
"Returns the complex value of the mode function for the specified value of \[Omega].  "<>
"The angular equation is solved based on the specified values of 'a' and \[Omega], "<>
"and the eigenvalue closest to Alm will be used.  This function is "<>
"useful for locating initial guesses for roots of the mode function.\n"<>
"When PolynomialMode (see SelectMode) is selected, values of n and Nrcf are ignored."


PlotModeFunctionL::usage=
"PlotModeFunctionL[n,s,m,a,L,\[Omega],Nrcf,Nm]\n"<>
"\t n : Overtone level (integer).\n"<>
"\t s : SpinWeight (integer).\n"<>
"\t m : Azimuthal index (integer).\n"<>
"\t a : Dimensionless spin parameter(rational or integer).\n"<>
"\t L : Index of the Angular separation constant to choose (integer).\n"<>
"\t \[Omega] : Angular Frequency (complex number).\n"<>
"\t Nrcf : Deepth into the continued fraction (integer).\n"<>
"\t Nm : Size of the angular spectral matrix (integer).\n\n"<>
"Returns the complex value of the mode function for the specified value of \[Omega].  "<>
"The angular equation is solved based on the specified values of 'a' and \[Omega], "<>
"and the \!\(\*SuperscriptBox[\(L\), \(th\)]\) eigenvalue will be used.  This function is "<>
"useful for locating initial guesses for roots of the mode function.\n"<>
"When PolynomialMode (see SelectMode) is selected, values of n and Nrcf are ignored."


SWSFLists::usage=
"SWSFLists[l,m,n,i]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n"<>
"\t i : index of specific mode in the sequence to be used\n"<>
"Returns the Spin-Weighted Spheroidal Function as functions of x = ArcCos[\[Theta]].  "<>
"By default, the function returns complex value of the Spin-Weighted Spheroidal "<>
"function.  If OutputType is set to ReIm, two lists are returned "<>
"{reallist,imaglisst} where each list is itself a list of {x,SWSH} pairs "<>
"suitable for plotting with ListPlot.  If OutputType is set to AbsArg, "<>
"two lists are returned {abslist,arglist} where each list is itself a list of "<>
"{x,SWSH} pairs suitable for plotting with ListPlot.\n\n"<>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden.\n"<>
"\t OutputType \[Rule] (Complex,ReIm,AbsArg), Defaults to Complex.\n"<>
"\t Chop \[Rule] False: If set data is Chopped at 10^(value)."


SWSFRealPath::usage=
"SWSFRealPath[l,m,n,i]\n"<>
"\t l : harmonic index\n"<>
"\t m : azimuthal index\n"<>
"\t n : overtone index (integer or overtone multiplet)\n"<>
"\t i : index of specific mode in the sequence to be used\n"<>
"Returns the complex path between the poles at z=\[PlusMinus]1 along which the Spin-Weighted "<> 
"Spheroidal Function is real, along with the the values of the Spin-Weighted "<>
"Spheroidal Function at those locations.  This information is returned as a list " <>
"of lists {path,values}.\n" <>
"Options:\n"<>
"\t ModeType \[Rule] (QNM,TTML,TTMR) Defaults according to which backage is run, \n"<>
"\t\t\t but can be overridden to plot any kind of sequence.\n"<>
"\t SpinWeight \[Rule] Defaults to values set by SetSpinWeight.  Can be overridden.\n" <>
"\t PathStart \[Rule] {+1,-1,Automatic} Defaults to Automatic.  Determines which pole \n" <>
"\t\t\t is fixed to be purely real.\n" <>
"\t MaxSteps \[Rule] Defaults to 5000.  Maximum number of path steps allowed.\n" <>
"\t StepSize \[Rule] Defaults to 1/50.  Maximum step size along path.\n" <>
"\t PlotStart \[Rule] Defaults to False.  Plot imaginary part of Spin-Weighted Spheroidal\n" <>
"\t\t\t Function around starting pole at distance of StepSize.\n"<>
"\t \[Phi]guess \[Rule] Defaults to Automatic.  Initial guess for starting direction.\n"<>
"\t PrintPoleValues \[Rule] Defaults to False."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[OTmultiple,QNM,TTML,TTMR,ModeType,ShortenBy];


Protect[PolynomialMode,ContinuedFractionMode,RunCFConvergence];


Protect[SpinWeight,ModePrecision,RadialCFDepth,RadialCFMinDepth,RadialDebug,RadialRelax,
        JacobianStep,Root\[Epsilon],SchDebug];


Protect[SolutionDebug,NoNeg\[Omega],ModePrecision,SolutionSlow,SolutionOscillate,SolutionIter,
		RadialCFDigits,NewtonRadius];


Protect[Minblevel,Maxblevel,CurvatureRatio,Max\[CapitalDelta]\[Omega],ExtrapolationOrder,Asymptote,LogLog];


Protect[ModeaStart,ModeGuess,SeqDirection,Maximala\[Epsilon],Minimala\[Epsilon],SolutionRelax,SolutionWindowl,SolutionWindowt];


Protect[Index,Refinement,RefinementAction,RefineAccuracy,RefinePrecision,RefineAdapt,FixAdapt,Update,RemoveLevels,TerminalHigh,TerminalLow,ForceRefinement,RefinementPlot,SeqLevel,RadialCFLevel,AccuracyLevel,PrecisionLevel,StepRatio,CurveRatio,LimitRefinement,RadialCFMaxGuess];


Protect[OutputType,ChopLevel,PathStart,StepSize,PlotStart,\[Phi]guess,PrintPoleValues];


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
Module[{\[Omega]r,Almr,sol,solpr,\[CapitalDelta]sol,pcount,pflag,\[CapitalDelta]\[Omega]r,RedRd\[Omega]r,RedRd\[Omega]i,ImdRd\[Omega]r,ImdRd\[Omega]i,\[Delta]\[Omega],WorkPrec,
		radialdebug=OptionValue[RadialDebug]},
	RadialLentzStep::increase="Excessive increase in Precision, Abort";
	RadialLentzStep::debug4="RadialLentzStep :`1` : `2`";
	RadialLentzStep::debug5Re="`1` \[Delta]\[Omega]r + `2` \[Delta]\[Omega]i==`3`";
	RadialLentzStep::debug5Im="`1` \[Delta]\[Omega]r + `2` \[Delta]\[Omega]i==`3`";
	RadialLentzStep::debug6="Precision - \[Delta]\[Omega] : `1`, `2`";
	RadialLentzStep::pflag="Set $MinPrecision =`1`";
	sol=ModeFunction[n,s,m,a,Alm,\[Omega],Nrcf];
	\[CapitalDelta]\[Omega]r=If[Chop[Re[\[Omega]]]==0,\[Omega]step,Re[\[Omega]]\[Omega]step,Re[\[Omega]]\[Omega]step];
	\[Omega]r=Re[\[Omega]]+\[CapitalDelta]\[Omega]r+I Im[\[Omega]];
	Almr=AngularSpectralRoot[s,m,a*\[Omega]r,Alm,Nm][[1]];
	pcount=0;pflag=False;
	While[True,
		solpr=ModeFunction[n,s,m,a,Almr,\[Omega]r,Nrcf];
		\[CapitalDelta]sol=solpr[[1]]-sol[[1]];
		If[Abs[\[CapitalDelta]sol]==0,Break[]];
		Off[Precision::mnprec,Accuracy::mnprec];
		If[MyPrecision[\[CapitalDelta]sol]==Precision[\[CapitalDelta]sol],
			On[Precision::mnprec,Accuracy::mnprec];Break[],
			If[Re[\[CapitalDelta]sol]==0 && Accuracy[Re[\[CapitalDelta]sol]]>=Precision[\[CapitalDelta]sol],
				On[Precision::mnprec,Accuracy::mnprec];Break[],
				$MinPrecision+=4,
				$MinPrecision+=4
			],
			If[Re[\[CapitalDelta]sol]==0 && Accuracy[Re[\[CapitalDelta]sol]]>=Precision[\[CapitalDelta]sol],
				On[Precision::mnprec,Accuracy::mnprec];Break[],
				$MinPrecision+=4,
				$MinPrecision+=4
			]
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


Options[TestRadialCFConvergence]=Union[Options[RadialLentzRoot],{NewtonRadius->10^(-3)}];


TestRadialCFConvergence[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
						Alm_?NumberQ,\[Omega]_?NumberQ,Nrcf_Integer,jacobian_,\[Epsilon]_Integer,Nrcfmin_Integer,
						Nm_Integer,\[Alpha]_Real|\[Alpha]_Rational|\[Alpha]_Integer,opts:OptionsPattern[]]:= 
Module[{N1,N2,Rem,CFval,CFval1,CFval2,cfpow,newNrcf,saveNrcf,diff,diffh,diffl,sol,cfpowcut=-2,
	newtonRadius=OptionValue[NewtonRadius]},
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
	sol=RadialLentzRoot[n,s,m,a,Alm,\[Omega],newNrcf,Nm,\[Epsilon],newtonRadius,
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


Options[ModeSolution]=Union[{SolutionDebug->0,NoNeg\[Omega]->False,RadialCFMinDepth->300,ModePrecision->24,
							SolutionSlow->10,SolutionOscillate->10,SolutionIter->50,
							RadialCFDigits->8,NewtonRadius->10^(-3)},Options[RadialLentzRoot]];


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
		err1,err2,N1,N2,\[CapitalDelta]N,
		newtonRadius=OptionValue[NewtonRadius],
		nonegfreq=OptionValue[NoNeg\[Omega]],
		slowval=OptionValue[SolutionSlow],oscval=OptionValue[SolutionOscillate],
		iterval=OptionValue[SolutionIter],jacobianstep=OptionValue[JacobianStep],
		solutiondebug=OptionValue[SolutionDebug],RCFmin=OptionValue[RadialCFMinDepth],
		precision=OptionValue[ModePrecision],RCFdigits=OptionValue[RadialCFDigits]},
	ModeSolution::soldebug1="\[Omega]g = `1` Almg = `2`";
	ModeSolution::soldebug2a="Initial Nradial : `1`";
	ModeSolution::soldebug2b="\[CapitalDelta]\[Omega] = `1`";
	ModeSolution::soldebug2c="Nradial : `1` , Nradialnew : `2` ";
	ModeSolution::soldebug3="a = `1` \[Omega] = `2` alm = `3` ";
	ModeSolution::soldebug4="\[Omega]0 = `1`, Alm0 = `2` ";
	ModeSolution::soldebug5a="Expconv = `1` , N = `2`";
	ModeSolution::soldebug5ad="Extrapolate angular matrix size: N = `1`";
	ModeSolution::soldebug5b="Increase Nradial to `1`";
	ModeSolution::soldebug5c="Decrease Nradial to `1`";
	ModeSolution::soldebug5d="Expconv =  `1`: count = `2`: N = `3`";
	ModeSolution::soldebug6="Jacobian matrix : `1` , rcferr : `2`";
	ModeSolution::underrelax="Reduce Under-relaxation parameter to \[Alpha] = `1`";
	ModeSolution::iterations="Too many iterations";
	ModeSolution::notnumber="Failure: Solution of RadialLentzRoot not a number";
	ModeSolution::slowconv="Persistance slow convergence \[Alpha]= `1`";
	ModeSolution::oscillations="Persistent oscillations: Abort";
	ModeSolution::failure="Persistent failure of RadialLentzRoot: Abort";
	ModeSolution::solwindow1="\[Omega] = `1`, is outside solution window";
	ModeSolution::solwindow2="Alm outside solution window";
	lmin = Max[Abs[m],Abs[s]];
	lmax = Max[l+Ceiling[Nm/2],lmin+Nm-1];
	Nmatrix=lmax-lmin+1;
	Nradialnew=Nradial=Nrcf;
	If[solutiondebug>1,Print[Style[StringForm[ModeSolution::soldebug2a,Nradial],{Medium,Darker[Blue,0.3]}]]];
	old\[Omega]=\[Omega]g;oldAlm=Almg;
	inversion=n;
	count=0;
	\[Alpha]=relax; (* under-relaxation parameter *)
	If[solutiondebug>0,Print[Style[StringForm[ModeSolution::soldebug1,\[Omega]g,Almg],{Medium,Darker[Green,0.3]}]]];
	If[solutiondebug>3,Print[Style[StringForm[ModeSolution::soldebug4,\[Omega]0,Alm0],{Medium,Darker[Orange,0.3]}]]];
	c=a \[Omega]g;
	angularsol=AngularSpectralRoot[s,m,c,oldAlm,Nmatrix];
	expconv=Max[Take[Abs[angularsol[[3]]],-2]];
	err1=expconv; N1=Nmatrix;
	While[expconv >= 10^\[Epsilon]2 && ++count<5,  (* Make sure Spectral resolution is good before doing lots of work *)
		err2=err1; N2=N1;
		err1=expconv; N1=Nmatrix;
		If[count>2,
			\[CapitalDelta]N=Max[1,Floor[-(N1-N2)Log[10^\[Epsilon]2/err1]/Log[err2/err1]]];
			If[\[CapitalDelta]N/N1>1/10,\[CapitalDelta]N=Max[Min[Floor[2\[CapitalDelta]N/3],10],1]];
			Nmatrix+=\[CapitalDelta]N-1;
			If[solutiondebug>4,Print[Style[StringForm[ModeSolution::soldebug5ad,Nmatrix+1],{Medium,Darker[Yellow,0.3]}]]]
		];
		angularsol=AngularSpectralRoot[s,m,c,oldAlm,++Nmatrix];
		expconv=Max[Take[Abs[angularsol[[3]]],-2]];
		If[solutiondebug>4,Print[Style[StringForm[ModeSolution::soldebug5a,expconv,Nmatrix],{Medium,Darker[Yellow,0.3]}]]];
	];
	count=0;
	While[Not[converged],
		If[(++iteration>iterval && \[Epsilon]2==\[Epsilon] && \[Alpha]==1) || (\[Alpha]<1 && iteration > (3/5)iterval/\[Alpha]),
			If[\[Alpha]>0.05,
				\[Alpha]*=3/5;iteration=0;
				Print[Style[StringForm[ModeSolution::underrelax,N[\[Alpha]]],{Medium,Darker[Red]}]],
				Message[ModeSolution::iterations];
				(*Print["a=",a," \[Omega]=",old\[Omega]," Alm=",oldAlm]*);
				Return[{False}];
			];
		];
		radialsol = RadialLentzRoot[inversion,s,m,a,oldAlm,old\[Omega],Nradial,Nmatrix,\[Epsilon]2,newtonRadius,
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
			radialsol = RadialLentzRoot[inversion,s,m,a,oldAlm,radialsol[[2,1]],Nradial,Nmatrix,\[Epsilon]2,newtonRadius,
											RadialRelax->\[Alpha],FilterRules[{opts}, Options[RadialLentzRoot]]];
			If[Not[radialsol[[1,1]]] && Length[radialsol]==1, (* Re[\[Omega]] flipping sign, solution fails *)
				Return[{False,\[Alpha],Nradialnew,{a,Null[],Null[]}}]
			];
			If[Head[radialsol[[1,3]]]==List,jacobianmatrix=radialsol[[1,3]]];
		];
		If[Not[radialsol[[1,1]]],++radialfail;(*Print["WARNING: RadialLentzRoot failed to converged"]*)];
		If[Not[NumberQ[radialsol[[2,1]]]],Message[ModeSolution::notnumber];Abort[]];
		c=a radialsol[[2,1]];
		angularsol=AngularSpectralRoot[s,m,c,oldAlm,Nmatrix];
		\[CapitalDelta]\[Omega]2=\[CapitalDelta]\[Omega];
		\[CapitalDelta]\[Omega]=radialsol[[2,5]];
		If[TrueQ[nonegfreq] && Re[radialsol[[2,1]]]<0,
				radialsol[[2,1]]=-Conjugate[radialsol[[2,1]]];
				angularsol[[1]]=Conjugate[angularsol[[1]]];
		];
		old\[Omega]=radialsol[[2,1]];oldAlm=angularsol[[1]];
		If[solutiondebug>2,Print[Style[StringForm[ModeSolution::soldebug3,a,old\[Omega],oldAlm],{Medium,Darker[Magenta,0.3]}]]];
		If[solutiondebug>1,Print[Style[StringForm[ModeSolution::soldebug2b,\[CapitalDelta]\[Omega]],{Medium,Darker[Blue,0.3]}]]];
		If[\[CapitalDelta]\[Omega]>0,
			If[Abs[\[CapitalDelta]\[Omega]2]/Abs[\[CapitalDelta]\[Omega]]<1,
				If[++slowcount2>slowval,
					(* Try under-relaxation if solution is slow *)
					slowcount2=0;radialfail=0;
					\[Alpha]*=3/5;iteration=0;
					Print[Style[StringForm[ModeSolution::slowconv,N[\[Alpha]]],{Medium,Darker[Red]}]];
					If[\[Alpha]<0.05,Return[{False}]]
				], 
				If[slowcount2 >0,slowcount2=0;
					(* Try lower precision if persisitent oscillations *)
					If[++oscillate>oscval,
						Message[ModeSolution::oscillations];Abort[]
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
			Message[ModeSolution::failure];Abort[]
		];
		If[\[CapitalDelta]\[Omega]<10^\[Epsilon]2,
			If[count==0,
				If[RunCFConvergence,
					rcfpower=0;
					rcferr=TestRadialCFConvergence[inversion,s,m,a,oldAlm(1+10^(\[Epsilon]+4)),old\[Omega](1-10^(\[Epsilon]+4)),
							Nradial,jacobianmatrix,\[Epsilon]2,RCFmin,Nmatrix,\[Alpha],FilterRules[{opts},Options[TestRadialCFConvergence]]];
					Nradialnew=rcferr[[1]];rcfpower=rcferr[[5]];
					If[solutiondebug>5,
						Print[Style[StringForm[ModeSolution::soldebug6,jacobianmatrix,rcferr],{Medium,Darker[Cyan,0.7]}]];
					];
					If[solutiondebug>1,Print[Style[StringForm[ModeSolution::soldebug2c,Nradial,Nradialnew],{Medium,Darker[Blue,0.3]}]]];
					If[Nradialnew>Nradial,
						If[Nradialnew>(11/10)Nradial,
							NradFlag=True;count=0;iteration/=2;slowcount2=0;oscillate=0;
							If[solutiondebug>4,Print[Style[StringForm[ModeSolution::soldebug5b,Nradialnew],{Medium,Darker[Yellow,0.3]}]]];
						];
						Nradial=Nradialnew;
					,
						If[solutiondebug>4 && Nradialnew<Nradial,Print[Style[StringForm[ModeSolution::soldebug5c,Nradialnew],{Medium,Darker[Yellow,0.3]}]]];
					];
				];
			];

			expconv=Max[Take[Abs[angularsol[[3]]],-2]];
			(*expconv=0;*)
			If[solutiondebug>4,Print[Style[StringForm[ModeSolution::soldebug5d,expconv,count,Nmatrix],{Medium,Darker[Yellow,0.3]}]]];
			If[expconv < 10^\[Epsilon]2,
				If[!NradFlag,++count];
				If[count>1&&radialsol[[1,1]],converged=True],
				err2=err1; N2=N1; err1=expconv; N1=Nmatrix;
				\[CapitalDelta]N=Max[1,Floor[-(N1-N2)Log[10^\[Epsilon]2/err1]/Log[err2/err1]]];
				If[\[CapitalDelta]N/N1>1/10,\[CapitalDelta]N=Max[Floor[2\[CapitalDelta]N/3],1]];
				Nmatrix+=\[CapitalDelta]N;
				If[solutiondebug>4,Print[Style[StringForm[ModeSolution::soldebug5ad,Nmatrix],{Medium,Darker[Yellow,0.3]}]]];
				count=0;iteration/=2;slowcount2=0;oscillate=0
			];
			NradFlag=False,
			count=0
		];
	];
	If[rl!=0 &&rt!=0,
		If[Not[SolutionWindow[\[Omega]0,\[Omega]g,old\[Omega],rl,rt,True]],
			Print[Style[StringForm[ModeSolution::solwindow1,old\[Omega]],{Medium,Darker[Red]}]];
			Return[{False,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix]}],angularsol}}]
		];
		If[Not[SolutionWindow[Alm0,Almg,oldAlm,rl,rt,True]],
			Print[Style[StringForm[ModeSolution::solwindow2],{Medium,Darker[Red]}]];
			Return[{False,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix]}],angularsol}}]
		];
	];
(*If[Nradialnew==Nrcf,Print["WARNING: Nradialnew not reset from Nrcf"]];*)
	{True,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix],$MinPrecision}],angularsol}}
]


(* ::Subsection::Closed:: *)
(*Adaptive Bisection sequencer*)


Options[AdaptCheck3]=Union[{Minblevel->0,Maxblevel->20,CurvatureRatio->1/2,Max\[CapitalDelta]\[Omega]->0.01,ExtrapolationOrder->2},Options[ModeSolution]];


Options[KerrModeSequence]=Union[{SpinWeight->Null[],ModeaStart->0,ModeGuess->0,
								SeqDirection->Forward,Maximala\[Epsilon]->10,Minimala\[Epsilon]->False,
								SolutionRelax->1,RadialCFDepth->1,
								SolutionWindowl->1/2,SolutionWindowt->1/3},
								Options[AdaptCheck3]];


KerrModeSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]max_Integer,
					opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],SpinWeightTable,KerrSEQ,KerrSEQret,AC3ret,SeqStatus,context,
		\[Omega]Guess,inversion,\[Omega],Alm,\[Omega]try,Almtry,ModeSol,a,index0=0,index0p,index0m,
		\[Epsilon]=\[Epsilon]max,Nrcf,Nm=4,rl=0,rt=0,NKMode=0,edat0,edat,ef,iv,afit,
		maximuma=1,minimuma=0,blevel=0,\[CapitalDelta]a=10^(-3),\[CapitalDelta]a2,\[CapitalDelta]a3,dir=1,\[CapitalDelta]aincflag=False,\[CapitalDelta]aincstep=0,blevelsave,
		\[Omega]0,\[Omega]p,\[Omega]m,Alm0,Almp,Almm,\[Omega]w=0,Almw=0,precisionsave,precisioncount=0,
		Minb=OptionValue[Minblevel],Maxb=OptionValue[Maxblevel],
		forward=If[OptionValue[SeqDirection]==Backward,False,True,True],
		modeastart=OptionValue[ModeaStart],guess=OptionValue[ModeGuess],
		precision=OptionValue[ModePrecision],extraporder=OptionValue[ExtrapolationOrder],
		solwinl=OptionValue[SolutionWindowl],solwint=OptionValue[SolutionWindowt],
		relax,srelax=Rationalize[OptionValue[SolutionRelax]],RCFmin=OptionValue[RadialCFMinDepth],
		rcfdepth=OptionValue[RadialCFDepth]},

	KerrModeSequence::minprecision="Set $MinPrecision = `1`";
	KerrModeSequence::invalidmax="Invalid value for Maximala\[Epsilon]";
	KerrModeSequence::invalidmin="Invalid value for Minimala\[Epsilon]";
	KerrModeSequence::entries="`1`[`2`,`3`,`4`] sequence exists with `5` entries";
	KerrModeSequence::untested1="Untested section of code! 1";
	KerrModeSequence::untested2="Untested section of code! 2";
	KerrModeSequence::guesses="Guesses set: `1` : `2` : `3` : `4`";
	KerrModeSequence::unusualstepsize="Sequence has unusual stepsizes, Aborting";
	KerrModeSequence::decblevel="Decreasing \[CapitalDelta]a, blevel = `1`";
	KerrModeSequence::decblevelmore="Further decreasing \[CapitalDelta]a, blevel = `1`";
	KerrModeSequence::incblevel="Increasing \[CapitalDelta]a, blevel = `1`";
	KerrModeSequence::missuse="Cannot use Accumulation extrapolation with backward sequencing";
	KerrModeSequence::startseq="Starting `1`[`2`,`3`,`4`] sequence";
	KerrModeSequence::status="Error determining status of `1`[`2`,`3`,`4`] sequence: Abort";
	KerrModeSequence::nosol="No solution found at a = `1`";
	KerrModeSequence::nosoltry="No solution found at a = `1`, try decreasing \[CapitalDelta]a, blevel = `2`";
	KerrModeSequence::nkmodeindex="NKMode <= 2 : index0 = `1`";
	KerrModeSequence::solfail="a+/- solution failed.";
	KerrModeSequence::invalidcall="Invalid call to ModeSolution"; 

	If[precision!=$MinPrecision,Print[Style[StringForm[KerrModeSequence::minprecision,precision],{Medium,Darker[Red]}]]];
	$MinPrecision=precisionsave=precision; (* Sets the minimum precision for entire calculation *)
	maximuma=If[OptionValue[Maximala\[Epsilon]] \[Element] Integers,
					1-(2^-OptionValue[Maximala\[Epsilon]])/1000,
					If[OptionValue[Maximala\[Epsilon]] \[Element] Booleans && OptionValue[Maximala\[Epsilon]]==False,
						1,
						Message[KerrModeSequence::invalidmax];Abort[]
					]
				];
	minimuma=If[OptionValue[Minimala\[Epsilon]] \[Element] Integers,
					(2^-OptionValue[Minimala\[Epsilon]])/1000,
					If[OptionValue[Minimala\[Epsilon]] \[Element] Booleans && OptionValue[Minimala\[Epsilon]]==False,
						0,
						Message[KerrModeSequence::invalidmin];Abort[]
					]
				];
	SpinWeightTable:=modeName;
	KerrSEQ:=modeName[l,m,n];
	SeqStatus=If[Head[KerrSEQ]==List,If[Length[KerrSEQ]>0,True,False,False],False,False];
	dir=If[forward,1,-1];
	relax=srelax;
	blevel=Minb;
	If[SeqStatus,
		(* Sequence exists, extend *)
		NKMode=Length[KerrSEQ];
		Print[Style[StringForm[KerrModeSequence::entries,SpinWeightTable,l,m,n,NKMode],{Medium,Darker[Green]}]];
		If[forward,
			\[Epsilon]=Min[\[Epsilon]max,KerrSEQ[[NKMode,2,4]]];
			If[Length[KerrSEQ[[NKMode,2]]]>=9,$MinPrecision=KerrSEQ[[NKMode,2,9]],Print[Style[StringForm[KerrModeSequence::minprecision,precision],{Medium,Darker[Red]}]]],
			\[Epsilon]=Min[\[Epsilon]max,KerrSEQ[[1,2,4]]];
			If[Length[KerrSEQ[[1,2]]]>=9,$MinPrecision=KerrSEQ[[1,2,9]],Print[Style[StringForm[KerrModeSequence::minprecision,precision],{Medium,Darker[Red]}]]];
		];
		a=If[forward,KerrSEQ[[NKMode,1]],KerrSEQ[[1,1]]];
		inversion=If[Head[n]==Integer,n,Null[],n[[1]]];
		If[forward,
			\[Omega]=SetPrecision[KerrSEQ[[NKMode,2,1]],Max[precision,$MinPrecision]];
			Alm=SetPrecision[KerrSEQ[[NKMode,3,1]],Max[precision,$MinPrecision]],
			\[Omega]=SetPrecision[KerrSEQ[[1,2,1]],Max[precision,$MinPrecision]];
			Alm=SetPrecision[KerrSEQ[[1,3,1]],Max[precision,$MinPrecision]]
		];
		Nm=If[forward,KerrSEQ[[NKMode,3,2]],KerrSEQ[[1,3,2]]];
		Nrcf=If[forward,
			If[Length[KerrSEQ[[NKMode,2]]]>=6,KerrSEQ[[NKMode,2,6]],KerrSEQ[[NKMode,2,3]]],
			If[Length[KerrSEQ[[1,2]]]>=6,KerrSEQ[[1,2,6]],KerrSEQ[[1,2,3]]]
		];
		If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
		If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
		Nrcf=Max[Nrcf,RCFmin];
		If[NKMode>=1,
			If[NKMode==2,
				\[CapitalDelta]a=KerrSEQ[[2,1]]-KerrSEQ[[1,1]];
				blevelsave=blevel=Round[-(3+Log10[\[CapitalDelta]a])/Log10[2]];
				If[forward,
					\[Omega]=2KerrSEQ[[2,2,1]]-KerrSEQ[[1,2,1]];Alm=2KerrSEQ[[2,3,1]]-KerrSEQ[[1,3,1]],
					\[Omega]=2KerrSEQ[[1,2,1]]-KerrSEQ[[2,2,1]];Alm=2KerrSEQ[[1,3,1]]-KerrSEQ[[2,3,1]]
				]
			];
			If[NKMode>2,
				index0=If[forward,NKMode-1,2];
				\[CapitalDelta]a=If[forward,KerrSEQ[[NKMode,1]]-KerrSEQ[[index0,1]],KerrSEQ[[index0,1]]-KerrSEQ[[1,1]]];
				blevelsave=blevel=Round[-(3+Log10[\[CapitalDelta]a])/Log10[2]];
				\[CapitalDelta]a2=If[forward,KerrSEQ[[index0,1]]-KerrSEQ[[index0-1,1]],KerrSEQ[[index0+1,1]]-KerrSEQ[[index0,1]]];
				If[\[CapitalDelta]a > 2\[CapitalDelta]a2,Message[KerrModeSequence::unusualstepsize];Abort[]];
				\[CapitalDelta]a3=2^(-Maxb)/1000;
				If[a+dir*\[CapitalDelta]a3>maximuma||a+dir*\[CapitalDelta]a3<minimuma,Return[]];
				If[\[CapitalDelta]a>=\[CapitalDelta]a2,
					If[\[CapitalDelta]a==2\[CapitalDelta]a2,
						{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
						AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,maximuma,minimuma,forward,True,False,FilterRules[{opts},Options[AdaptCheck3]]],
						{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
						AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,maximuma,minimuma,forward,False,False,FilterRules[{opts},Options[AdaptCheck3]]]
					];
					modeName[l,m,n]=KerrSEQret;
				];
				NKMode=Length[KerrSEQ];
				index0=If[forward,NKMode-1,2];
				\[CapitalDelta]a=2^(-blevel)/1000;
				If[blevel>blevelsave,Print[Style[StringForm[KerrModeSequence::decblevel,blevel],{Medium,Darker[Green]}]]];
				If[blevel<blevelsave,Print[Style[StringForm[KerrModeSequence::incblevel,blevel],{Medium,Darker[Green]}]]];
				index0p=If[forward,index0+1,index0+\[CapitalDelta]aincstep];
				index0m=If[forward,index0-\[CapitalDelta]aincstep,index0-1];
				\[Omega]0=SetPrecision[KerrSEQ[[index0,2,1]],Max[precision,$MinPrecision]];
				\[Omega]p=SetPrecision[KerrSEQ[[index0p,2,1]],Max[precision,$MinPrecision]];
				\[Omega]m=SetPrecision[KerrSEQ[[index0m,2,1]],Max[precision,$MinPrecision]];
				Alm0=SetPrecision[KerrSEQ[[index0,3,1]],Max[precision,$MinPrecision]];
				Almp=SetPrecision[KerrSEQ[[index0p,3,1]],Max[precision,$MinPrecision]];
				Almm=SetPrecision[KerrSEQ[[index0m,3,1]],Max[precision,$MinPrecision]];
				\[Omega]w=If[forward,\[Omega]p,\[Omega]m];
				Almw=If[forward,Almp,Almm];
				rl=solwinl;rt=solwint;
				If[\[CapitalDelta]aincflag,
					(* Do not reset \[CapitalDelta]aincflag, it is needed for next call to AdaptCheck3 *)
					If[forward,
						\[Omega]=6\[Omega]p-8\[Omega]0+3\[Omega]m;Alm=6Almp-8Alm0+3Almm,
						\[Omega]=6\[Omega]m-8\[Omega]0+3\[Omega]p;Alm=6Almm-8Alm0+3Almp
					],
					If[forward,
						\[Omega]=3(\[Omega]p-\[Omega]0)+\[Omega]m;Alm=3(Almp-Alm0)+Almm,
						\[Omega]=3(\[Omega]m-\[Omega]0)+\[Omega]p;Alm=3(Almm-Alm0)+Almp
					];
					Switch[extraporder,
						Accumulate,
						If[!forward,
							Message[KerrModeSequence::missuse];
							Abort[]
						];
						edat0=SetPrecision[Take[KerrSEQ,-10],Max[precision,$MinPrecision]];
						edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
						\[Omega]=afit[1-(KerrSEQ[[NKMode,1]]+\[CapitalDelta]a)];
						edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						\[Omega]+=I afit[1-(KerrSEQ[[NKMode,1]]+\[CapitalDelta]a)],
						LogLog,
						edat0=Take[KerrSEQ,If[forward,-2,2]];
						edat=Table[{Log[edat0[[i,1]]],Log[Abs[Re[edat0[[i,2,1]]]]]},{i,1,2}];
						ef=LinearModelFit[edat,{1,af},af];
						\[Omega]=Sign[Re[edat0[[1,2,1]]]]Exp[ef[Log[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]]]];
						edat=Table[{Log[edat0[[i,1]]],Log[Abs[Im[edat0[[i,2,1]]]]]},{i,1,2}];
						ef=LinearModelFit[edat,{1,af},af];
						\[Omega]+=I Sign[Im[edat0[[1,2,1]]]]Exp[ef[Log[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]]]],
						Asymptote,
						{\[Omega],Alm}=AsymptoteFunction[s,l,m,If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]],
						_,
						If[extraporder>2,
							edat0=Take[KerrSEQ,If[forward,-(extraporder+1),extraporder+1]];
							edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
							ef[iv_]=InterpolatingPolynomial[edat,iv];
							\[Omega]=ef[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]];
						];
					];
				];
			];
			If[Head[guess]==List,
Print[KerrModeSequence::untested1];
				\[Omega]=guess[[1]];
				Alm=guess[[2]];
				If[Length[guess]>=3,Nrcf=guess[[3]]];
				If[Length[guess]==4,Nm=guess[[4]]];
				Print[Style[StringForm[KerrModeSequence::guesses,\[Omega],Alm,Nrcf,Nm],{Medium,Darker[Red]}]];
			];
		],
		(* Sequence does not exist, start it *)
		If[Head[KerrSEQ]==SpinWeightTable || Length[KerrSEQ]==0,
			Print[Style[StringForm[KerrModeSequence::startseq,SpinWeightTable,l,m,n],{Medium,Darker[Green]}]];
			a=0;
			inversion=If[Head[n]==Integer,n,Null[],n[[1]]];
			If[Head[modeastart]==List,
			(* For cases that cannot start at a=0 *)
				a=modeastart[[1]];
				\[Omega]=SetPrecision[modeastart[[2]],Max[precision,$MinPrecision]];
				Alm=SetPrecision[modeastart[[3]],Max[precision,$MinPrecision]];
				If[Length[modeastart]==4,Nm=modeastart[[4]]],
				Null[],
				\[Omega]Guess=If[Head[n]==Integer,
							SetPrecision[SchGuess[l,n,FilterRules[{opts},Options[SchGuess]]],Max[precision,$MinPrecision]],
							Null[],
							SetPrecision[SchGuess[l,n[[1]],FilterRules[{opts},Options[SchGuess]]],Max[precision,$MinPrecision]]
							];
				\[Omega]=SetPrecision[\[Omega]Guess[[1]],Max[precision,$MinPrecision]];
				Alm = l(l+1)-s(s+1);
			];
			a=a-dir*\[CapitalDelta]a; (* offset to "previous" a *)
			Nrcf=RCFmin;
			If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
			modeName[l,m,n]={},
			Message[KerrModeSequence::status,SpinWeightTable,l,m,n];
			Abort[],
			Message[KerrModeSequence::status,SpinWeightTable,l,m,n];
			Abort[]
		];
	];
	While[minimuma <= a+dir*\[CapitalDelta]a <= maximuma,
		(* Try lowering precision every so often *)
		If[$MinPrecision<precisionsave,precisionsave=$MinPrecision;precisioncount=0];
		If[$MinPrecision>precisionsave,
			precisioncount=0;precisionsave=$MinPrecision,
			If[$MinPrecision==precisionsave && ++precisioncount>10,
				$MinPrecision=Max[precision,precisionsave-=4];precisioncount=0]
		];
		(* Print["lasta = ",Block[{$MinPrecision=0},N[a,{Infinity,20}]]," \[CapitalDelta]a = ",Block[{$MinPrecision=0},N[dir*\[CapitalDelta]a,{Infinity,20}]]]; *)
		ModeSol = ModeSolution[inversion,s,l,m,a+dir*\[CapitalDelta]a,SetPrecision[\[Omega],Max[precision,$MinPrecision]],SetPrecision[Alm,Max[precision,$MinPrecision]],\[Epsilon],relax,Nrcf,Nm,\[Omega]w,Almw,rl,rt,FilterRules[{opts},Options[ModeSolution]]];
		If[ModeSol[[1]],
			(* Solution found, save solution check sequence smoothness *)
			a=a+dir*\[CapitalDelta]a;
			Print["ModeSol a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision]];
			(* Print["levelcount = ",levelcount]; *)
			relax=Min[srelax,5/3ModeSol[[2]]];
			Nm=ModeSol[[4,3,2]];
			Nrcf=ModeSol[[3]];
			Nrcf=Max[Nrcf,RCFmin];
			(* Nrcf=Max[RCFmin,Nrcf 4/5]; speed up solution? *)
			If[forward, 
				AppendTo[modeName[l,m,n],ModeSol[[4]]],
				PrependTo[modeName[l,m,n],ModeSol[[4]]];
					   
			];
			NKMode=Length[KerrSEQ];
			If[NKMode==1,\[Omega]=KerrSEQ[[1,2,1]];Alm=KerrSEQ[[1,3,1]]];
			If[NKMode>1,index0=If[forward,NKMode-1,2]];
			If[NKMode==2,
				If[forward,
					\[Omega]=2KerrSEQ[[2,2,1]]-KerrSEQ[[1,2,1]];Alm=2KerrSEQ[[2,3,1]]-KerrSEQ[[1,3,1]],
					\[Omega]=2KerrSEQ[[1,2,1]]-KerrSEQ[[2,2,1]];Alm=2KerrSEQ[[1,3,1]]-KerrSEQ[[2,3,1]]
				]
			];
			If[NKMode>2,
				blevelsave=blevel;
				{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
					AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,maximuma,minimuma,forward,\[CapitalDelta]aincflag,False,FilterRules[{opts},Options[AdaptCheck3]]];
				modeName[l,m,n]=KerrSEQret;
				NKMode=Length[KerrSEQ];
				index0=If[forward,NKMode-1,2];
				\[CapitalDelta]a=2^(-blevel)/1000;
				index0p=If[forward,index0+1,index0+\[CapitalDelta]aincstep];
				index0m=If[forward,index0-\[CapitalDelta]aincstep,index0-1];
				\[Omega]0=KerrSEQ[[index0,2,1]];
				\[Omega]p=KerrSEQ[[index0p,2,1]];
				\[Omega]m=KerrSEQ[[index0m,2,1]];
				Alm0=KerrSEQ[[index0,3,1]];
				Almp=KerrSEQ[[index0p,3,1]];
				Almm=KerrSEQ[[index0m,3,1]];
				\[Omega]w=If[forward,\[Omega]p,\[Omega]m];
				Almw=If[forward,Almp,Almm];
				rl=solwinl;rt=solwint;
				If[\[CapitalDelta]aincflag,
					(* Do not reset \[CapitalDelta]aincflag, it is needed for next call to AdaptCheck3 *)
					Print[Style[StringForm[KerrModeSequence::incblevel,blevel],{Medium,Darker[Green]}]];
					If[forward,
						\[Omega]=6\[Omega]p-8\[Omega]0+3\[Omega]m;Alm=6Almp-8Alm0+3Almm,
						\[Omega]=6\[Omega]m-8\[Omega]0+3\[Omega]p;Alm=6Almm-8Alm0+3Almp
					];
					Switch[extraporder,
						Accumulate,
						If[!forward,
							Print[KerrModeSequence::missuse];
							Abort[]
						];
						edat0=SetPrecision[Take[KerrSEQ,-10],Max[precision,$MinPrecision]];
						edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
						\[Omega]=afit[1-(KerrSEQ[[NKMode,1]]+\[CapitalDelta]a)];
						edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						\[Omega]+=I afit[1-(KerrSEQ[[NKMode,1]]+\[CapitalDelta]a)],
						LogLog,
						edat0=Take[KerrSEQ,If[forward,-2,2]];
						edat=Table[{Log[edat0[[i,1]]],Log[Abs[Re[edat0[[i,2,1]]]]]},{i,1,2}];
						ef=LinearModelFit[edat,{1,af},af];
						\[Omega]=Sign[Re[edat0[[1,2,1]]]]Exp[ef[Log[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]]]];
						edat=Table[{Log[edat0[[i,1]]],Log[Abs[Im[edat0[[i,2,1]]]]]},{i,1,2}];
						ef=LinearModelFit[edat,{1,af},af];
						\[Omega]+=I Sign[Im[edat0[[1,2,1]]]]Exp[ef[Log[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]]]],
						Asymptote,
						{\[Omega],Alm}=AsymptoteFunction[s,l,m,If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]],
						_,
						If[extraporder>2,
							edat0=Take[KerrSEQ,If[forward,-(extraporder+1),extraporder+1]];
							edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
							ef[iv_]=InterpolatingPolynomial[edat,iv];
							\[Omega]=ef[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]];
						];
					];
					,
					If[blevel>blevelsave,Print[Style[StringForm[KerrModeSequence::decblevel,blevel],{Medium,Darker[Green]}]]];
					If[forward,
						\[Omega]=3(\[Omega]p-\[Omega]0)+\[Omega]m;Alm=3(Almp-Alm0)+Almm,
						\[Omega]=3(\[Omega]m-\[Omega]0)+\[Omega]p;Alm=3(Almm-Alm0)+Almp
					];
					Switch[extraporder,
						Accumulate,
						If[!forward,
							Print[KerrModeSequence::missuse];
							Abort[]
						];
						edat0=SetPrecision[Take[KerrSEQ,-10],Max[precision,$MinPrecision]];
						edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
						\[Omega]=afit[1-(KerrSEQ[[NKMode,1]]+\[CapitalDelta]a)];
						edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						\[Omega]+=I afit[1-(KerrSEQ[[NKMode,1]]+\[CapitalDelta]a)],
						LogLog,
						edat0=Take[KerrSEQ,If[forward,-2,2]];
						edat=Table[{Log[edat0[[i,1]]],Log[Abs[Re[edat0[[i,2,1]]]]]},{i,1,2}];
						ef=LinearModelFit[edat,{1,af},af];
						\[Omega]=Sign[Re[edat0[[1,2,1]]]]Exp[ef[Log[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]]]];
						edat=Table[{Log[edat0[[i,1]]],Log[Abs[Im[edat0[[i,2,1]]]]]},{i,1,2}];
						ef=LinearModelFit[edat,{1,af},af];
						\[Omega]+=I Sign[Im[edat0[[1,2,1]]]]Exp[ef[Log[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]]]],
						Asymptote,
						{\[Omega],Alm}=AsymptoteFunction[s,l,m,If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]],
						_,
						If[extraporder>2,
							edat0=Take[KerrSEQ,If[forward,-(extraporder+1),extraporder+1]];
							edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
							ef[iv_]=InterpolatingPolynomial[edat,iv];
							\[Omega]=ef[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]];
						];
					];
				];
			],
			(* No solution found, try decreasing a *)
			blevelsave= ++blevel;
			If[blevel > Maxb,Return[]];
			If[NKMode==0,Message[KerrModeSequence::nosol,Block[{$MinPrecision=0},N[a+dir*\[CapitalDelta]a,{Infinity,20}]]];Abort[]];
			Print[Style[StringForm[KerrModeSequence::nosoltry,Block[{$MinPrecision=0},N[a+dir*\[CapitalDelta]a,{Infinity,20}]],blevel],{Medium,Darker[Red]}]];
			If[NKMode>2,
				index0p=If[forward,index0+1,index0+\[CapitalDelta]aincstep];
				index0m=If[forward,index0-\[CapitalDelta]aincstep,index0-1],
				Print[Style[StringForm[KerrModeSequence::nkmodeindex,index0],{Medium,Darker[Green]}]];
				If[NKMode==2,
					index0=If[forward,1,2];
					index0p=If[forward,2,1];
					index0m=If[forward,1,2],
					index0p=index0m=index0=1;
				];
			];
			\[Omega]0=KerrSEQ[[index0,2,1]];
			\[Omega]p=KerrSEQ[[index0p,2,1]];
			\[Omega]m=KerrSEQ[[index0m,2,1]];
			Alm0=KerrSEQ[[index0,3,1]];
			Almp=KerrSEQ[[index0p,3,1]];
			Almm=KerrSEQ[[index0m,3,1]];
			\[Omega]w=If[forward,\[Omega]p,\[Omega]m];
			Almw=If[forward,Almp,Almm];
			rl=solwinl;rt=solwint;
			If[\[CapitalDelta]aincflag,
			(* \[CapitalDelta]a increase caused failure, reset and try again *)
				\[CapitalDelta]a=2^(-blevel)/1000;
				\[CapitalDelta]aincflag=False; (* \[CapitalDelta]a caused failure, reset and try again *)
				If[forward,
					\[Omega]=3(\[Omega]p-\[Omega]0)+\[Omega]m;Alm=3(Almp-Alm0)+Almm,
					\[Omega]=3(\[Omega]m-\[Omega]0)+\[Omega]p;Alm=3(Almm-Alm0)+Almp
				],
			(* failure not due to increase in \[CapitalDelta]a, refine previous step *)
				If[forward,
					\[Omega] = (6\[Omega]0+3\[Omega]p-\[Omega]m)/8; Alm = (6Alm0+3Almp-Almm)/8,
					\[Omega] = (6\[Omega]0+3\[Omega]m-\[Omega]p)/8; Alm = (6Alm0+3Almm-Almp)/8;
				];
				(*Print["Fail call ModeSol, a= ",Block[{$MinPrecision=0},N[KerrSEQ[[index0,1]]+dir \[CapitalDelta]a/2,{Infinity,20}]]];*)
				If[forward,
					Nrcf=Max[If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]],
								If[Length[KerrSEQ[[index0p,2]]]>=6,KerrSEQ[[index0p,2,6]],KerrSEQ[[index0p,2,3]]]],
					Nrcf=Max[If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]],
								If[Length[KerrSEQ[[index0m,2]]]>=6,KerrSEQ[[index0m,2,6]],KerrSEQ[[index0m,2,3]]]]
				];
				Nrcf=Max[Nrcf,RCFmin];
				Nm=KerrSEQ[[index0,3,2]];
				ModeSol=ModeSolution[inversion,s,l,m,
									KerrSEQ[[index0,1]]+dir \[CapitalDelta]a/2,
									SetPrecision[\[Omega],Max[precision,$MinPrecision]],
									SetPrecision[Alm,Max[precision,$MinPrecision]],\[Epsilon],relax,
									Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
				If[Not[ModeSol[[1]]],Message[KerrModeSequence::solfail];Abort[]];
			    Print["ModeSol+/- a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision]];
				If[forward,index0=index0+1];
				modeName[l,m,n]=Insert[KerrSEQ,ModeSol[[4]],index0];
				NKMode=Length[KerrSEQ];
				blevelsave=blevel;
				If[NKMode>2,
					{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
						AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,maximuma,minimuma,forward,False,False,Minblevel->Max[blevel,OptionValue[Minblevel]],FilterRules[{opts},Options[AdaptCheck3]]];
					modeName[l,m,n]=KerrSEQret;
					NKMode=Length[KerrSEQ];
					index0=If[forward,NKMode-1,2];
					\[CapitalDelta]a=2^(-blevel)/1000;
					a=KerrSEQ[[index0,1]]+dir*\[CapitalDelta]a;
					If[\[CapitalDelta]aincflag,\[CapitalDelta]aincflag=False,blevel=blevelsave]; (* Don't allow \[CapitalDelta]a to increase *)
					index0p=If[forward,index0+1,index0+\[CapitalDelta]aincstep];
					index0m=If[forward,index0-\[CapitalDelta]aincstep,index0-1];
					\[Omega]0=KerrSEQ[[index0,2,1]];
					\[Omega]p=KerrSEQ[[index0p,2,1]];
					\[Omega]m=KerrSEQ[[index0m,2,1]];
					Alm0=KerrSEQ[[index0,3,1]];
					Almp=KerrSEQ[[index0p,3,1]];
					Almm=KerrSEQ[[index0m,3,1]];
					\[Omega]w=If[forward,\[Omega]p,\[Omega]m];
					Almw=If[forward,Almp,Almm];
					rl=solwinl;rt=solwint;
					If[blevel>blevelsave,Print[Style[StringForm[KerrModeSequence::decblevelmore,blevel],{Medium,Darker[Green]}]]];
					If[forward,
						\[Omega]=3(\[Omega]p-\[Omega]0)+\[Omega]m;Alm=3(Almp-Alm0)+Almm,
						\[Omega]=3(\[Omega]m-\[Omega]0)+\[Omega]p;Alm=3(Almm-Alm0)+Almp
					];
					Switch[extraporder,
						Accumulate,
						If[!forward,
							Print[KerrModeSequence::missuse];
							Abort[]
						];
						edat0=SetPrecision[Take[KerrSEQ,-10],Max[precision,$MinPrecision]];
						edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
						\[Omega]=afit[1-(KerrSEQ[[NKMode,1]]+\[CapitalDelta]a)];
						edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						\[Omega]+=I afit[1-(KerrSEQ[[NKMode,1]]+\[CapitalDelta]a)],
						LogLog,
						edat0=Take[KerrSEQ,If[forward,-2,2]];
						edat=Table[{Log[edat0[[i,1]]],Log[Abs[Re[edat0[[i,2,1]]]]]},{i,1,2}];
						ef=LinearModelFit[edat,{1,af},af];
						\[Omega]=Sign[Re[edat0[[1,2,1]]]]Exp[ef[Log[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]]]];
						edat=Table[{Log[edat0[[i,1]]],Log[Abs[Im[edat0[[i,2,1]]]]]},{i,1,2}];
						ef=LinearModelFit[edat,{1,af},af];
						\[Omega]+=I Sign[Im[edat0[[1,2,1]]]]Exp[ef[Log[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]]]],
						Asymptote,
						{\[Omega],Alm}=AsymptoteFunction[s,l,m,If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]],
						_,
						If[extraporder>2,
							edat0=Take[KerrSEQ,If[forward,-(extraporder+1),extraporder+1]];
							edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
							ef[iv_]=InterpolatingPolynomial[edat,iv];
							\[Omega]=ef[If[forward,KerrSEQ[[NKMode,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]];
						];
					],
					a=If[forward,KerrSEQ[[2,1]],KerrSEQ[[1,1]]];
					\[CapitalDelta]a=KerrSEQ[[2,1]]-KerrSEQ[[1,1]];
					blevelsave=blevel=Round[-(3+Log10[\[CapitalDelta]a])/Log10[2]];
					If[forward,
						\[Omega]=2KerrSEQ[[2,2,1]]-KerrSEQ[[1,2,1]];Alm=2KerrSEQ[[2,3,1]]-KerrSEQ[[1,3,1]],
						\[Omega]=2KerrSEQ[[1,2,1]]-KerrSEQ[[2,2,1]];Alm=2KerrSEQ[[1,3,1]]-KerrSEQ[[2,3,1]]
					]
				]
			],
			(* Unknown Failure *)
			Message[KerrModeSequence::invalidcall];Abort[];
		];
	];
]


AdaptCheck3[KerrTMP_List,inversion_Integer,s_Integer,l_Integer,m_Integer,\[Epsilon]max_Integer,
			relax_Real|relax_Rational|relax_Integer,
			index0_Integer,blevel_Integer,
			maximuma_Rational|maximuma_Integer,minimuma_Rational|minimuma_Integer,
			forward_/;forward \[Element] Booleans,
			incflag_/;incflag \[Element] Booleans,
			recursflag_/;recursflag \[Element] Booleans,
			opts:OptionsPattern[]]:=
Module[{KerrSEQ=KerrTMP,AC3ret,ind0,index0p=index0+1,index0m=index0-1,blevelp=blevel,blevelm=blevel,
		a0,\[Epsilon]=\[Epsilon]max,\[Epsilon]2,\[Epsilon]p,\[Epsilon]m,\[Epsilon]min=-13,Nrcf,Nm,\[CapitalDelta]a,\[CapitalDelta]ap,\[CapitalDelta]am,\[Omega]m,\[Omega]0,\[Omega]p,curvrat,\[CapitalDelta]\[Omega],d\[Omega],dd\[Omega],Almm,Alm0,Almp,\[Omega]g,Almg,
		ModeSol,incres=False,incstep=1,inctmp,
		edat,edat0,ef,iv,afit,extraporder=OptionValue[ExtrapolationOrder],
		Minb=OptionValue[Minblevel],Maxb=OptionValue[Maxblevel],RCFmin=OptionValue[RadialCFMinDepth],
		maxcurvrat=OptionValue[CurvatureRatio],max\[CapitalDelta]\[Omega]=Abs[OptionValue[Max\[CapitalDelta]\[Omega]]],precision=OptionValue[ModePrecision]},
	AdaptCheck3::incorrectam="Incorrect \[CapitalDelta]a(-)";
	AdaptCheck3::incorrectap="Incorrect \[CapitalDelta]a(+)";
	AdaptCheck3::missuse="Cannot use Accumulation extrapolation with backward sequencing";
	AdaptCheck3::apsolfail="a+ solution failed.";
	AdaptCheck3::amsolfail="a- solution failed.";
	AdaptCheck3::appsolfail="a++ solution failed.";
	AdaptCheck3::ammsolfail="a-- solution failed.";

	(*If[blevel>Maxb || blevel<Minb,Print["\[CapitalDelta]a level out of bounds."];Abort[]];*)
	a0=KerrSEQ[[index0,1]];
	\[CapitalDelta]a=2^(-blevel)/1000;
	If[Not[forward] && a0==\[CapitalDelta]a,Return[{KerrSEQ,blevel,incres,incstep,\[Epsilon]}]]; (* Don't step past a=0 *)
	Nm=KerrSEQ[[index0,3,2]];
	If[incflag,
		If[forward,
			While[KerrSEQ[[index0,1]]-KerrSEQ[[index0m,1]]<\[CapitalDelta]a,--index0m;++incstep],
			While[KerrSEQ[[index0p,1]]-KerrSEQ[[index0,1]]<\[CapitalDelta]a,++index0p;++incstep]
		]
	];
	If[KerrSEQ[[index0,1]]-KerrSEQ[[index0m,1]]!=\[CapitalDelta]a,Message[AdaptCheck3::incorrectam];Abort[]];
	If[KerrSEQ[[index0p,1]]-KerrSEQ[[index0,1]]!=\[CapitalDelta]a,Message[AdaptCheck3::incorrectap];Abort[]];
	\[Omega]0=KerrSEQ[[index0,2,1]];
	\[Omega]p=KerrSEQ[[index0p,2,1]];
	\[Omega]m=KerrSEQ[[index0m,2,1]];
	Alm0=KerrSEQ[[index0,3,1]];
	Almp=KerrSEQ[[index0p,3,1]];
	Almm=KerrSEQ[[index0m,3,1]];
	index0p=index0m=index0; (* reset to central index *)
	blevelp=blevelm=blevel;
	d\[Omega]=\[Omega]p-\[Omega]m; dd\[Omega]=\[Omega]p-2\[Omega]0+\[Omega]m;
	curvrat=4Sqrt[Chop[Abs[d\[Omega]]^2Abs[dd\[Omega]]^2-(Re[d\[Omega]]Re[dd\[Omega]]+Im[d\[Omega]]Im[dd\[Omega]])^2]]/(Abs[d\[Omega]]^2);
	\[CapitalDelta]\[Omega]=Abs[\[Omega]0-If[forward,\[Omega]p,\[Omega]m]];
	\[Epsilon]p=\[Epsilon]m=\[Epsilon]=Max[Min[\[Epsilon]max,Floor[Log10[Abs[\[CapitalDelta]\[Omega]]]-2.5]],\[Epsilon]min];
	\[Epsilon]2=Max[Min[\[Epsilon]max,Floor[Log10[Abs[\[CapitalDelta]\[Omega]/2]]-2.5]],\[Epsilon]min];
	If[blevel < Minb || 
		(Maxb > blevel && (curvrat > maxcurvrat || \[CapitalDelta]\[Omega] > max\[CapitalDelta]\[Omega] || 
		(forward && a0+\[CapitalDelta]a<maximuma && a0+2\[CapitalDelta]a>maximuma) || (Not[forward] && a0-\[CapitalDelta]a>minimuma && a0-2\[CapitalDelta]a<minimuma))),
		If[incflag,incstep=1]; (* shouldn't have increased step size *)
		(* Increase resolution *)
		If[forward || (!forward && !incflag),
			(* Compute solution at a+ *)
			\[Omega]g = (6\[Omega]0+3\[Omega]p-\[Omega]m)/8;
			Almg = (6Alm0+3Almp-Almm)/8;
			Switch[extraporder,
				Accumulate,
				If[!forward,
					Message[AdaptCheck3::missuse];
					Abort[]
				];
				edat0=SetPrecision[Take[KerrSEQ,{index0p-9,index0p}],Max[precision,$MinPrecision]];
				edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
				afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
											{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
				(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
				\[Omega]g=afit[1-(KerrSEQ[[index0,1]]+\[CapitalDelta]a/2)];
				edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
				afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
											{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
				\[Omega]g+=I afit[1-(KerrSEQ[[index0,1]]+\[CapitalDelta]a/2)],
				LogLog,
				edat0=Take[KerrSEQ,If[forward,-2,2]];
				edat=Table[{Log[edat0[[i,1]]],Log[Abs[Re[edat0[[i,2,1]]]]]},{i,1,2}];
				ef=LinearModelFit[edat,{1,a},a];
				\[Omega]g=Sign[Re[edat0[[1,2,1]]]]Exp[ef[Log[KerrSEQ[[index0,1]]+\[CapitalDelta]a/2]]];
				edat=Table[{Log[edat0[[i,1]]],Log[Abs[Im[edat0[[i,2,1]]]]]},{i,1,2}];
				ef=LinearModelFit[edat,{1,a},a];
				\[Omega]g+=I Sign[Im[edat0[[1,2,1]]]]Exp[ef[Log[KerrSEQ[[index0,1]]+\[CapitalDelta]a/2]]],
				Asymptote,
				{\[Omega]g,Almg}=AsymptoteFunction[s,l,m,KerrSEQ[[index0,1]]+\[CapitalDelta]a/2],
				_,
				If[extraporder>2,
						edat0=Take[KerrSEQ,If[forward,{index0p-extraporder,index0p},{index0m,index0m+extraporder}]];
						edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
						ef[iv_]=InterpolatingPolynomial[edat,iv];
						\[Omega]g=ef[KerrSEQ[[index0,1]]+\[CapitalDelta]a/2];
				];
			];
			Nrcf=Max[If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]],
						If[Length[KerrSEQ[[index0p,2]]]>=6,KerrSEQ[[index0p,2,6]],KerrSEQ[[index0p,2,3]]]];
			Nrcf=Max[Nrcf,RCFmin];
			ModeSol=ModeSolution[inversion,s,l,m,a0+\[CapitalDelta]a/2,
								SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
								SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon]2,relax,
								Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
			If[ModeSol[[1]],(* valid solution *)
				Print["ModeSol+ a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision]];
				index0p=index0+1;
				blevelp=blevel+1;
				KerrSEQ=Insert[KerrSEQ,ModeSol[[4]],index0p];
				AC3ret=AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon]2,relax,index0p,blevelp,maximuma,minimuma,forward,False,True,FilterRules[{opts},Options[AdaptCheck3]]];
				KerrSEQ=AC3ret[[1]];
				blevelp=AC3ret[[2]];
				\[Epsilon]p=AC3ret[[5]];
				,(* invalid solution *)
				Message[AdaptCheck3::apsolfail];
				Abort[];
			];
		];
		If[!forward || (forward && !incflag),
			(* Compute solution at a- *)
			\[Omega]g = (6\[Omega]0+3\[Omega]m-\[Omega]p)/8;
			Almg = (6Alm0+3Almm-Almp)/8;
			Switch[extraporder,
				Accumulate,
				If[!forward,
					Message[AdaptCheck3::missuse];
					Abort[]
				];
				edat0=SetPrecision[Take[KerrSEQ,{index0p-9,index0p}],Max[precision,$MinPrecision]];
				edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
				afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
											{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
				(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
				\[Omega]g=afit[1-(KerrSEQ[[index0,1]]-\[CapitalDelta]a/2)];
				edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
				afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
											{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
				\[Omega]g+=I afit[1-(KerrSEQ[[index0,1]]-\[CapitalDelta]a/2)],
				LogLog,
				edat0=Take[KerrSEQ,If[forward,-2,2]];
				edat=Table[{Log[edat0[[i,1]]],Log[Abs[Re[edat0[[i,2,1]]]]]},{i,1,2}];
				ef=LinearModelFit[edat,{1,a},a];
				\[Omega]g=Sign[Re[edat0[[1,2,1]]]]Exp[ef[Log[KerrSEQ[[index0,1]]-\[CapitalDelta]a/2]]];
				edat=Table[{Log[edat0[[i,1]]],Log[Abs[Im[edat0[[i,2,1]]]]]},{i,1,2}];
				ef=LinearModelFit[edat,{1,a},a];
				\[Omega]g+=I Sign[Im[edat0[[1,2,1]]]]Exp[ef[Log[KerrSEQ[[index0,1]]-\[CapitalDelta]a/2]]],
				Asymptote,
				{\[Omega]g,Almg}=AsymptoteFunction[s,l,m,KerrSEQ[[index0,1]]-\[CapitalDelta]a/2],
				_,
				If[extraporder>2,
						edat0=Take[KerrSEQ,If[forward,{index0p-extraporder,index0p},{index0m,index0m+extraporder}]];
						edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
						ef[iv_]=InterpolatingPolynomial[edat,iv];
						\[Omega]g=ef[KerrSEQ[[index0,1]]-\[CapitalDelta]a/2];
				];
			];
			Nrcf=Max[If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]],
						If[Length[KerrSEQ[[index0m,2]]]>=6,KerrSEQ[[index0m,2,6]],KerrSEQ[[index0m,2,3]]]];
			Nrcf=Max[Nrcf,RCFmin];
			ModeSol=ModeSolution[inversion,s,l,m,a0-\[CapitalDelta]a/2,
								SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
								SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon]2,relax,
								Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
			If[ModeSol[[1]],(* valid solution *)
				Print["ModeSol- a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision]];
				blevelm=blevel+1;
				KerrSEQ=Insert[KerrSEQ,ModeSol[[4]],index0];
				AC3ret=AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon]2,relax,index0m,blevelm,maximuma,minimuma,forward,False,True,FilterRules[{opts},Options[AdaptCheck3]]];
				KerrSEQ=AC3ret[[1]];
				blevelm=AC3ret[[2]];
				\[Epsilon]m=AC3ret[[5]];
				,(* invalid solution *)
				Message[AdaptCheck3::amsolfail];
				Abort[];
			];
		];
		(* get \[CapitalDelta]a left and right of index0.  If off by more than 2x, refine once to help smooth near rapid changes *)
		ind0=index0;
		While[KerrSEQ[[ind0,1]]>a0,--ind0];
		While[KerrSEQ[[ind0,1]]<a0,++ind0];
		\[CapitalDelta]ap = KerrSEQ[[ind0+1,1]]-KerrSEQ[[ind0,1]];
		\[CapitalDelta]am = KerrSEQ[[ind0,1]]-KerrSEQ[[ind0-1,1]];
		While[\[CapitalDelta]ap>2\[CapitalDelta]am || \[CapitalDelta]am>2\[CapitalDelta]ap,
			If[\[CapitalDelta]ap>2\[CapitalDelta]am,
				\[Omega]g = (KerrSEQ[[ind0+1,2,1]]+KerrSEQ[[ind0,2,1]])/2;
				Almg = (KerrSEQ[[ind0+1,3,1]]+KerrSEQ[[ind0,3,1]])/2;
				Nrcf=Max[If[Length[KerrSEQ[[ind0,2]]]>=6,KerrSEQ[[ind0,2,6]],KerrSEQ[[ind0,2,3]]],
							If[Length[KerrSEQ[[ind0+1,2]]]>=6,KerrSEQ[[ind0+1,2,6]],KerrSEQ[[ind0+1,2,3]]]];
				Nrcf=Max[Nrcf,RCFmin];
				\[Epsilon]2=Max[Min[\[Epsilon]max,Floor[Log10[Abs[(KerrSEQ[[ind0+1,2,1]]+KerrSEQ[[ind0,2,1]])/2]]-2.5]],\[Epsilon]min];
				ModeSol=ModeSolution[inversion,s,l,m,a0+\[CapitalDelta]ap/2,
									SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
									SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon]2,relax,
									Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
				If[ModeSol[[1]],
					Print["ModeSol++ a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision]];
					KerrSEQ=Insert[KerrSEQ,ModeSol[[4]],ind0+1];
					,(* invalid solution *)
					Message[AdaptCheck3::appsolfail];
					Abort[];
				];
			,
			If[\[CapitalDelta]am>2\[CapitalDelta]ap,
				\[Omega]g = (KerrSEQ[[ind0,2,1]]+KerrSEQ[[ind0-1,2,1]])/2;
				Almg = (KerrSEQ[[ind0,3,1]]+KerrSEQ[[ind0-1,3,1]])/2;
				Nrcf=Max[If[Length[KerrSEQ[[ind0,2]]]>=6,KerrSEQ[[ind0,2,6]],KerrSEQ[[ind0,2,3]]],
							If[Length[KerrSEQ[[ind0-1,2]]]>=6,KerrSEQ[[ind0-1,2,6]],KerrSEQ[[ind0-1,2,3]]]];
				Nrcf=Max[Nrcf,RCFmin];
				\[Epsilon]2=Max[Min[\[Epsilon]max,Floor[Log10[Abs[(KerrSEQ[[ind0,2,1]]-KerrSEQ[[ind0-1,2,1]])/2]]-2.5]],\[Epsilon]min];
				ModeSol=ModeSolution[inversion,s,l,m,a0-\[CapitalDelta]am/2,
									SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
									SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon]2,relax,
									Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
				If[ModeSol[[1]],
					Print["ModeSol-- a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision]];
					KerrSEQ=Insert[KerrSEQ,ModeSol[[4]],ind0++]; (* must increment to keep ind0 at same a *)
					,(* invalid solution *)
					Message[AdaptCheck3::ammsolfail];
					Abort[];
				];
			]];
			\[CapitalDelta]ap = KerrSEQ[[ind0+1,1]]-KerrSEQ[[ind0,1]];
			\[CapitalDelta]am = KerrSEQ[[ind0,1]]-KerrSEQ[[ind0-1,1]];
		];
		,
		If[Not[recursflag] && 
			(blevel > Maxb || (Minb < blevel && curvrat < maxcurvrat/2 && \[CapitalDelta]\[Omega] < max\[CapitalDelta]\[Omega]/2
								&& ((forward && a0+3\[CapitalDelta]a<=maximuma) || (Not[forward] && a0-3\[CapitalDelta]a>=minimuma)))),
			(* Reduce resoltuion *)
			If[Mod[1000(a0+\[CapitalDelta]a),2000\[CapitalDelta]a]==0, (* Make sure we can end up one sensible values of a *)
				blevelm=blevelp=blevel-1;
				incres=True;
			]
		]
	];
(*
	KerrSEQ: full sequence list [must be saved at some point]
	blevel[p,m]: blevel to use for next extrapolated step
	incres: true if next step will have larger step size
	incstep: for a successful increas in step size, offset to [+/-] data
*)
	If[forward,
		{KerrSEQ,blevelp,incres,incstep,\[Epsilon]p},
		{KerrSEQ,blevelm,incres,incstep,\[Epsilon]m}
	]
]


Options[KerrModeRefineSequence]=Union[{SpinWeight->Null[],Index->False,
								Refinement->All,RefinementAction->None,ForceRefinement->False,
								RefinementPlot->SeqLevel,LimitRefinement->None,
								SolutionRelax->1,RadialCFDepth->1,RadialCFMaxGuess->20000000},
								Options[AdaptCheck3],Options[ListLinePlot]];


KerrModeRefineSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]max_Integer,
				opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],SpinWeightTable,KerrSEQ,
		NKMode,index0,index0m,index0p,alow,ahigh,refinechange=False,
		i,plotdata,plotdata1,\[Omega]dat,d\[Omega],dd\[Omega],ind0,a0,\[CapitalDelta]ap,\[CapitalDelta]am,\[Omega]g,Almg,ModeSol,inversion,\[Epsilon]=\[Epsilon]max,Nrcf,Nm,
		KerrSEQret,dummy,blevel,forward,incflag,limitlist={},ll,inc,dec,last,re\[Omega],width,
		indexmin,indexmax,offset,oldNrcf,newNrcf,oldCf,newCf,rcfmin,ref\[Epsilon],
		saveSeq,terminalSeq,oldbstruct,newbstruct,addlist,removelist,nearesta,
		takerange,edat0,edat,ef,a,afit,alist,lev,
		useindex=OptionValue[Index],
		precision=OptionValue[ModePrecision],refinement=OptionValue[Refinement],
		action=OptionValue[RefinementAction],forcerefine=OptionValue[ForceRefinement],
		plottype=OptionValue[RefinementPlot],limitrefine=OptionValue[LimitRefinement],
		Minb=OptionValue[Minblevel],Maxb=OptionValue[Maxblevel],
		relax=Rationalize[OptionValue[SolutionRelax]],
		RCFmin=OptionValue[RadialCFMinDepth],RCFmax=OptionValue[RadialCFMaxGuess],
		rcfdepth=OptionValue[RadialCFDepth],extraporder=OptionValue[ExtrapolationOrder]},
	KerrModeRefineSequence::Refinement="The value of Refinement (`1`) is not an integer index, real value for a, a list of either specifying a range, or ALL.";
	KerrModeRefineSequence::index="using index `1` instead of `2`";
	KerrModeRefineSequence::value="using a = `1` instead of `2`";
	KerrModeRefineSequence::range="using range `1` instead of `2`";
	KerrModeRefineSequence::list="`1` must be a 2 element list of either integers or reals";
	KerrModeRefineSequence::invalidlist="`1` is an invalid range of elements";
	KerrModeRefineSequence::sequence="Sequence has `1` elements; must have at least 3 to refine";
	KerrModeRefineSequence::badplot="Invalid RefinementPlot `1` given";
	KerrModeRefineSequence::badaction="Invalid RefinementAction `1` given";
	KerrModeRefineSequence::limits="Invalid LimitRefinement `1` given";
	KerrModeRefineSequence::precision="Set $MinPrecision = `1`";
	KerrModeRefineSequence::largedepth="Warning: Computed radial CF depth too large.";
	KerrModeRefineSequence::setdepth="         Setting CF depth to `1`";
	KerrModeRefineSequence::indexfail="Solution failed at index `1`";
	KerrModeRefineSequence::TerminalLoRange="Refinement `1` does not start at beginning of sequence";
	KerrModeRefineSequence::TerminalHiRange="Refinement `1` does not end at end of sequence";
	KerrModeRefineSequence::TerminalLowFailure="Interpolated solution at `1` failed";
	KerrModeRefineSequence::TerminalHiFailure="Interpolated solution at `1` failed";
	SpinWeightTable:=modeName;
	KerrSEQ:=modeName[l,m,n];
	NKMode=Length[KerrSEQ];
	If[NKMode<3,Message[KerrModeRefineSequence::sequence,NKMode];Return[]];
	If[action==RefineAccuracy || action==RefinePrecision || action==Update || action==None || action==TerminalDescent || action==TerminalAscent,
		indexmin=1;indexmax=NKMode;offset=0,
		indexmin=2;indexmax=NKMode-1;offset=1,
		indexmin=2;indexmax=NKMode-1;offset=1
	];
(* Parse Refinement option to set range of values to refine *)
	Switch[refinement
		,_Symbol,
			Switch[refinement
				,All,
					index0m=1;
					index0p=NKMode;
					alow=KerrSEQ[[1,1]];
					ahigh=KerrSEQ[[NKMode,1]];
					index0=indexmax;
				,_,Message[KerrModeRefineSequence::Refinement,refinement];Return[]
			];
		,_Integer,
			index0=refinement;
			If[index0<0,index0+=NKMode+1];
			If[index0<indexmin,refinechange=True;index0=indexmin];
			If[index0>indexmax,refinechange=True;index0=indexmax];
			If[refinechange,If[refinement>=0,
				Print[Style[StringForm[KerrModeRefineSequence::index,index0,refinement],{Medium,Darker[Green]}]],
				Print[Style[StringForm[KerrModeRefineSequence::index,index0-NKMode-1,refinement],{Medium,Darker[Green]}]]
			]];
			index0m=index0-offset;
			index0p=index0+offset;
			alow=KerrSEQ[[index0m,1]];ahigh=KerrSEQ[[index0p,1]];
		,_Rational|_Real,
			index0=1;
			While[KerrSEQ[[index0,1]]<refinement && index0<NKMode,++index0];
			If[index0<indexmin,refinechange=True;index0=indexmin];
			If[index0>indexmax,refinechange=True;index0=indexmax];
			If[KerrSEQ[[index0,1]]!=refinement,
				refinechange=True;
				If[index0<NKMode-1 &&
					(Abs[KerrSEQ[[index0,1]]-refinement]>Abs[KerrSEQ[[index0+1,1]]-refinement]),
						++index0];
			];
			If[refinechange,Print[Style[StringForm[KerrModeRefineSequence::value,Block[{$MinPrecision=0},N[KerrSEQ[[index0,1]],{Infinity,20}]],refinement]],{Medium,Darker[Green]}]];
			index0m=index0-offset;
			index0p=index0+offset;
			alow=KerrSEQ[[index0m,1]];ahigh=KerrSEQ[[index0p,1]];
		,_List,
			If[Length[refinement]!=2,Message[KerrModeRefineSequence::list,refinement];Return[]];
			Switch[refinement[[1]]
				,_Integer,
					index0m=refinement[[1]];index0p=refinement[[2]];
					If[Head[index0p]==Integer,Null,Null,Message[KerrModeRefineSequence::list,refinement];Return[]]; 
					If[index0m<0,index0m+=NKMode+1];
					If[index0p<0,index0p+=NKMode+1];
					If[index0m>index0p,Message[KerrModeRefineSequence::invalidlist,refinement];Return[]];
					If[index0m<1,refinechange=True;index0m=1];
					If[index0p>NKMode,refinechange=True;index0p=NKMode];
					If[index0m>index0p-offset,refinechange=True;index0m=index0p-offset];
					If[index0p<index0m+offset,refinechange=True;index0p=index0m+offset];
					If[index0m==NKMode,index0m=indexmax-offset];
					If[refinechange,Print[Style[StringForm[KerrModeRefineSequence::range,
						{If[refinement[[1]]<0,index0m-NKMode-1,index0m],
						 If[refinement[[2]]<0,index0p-NKMode-1,index0p]},refinement]],{Medium,Darker[Green]}]];

					alow=KerrSEQ[[index0m,1]];
					ahigh=KerrSEQ[[index0p,1]];
					index0=index0p-1;
				,_Rational|_Real,
					If[Head[refinement[[2]]]==Real||Head[refinement[[2]]]==Rational,Null,Null,Message[KerrModeRefineSequence::list,refinement];Return[]]; 
					If[refinement[[1]]>refinement[[2]],Message[KerrModeRefineSequence::invalidlist,refinement];Return[]];

					index0m=1;
					While[KerrSEQ[[index0m,1]]<refinement[[1]] && index0m<NKMode,++index0m];
					If[KerrSEQ[[index0m,1]]!=refinement[[1]],
						refinechange=True;
						If[index0m<NKMode &&
							(Abs[KerrSEQ[[index0m,1]]-refinement[[1]]]>Abs[KerrSEQ[[index0m+1,1]]-refinement[[1]]]),
								++index0m];
					];
					index0p=1;
					While[KerrSEQ[[index0p,1]]<refinement[[2]] && index0p<NKMode,++index0p];
					If[KerrSEQ[[index0p,1]]!=refinement[[2]],
						refinechange=True;
						If[index0p<NKMode &&
							(Abs[KerrSEQ[[index0p,1]]-refinement[[2]]]>Abs[KerrSEQ[[index0p+1,1]]-refinement[[2]]]),
								++index0p];
					];
					If[index0p>NKMode,refinechange=True;index0p=NKMode];
					If[index0m>index0p-1-offset,refinechange=True;index0m=index0p-1-offset;];
					If[index0p<index0m+1+offset,refinechange=True;index0p=index0m+1+offset;];
					alow=KerrSEQ[[index0m,1]];
					ahigh=KerrSEQ[[index0p,1]];
					If[refinechange,Print[Style[StringForm[KerrModeRefineSequence::range,Block[{$MinPrecision=0},N[{alow,ahigh},{Infinity,20}]],refinement],{Medium,Darker[Green]}]]];
					index0=index0p-1;
				,_,
				Message[KerrModeRefineSequence::list,refinement];Return[]
			];
		,_,Message[KerrModeRefineSequence::Refinement,refinement];Return[] 

	];
	Switch[limitrefine
		,None,
			limitlist={{index0m,index0p}};
		,Minima,
			inc=False;dec=False;
			last=Re[KerrSEQ[[index0m,2,1]]];
			For[i=index0m+1,i<=index0p,++i,
				re\[Omega]=Re[KerrSEQ[[i,2,1]]];
				If[re\[Omega]>=last,
					If[dec && Not[inc],PrependTo[limitlist,{i-1,i-1}]];
					If[Not[inc],inc=True,dec=False],
					(*re\[Omega]<last*)
					If[Not[dec],dec=True,inc=False]
				];
				last=re\[Omega];
			];
		,_List,
			Switch[limitrefine[[1]]
				,Minima,
					Switch[limitrefine[[2]]
						,_Integer,
							width=limitrefine[[2]];
						,_,
							Message[KerrModeRefineSequence::limits,limitrefine];Return[] 

					];
					inc=False;dec=False;
					last=Re[KerrSEQ[[index0m,2,1]]];
					For[i=index0m+1,i<=index0p,++i,
						re\[Omega]=Re[KerrSEQ[[i,2,1]]];
						If[re\[Omega]>=last,
							If[dec && Not[inc],PrependTo[limitlist,{i-1-width,i-1+width}]];
							If[Not[inc],inc=True,dec=False],
							(*re\[Omega]<last*)
							If[Not[dec],dec=True,inc=False]
						];
						last=re\[Omega];
					];
					If[Length[limitrefine]==3,
						rcfmin=limitrefine[[3]];
						ll=limitlist;limitlist={};
						For[i=1,i<=Length[ll],++i,
							PrependTo[limitlist,Table[j,{j,ll[[i,1]],ll[[i,2]]}]];
						];
						ll=Flatten[limitlist];limitlist={};
						For[i=1,i<=Length[ll],++i,
							If[Length[KerrSEQ[[ll[[i]],2]]]>=6 && KerrSEQ[[ll[[i]],2,6]]<=rcfmin,
								PrependTo[limitlist,{ll[[i]],ll[[i]]}];
							];
						];
					];
				,RadialCFMinDepth,
					Switch[limitrefine[[2]]
						,_Integer,
							rcfmin=limitrefine[[2]];
						,_,
							Message[KerrModeRefineSequence::limits,limitrefine];Return[] 

					];
					For[i=index0m,i<=index0p,++i,
						If[Length[KerrSEQ[[i,2]]]>=6 && KerrSEQ[[i,2,6]]<=rcfmin,
							PrependTo[limitlist,{i,i}];
						];
					];
				,_,
					Message[KerrModeRefineSequence::limits,limitrefine];Return[] 
			];
		,_,
			Message[KerrModeRefineSequence::limits,limitrefine];Return[]  
	];
	Switch[plottype
		,None,
			plotdata=None;
		,SeqLevel,
			alist=Flatten[Take[KerrSEQ,All,1]];
			lev=Round[-(3+Log10[Take[Min[#[[1]],#[[2]]]&/@
									Transpose[{alist-Prepend[Drop[alist,-1],-\[Infinity]],
											   Append[Drop[alist,1],\[Infinity]]-Take[alist]}],{index0m,index0p}]])/Log10[2]];
			alist=Take[alist,{index0m,index0p}];
			plotdata=Transpose[{If[useindex,Range[Length[alist]],alist],lev}];
			(*plotdata=Table[{If[useindex,i,KerrSEQ[[i,1]]],Round[-(3+Log10[KerrSEQ[[i+1,1]]-KerrSEQ[[i,1]]])/Log10[2]]},{i,Min[NKMode-1,index0m],Min[NKMode-1,index0p]}];*)
		,RadialCFLevel,
			plotdata=Table[{If[useindex,i,KerrSEQ[[i,1]]],If[Length[KerrSEQ[[i,2]]]>=6,KerrSEQ[[i,2,6]],KerrSEQ[[i,2,3]],0]},{i,Min[NKMode-1,index0m],Min[NKMode-1,index0p]}];
		,AccuracyLevel,
			plotdata=Table[{If[useindex,i,KerrSEQ[[i,1]]],KerrSEQ[[i,2,4]]},{i,Max[1,index0m],Min[NKMode,index0p]}];
		,PrecisionLevel,
			plotdata=Table[{If[useindex,i,KerrSEQ[[i,1]]],MyPrecision[KerrSEQ[[i,2,1]]]},{i,Max[1,index0m],Min[NKMode,index0p]}];
		,StepRatio,
			plotdata1=Table[{KerrSEQ[[i,1]],i},{i,index0m,index0p}];
			plotdata=RotateLeft[plotdata1,1]-plotdata1;
			plotdata=RotateRight[plotdata,1]/plotdata;
			plotdata=Table[{If[useindex,plotdata1[[i,2]],plotdata1[[i,1]]],plotdata[[i,1]]},{i,2,Length[plotdata1]-1}];
		,CurveRatio,
			plotdata1=Table[{KerrSEQ[[i,1]],i},{i,index0m,index0p}];
			plotdata=\[Omega]dat=Table[KerrSEQ[[i,2,1]],{i,index0m,index0p}];
			For[i=2,i<=Length[plotdata1]-1,++i,
				If[plotdata1[[i+1,1]]-plotdata1[[i,1]]==plotdata1[[i,1]]-plotdata1[[i-1,1]] || i==2 || i==Length[plotdata1]-1,
					d\[Omega]=\[Omega]dat[[i+1]]-\[Omega]dat[[i-1]];
					dd\[Omega]=\[Omega]dat[[i+1]]-2\[Omega]dat[[i]]+\[Omega]dat[[i-1]];
				,If[plotdata1[[i+1,1]]-plotdata1[[i,1]]>plotdata1[[i,1]]-plotdata1[[i-1,1]],
					d\[Omega]=\[Omega]dat[[i+1]]-\[Omega]dat[[i-2]];
					dd\[Omega]=\[Omega]dat[[i+1]]-2\[Omega]dat[[i]]+\[Omega]dat[[i-2]];
				,
					d\[Omega]=\[Omega]dat[[i+2]]-\[Omega]dat[[i-1]];
					dd\[Omega]=\[Omega]dat[[i+2]]-2\[Omega]dat[[i]]+\[Omega]dat[[i-1]];
				]];
				plotdata[[i]]=4Sqrt[Abs[d\[Omega]]^2Abs[dd\[Omega]]^2-(Re[d\[Omega]]Re[dd\[Omega]]+Im[d\[Omega]]Im[dd\[Omega]])^2]/(Abs[d\[Omega]]^2);
			];
			plotdata=Table[{If[useindex,plotdata1[[i,2]],plotdata1[[i,1]]],plotdata[[i]]},{i,2,Length[plotdata1]-1}];
		,_,Message[KerrModeRefineSequence::badplot,plottype];Return[]
	];
	If[plottype==SeqLevel || plottype==RadialCFLevel || plottype==AccuracyLevel || plottype==PrecisionLevel || plottype==StepRatio || plottype==CurveRatio,
		Print[ListLinePlot[plotdata,FilterRules[{opts},Options[ListLinePlot]]]];
	];
	Switch[action
		,None,Null
		,RefineAccuracy,
			If[forcerefine && precision!=$MinPrecision,Print[Style[StringForm[KerrModeRefineSequence::precision,precision],{Medium,Darker[Red]}]]];
			$MinPrecision=precision;
			For[index0=index0m,index0<=index0p,++index0,
				While[index0>limitlist[[-1,2]],
					limitlist=Drop[limitlist,-1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0p+1;Continue[]];
				If[(\[Epsilon]<KerrSEQ[[index0,2,4]] || forcerefine) && (limitlist[[-1,1]]<=index0<=limitlist[[-1,2]]),
Print["RefineAcc at : ",index0];
					inversion=KerrSEQ[[index0,2,2]];
					oldNrcf=If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]];
(*Print["Old Nrcf : ",oldNrcf];*)
					If[RunCFConvergence, (* only extrapolate Nrcf is using CF *)
						Nrcf=If[Length[KerrSEQ[[index0,2]]]>=8,
							If[NumberQ[KerrSEQ[[index0,2,7]]],
								IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/KerrSEQ[[index0,2,7]])],
								oldNrcf],
							IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/3)]];
(*Print["New Nrcf : ",Nrcf];*)
						If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
						If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
						Nrcf=Max[Nrcf,RCFmin];
						If[Nrcf>RCFmax,
							Print[Style[StringForm[KerrModeRefineSequence::largedepth,Nrcf],{Medium,Darker[Red]}]];
							Print[Style[StringForm[KerrModeRefineSequence::setdepth,RCFmax],{Medium,Darker[Red]}]];
							Nrcf=RCFmax;
						],
						Nrcf=oldNrcf;
					];
					Nm=KerrSEQ[[index0,3,2]];
					a0=KerrSEQ[[index0,1]];
					\[Omega]g=KerrSEQ[[index0,2,1]];
					Almg=KerrSEQ[[index0,3,1]];
					If[!forcerefine,
						$MinPrecision=Max[precision,
							If[Length[KerrSEQ[[index0,2]]]>=9,KerrSEQ[[index0,2,9]],IntegerPart[MyPrecision[\[Omega]g]]]]
					];
					ModeSol=ModeSolution[inversion,s,l,m,a0,
										SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
										SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon],relax,
										Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
					If[ModeSol[[1]],
						Print["ModeSol a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision],
							" Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision],
							"  |\[CapitalDelta]\[Omega]| = ",SetPrecision[Abs[\[Omega]g-ModeSol[[4,2,1]]],MachinePrecision]];
						(*
						oldCf=RadialCFRem[inversion,s,m,a0,Almg,\[Omega]g,oldNrcf];
						newNrcf=ModeSol[[4,2,3]];(* Must be the value used, not the "best" *)
						newCf=RadialCFRem[inversion,s,m,a0,Almg,\[Omega]g,newNrcf];
						Print["Prior Accuracy : ",KerrSEQ[[index0,2,4]]," Pred|\[CapitalDelta]\[Omega]| = ",1/Sqrt[Det[ModeSol[[5]]]]Abs[newCf[[1]]-oldCf[[1]]] ];
						*)
						modeName[l,m,n]=ReplacePart[KerrSEQ,index0->ModeSol[[4]]];
						(*Switch[s,
							   -2,Global`KerrQNM[l,m,n]=ReplacePart[KerrSEQ,index0\[Rule]ModeSol[[4]]],
							   -1,Global`KerrQNMe[l,m,n]=ReplacePart[KerrSEQ,index0\[Rule]ModeSol[[4]]],
								0,Global`KerrQNMs[l,m,n]=ReplacePart[KerrSEQ,index0\[Rule]ModeSol[[4]]]
							  ];*)
						,(* invalid solution *)
						Print[Style[StringForm[KerrModeRefineSequence::indexfail,index0],{Medium,Darker[Red]}]];
					];
				];
			];
		,RefinePrecision,
			If[precision!=$MinPrecision,Print[Style[StringForm[KerrModeRefineSequence::precision,precision],{Medium,Darker[Red]}]]];
			$MinPrecision=precision;
			For[index0=index0m,index0<=index0p,++index0,
				While[index0>limitlist[[-1,2]],
					limitlist=Drop[limitlist,-1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0p+1;Continue[]];
				If[(precision>MyPrecision[KerrSEQ[[index0,2,1]]] || forcerefine) && (limitlist[[-1,1]]<=index0<=limitlist[[-1,2]]),
					inversion=KerrSEQ[[index0,2,2]];
					oldNrcf=If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]];
					If[RunCFConvergence, (* only extrapolate Nrcf is using CF *)
						Nrcf=If[Length[KerrSEQ[[index0,2]]]>=8,
							If[NumberQ[KerrSEQ[[index0,2,7]]],
								IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/KerrSEQ[[index0,2,7]])],
								oldNrcf],
							IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/3)]];
						If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
						If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
						Nrcf=Max[Nrcf,RCFmin];
						If[Nrcf>RCFmax,
							Print[Style[StringForm[KerrModeRefineSequence::largedepth,Nrcf],{Medium,Darker[Red]}]];
							Print[Style[StringForm[KerrModeRefineSequence::setdepth,RCFmax],{Medium,Darker[Red]}]];
							Nrcf=RCFmax;
						],
						Nrcf=oldNrcf;
					];
					Nm=KerrSEQ[[index0,3,2]];
					a0=KerrSEQ[[index0,1]];
					\[Omega]g=KerrSEQ[[index0,2,1]];
					Almg=KerrSEQ[[index0,3,1]];
					ref\[Epsilon]=If[forcerefine,\[Epsilon],Min[\[Epsilon],KerrSEQ[[index0,2,4]]]];
					If[!forcerefine,
						$MinPrecision=Max[precision,
							If[Length[KerrSEQ[[index0,2]]]>=9,KerrSEQ[[index0,2,9]],IntegerPart[MyPrecision[\[Omega]g]]]]
					];
					ModeSol=ModeSolution[inversion,s,l,m,a0,
										SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
										SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
										Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
					If[ModeSol[[1]],
						Print["ModeSol a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision],
							" Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision],
							"  |\[CapitalDelta]\[Omega]| = ",SetPrecision[Abs[\[Omega]g-ModeSol[[4,2,1]]],MachinePrecision]];
						modeName[l,m,n]=ReplacePart[KerrSEQ,index0->ModeSol[[4]]];
						(*Switch[s,
							   -2,Global`KerrQNM[l,m,n]=ReplacePart[KerrSEQ,index0\[Rule]ModeSol[[4]]],
							   -1,Global`KerrQNMe[l,m,n]=ReplacePart[KerrSEQ,index0\[Rule]ModeSol[[4]]],
								0,Global`KerrQNMs[l,m,n]=ReplacePart[KerrSEQ,index0\[Rule]ModeSol[[4]]]
							  ];*)
						,(* invalid solution *)
						Print[Style[StringForm[KerrModeRefineSequence::indexfail,index0],{Medium,Darker[Red]}]];
					];
				];
			];
		,RefineAdapt,
			If[precision!=$MinPrecision,Print[Style[StringForm[KerrModeRefineSequence::precision,precision],{Medium,Darker[Red]}]]];
			$MinPrecision=precision;
			incflag=False;
			While[index0>index0m,
				While[index0<limitlist[[1,1]],
					limitlist=Drop[limitlist,1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0m;Continue[]];
				If[index0>limitlist[[1,2]],--index0;Continue[]];
				inversion=KerrSEQ[[index0,2,2]];
				Nrcf=If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]];
				If[RunCFConvergence, (* only update Nrcf is using CF *)
					If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
					If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
					Nrcf=Max[Nrcf,RCFmin];
					If[Nrcf>RCFmax,
						Print[Style[StringForm[KerrModeRefineSequence::largedepth,Nrcf],{Medium,Darker[Red]}]];
						Print[Style[StringForm[KerrModeRefineSequence::setdepth,RCFmax],{Medium,Darker[Red]}]];
						Nrcf=RCFmax;
					];
				];
				Nm=KerrSEQ[[index0,3,2]];
				blevel=Round[-(3+Log10[KerrSEQ[[index0+1,1]]-KerrSEQ[[index0,1]]])/Log10[2]];
				forward=True;
				incflag=False;
				If[KerrSEQ[[index0+1,1]]-KerrSEQ[[index0,1]]>KerrSEQ[[index0,1]]-KerrSEQ[[index0-1,1]],
					incflag=True,
					forward=False;
					blevel=Round[-(3+Log10[KerrSEQ[[index0,1]]-KerrSEQ[[index0-1,1]]])/Log10[2]];
					If[KerrSEQ[[index0+1,1]]-KerrSEQ[[index0,1]]<KerrSEQ[[index0,1]]-KerrSEQ[[index0-1,1]],
						incflag=True;
					];
				];
				\[Epsilon]=Min[\[Epsilon]max,KerrSEQ[[index0,2,4]]];
				$MinPrecision=If[Length[KerrSEQ[[index0,2]]]>=9,KerrSEQ[[index0,2,9]],IntegerPart[MyPrecision[KerrSEQ[[index0,2,1]]]]];
(*Print["At a = ",N[KerrSEQ[[index0,1]]]," index0 = ",index0," forward : ",forward," incflag : ",incflag];*)
				{KerrSEQret,blevel,dummy,dummy,\[Epsilon]}=
					AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,1,0,forward,incflag,False,FilterRules[{opts},Options[AdaptCheck3]]];
				modeName[l,m,n]=KerrSEQret;
				(*Switch[s,
					   -2,Global`KerrQNM[l,m,n]=KerrSEQret,
					   -1,Global`KerrQNMe[l,m,n]=KerrSEQret,
						0,Global`KerrQNMs[l,m,n]=KerrSEQret
					  ];*)
				incflag=True;
				--index0;
			];
		,FixAdapt,
			If[precision!=$MinPrecision,Print[Style[StringForm[KerrModeRefineSequence::precision,precision],{Medium,Darker[Red]}]]];
			$MinPrecision=precision;
			plotdata1=Table[KerrSEQ[[i,1]],{i,index0m,index0p}];
			plotdata=RotateLeft[plotdata1,1]-plotdata1;
			plotdata=RotateRight[plotdata,1]/plotdata;
			For[i=Length[plotdata1]-1,i>1,--i,
				ind0=index0m+i-1;
				While[ind0<limitlist[[1,1]],
					limitlist=Drop[limitlist,1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,i=1;Continue[]];
				If[ind0>limitlist[[1,2]],Continue[]];
				If[Not[1/2<=plotdata[[i]]<=2],
					a0=plotdata1[[i]];
					inversion=KerrSEQ[[ind0,2,2]];
					Nrcf=If[Length[KerrSEQ[[ind0,2]]]>=6,KerrSEQ[[ind0,2,6]],KerrSEQ[[ind0,2,3]]];
					If[RunCFConvergence, (* only update Nrcf is using CF *)
						If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
						If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
						Nrcf=Max[Nrcf,RCFmin];
						If[Nrcf>RCFmax,
							Print[Style[StringForm[KerrModeRefineSequence::largedepth,Nrcf],{Medium,Darker[Red]}]];
							Print[Style[StringForm[KerrModeRefineSequence::setdepth,RCFmax],{Medium,Darker[Red]}]];
							Nrcf=RCFmax;
						];
					];
					Nm=KerrSEQ[[ind0,3,2]];
					Print["Adaptation error at a = ",N[plotdata1[[i]]]," : ratio = ",plotdata[[i]]];
					\[CapitalDelta]ap=plotdata1[[i+1]]-plotdata1[[i]];
					\[CapitalDelta]am=plotdata1[[i]]-plotdata1[[i-1]];
					While[\[CapitalDelta]ap>2\[CapitalDelta]am || \[CapitalDelta]am>2\[CapitalDelta]ap,
						If[\[CapitalDelta]ap>2\[CapitalDelta]am,
							\[Omega]g=(KerrSEQ[[ind0+1,2,1]]+KerrSEQ[[ind0,2,1]])/2;
							Almg = (KerrSEQ[[ind0+1,3,1]]+KerrSEQ[[ind0,3,1]])/2;
							ref\[Epsilon]=If[forcerefine,\[Epsilon],Min[\[Epsilon],Max[KerrSEQ[[index0,2,4]],KerrSEQ[[index0+1,2,4]]]]];
							$MinPrecision=If[Length[KerrSEQ[[ind0,2]]]>=9,KerrSEQ[[ind0,2,9]],IntegerPart[MyPrecision[KerrSEQ[[ind0,2,1]]]]];
							ModeSol=ModeSolution[inversion,s,l,m,a0+\[CapitalDelta]ap/2,
												SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
												SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
											Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
							If[ModeSol[[1]],
								Print["ModeSol++ a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision]];
								modeName[l,m,n]=Insert[KerrSEQ,ModeSol[[4]],ind0+1];
								(*Switch[s,
								   -2,Global`KerrQNM[l,m,n] =Insert[KerrSEQ,ModeSol[[4]],ind0+1],
								   -1,Global`KerrQNMe[l,m,n]=Insert[KerrSEQ,ModeSol[[4]],ind0+1],
									0,Global`KerrQNMs[l,m,n]=Insert[KerrSEQ,ModeSol[[4]],ind0+1]
								  ];*)
								,(* invalid solution *)
								Print["a++ solution failed."];
								Abort[];
							];
						,If[\[CapitalDelta]am>2\[CapitalDelta]ap,
							\[Omega]g = (KerrSEQ[[ind0,2,1]]+KerrSEQ[[ind0-1,2,1]])/2;
							Almg = (KerrSEQ[[ind0,3,1]]+KerrSEQ[[ind0-1,3,1]])/2;
							ref\[Epsilon]=If[forcerefine,\[Epsilon],Min[\[Epsilon],Max[KerrSEQ[[index0,2,4]],KerrSEQ[[index0-1,2,4]]]]];
							$MinPrecision=If[Length[KerrSEQ[[ind0,2]]]>=9,KerrSEQ[[ind0,2,9]],IntegerPart[MyPrecision[KerrSEQ[[ind0,2,1]]]]];
							ModeSol=ModeSolution[inversion,s,l,m,a0-\[CapitalDelta]am/2,
												SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
												SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
												Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
							If[ModeSol[[1]],
								Print["ModeSol-- a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision]];
								modeName[l,m,n]=Insert[KerrSEQ,ModeSol[[4]],ind0];
								(*Switch[s,
								   -2,Global`KerrQNM[l,m,n] =Insert[KerrSEQ,ModeSol[[4]],ind0],
								   -1,Global`KerrQNMe[l,m,n]=Insert[KerrSEQ,ModeSol[[4]],ind0],
									0,Global`KerrQNMs[l,m,n]=Insert[KerrSEQ,ModeSol[[4]],ind0]
								  ];*)
								++ind0; (* must increment to keep ind0 at same a *)
								,(* invalid solution *)
								Print["a-- solution failed."];
								Abort[];
							];
						]];
						\[CapitalDelta]ap = KerrSEQ[[ind0+1,1]]-KerrSEQ[[ind0,1]];
						\[CapitalDelta]am = KerrSEQ[[ind0,1]]-KerrSEQ[[ind0-1,1]];
					];
				];
			];
		,RemoveLevels,
			While[index0>index0m,
				While[index0<limitlist[[1,1]],
					limitlist=Drop[limitlist,1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0m;Continue[]];
				If[index0>limitlist[[1,2]],--index0;Continue[]];
				If[Mod[1000KerrSEQ[[index0,1]],2^(-Maxb)]!=0,
					modeName[l,m,n]=Drop[KerrSEQ,{index0}];
					(*Switch[s,
					   -2,Global`KerrQNM[l,m,n] =Drop[KerrSEQ,{index0}],
					   -1,Global`KerrQNMe[l,m,n]=Drop[KerrSEQ,{index0}],
						0,Global`KerrQNMs[l,m,n]=Drop[KerrSEQ,{index0}]
					  ];*)
				];
				--index0;
			];
		,Update,
			If[precision!=$MinPrecision,Print[Style[StringForm[KerrModeRefineSequence::precision,precision],{Medium,Darker[Red]}]]];
			$MinPrecision=precision;
			For[index0=index0m,index0<=index0p,++index0,
				While[index0>limitlist[[-1,2]],
					limitlist=Drop[limitlist,-1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0p+1;Continue[]];
				If[(Length[KerrSEQ[[index0,2]]]<9 || forcerefine) && (limitlist[[-1,1]]<=index0<=limitlist[[-1,2]]),
					inversion=KerrSEQ[[index0,2,2]];
					oldNrcf=If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]];
					If[RunCFConvergence, (* only extrapolate Nrcf is using CF *)
						Nrcf=If[Length[KerrSEQ[[index0,2]]]>=8,
							If[NumberQ[KerrSEQ[[index0,2,7]]],
								IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/KerrSEQ[[index0,2,7]])],
								oldNrcf],
							IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/3)]];
						If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
						If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
						Nrcf=Max[Nrcf,RCFmin];
						If[Nrcf>RCFmax,
							Print[Style[StringForm[KerrModeRefineSequence::largedepth,Nrcf],{Medium,Darker[Red]}]];
							Print[Style[StringForm[KerrModeRefineSequence::setdepth,RCFmax],{Medium,Darker[Red]}]];
							Nrcf=RCFmax;
						],
						Nrcf=oldNrcf;
					];
					Nm=KerrSEQ[[index0,3,2]];
					a0=KerrSEQ[[index0,1]];
					\[Omega]g=KerrSEQ[[index0,2,1]];
					Almg=KerrSEQ[[index0,3,1]];
					ref\[Epsilon]=If[forcerefine,\[Epsilon],Min[\[Epsilon],KerrSEQ[[index0,2,4]]]];
					If[!forcerefine,
						$MinPrecision=Max[precision,
							If[Length[KerrSEQ[[index0,2]]]>=9,KerrSEQ[[index0,2,9]],IntegerPart[MyPrecision[\[Omega]g]]]]
					];
					ModeSol=ModeSolution[inversion,s,l,m,a0,
										SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
										SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
										Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
					If[ModeSol[[1]],
						Print["ModeSol a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision],
							" Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision],
							"  |\[CapitalDelta]\[Omega]| = ",SetPrecision[Abs[\[Omega]g-ModeSol[[4,2,1]]],MachinePrecision]];
						modeName[l,m,n]=ReplacePart[KerrSEQ,index0->ModeSol[[4]]];
						(*Switch[s,
							   -2,Global`KerrQNM[l,m,n]=ReplacePart[KerrSEQ,index0->ModeSol[[4]]],
							   -1,Global`KerrQNMe[l,m,n]=ReplacePart[KerrSEQ,index0->ModeSol[[4]]],
								0,Global`KerrQNMs[l,m,n]=ReplacePart[KerrSEQ,index0->ModeSol[[4]]]
							  ];*)
						,(* invalid solution *)
						Print[Style[StringForm[KerrModeRefineSequence::indexfail,index0],{Medium,Darker[Red]}]];
					];
				];
			];
		,TerminalLow,
			If[precision!=$MinPrecision,Print[Style[StringForm[KerrModeRefineSequence::precision,precision],{Medium,Darker[Red]}]]];
			$MinPrecision=precision;
			index0m=1;
			index0p=0;
			While[Length[limitlist]>0,
				If[limitlist[[-1,1]]==1,index0p=limitlist[[-1,2]];Break[]];
				limitlist=Drop[limitlist,-1];
			];
			If[index0p==0,Message[KerrModeRefineSequence::TerminalLoRange,refinement];Return[]];
			{terminalSeq,saveSeq}=TakeDrop[KerrSEQ,index0p];
			oldbstruct=bStructure[Take[KerrSEQ,All,1],{1,index0p}];
			newbstruct=DescendFast[oldbstruct];
			addlist=Complement[afrombstruct[#]&/@newbstruct,afrombstruct[#]&/@oldbstruct];
			removelist=Complement[afrombstruct[#]&/@oldbstruct,afrombstruct[#]&/@newbstruct];
			While[Length[addlist]>0,
				nearesta=Nearest[#[[1]]&/@terminalSeq,addlist[[-1]]][[1]];
				ind0=Flatten[Position[#[[1]]&/@terminalSeq,nearesta]][[1]];
				(*Print["TerminalLow (nearest) \[Omega]=",terminalSeq[[ind0,2,1]]," : Alm=",terminalSeq[[ind0,3,1]]];*)
				takerange=If[ind0==1||(addlist[[-1]]>nearesta&&ind0<Length[terminalSeq]),{ind0,ind0+1},{ind0-1,ind0}];
				Switch[extraporder,
					LogLog,
					edat0=Take[terminalSeq,takerange];
					edat=Table[{Log[edat0[[i,1]]],Log[Abs[Re[edat0[[i,2,1]]]]]},{i,1,2}];
					ef=LinearModelFit[edat,{1,a},a];
					\[Omega]g=Sign[Re[edat0[[1,2,1]]]]Exp[ef[Log[addlist[[-1]]]]];
					edat=Table[{Log[edat0[[i,1]]],Log[Abs[Im[edat0[[i,2,1]]]]]},{i,1,2}];
					ef=LinearModelFit[edat,{1,a},a];
					\[Omega]g+=I Sign[Im[edat0[[1,2,1]]]]Exp[ef[Log[addlist[[-1]]]]];
					takerange={ind0-1,ind0+1};
					takerange[[2]]=Min[Length[terminalSeq],takerange[[1]]+2];
					takerange[[1]]=Max[1,takerange[[2]]-2];
					edat0=Take[terminalSeq,takerange];
					edat=Table[{edat0[[i,1]],edat0[[i,3,1]]},{i,1,Length[edat0]}];
					Almg=InterpolatingPolynomial[edat,addlist[[-1]]],
					Asymptote,
					{\[Omega]g,Almg}=AsymptoteFunction[s,l,m,addlist[[-1]]],
					_,
					If[!IntegerQ[extraporder]||extraporder<2,extraporder=2];
					takerange={ind0-1,ind0+1};
					If[extraporder>2,
						takerange[[1]]=Max[1,If[addlist[[-1]]>nearesta,takerange[[1]]-Floor[extraporder/2]-1,takerange[[1]]-Floor[extraporder/2]]];
					];
					takerange[[2]]=Min[Length[terminalSeq],takerange[[1]]+extraporder];
					takerange[[1]]=Max[1,takerange[[2]]-extraporder];
					edat0=Take[terminalSeq,takerange];
					edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
					\[Omega]g=InterpolatingPolynomial[edat,addlist[[-1]]];
					edat=Table[{edat0[[i,1]],edat0[[i,3,1]]},{i,1,Length[edat0]}];
					Almg=InterpolatingPolynomial[edat,addlist[[-1]]];
				];	
				(*Print["TerminalLow (guess) \[Omega]=",\[Omega]g," : Alm=",Almg];*)
				inversion=terminalSeq[[ind0,2,2]];
				ref\[Epsilon]=Min[\[Epsilon],terminalSeq[[ind0,2,4]]];
				Nrcf=If[Length[terminalSeq[[ind0,2]]]>=6,terminalSeq[[ind0,2,6]],terminalSeq[[ind0,2,3]]];
				Nm=terminalSeq[[ind0,3,2]];
				ModeSol=ModeSolution[inversion,s,l,m,addlist[[-1]],
									SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
									SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
									Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
				If[ModeSol[[1]],
					Print["ModeSol a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision],
							" Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision],
							"  |\[CapitalDelta]\[Omega]| = ",SetPrecision[Abs[\[Omega]g-ModeSol[[4,2,1]]],MachinePrecision]];
					terminalSeq=Insert[terminalSeq,ModeSol[[4]],If[addlist[[-1]]>nearesta,ind0+1,ind0]],
					Message[KerrModeRefineSequence::TerminalLowFailure,N[addlist[[1]],{Infinity,20}]];Return[]		
				];
				addlist=Drop[addlist,-1];
			];
			While[Length[removelist]>0,
				terminalSeq=Drop[terminalSeq,Position[#[[1]]&/@terminalSeq,removelist[[1]]][[1]]];
				removelist=Drop[removelist,1];
			];
			modeName[l,m,n]=Join[terminalSeq,saveSeq];
		,TerminalHigh,
			If[precision!=$MinPrecision,Print[Style[StringForm[KerrModeRefineSequence::precision,precision],{Medium,Darker[Red]}]]];
			$MinPrecision=precision;
			index0m=0;
			index0p=NKMode;
			While[Length[limitlist]>0,
				If[limitlist[[-1,2]]==NKMode,index0m=limitlist[[-1,1]];Break[]];
				limitlist=Drop[limitlist,-1];
			];
			If[index0m==0,Message[KerrModeRefineSequence::TerminalHiRange,refinement];Return[]];
			{saveSeq,terminalSeq}=TakeDrop[KerrSEQ,index0m-1];
			oldbstruct=bStructure[Take[KerrSEQ,All,1],{index0m,-1}];
			newbstruct=AscendFast[oldbstruct];
			addlist=Complement[afrombstruct[#]&/@newbstruct,afrombstruct[#]&/@oldbstruct];
			removelist=Complement[afrombstruct[#]&/@oldbstruct,afrombstruct[#]&/@newbstruct];
			While[Length[addlist]>0,
				nearesta=Nearest[#[[1]]&/@terminalSeq,addlist[[1]]][[1]];
				ind0=Flatten[Position[#[[1]]&/@terminalSeq,nearesta]][[1]];
				(*Print["TerminalHigh (nearest) \[Omega]=",terminalSeq[[ind0,2,1]]," : Alm=",terminalSeq[[ind0,3,1]]];*)
				takerange=If[ind0==1||(addlist[[-1]]>nearesta&&ind0<Length[terminalSeq]),{ind0,ind0+1},{ind0-1,ind0}];
				Switch[extraporder,
					Accumulate,
					takerange[[2]]=Min[Length[terminalSeq],takerange[[1]]+2];
					takerange[[1]]=Max[1,takerange[[2]]-2];
					edat0=Take[terminalSeq,takerange];
					edat=Table[{edat0[[i,1]],edat0[[i,3,1]]},{i,1,Length[edat0]}];
					Almg=InterpolatingPolynomial[edat,addlist[[1]]];
					takerange={ind0-5,ind0+5};
					takerange[[2]]=Min[Length[terminalSeq],takerange[[1]]+9];
					takerange[[1]]=Max[1,takerange[[2]]-9];
					edat0=SetPrecision[Take[terminalSeq,takerange],Max[precision,$MinPrecision]];
					edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
					afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
												{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
					\[Omega]g=afit[1-addlist[[1]]];
					edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
					afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
												{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
					\[Omega]g+=I afit[1-addlist[[1]]],
					LogLog,
					edat0=Take[terminalSeq,takerange];
					edat=Table[{Log[edat0[[i,1]]],Log[Abs[Re[edat0[[i,2,1]]]]]},{i,1,2}];
					ef=LinearModelFit[edat,{1,a},a];
					\[Omega]g=Sign[Re[edat0[[1,2,1]]]]Exp[ef[Log[addlist[[-1]]]]];
					edat=Table[{Log[edat0[[i,1]]],Log[Abs[Im[edat0[[i,2,1]]]]]},{i,1,2}];
					ef=LinearModelFit[edat,{1,a},a];
					\[Omega]g+=I Sign[Im[edat0[[1,2,1]]]]Exp[ef[Log[addlist[[1]]]]];
					takerange={ind0-1,ind0+1};
					takerange[[2]]=Min[Length[terminalSeq],takerange[[1]]+2];
					takerange[[1]]=Max[1,takerange[[2]]-2];
					edat0=Take[terminalSeq,takerange];
					edat=Table[{edat0[[i,1]],edat0[[i,3,1]]},{i,1,Length[edat0]}];
					Almg=InterpolatingPolynomial[edat,addlist[[1]]],
					Asymptote,
Print["Shouldn't Call this! : TerminalHigh : Asymptote"];
					{\[Omega]g,Almg}=AsymptoteFunction[s,l,m,addlist[[1]]],
					_,
					If[!IntegerQ[extraporder]||extraporder<2,extraporder=2];
					takerange={ind0-1,ind0+1};
					If[extraporder>2,
						takerange[[1]]=Max[1,If[addlist[[-1]]>nearesta,takerange[[1]]-Floor[extraporder/2]-1,takerange[[1]]-Floor[extraporder/2]]];
					];
					takerange[[2]]=Min[Length[terminalSeq],takerange[[1]]+extraporder];
					takerange[[1]]=Max[1,takerange[[2]]-extraporder];
					edat0=Take[terminalSeq,takerange];
					edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
					\[Omega]g=InterpolatingPolynomial[edat,addlist[[-1]]];
					edat=Table[{edat0[[i,1]],edat0[[i,3,1]]},{i,1,Length[edat0]}];
					Almg=InterpolatingPolynomial[edat,addlist[[1]]];
				];
				(*Print["TerminalHigh (guess) \[Omega]=",\[Omega]g," : Alm=",Almg];*)
				inversion=terminalSeq[[ind0,2,2]];
				ref\[Epsilon]=Min[\[Epsilon],terminalSeq[[ind0,2,4]]];
				Nrcf=If[Length[terminalSeq[[ind0,2]]]>=6,terminalSeq[[ind0,2,6]],terminalSeq[[ind0,2,3]]];
				Nm=terminalSeq[[ind0,3,2]];
				ModeSol=ModeSolution[inversion,s,l,m,addlist[[1]],
									SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
									SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
									Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[ModeSolution]]];
				If[ModeSol[[1]],
					Print["ModeSol a=",Block[{$MinPrecision=0},N[ModeSol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[ModeSol[[4,2,1]],MachinePrecision],
							" Alm=",SetPrecision[ModeSol[[4,3,1]],MachinePrecision],
							"  |\[CapitalDelta]\[Omega]| = ",SetPrecision[Abs[\[Omega]g-ModeSol[[4,2,1]]],MachinePrecision]];
					terminalSeq=Insert[terminalSeq,ModeSol[[4]],If[addlist[[1]]>nearesta,ind0+1,ind0]],
					Message[KerrModeRefineSequence::TerminalHiFailure,N[addlist[[1]],{Infinity,20}]];Return[]	
				];
				addlist=Drop[addlist,1];
			];
			While[Length[removelist]>0,
				terminalSeq=Drop[terminalSeq,Position[#[[1]]&/@terminalSeq,removelist[[1]]][[1]]];
				removelist=Drop[removelist,1];
			];
			modeName[l,m,n]=Join[saveSeq,terminalSeq];
		,_,Message[KerrModeRefineSequence::badaction,action];Return[]
	];
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


GetKerrName[modeList_,s_Integer]:=
Module[{},
	GetKerrName::spinweight="Invalid Spin Weight : `1`";
	GetKerrName::invalid="Invalid Mode, select QNM, TTML, or TTMR";
		Switch[modeList,
			QNM,
				Switch[s,
					-2,Global`KerrQNM,
					-1,Global`KerrQNMe,
					 0,Global`KerrQNMs,
					_,Message[GetKerrName::spinweight,s];Abort[]
					],
			TTML,
				Switch[s,
		             -2,Global`KerrTTML,
		             -1,Global`KerrTTMLe,
		              0,Global`KerrTTMLs,
			         _,Message[GetKerrName::spinweight,s];Abort[]
		            ],
			TTMR,
				Switch[s,
		             2,Global`KerrTTMR,
		             1,Global`KerrTTMRe,
		              0,Global`KerrTTMRs,
			         _,Message[GetKerrName::spinweight,s];Abort[]
		            ],
				_,
					Message[GetKerrName::invalid];
					Abort[]
			]
	]


GetSchName[modeList_,s_Integer]:=
Module[{},
	GetSchName::spinweight="Invalid Spin Weight : `1`";
	GetSchName::invalid="Invalid Mode, select QNM, TTML, or TTMR";
		Switch[modeList,
			QNM,
				Switch[s,
					-2,Global`SchQNMTable,
					-1,Global`SchQNMeTable,
					 0,Global`SchQNMsTable,
					_,Message[GetSchName::spinweight,s];Abort[]
					],
			TTML,
				Switch[s,
		             -2,Global`SchTTMLTable,
		             -1,Global`SchTTMLeTable,
		              0,Global`SchTTMLsTable,
			         _,Message[GetSchName::spinweight,s];Abort[]
		            ],
			TTMR,
				Switch[s,
		             2,Global`SchTTMRTable,
		             1,Global`SchTTMReTable,
		              0,Global`SchTTMRsTable,
			         _,Message[GetSchName::spinweight,s];Abort[]
		            ],
				_,
					Message[GetSchName::invalid];
					Abort[]
			]
	]


Options[KerrModeMakeMultiplet]={OTmultiple->0};


KerrModeMakeMultiplet[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:=
Module[{nm=OptionValue[OTmultiple],KerrSEQ},
	KerrSEQ=modeName[l,m,n];
	If[Head[n]==Integer,
		modeName[l,m,{n,nm}]=KerrSEQ;modeName[l,m,n]=.,
		Null[],
		If[n[[2]]!=nm,
			modeName[l,m,{n[[1]],nm}]=KerrSEQ;modeName[l,m,n]=.
		]
	];
]


Options[ShortenModeSequence]={ShortenBy->Drop};


ShortenModeSequence[l_Integer,m_Integer,n_Integer|n_List,Ns_Integer,OptionsPattern[]]:=
Module[{shorten=OptionValue[ShortenBy],KerrSEQ,SeqStatus,na},
	ShortenModeSequence::seqstatus="`1`[`2`,`3`,`4`] does not exist.";	
	ShortenModeSequence::origlength="Original Length of `1`[`2`,`3`,`4`] = `5`";
	ShortenModeSequence::removelast="Removing last `1` elements";
	ShortenModeSequence::removefirst="Removing first `1` elements";
	ShortenModeSequence::takelast="Taking last `1` elements";
	ShortenModeSequence::takefirst="Taking last `1` elements";
	ShortenModeSequence::invalid="Invalid Shortening Procedure";
	KerrSEQ:=modeName[l,m,n];
	SeqStatus=If[Head[KerrSEQ]==List,If[Length[KerrSEQ]>0,True,False,False],False,False];
	If[!SeqStatus,Print[Style[StringForm[ShortenModeSequence::seqstatus,modeName,l,m,n],{Medium,Darker[Red]}]];Return[]];
	na=Length[KerrSEQ];
	Print[StringForm[ShortenModeSequence::origlength,modeName,l,m,n,na]];
	Switch[shorten,
		Drop,
			If[Ns<0,Print[Style[StringForm[ShortenModeSequence::removelast,-Ns],{Medium,Darker[Green]}]],
					Print[Style[StringForm[ShortenModeSequence::removefirst,Ns],{Medium,Darker[Green]}]]],
		Take,
			If[Ns<0,Print[Style[StringForm[ShortenModeSequence::takelast,-Ns],{Medium,Darker[Green]}]],
					Print[Style[StringForm[ShortenModeSequence::takefirst,Ns],{Medium,Darker[Green]}]]],
		_,Message[ShortenModeSequence::invalid];Abort[];
	];
	modeName[l,m,n]=shorten[KerrSEQ,Ns];
]


bStructure[a_,range_]:=Module[{Na=Length[a],a0=Take[a,range],ap,am},
	am=Prepend[Drop[a0,-1],If[range[[1]]==1 || Na+range[[1]]==0,-\[Infinity],a[[range[[1]]-1]]]];
	ap=Append[Drop[a0,1],If[range[[2]]==-1 || Na-range[[2]]==0,\[Infinity],a[[range[[2]]+1]]]];
	Join[{#[[2]]},
		Flatten[Reap[NestWhileList[
						Function[w,Function[u,Sow[u[[1]]];Drop[u,1]][Join[QuotientRemainder[w[[1]],10^(-3)2^(-w[[2]])],{w[[2]]+1}]]],
						{#[[1]],0},
						#[[1]]!=0&]
					][[2]]]]&/@
		Transpose[{Flatten[a0],
				   (-FullSimplify[Log[1000#]/Log[2]]&/@(Min[#[[1]],#[[2]]]&/@Transpose[{Flatten[a0-am],Flatten[ap-a0]}]))}
				 ]
]


afrombstruct[bstruct_]:=10^(-3)Sum[bstruct[[i+2]]2^(-i),{i,0,Length[bstruct]-2}]


DescendFast[bstruct_]:=Module[{Posmax,count,b,bmax=bstruct[[1,1]],bmin=bstruct[[-1,1]],newa},
	Posmax=Flatten[Position[#[[1]]&/@bstruct,bmax]];
	While[Posmax[[-1]]!=Length[Posmax],Posmax=Drop[Posmax,-1]];
	newa={afrombstruct[bstruct[[1]]]};
	count=4;While[count-->0,AppendTo[newa,newa[[-1]]+ 10^(-3)2^(-bmax)]];
	While[Length[bStructure[newa,{1,-1}][[-1]]]-2>bmax-1,AppendTo[newa,newa[[-1]]+ 10^(-3)2^(-bmax)]];
	For[b=bmax-1,b>bmin,--b,
		count=3;While[count-->0,AppendTo[newa,newa[[-1]]+ 10^(-3)2^(-b)]]
		While[Length[bStructure[newa,{1,-1}][[-1]]]-2>b-1,AppendTo[newa,newa[[-1]]+ 10^(-3)2^(-b)]];
	];
	While[newa[[-1]]<afrombstruct[bstruct[[-1]]],AppendTo[newa,newa[[-1]]+ 10^(-3)2^(-bmin)]];
	While[newa[[-1]]>afrombstruct[bstruct[[-1]]],newa=Drop[newa,-1]];
	bStructure[newa,{1,-1}]
]


AscendFast[bstruct_]:=Module[{Posmax,count,b,bmax=bstruct[[-1,1]],bmin=bstruct[[1,1]],newa},
	Posmax=Flatten[Position[#[[1]]&/@bstruct,bmax]];
	While[Posmax[[1]]!=Length[bstruct]-Length[Posmax]+1,Posmax=Drop[Posmax,1]];
	newa={afrombstruct[bstruct[[-1]]]};
	count=4;While[count-->0,PrependTo[newa,newa[[1]]- 10^(-3)2^(-bmax)]];
	While[Length[bStructure[newa,{1,-1}][[1]]]-2>bmax-1,PrependTo[newa,newa[[1]]- 10^(-3)2^(-bmax)]];
	For[b=bmax-1,b>bmin,--b,
		count=3;While[count-->0,PrependTo[newa,newa[[1]]- 10^(-3)2^(-b)]];
		While[Length[bStructure[newa,{1,-1}][[1]]]-2>b-1,PrependTo[newa,newa[[1]]- 10^(-3)2^(-b)]];
	];
	While[newa[[1]]>afrombstruct[bstruct[[1]]],PrependTo[newa,newa[[1]]- 10^(-3)2^(-bmin)]];
	While[newa[[1]]<afrombstruct[bstruct[[1]]],newa=Drop[newa,1]];
	bStructure[newa,{1,-1}]
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
								RadialCFMinDepth->300,RadialCFDepth->1,ModePrecision->24,NewtonRadius->10^(-3)},
								Options[TestRadialCFConvergence]];


SchwarzschildMode[l_Integer,n_Integer,opts:OptionsPattern[]] :=
Module[{s=OptionValue[SpinWeight],debug=OptionValue[SchDebug],
		rcfdepth=OptionValue[RadialCFDepth],
		Nrcf=300,testinv0,testinv1,Ninv,rcferr,rcfpow,nrpow,jacobianmatrix,Nradialnew,
		newtonRadius=OptionValue[NewtonRadius],
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
		sol1=RadialLentzRoot[Ninv,s,0,0,l(l+1)-s(s+1),SetPrecision[sol1[[1]],precision],Nrcf,l+2,\[Epsilon],newtonRadius,FilterRules[{opts},Options[RadialLentzRoot]]];
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
				],
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
(*Graphics*)


PlotModeFunction[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			Alm_?NumberQ,\[Omega]_?NumberQ,Nrcf_Integer,Nm_Integer]:= 
Module[{Alm\[Omega]},
	Alm\[Omega]=AngularSpectralRoot[s,m,a*\[Omega],Alm,Nm][[1]];
     ModeFunction[n,s,m,a,Alm\[Omega],\[Omega],Nrcf][[1]]
]


PlotModeFunctionL[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			lindex_Integer,\[Omega]_?NumberQ,Nrcf_Integer,Nm_Integer]:= 
Module[{angularsol,expconv,count=0,Nmatrix=Nm,Alm\[Omega]},	
	PlotModeFunctionL::inaccurateAlm="Angular solution matrix too small; expconv = `1`";
	angularsol=AngularSpectralRootIndex[s,m,a*\[Omega],lindex,Nmatrix];
	expconv=Max[Take[Abs[angularsol[[3]]],-2]];
	While[expconv >= 10^-8 && ++count<5,  (* Make sure Spectral resolution is good before doing lots of work *)
		angularsol=AngularSpectralRootIndex[s,m,a*\[Omega],lindex,++Nmatrix];
		expconv=Max[Take[Abs[angularsol[[3]]],-2]];
	];
	If[count>=5,Message[PlotModeFunctionL::inaccurateAlm,expconv]];
	Alm\[Omega]=angularsol[[1]];
	ModeFunction[n,s,m,a,Alm\[Omega],\[Omega],Nrcf][[1]]
]


Options[PlotSchModes]=Union[{PlotTable->Null[]},Options[ListPlot]];


PlotSchModes[l_Integer,opts:OptionsPattern[]]:=
Module[{ptable=OptionValue[PlotTable]},
	PlotSchModes::plottable="Invalid PlotTable : `1`";
	If[ptable==Null[],Message[PlotSchModes::plottable,ptable];Abort[]];
	plist=Table[{Re[ptable[l,n][[1]]],-Im[ptable[l,n][[1]]]},{n,0,ptable[l]-1}];
	mlist=Table[{-Re[ptable[l,n][[1]]],-Im[ptable[l,n][[1]]]},{n,0,ptable[l]-1}];
	ListLinePlot[{plist,mlist},PlotMarkers->Automatic,FilterRules[{opts},Options[ListLinePlot]]]
]


Options[SchwarzschildOmega]={ModeType->Null[],SpinWeight->Null[]};


SchwarzschildOmega[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,not},
	SchwarzschildOmega::noschwarzt="Schwarzschild limit not found, a=`1`.";
	If[Head[n]==Integer,not=n,Null[],not=n[[1]]];
	If[l==2 && not==8,Return[{0,2}]]; (* Special Cases *)
	If[l==3 && not==40,Return[{0,10}]]; (* Special Cases *)
	KerrSEQ:= modeName[l,m,n];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	If[Chop[KerrSEQ[[1,1]]]!=0,Message[SchwarzschildOmega::noschwarzt,KerrSEQ[[1,1]]];Abort[]];
	{Re[KerrSEQ[[1,2,1]]],-Im[KerrSEQ[[1,2,1]]]}
]


Options[KerrOmegaList]={ModeType->Null[],SpinWeight->Null[]};


KerrOmegaList[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:= 
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,Na},
	KerrSEQ:= modeName[l,m,n];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	Na = Length[KerrSEQ];
	Table[{Re[KerrSEQ[[i,2,1]]],-Im[KerrSEQ[[i,2,1]]]},{i,1,Na}]
]


Options[KerrOmegaListS]={ModeType->Null[],SpinWeight->Null[]};


KerrOmegaListS[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:= 
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,Na,Nend,i,Slist={}},
	KerrSEQ:= modeName[l,m,n];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	Na = Length[KerrSEQ];
	Nend = If[KerrSEQ[[Na,1]]<999999/1000000,Na,Na-1];
	For[i=1,i<=Nend,++i,
		If[Mod[KerrSEQ[[i,1]],1/20]==0,
			AppendTo[Slist,{Re[KerrSEQ[[i,2,1]]],-Im[KerrSEQ[[i,2,1]]]}]
		];
	];
	If[KerrSEQ[[Na,1]]>=999999/1000000,
		Append[Slist,{Re[KerrSEQ[[Na,2,1]]],-Im[KerrSEQ[[Na,2,1]]]}],
		Slist
	]
]


Options[KerrAList]={ModeType->Null[],SpinWeight->Null[]};


KerrAList[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:= 
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,Na},
	KerrSEQ:= modeName[l,m,n];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	Na = Length[KerrSEQ];
	Table[{Re[KerrSEQ[[i,3,1]]],Im[KerrSEQ[[i,3,1]]]},{i,1,Na}]
]


Options[KerrAListS]={ModeType->Null[],SpinWeight->Null[]};


KerrAListS[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:= 
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,Na,Nend,i,Slist={}},
	KerrSEQ:= modeName[l,m,n];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	Na = Length[KerrSEQ];
	Nend = If[KerrSEQ[[Na,1]]<999999/1000000,Na,Na-1];
	For[i=1,i<=Nend,++i,
		If[Mod[KerrSEQ[[i,1]],1/20]==0,
			AppendTo[Slist,{Re[KerrSEQ[[i,3,1]]],Im[KerrSEQ[[i,3,1]]]}]
		];
	];
	If[KerrSEQ[[Na,1]]>=999999/1000000,
		Append[Slist,{Re[KerrSEQ[[Na,3,1]]],Im[KerrSEQ[[Na,3,1]]]}],
		Slist
	]
]


Options[KerraOmegaList]={ModeType->Null[],SpinWeight->Null[]};


KerraOmegaList[l_Integer,m_Integer,n_Integer|n_List,ReIm_Symbol,OptionsPattern[]]:= 
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,Na},
	KerraOmegaList::ReIm="ReIm must be Re, Im, or Abs, set to `1`";
	KerrSEQ:= modeName[l,m,n];
	If[ReIm==Re || ReIm==Im || ReIm==Abs,Null[],Null[],Message[KerraOmegaList::ReIm,ReIm];Abort[]];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	Na = Length[KerrSEQ];
	Table[{KerrSEQ[[i,1]],If[ReIm==Im,-1,1,1]ReIm[KerrSEQ[[i,2,1]]]},{i,1,Na}]
]


Options[KerraOmegaListS]={ModeType->Null[],SpinWeight->Null[]};


KerraOmegaListS[l_Integer,m_Integer,n_Integer|n_List,ReIm_Symbol,OptionsPattern[]]:= 
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,Na,Nend,i,Slist={}},
	KerraOmegaListS::ReIm="ReIm must be Re, Im, or Abs, set to `1`";
	KerrSEQ:= modeName[l,m,n];
	If[ReIm==Re || ReIm==Im || ReIm==Abs,Null[],Null[],Message[KerraOmegaListS::ReIm,ReIm];Abort[]];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	Na = Length[KerrSEQ];
	Nend = If[KerrSEQ[[Na,1]]<999999/1000000,Na,Na-1];
	For[i=1,i<=Nend,++i,
		If[Mod[KerrSEQ[[i,1]],1/20]==0,
			AppendTo[Slist,{KerrSEQ[[i,1]],If[ReIm==Im,-1,1,1]ReIm[KerrSEQ[[i,2,1]]]}]
		];
	];
	If[KerrSEQ[[Na,1]]>=999999/1000000,
		Append[Slist,{KerrSEQ[[Na,1]],If[ReIm==Im,-1,1,1]ReIm[KerrSEQ[[Na,2,1]]]}],
		Slist
	]
]


Options[KerraAList]={ModeType->Null[],SpinWeight->Null[]};


KerraAList[l_Integer,m_Integer,n_Integer|n_List,ReIm_Symbol,OptionsPattern[]]:= 
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,Na},
	KerraOmegaList::ReIm="ReIm must be RRe, Im, or Abs, set to `1`";
	KerrSEQ:= modeName[l,m,n];
	If[ReIm==Re || ReIm==Im || ReIm==Abs,Null[],Null[],Message[KerraOmegaList::ReIm,ReIm];Abort[]];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	Na = Length[KerrSEQ];
	Table[{KerrSEQ[[i,1]],ReIm[KerrSEQ[[i,3,1]]]},{i,1,Na}]
]


Options[KerraAListS]={ModeType->Null[],SpinWeight->Null[]};


KerraAListS[l_Integer,m_Integer,n_Integer|n_List,ReIm_Symbol,OptionsPattern[]]:= 
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],KerrSEQ,Na,Nend,i,Slist={}},
	KerraOmegaListS::ReIm="ReIm must be Re, Im, or Abs, set to `1`";
	KerrSEQ:= modeName[l,m,n];
	If[ReIm==Re || ReIm==Im || ReIm==Abs,Null[],Null[],Message[KerraOmegaListS::ReIm,ReIm];Abort[]];
	If[MemberQ[{QNM,TTML,TTMR},modetype],KerrSEQ:=GetKerrName[modetype,s][l,m,n]];
	Na = Length[KerrSEQ];
	Nend = If[KerrSEQ[[Na,1]]<999999/1000000,Na,Na-1];
	For[i=1,i<=Nend,++i,
		If[Mod[KerrSEQ[[i,1]],1/20]==0,
			AppendTo[Slist,{KerrSEQ[[i,1]],ReIm[KerrSEQ[[i,3,1]]]}]
		];
	];
	If[KerrSEQ[[Na,1]]>=999999/1000000,
		Append[Slist,{KerrSEQ[[Na,1]],ReIm[KerrSEQ[[Na,3,1]]]}],
		Slist
	]
]


Options[ModePlotOmega]=Union[{ModeType->Null[],SpinWeight->Null[],OTmultiple->{}},
							Options[ListLinePlot],Options[ListPlot]];


ModePlotOmega[l_Integer,n_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],multiple=OptionValue[OTmultiple],
		legend=OptionValue[PlotLegends],modetype=OptionValue[ModeType],autolegend,
		SpinWeightTable,KerrSEQ,
		mmodes={},multints,i,pos,m,linelist,pointlist,mainplot},
	SpinWeightTable:=modeName;
	If[MemberQ[{QNM,TTML,TTMR},modetype],SpinWeightTable:=GetKerrName[modetype,s]];
	multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
	For[m=-l,m<=l,++m,
		If[MemberQ[multints,m],
			(* m,n is a multplet *)
			pos=Flatten[Position[multiple,{m,_}]][[1]];
			For[i=0,i<multiple[[pos,2]],++i,
				KerrSEQ:=SpinWeightTable[l,m,{n,i}];
				If[Head[KerrSEQ]==List,AppendTo[mmodes,{m,{n,i}}]]
			],
			(* Not a multplet *)
			KerrSEQ:=SpinWeightTable[l,m,n];
			If[Head[KerrSEQ]==List,AppendTo[mmodes,{m,n}]]
		]
	];
	autolegend=Table[If[Head[mmodes[[i,2]]]==List,Subscript[mmodes[[i,1]], mmodes[[i,2,2]]],Null,mmodes[[i,1]]],{i,1,Length[mmodes]}];
	If[legend==Automatic,legend=autolegend];
	If[Head[legend]==Placed,If[legend[[1]]==Automatic,legend=Placed[autolegend,legend[[2]]]]];
	linelist=KerrOmegaList[l,#[[1]],#[[2]],FilterRules[{opts},Options[KerrOmegaList]]]&/@  mmodes;
	pointlist=KerrOmegaListS[l,#[[1]],#[[2]],FilterRules[{opts},Options[KerrOmegaListS]]]&/@  mmodes;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,PlotLegends->legend,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic]],
		Show[mainplot]
	]
]


ModePlotOmega[l_Integer,m_Integer,n_Integer|n_List,opts:OptionsPattern[]]:=
Module[{mmodes={},linelist,pointlist,mainplot},
	mmodes={m};
	linelist=KerrOmegaList[l,#,n,FilterRules[{opts},Options[KerrOmegaList]]]&/@  mmodes;
	pointlist=KerrOmegaListS[l,#,n,FilterRules[{opts},Options[KerrOmegaListS]]]&/@  mmodes;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic]],
		Show[mainplot]
	]
]


Options[ModePlotA]=Union[{ModeType->Null[],SpinWeight->Null[],OTmultiple->{}},
							Options[ListLinePlot],Options[ListPlot]];


ModePlotA[l_Integer,n_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],multiple=OptionValue[OTmultiple],
		legend=OptionValue[PlotLegends],modetype=OptionValue[ModeType],autolegend,
		SpinWeightTable,KerrSEQ,
		mmodes={},multints,i,pos,m,linelist,pointlist,mainplot},
	SpinWeightTable:=modeName;
	If[MemberQ[{QNM,TTML,TTMR},modetype],SpinWeightTable:=GetKerrName[modetype,s]];
	multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
	For[m=-l,m<=l,++m,
		If[MemberQ[multints,m],
			(* m,n is a multplet *)
			pos=Flatten[Position[multiple,{m,_}]][[1]];
			For[i=0,i<multiple[[pos,2]],++i,
				KerrSEQ:=SpinWeightTable[l,m,{n,i}];
				If[Head[KerrSEQ]==List,AppendTo[mmodes,{m,{n,i}}]]
			],
			(* Not a multplet *)
			KerrSEQ:=SpinWeightTable[l,m,n];
			If[Head[KerrSEQ]==List,AppendTo[mmodes,{m,n}]]
		]
	];
	autolegend=Table[If[Head[mmodes[[i,2]]]==List,Subscript[mmodes[[i,1]], mmodes[[i,2,2]]],Null,mmodes[[i,1]]],{i,1,Length[mmodes]}];
	If[legend==Automatic,legend=autolegend];
	If[Head[legend]==Placed,If[legend[[1]]==Automatic,legend=Placed[autolegend,legend[[2]]]]];
	linelist=KerrAList[l,#[[1]],#[[2]],FilterRules[{opts},Options[KerrAList]]]&/@  mmodes;
	pointlist=KerrAListS[l,#[[1]],#[[2]],FilterRules[{opts},Options[KerrAListS]]]&/@  mmodes;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,PlotLegends->legend,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic]],
		Show[mainplot]
	]
]


ModePlotA[l_Integer,m_Integer,n_Integer|n_List,opts:OptionsPattern[]]:=
Module[{mmodes={},linelist,pointlist,mainplot},
	mmodes={m};
	linelist=KerrAList[l,#,n,FilterRules[{opts},Options[KerrAList]]]&/@  mmodes;
	pointlist=KerrAListS[l,#,n,FilterRules[{opts},Options[KerrAListS]]]&/@  mmodes;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic]],
		Show[mainplot]
	]
]


Options[ModePlotOmegaTones]=Union[{ModeType->Null[],SpinWeight->Null[],OTmultiple->{}},
							Options[ListLinePlot],Options[ListPlot]];


ModePlotOmegaTones[l_Integer,m_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],multiple=OptionValue[OTmultiple],
		legend=OptionValue[PlotLegends],modetype=OptionValue[ModeType],autolegend,
		SpinWeightTable,KerrSEQ,
		ntones={},schtones={},n,multints,pos,i,j,linelist,pointlist,mainplot,Schlist,amin,amini},
	SpinWeightTable:=modeName;
	If[MemberQ[{QNM,TTML,TTMR},modetype],SpinWeightTable:=GetKerrName[modetype,s]];
	multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
	For[n=0,n<=300,++n,
		If[MemberQ[multints,n],
			pos=Flatten[Position[multiple,{n,_}]][[1]];
			For[i=0,i<multiple[[pos,2]],++i,
				KerrSEQ:=SpinWeightTable[l,m,{n,i}];
				If[Head[KerrSEQ]==List,AppendTo[ntones,{n,i}]]
			]
			,
			KerrSEQ:=SpinWeightTable[l,m,n];
			If[Head[KerrSEQ]==List,AppendTo[ntones,n]]
		]
	];
	schtones=ntones;
	For[i=1,i<=Length[multints],++i,
		pos=Position[schtones,{multints[[i]],_}];
		amini=1;
		amin=SpinWeightTable[l,m,schtones[[pos[[1,1]]]]][[1,1]];
		For[j=2,j<=Length[pos],++j,
			If[SpinWeightTable[l,m,schtones[[pos[[j,1]]]]][[1,1]]<amin,
				amini=j;
				amin=SpinWeightTable[l,m,schtones[[pos[[j,1]]]]][[1,1]];
			];
		];
		pos=Drop[pos,amini];
		schtones=Delete[schtones,pos];
	];
	autolegend=Table[If[Head[ntones[[i]]]==List,Subscript[ntones[[i,1]], ntones[[i,2]]],Null,ntones[[i]]],{i,1,Length[ntones]}];
	If[legend==Automatic,legend=autolegend];
	If[Head[legend]==Placed,If[legend[[1]]==Automatic,legend=Placed[autolegend,legend[[2]]]]];
	linelist=KerrOmegaList[l,m,#,FilterRules[{opts},Options[KerrOmegaList]]]&/@ ntones;
	pointlist=KerrOmegaListS[l,m,#,FilterRules[{opts},Options[KerrOmegaListS]]]&/@ ntones;
	Schlist=SchwarzschildOmega[l,m,#,FilterRules[{opts},Options[SchwarzschildOmega]]]&/@ schtones;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,PlotLegends->legend,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic],
			ListLinePlot[Schlist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers,PlotStyle}]],PlotStyle->{Gray,Dashed},PlotRange->All],ImageSize->800],
		Show[mainplot,
			ListLinePlot[Schlist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers,PlotStyle}]],PlotStyle->{Gray,Dashed},PlotRange->All],ImageSize->800]
	]
]


Options[SWSFLists]=Union[{ModeType->Null[],SpinWeight->Null[],OutputType->Complex,ChopLevel->False}];


SWSFLists[l_Integer,m_Integer,n_Integer|n_List,index_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],
		SpinWeightTable,KerrSEQ,SWdat,NC,x,theta,Ntheta,lmin,Matdlx,SWSF,choplev},
	SpinWeightTable:=modeName;
	If[MemberQ[{QNM,TTML,TTMR},modetype],SpinWeightTable:=GetKerrName[modetype,s]];
	KerrSEQ:=SpinWeightTable[l,m,n];
	SWdat=KerrSEQ[[index,3,3]];
	NC=KerrSEQ[[index,3,2]];
	x=Table[x,{x,-1,1,1/100}];
	theta=ArcCos[#]&/@x;
	Ntheta=Length[theta];
	lmin=Max[Abs[m],Abs[s]];
	Matdlx = Table[N[WignerD[{j-1+lmin,m,-s},0,theta[[k]],0]],{k,1,Ntheta},{j,1,NC}];
	SWSF = N[(-1)^m Sqrt[\[Pi]] Matdlx . Table[N[Sqrt[2(j-1+lmin)+1]SWdat[[j]]],{j,1,NC}]];
	If[OptionValue[ChopLevel]!=False,
		choplev=OptionValue[ChopLevel];
		If[IntegerQ[choplev],
			SWSF = Chop[SWSF,10^(choplev)]]
	];
	Switch[OptionValue[OutputType],
		Complex,Transpose[{x,SWSF}],
		ReIm,{Transpose[{x,Re[SWSF]}],Transpose[{x,Im[SWSF]}]},
		AbsArg,{Transpose[{x,Abs[SWSF]}],Transpose[{x,Arg[SWSF]}]}
	]
]


Options[SWSFRealPath]=Union[{ModeType->Null[],SpinWeight->Null[],PathStart->Automatic,MaxSteps->5000,StepSize->1/50,PlotStart->False,\[Phi]guess->Automatic,PrintPoleValues->False}]


SWSFRealPath[l_Integer,m_Integer,n_Integer|n_List,index_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],modetype=OptionValue[ModeType],SpinWeightTable,KerrSEQ,
		SWdat,NC,lmin,Matdlx,SWSF,z,z0,zf,\[Delta]z,\[Phi],\[Phi]g,\[Phi]0,phase,zlist,SWSFlist,count=0,
		max\[Delta]r=OptionValue[StepSize],maxcount=OptionValue[MaxSteps],
		plot=OptionValue[PlotStart]},
	SWSFRealPath::PathStart="Invalid PathStart : `1`";
	SpinWeightTable:=modeName;
	If[MemberQ[{QNM,TTML,TTMR},modetype],SpinWeightTable:=GetKerrName[modetype,s]];
	KerrSEQ:=SpinWeightTable[l,m,n];
	SWdat=KerrSEQ[[index,3,3]];
	NC=KerrSEQ[[index,3,2]];
	lmin=Max[Abs[m],Abs[s]];
	Matdlx = Table[WignerD[{j-1+lmin,m,-s},0,ArcCos[z],0],{j,1,NC}];
	SWSF = (-1)^m Sqrt[\[Pi]] Matdlx . Table[Sqrt[2(j-1+lmin)+1]SWdat[[j]],{j,1,NC}];
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
(*End of KerrModes Package*)


End[] (* `Private` *)


EndPackage[]
