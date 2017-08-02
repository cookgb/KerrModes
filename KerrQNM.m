(* ::Package:: *)

(* ::Title:: *)
(*QuasiNormal Modes of Kerr*)


(* ::Section:: *)
(*Begin KerrQNM Package*)


BeginPackage["KerrQNM`",{"KerrModes`"}]


Unprotect[KerrQNMDebug];
KerrQNMDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[KerrQNMDebug,Unprotect["KerrQNM`*"];Unprotect["KerrQNM`Private`*"]];
Protect[KerrQNMDebug];


(* ::Section:: *)
(*Documentation of External Functions in KerrModes Namespace*)


SetSpinWeight::usage=
	"SetSpinWeight[s] sets the value of the spin-weight used in all subsequent "<>
	"QNM computations:\n"<>
	"\t s=-2 : Gravitational perturbations\n"<>
	"\t s=-1 : Electro-Magnetic perturbations\n"<>
	"\t s= 0 : Scalar perturbations."


(* ::Section:: *)
(*Definitions for KerrModes Namespace*)


Begin["KerrModes`Private`"]


(* ::Subsection::Closed:: *)
(*Define the QNM Radial Equation recurrence relation coefficients*)


rp=1+Sqrt[1-a^2]; rm=1-Sqrt[1-a^2]; \[Sigma]p=(2\[Omega]*rp-m a)/(rp-rm); \[Sigma]m=(2\[Omega]*rm-m a)/(rp-rm);


\[Zeta]=I*\[Omega]; \[Xi]=-s-I*\[Sigma]p; \[Eta]=-I*\[Sigma]m;


p=(rp-rm)\[Zeta]/2; \[Alpha]=1+s+\[Xi]+\[Eta]-2\[Zeta]+s*I*\[Omega]/\[Zeta]; \[Gamma]=1+s+2\[Eta]; \[Delta]=1+s+2\[Xi]; \[Sigma]=Alm+a^2*\[Omega]^2-8\[Omega]^2+p(2\[Alpha]+\[Gamma]-\[Delta])+(1+s-(\[Gamma]+\[Delta])/2)(s+(\[Gamma]+\[Delta])/2);


D0=Simplify[\[Delta]]; D1=Simplify[-2+4p-2\[Alpha]+\[Gamma]-\[Delta]]; D2=Simplify[2+2\[Alpha]-\[Gamma]];
D3=Simplify[4p*\[Alpha]-\[Alpha]*\[Delta]-\[Sigma]]; D4=Simplify[\[Alpha](1+\[Alpha]-\[Gamma])];


\[Alpha]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D0+1)n+D0;
\[Beta]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=-2n^2+(D1+2)n+D3;
\[Gamma]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D2-3)n+D4-D2+2;


Clear[rp,rm,\[Sigma]p,\[Sigma]m,\[Zeta],\[Xi],\[Eta],p,\[Alpha],\[Gamma],\[Delta],\[Sigma],D0,D1,D2,D3,D4];
If[!KerrQNMDebug,Protect[\[Alpha]r,\[Beta]r,\[Gamma]r]];


(* ::Subsection::Closed:: *)
(*Approximation of remainder at the Nmax element of the continued fraction*)


RadialCFRemainder[s_Integer,m_Integer,a_Rational|a_Integer,
				  Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:=
Module[{C12tmp,C1,C2,C3,C4,C5,Rem,Err},
	C12tmp=-4I Sqrt[1-a^2]\[Omega];
	C1 = Sqrt[C12tmp];
	If[Re[C1]>0,C1=-C1];
	C2=-(3-4s+8I(\[Omega]+Sqrt[1-a^2]\[Omega]))/4;
	C3=(3+16Alm+16s+16s^2-32a m \[Omega]-128\[Omega]^2+80a^2\[Omega]^2-128Sqrt[1-a^2]\[Omega]^2-32I(-5Sqrt[1-a^2]\[Omega]+4Sqrt[1-a^2]s \[Omega]))/(32C1);
	C4=(32\[Omega](9+4Alm+2s(-5+4s)+4a \[Omega](-2m+a \[Omega]))+(I(3+16Alm+96a m \[Omega]+16(s^2+(32-51a^2+32Sqrt[1-a^2])\[Omega]^2+s(1+16a \[Omega](-m+2a \[Omega])))))/Sqrt[1-a^2])/(256\[Omega]);
	C5=-(1/(8192(-1+a^2)C1 \[Omega]))I(-256Sqrt[1-a^2]Alm^2+Sqrt[1-a^2](63+288s+32s^2-512s^3-256s^4)+64I(21+100s+48s^2-64s^3+I a Sqrt[1-a^2]m(9+48s+112s^2)+a^2(-21-100s-48s^2+64s^3))\[Omega]+32(320I a m(-3+4s)-320I a^3m(-3+4s)-8(-3+225Sqrt[1-a^2]-16(1+19Sqrt[1-a^2])s+16(-1+3Sqrt[1-a^2])s^2)+a^2(3(-8+591Sqrt[1-a^2])+224Sqrt[1-a^2]m^2-16(8+133Sqrt[1-a^2])s+16(-8+59Sqrt[1-a^2])s^2))\[Omega]^2-1024I(-8I a(1+Sqrt[1-a^2])m-I a^3(-8+15Sqrt[1-a^2])m-5a^4(-7+4s)+8(1+Sqrt[1-a^2])(3+4s)-a^2(59+24Sqrt[1-a^2]+4(3+8Sqrt[1-a^2])s))\[Omega]^3+256(-128(1+Sqrt[1-a^2])+a^4(-80+7Sqrt[1-a^2])+16a^2(13+9Sqrt[1-a^2]))\[Omega]^4-32Alm(Sqrt[1-a^2](-9+16s+16s^2)-32I(7+3I a Sqrt[1-a^2]m-4s+a^2(-7+4s))\[Omega]-16(8(1+Sqrt[1-a^2])+a^2(-8+11Sqrt[1-a^2]))\[Omega]^2));
	Rem=1+C1/Sqrt[Nmax]+C2/Nmax+C3/(Nmax^(3/2))+C4/Nmax^2+C5/(Nmax^(5/2));
	Err=Abs[(C5/(Nmax^(5/2)))/Rem];
	{Rem,Err}
]


If[!KerrQNMDebug,Protect[RadialCFRemainder]];


(* ::Subsection:: *)
(*Set SpinWeight, and Data-Variable Names*)


SetSpinWeight[s_Integer]:=
Module[{},
(* Mode selection preset for QNM *)
	Unprotect[RunCFConvergence];RunCFConvergence=True;Protect[RunCFConvergence];
	ModeFunction[n_,s1_,m_,a_,Alm_,\[Omega]_,Nrcf_]=RadialCFRem[n,s1,m,a,Alm,\[Omega],Nrcf];
	SetSpinWeight::spinweight="Invalid QNM Spin Weight : `1`";
	SetSpinWeight::confirm="All KerrMode routines (QNM) set for Spin-Weight s = `1`";
	Switch[s,
		   -2,modeName:=Global`KerrQNM; SchTable:=Global`SchQNMTable,
		   -1,modeName:=Global`KerrQNMe; SchTable:=Global`SchQNMeTable,
		    0,modeName:=Global`KerrQNMs; SchTable:=Global`SchQNMsTable,
			_,Message[SetSpinWeight::spinweight,s];Abort[]
		  ];
	SetOptions[KerrQNM`SchwarzschildQNM,SpinWeight->s];
	SetOptions[KerrQNM`KerrQNMSequence,SpinWeight->s];
	SetOptions[KerrModes`KerrOmegaList,SpinWeight->s];
	SetOptions[KerrModes`KerrOmegaListS,SpinWeight->s];
	Print[Style[StringForm[SetSpinWeight::confirm,s],{Medium,Darker[Green]}]];
]


If[!KerrModeDebug,Protect[SetSpinWeight]];


End[] (* KerrModes`Private` *)


(* ::Section:: *)
(*Documentation of External Functions in KerrQNM Namespace*)





KerrQNMSequence::usage=
	"KerrQNMSequence[l,m,n,\[Epsilon]] computes a sequence of Quasi-Normal Mode solutions "<>
	"for overtone n of mode (l,m).  The solutions are computed to an absolute "<>
	"accuracty of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  The sequence is "<>
	"parameterized by the dimensionless angular momentum 'a', with 0\[LessEqual]a\[LessEqual]1.  "<>
	"Steps along the sequence take sizes \[CapitalDelta]a=\!\(\*SuperscriptBox[\(2\), \(-b\)]\)/1000 "<>
    "where \!\(\*SubscriptBox[\(b\),\(min\)]\)\[LessEqual]b\[LessEqual]\!\(\*SubscriptBox[\(b\),\(max\)]\).\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  'n' "<>
	"can be either an Integer or an overtone multiplet index.  An overtone multilplet index "<>
	"is a 2 element list {n,mult}, where 'n' is the Integer overtone number, and 'mult' is "<>
	"in the range 0,1,...,(Nmult-1), with 'Nmult' the number of sequences with the same overtone "<>
	"index.\n\n"<>
	"Options:\n"<>
	"\t Root\[Epsilon] \[Rule] Null[ ]\n"<>
	"\t\t Solutions are not considered valid until the value of the continued fraction\n"<>
    "\t\t is less than \!\(\*SuperscriptBox[\(10\), \(Root\[Epsilon]\)]\).  "<>
    " By default, Root\[Epsilon] is set to \[Epsilon].\n"<>
	"\t SeqDirection \[Rule] Forward : Forward, Backward\n"<>
	"\t\t Direction for new elements of the sequence.  With Forward, new elements \n"<>
    "\t\t are added after the maximum value of 'a'.  With Backward, new elements \n"<>
    "\t\t are added before the minimum value of 'a'.\n"<>
	"\t Minblevel \[Rule] 0\n"<>
	"\t\t Value for \!\(\*SubscriptBox[\(b\), \(min\)]\)\n"<>
	"\t Maxblevel \[Rule] 20\n"<>
	"\t\t Value for \!\(\*SubscriptBox[\(b\), \(max\)]\)\n"<>
	"\t Maximala\[Epsilon] \[Rule] 10\n"<>
	"\t\t Sequencer terminates at a=1 if Maximala\[Epsilon]=False.\n"<>
    "\t\t If Maximala\[Epsilon] is an integer (n), it terminates at 1-\!\(\*SuperscriptBox[\(2\),\(-n\)]\)/1000.\n"<>
	"\t Max\[CapitalDelta]\[Omega] \[Rule] 0.01\n"<>
	"\t\t Maximum distance in \[Omega] space between solutions on a sequence\n"<>
	"\t ModeaStart\[Rule]0 : {a,\[Omega],\!\(\*SubscriptBox[\(A\), \(lm\)]\)}\n"<>
	"\t\t This option is only used when starting a new sequence.\n"<>
	"\t\t The List contains the initial value of 'a', and initial guesses for \[Omega] and \!\(\*SubscriptBox[\(A\), \(lm\)]\).\n"<>
	"\t\t An option 4th argument, specifying the initial Integer size of the spectral\n"<>
	"\t\t matrix used to solve the angular Teukolsky equation.\n"<>
	"\t ModeGuess \[Rule] 0 : {a,\[Omega]}\n"<>
	"\t\t This option is used to override the initial guess (only the first solution).\n"<>
	"\t\t The list contains initial guesses for \[Omega] and \!\(\*SubscriptBox[\(A\), \(lm\)]\).\n"<>
	"\t\t An option 3rd argument contains the initial depth of the continued fraction\n"<>
	"\t\t An option 4th argument, specifying the initial Integer size of the spectral\n"<>
	"\t\t matrix used to solve the angular Teukolsky equation.\n"<>
	"\t NoNeg\[Omega] \[Rule] False : True, False\n"<>
	"\t\t If True, solutions are forced to have a positive real part for \[Omega]\n"<>
	"\t SolutionWindowl\[Rule]1/2\n"<>
	"\t SolutionWindowt\[Rule]1/3\n"<>
	"\t\t Set the size of the solution window.  Solutions must fall within this window\n"<>
	"\t\t to be accepted.  Setting either to 0 causes any solutions to be accepted.\n"<>
    "\t ExtrapolationOrder \[Rule] 2 : An integer or 'Accumulate\n"<>
	"\t\t An integer specifies the  order of polynomial extrapolation for the next point.\n"<> 
	"\t\t Otherwise for sequences where \[Omega]\[Rule]m/2 as a\[Rule]1, set extraporder to \n"<>
    "\t\t 'Accumulate'to provide a more accurate extrapolation\n"<>
	"\t CurvatureRatio \[Rule] 1/2\n"<>
	"\t\t In deciding the size of \[CapitalDelta]a, the local ratio of stepsize \[CapitalDelta]\[Omega] betweeen points\n"<>
	"\t\t to the ratius of curvature is compared to CurvatureRatio.\n"<>
    "\t\t Smaller numbers give higher resolution.\n"<>
	"\t RadialCFMinDepth \[Rule] 300\n"<>
	"\t\t Minumum depth of the continued fraction\n"<>
	"\t RadialCFDepth\[Rule]1\n"<>
	"\t\t The initial Radial Continued Fraction Depth is usually taken from the prior \n"<>
	"\t\t solution in the sequence.  Fractional values reduce this initial vaue by that\n"<>
	"\t\t fraction.  Integer values larger than RadialCFMinDepth specify the next value.\n"<>
	"\t RadialCFDigits\[Rule]8\n"<>
	"\t\t Approximate number of digits of agreement when testing the continued fraction\n"<>
	"\t\t approximation.  The continued fraction is evaluated in the proximity of the \n"<>
	"\t\t guessed root, then evaluated again with a depth half-again deeper.  The first\n"<>
	"\t\t RadialCFDigits non-vanishing digits must agree for the solution to be accepted.\n"<>
	"\t ModePrecision\[Rule]24\n"<>
    "\t\t The minimum number of digits used in the computation.\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t The radial solver make use of numerical derivatives.  The relative step size\n"<>
	"\t\t is set to \!\(\*SuperscriptBox[\(10\), \(JacobianStep\)]\).\n"<>
	"\t SolutionRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial Newton iterations.\n"<>
    "\t\t (Note this overrides RadialRelax).\n"<>
	"\t SolutionIter \[Rule] 50\n"<>
    "\t\t Controls number of 'solution' iterations before under-relaxation is tried,\n"<>
    "\t\t and ultimately the total number of iterations allowed.\n"<>
	"\t SolutionOscillate \[Rule] 10\n"<>
    "\t\t Controls number of times the iteration is allowed to become non-convergent.\n"<>
	"\t SolutionSlow \[Rule] 10\n"<>
    "\t\t Controls number of slow Newton iterations are allowed before under-relaxation\n"<>
    "\t\t is tried.\n"<>
	"\t SolutionDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial/Angular iterations.\n"<>
	"\t RadialDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial Newton iterations.\n"<>
	"\t SpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight is set by default by SetSpinWeight[]."
	


SchwarzschildQNM::usage=
	"SchwarzschildQNM[l,n] computes the Quasi-Normal Mode solutions for overtone n of mode l.  "<>
	"The mode is computed to an accuracy of \!\(\*SuperscriptBox[\(10\), \(-14\)]\).  For given 'l', if solutions with overtones "<>
	"(n-1) and (n-2) have not been computed, then the routine is recersively called for overtone "<>
	"'(n-1).  If no solutions exist for mode l, then the first two overtones of (l-1) and (l-2) "<>
	"are used to extrapolate initial guesses for these modes.\n\n"<>
	"Options:\n"<>
	"\t SpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set via a call to SetSpinWeight before any KerrQNM\n"<>
	"\t\t function call.\n"<>
	"\t ModePrecision\[Rule]24\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t The radial solver make use of numerical derivatives.  The relative step size\n"<>
	"\t\t is set to \!\(\*SuperscriptBox[\(10\), \(JacobianStep\)]\).\n"<>
	"\t RadialCFMinDepth\[Rule]300\n"<>
	"\t\t The minimum value for the Radial Continued Fraction Depth.\n"<>
	"\t RadialCFDepth\[Rule]1\n"<>
	"\t\t The initial Radial Continued Fraction Depth is usually taken from the prior \n"<>
	"\t\t solution in the sequence.  Fractional values reduce this initial value by that\n"<>
	"\t\t fraction.  Integer values larger than RadialCFMinDepth are used as the new\n"<>
	"\t\t  initial value.\n"<>
	"\t RadialRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial Newton iterations.\n"<>
	"\t Root\[Epsilon]\[Rule]-14\n"<>
	"\t\t \!\(\*SuperscriptBox[\(10\), \(Root\[Epsilon]\)]\) is the maximum allowed magnitude for the Radial Continued Fraction.\n"<>
	"\t SchDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Schwarzschile iterations.\n"<>
	"\t RadialDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial Newton iterations.\n"


PlotSchQNM::usage=
	"PlotSchQNM[l] plots both the \"positive\" and \"negative\" frequency QNMs.  "<>
	"By default, the gravitational modes are plotted, but the PlotSpinWeight option "<>
	"can be set to change this.\n\n"<>
	"Options:\n"<>
	"\t PlotSpinWeight->-2 : -2,-1,0\n\n"<>
	"PlotSchQNM also take all of the options available to ListPlot.\n"


(* ::Subsection:: *)
(*Reserved Globals*)


Protect[PlotSpinWeight];


Begin["`Private`"]


(* ::Section:: *)
(*Kerr QNM methods*)


(* ::Subsection:: *)
(*Adaptive Bisection sequencer*)


Options[KerrQNMSequence]=Options[KerrModes`Private`KerrModeSequence];


KerrQNMSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{ModeSavePrecision=$MinPrecision,saneopts},
	KerrQNMSequence::spinweight="SpinWeight has not been set.  You must call SetSpinWeight[]";
	KerrQNMSequence::argl="The order l is set to `1`, but must be \[GreaterEqual] |`2`|";
	KerrQNMSequence::argm="The index m is set to `1`, but must be between -`2` and `2`";
	KerrQNMSequence::argn="The overtone n is set to `1`, but must be \[GreaterEqual] 0";
	If[OptionValue[SpinWeight]==Null[],Message[SchwarzschildQNM::spinweight];Abort[]];
	If[l<Abs[OptionValue[SpinWeight]],
			Message[SchwarzschildQNM::argl,l,OptionValue[SpinWeight]];Abort[]];
	If[Abs[m]>l,
			Message[SchwarzschildQNM::argm,m,l];Abort[]];
	If[n<0,Message[SchwarzschildQNM::argn,n];Abort[]];
	(* saneopts ensures options set via SetOptions[KerQNMSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrQNMSequence],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`KerrModeSequence[l,m,n,\[Epsilon],FilterRules[saneopts,Options[KerrQNMSequence]]],
				$MinPrecision=ModeSavePrecision;Abort[]];
	$MinPrecision=ModeSavePrecision;
]


(* ::Section:: *)
(*Initial Guesses*)


Options[SchwarzschildQNM]=Options[KerrModes`Private`SchwarzschildMode];


SchwarzschildQNM[l_Integer,n_Integer,
				opts:OptionsPattern[]]:=
Module[{SavePrecision=$MinPrecision,saneopts},
	SchwarzschildQNM::spinweight="SpinWeight has not been set.  You must call SetSpinWeight[]";
	SchwarzschildQNM::argl="The order l is set to `1`, but must be \[GreaterEqual] |`2`|";
	SchwarzschildQNM::argn="The overtone n is set to `1`, but must be \[GreaterEqual] 0";
	If[OptionValue[SpinWeight]==Null[],Message[SchwarzschildQNM::spinweight];Abort[]];
	If[l<Abs[OptionValue[SpinWeight]],
			Message[SchwarzschildQNM::argl,l,OptionValue[SpinWeight]];Abort[]];
	If[n<0,Message[SchwarzschildQNM::argn,n];Abort[]];
	(* saneopts ensures options set via SetOptions[SchwarzschildQNM,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[SchwarzschildQNM],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`SchwarzschildMode[l,n,FilterRules[saneopts,Options[SchwarzschildQNM]]],
				$MinPrecision=SavePrecision;Abort[]];
	$MinPrecision=SavePrecision;
]


(* ::Section:: *)
(*Graphics*)


Options[PlotSchQNM]=Union[{PlotSpinWeight->-2},Options[KerrModes`Private`PlotSchModes]];


PlotSchQNM[l_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[PlotSpinWeight],ptable},
	PlotSchQNM::spinweight="Invalid QNM Spin Weight : `1`";
	Switch[s,
		   -2,ptable:=Global`SchQNMTable,
		   -1,ptable:=Global`SchQNMeTable,
		    0,ptable:=Global`SchQNMsTable,
			_,Message[PlotSchQNM::spinweight,s];Abort[]
		  ];
	KerrModes`Private`PlotSchModes[l,PlotTable->ptable,
									FilterRules[{opts},Options[KerrModes`Private`PlotSchModes]]]
]


(* ::Section::Closed:: *)
(*End of KerrQNM Package*)


End[] (* `Private` *)


EndPackage[]
