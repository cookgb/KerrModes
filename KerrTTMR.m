(* ::Package:: *)

(* ::Title:: *)
(*Right Total Transmission Modes of Kerr*)


(* ::Section:: *)
(*Begin KerrTTMR Package*)


BeginPackage["KerrTTMR`",{"KerrModes`"}]


Unprotect[KerrTTMRDebug];
KerrTTMRDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[KerrTTMRDebug,Unprotect["KerrTTMR`*"];Unprotect["KerrTTMR`Private`*"]];
Protect[KerrTTMRDebug];


(* ::Section:: *)
(*Documentation of External Functions in KerrModes Namespace*)


SetSpinWeight::usage=
	"SetSpinWeight[s] sets the value of the spin-weight used in all subsequent "<>
	"TTMR computations:\n"<>
	"\t s=2 : Gravitational perturbations\n"<>
	"\t s=1 : Electro-Magnetic perturbations\n"<>
	"\t s=0 : Scalar perturbations."


SelectMode::usage=
	"Must set desired mode before beginning package"<>
	"Either PolynomialMode or ContinuedFractionMode"


(* ::Section:: *)
(*Definitions for KerrModes Namespace*)


Begin["KerrModes`Private`"]


(* ::Subsection::Closed:: *)
(*Define the TTMR Radial Equation recurrence relation coefficients*)


rp=1+Sqrt[1-a^2]; rm=1-Sqrt[1-a^2]; \[Sigma]p=(2\[Omega]*rp-m a)/(rp-rm); \[Sigma]m=(2\[Omega]*rm-m a)/(rp-rm);


\[Zeta]=(-I)*\[Omega]; \[Xi]=-s-I*\[Sigma]p; \[Eta]=-I*\[Sigma]m;


p=(rp-rm)\[Zeta]/2; \[Alpha]=1+s+\[Xi]+\[Eta]-2\[Zeta]+s*I*\[Omega]/\[Zeta]; \[Gamma]=1+s+2\[Eta]; \[Delta]=1+s+2\[Xi]; \[Sigma]=Alm+a^2*\[Omega]^2-8\[Omega]^2+p(2\[Alpha]+\[Gamma]-\[Delta])+(1+s-(\[Gamma]+\[Delta])/2)(s+(\[Gamma]+\[Delta])/2);


D0=Simplify[\[Delta]]; D1=Simplify[-2+4p-2\[Alpha]+\[Gamma]-\[Delta]]; D2=Simplify[2+2\[Alpha]-\[Gamma]];
D3=Simplify[4p*\[Alpha]-\[Alpha]*\[Delta]-\[Sigma]]; D4=Simplify[\[Alpha](1+\[Alpha]-\[Gamma])];


\[Alpha]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D0+1)n+D0;
\[Beta]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=-2n^2+(D1+2)n+D3;
\[Gamma]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D2-3)n+D4-D2+2;


Clear[rp,rm,\[Sigma]p,\[Sigma]m,\[Zeta],\[Xi],\[Eta],p,\[Alpha],\[Gamma],\[Delta],\[Sigma],D0,D1,D2,D3,D4];
If[!KerrTTMRDebug,Protect[\[Alpha]r,\[Beta]r,\[Gamma]r]];


(* ::Subsection::Closed:: *)
(*Approximation of remainder at the Nmax element of the continued fraction*)


RadialCFRemainder[s_Integer,m_Integer,a_Rational|a_Integer,
				  Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:=
Module[{C12tmp,C1,C2,C3,C4,C5,Rem,Err},
	C12tmp=4I Sqrt[1-a^2]\[Omega];
	C1 = Sqrt[C12tmp];
	If[Re[C1]>0,C1=-C1];
	C2=(-3-4s+(8I)(\[Omega]+Sqrt[1-a^2]\[Omega]))/4;
	C3=(3+16Alm+16s+16s^2+32a*m*\[Omega]-256\[Omega]^2+80a^2\[Omega]^2-256Sqrt[1-a^2]\[Omega]^2-(32I)(5Sqrt[1-a^2]\[Omega]+2Sqrt[1-a^2]s*\[Omega]))/(32C1);
	C4=(-3I-(16I)s(1+s)+288Sqrt[1-a^2]\[Omega]+16(Alm(-I+8Sqrt[1-a^2]\[Omega])+\[Omega](2a*m(3I+8(2+Sqrt[1-a^2])\[Omega])+a^2\[Omega](51I+8(8+Sqrt[1-a^2])\[Omega])+8(3Sqrt[1-a^2]s-(1+Sqrt[1-a^2])\[Omega](7I+16\[Omega])))))/(256Sqrt[1-a^2]\[Omega]);
	C5=(-63I+(256I)Alm^2+(256I)s^4-1344Sqrt[1-a^2]\[Omega]-(576I)a*m*\[Omega]+(59904I)\[Omega]^2-(56736I)a^2\[Omega]^2+(1536I)Sqrt[1-a^2]\[Omega]^2-30720a*Sqrt[1-a^2]m*\[Omega]^2-(7168I)a^2m^2\[Omega]^2+147456\[Omega]^3-147456a^2\[Omega]^3+147456Sqrt[1-a^2]\[Omega]^3-35840a^2Sqrt[1-a^2]\[Omega]^3+(81920I)a*m*\[Omega]^3-(15360I)a^3m*\[Omega]^3+(81920I)a*Sqrt[1-a^2]m*\[Omega]^3-(262144I)\[Omega]^4+(172032I)a^2\[Omega]^4-(1792I)a^4\[Omega]^4-(262144I)Sqrt[1-a^2]\[Omega]^4+(40960I)a^2Sqrt[1-a^2]\[Omega]^4+512s^3(I+4Sqrt[1-a^2]\[Omega])+(32I)s^2(-1+32((5I)Sqrt[1-a^2]+a*m)\[Omega]+16(-24+13a^2+16Sqrt[1-a^2])\[Omega]^2)+32s(-9I-4(53Sqrt[1-a^2]+(24I)a*m)\[Omega]+16((-35I)a^2+(8I)(5+2Sqrt[1-a^2])+8a*Sqrt[1-a^2]m)\[Omega]^2+64(-16(1+Sqrt[1-a^2])+a^2(16+5Sqrt[1-a^2]))\[Omega]^3)+32Alm((16I)s^2+16s(I+4Sqrt[1 - a^2]\[Omega])+I(-9+(32*I)(7Sqrt[1-a^2]+(3I)a*m)\[Omega]+16(16-11a^2+16Sqrt[1-a^2])\[Omega]^2)))/(8192Sqrt[1-a^2]C1*\[Omega]);
	Rem=1+C1/Sqrt[Nmax]+C2/Nmax+C3/(Nmax^(3/2))+C4/Nmax^2+C5/(Nmax^(5/2));
	Err=Abs[(C5/(Nmax^(5/2)))/Rem];
	{Rem,Err}
]


If[!KerrTTMRDebug,Protect[RadialCFRemainder]];


(* ::Subsection:: *)
(*Set SpinWeight, SelectMode, and Data-Variable Names*)


SetSpinWeight[s_Integer]:=
Module[{},
	SetSpinWeight::spinweight="Invalid TTMR Spin Weight : `1`";
	SetSpinWeight::confirm="All KerrMode routines (TTMR) set for Spin-Weight s = `1`";
	Switch[s,
		   2,modeName:=Global`KerrTTMR; SchTable:=Global`SchTTMRTable,
		   1,modeName:=Global`KerrTTMRe; SchTable:=Global`SchTTMReTable,
		   0,modeName:=Global`KerrTTMRs; SchTable:=Global`SchTTMRsTable,
		   _,Message[SetSpinWeight::spinweight,s];Abort[]
		  ];
	SetOptions[KerrTTMR`SchwarzschildTTMR,SpinWeight->s];
	SetOptions[KerrTTMR`KerrTTMRSequence,SpinWeight->s];
	SetOptions[KerrModes`KerrOmegaList,SpinWeight->s];
	SetOptions[KerrModes`KerrOmegaListS,SpinWeight->s];
	Print[Style[StringForm[SetSpinWeight::confirm,s],{Medium,Darker[Green]}]];
]


SelectMode[funcmodechoice_]:=
Module[{},
	SelectMode::polynomial="Mode set to find Polynomial solutions";
	SelectMode::continuedfraction="Mode set to find Continued Fraction solutions";
	SelectMode::choose="Must Choose PolynomialMode or ContinuedFractionMode for respective solutions";
	Switch[funcmodechoice,
			PolynomialMode,
				Print[Style[StringForm[SelectMode::polynomial],{Medium,Darker[Green]}]];
				Unprotect[RunCFConvergence];RunCFConvergence=False;Protect[RunCFConvergence];
				ModeFunction[n_,s_,m_,a_,Alm_,\[Omega]_,Nrcf_]=Starobinsky[s,m,a,Alm,\[Omega]],
			ContinuedFractionMode,
				Print[Style[StringForm[SelectMode::continuedfraction],{Medium,Darker[Green]}]];
				Unprotect[RunCFConvergence];RunCFConvergence=True;Protect[RunCFConvergence];
				ModeFunction[n_,s_,m_,a_,Alm_,\[Omega]_,Nrcf_]=RadialCFRem[n,s,m,a,Alm,\[Omega],Nrcf],
			_,
				Message[SelectMode::choose];
				Unprotect[RunCFConvergence];RunCFConvergence=Null[];Protect[RunCFConvergence];
				Abort[]
		];
]


If[!KerrModeDebug,Protect[SetSpinWeight]];


End[] (* KerrModes`Private` *)


(* ::Section::Closed:: *)
(*Documentation of External Functions in KerrTTMR Namespace*)


KerrTTMRSequence::usage=""


SchwarzschildTTMR::usage=
	"SchwarzschildTTMR[l,n] computes the Total Transmission Mode solutions for overtone n of mode l.  "<>
	"The mode is computed to an accuracy of \!\(\*SuperscriptBox[\(10\), \(-14\)]\).  For given 'l', if solutions with overtones "<>
	"(n-1) and (n-2) have not been computed, then the routine is recersively called for overtone "<>
	"'(n-1).  If no solutions exist for mode l, then the first two overtones of (l-1) and (l-2) "<>
	"are used to extrapolate initial guesses for these modes.\n\n"<>
	"Options:\n"<>
	"\t SpinWeight\[Rule]Null : 2,1,0\n"<>
	"\t\t The spin weight must be set via a call to SetSpinWeight before any KerrTTMR\n"<>
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


PlotSchTTMR::usage=
	"PlotSchTTMR[l] plots both the \"positive\" and \"negative\" frequency TTMRs.  "<>
	"By default, the gravitational modes are plotted, but the PlotSpinWeight option "<>
	"can be set to change this.\n\n"<>
	"Options:\n"<>
	"\t PlotSpinWeight->2 : 2,1,0\n\n"<>
	"PlotSchTTMR also take all of the options available to ListPlot.\n"


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[PlotSpinWeight];


Begin["`Private`"]


(* ::Section:: *)
(*Kerr TTMR methods*)


(* ::Subsection:: *)
(*Adaptive Bisection sequencer*)


Options[KerrTTMRSequence]=Options[KerrModes`Private`KerrModeSequence];


KerrTTMRSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{ModeSavePrecision=$MinPrecision,saneopts},
	(* saneopts ensures options set via SetOptions[KerQNMSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrQNMSequence],Except[Flatten[{opts}]]]]];
	CheckAbort[Print["Debug 0"];KerrModes`Private`KerrModeSequence[l,m,n,\[Epsilon],FilterRules[saneopts,Options[KerrQNMSequence]]],
				$MinPrecision=ModeSavePrecision;Abort[]];
	$MinPrecision=ModeSavePrecision;
]


(* ::Section::Closed:: *)
(*Initial Guesses*)


Options[SchwarzschildTTMR]=Options[KerrModes`Private`SchwarzschildMode];


SchwarzschildTTMR[l_Integer,n_Integer,
				opts:OptionsPattern[]]:=
Module[{SavePrecision=$MinPrecision,saneopts},
	SchwarzschildTTMR::spinweight="SpinWeight has not been set.  You must call SetSpinWeight[]";
	SchwarzschildTTMR::argl="The order l is set to `1`, but must be \[GreaterEqual] |`2`|";
	SchwarzschildTTMR::argn="The overtone n is set to `1`, but must be \[GreaterEqual] 0";
	If[OptionValue[SpinWeight]==Null[],Message[SchwarzschildTTMR::spinweight];Abort[]];
	If[l<Abs[OptionValue[SpinWeight]],
			Message[SchwarzschildTTMR::argl,l,OptionValue[SpinWeight]];Abort[]];
	If[n<0,Message[SchwarzschildTTMR::argn,n];Abort[]];
	(* saneopts ensures options set via SetOptions[SchwarzschildTTMR,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[SchwarzschildTTMR],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`SchwarzschildMode[l,n,FilterRules[saneopts,Options[SchwarzschildTTMR]]],
				$MinPrecision=SavePrecision;Abort[]];
	$MinPrecision=SavePrecision;
]


(* ::Section::Closed:: *)
(*Graphics*)


Options[PlotSchTTMR]=Union[{PlotSpinWeight->2},Options[KerrModes`Private`PlotSchModes]];


PlotSchTTMR[l_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[PlotSpinWeight],ptable},
	PlotSchTTMR::spinweight="Invalid TTMR Spin Weight : `1`";
	Switch[s,
		   2,ptable:=Global`SchTTMRTable,
		   1,ptable:=Global`SchTTMReTable,
		   0,ptable:=Global`SchTTMRsTable,
		   _,Message[PlotSchTTMR::spinweight,s];Abort[]
		  ];
	KerrModes`Private`PlotSchModes[l,PlotTable->ptable,
									FilterRules[{opts},Options[KerrModes`Private`PlotSchModes]]]
]


(* ::Section::Closed:: *)
(*End of KerrTTMR Package*)


End[] (* `Private` *)


EndPackage[]
