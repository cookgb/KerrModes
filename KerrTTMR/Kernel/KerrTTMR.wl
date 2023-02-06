(* ::Package:: *)

(* ::Title:: *)
(*Right Total Transmission Modes of Kerr*)


(* ::Section::Closed:: *)
(*Begin KerrTTMR Package*)


BeginPackage["KerrTTMR`",{"KerrModes`"}]


Unprotect[KerrTTMRDebug];
KerrTTMRDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[KerrTTMRDebug,Unprotect["KerrTTMR`*"];Unprotect["KerrTTMR`Private`*"]];
Protect[KerrTTMRDebug];


(* ::Section::Closed:: *)
(*Documentation of External Functions in KerrModes Namespace*)


SelectMode::usage=
"SelectMode[type] "<>
"specifies the type of mode function when solving for TTMs.  This function must be "<>
"called to set the type of the mode function before any calls to functions in the TTMR` "<>
"package."


(* ::Section::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*AsymptoteFunction*)


AsymptoteFunction[s_Integer,l_Integer,m_Integer,a_Rational|a_Integer]:=
Module[{\[Omega],Alm},
	\[Omega]=SetPrecision[-(2^(2/3)3^(1/3)a^(-4/3)+(4/3l-2)/a
		-(39+6l-2l^2)/(18 2^(2/3) 3^(1/3))a^(-2/3)
		+(6561-4545l+171l^2-38l^3)/(1296 2^(1/3) 3^(2/3))a^(-1/3)
		+(-45+27l-6l^3+l^4)/384)I,$MinPrecision];
	Alm=SetPrecision[((2l-3)(I \[Omega] a)
		+ 1/4+3/2*l-1/2l^2
		+ (-24 2^(1/3) 3^(2/3) a^(1/3)(177-109l-9l^2+2l^3)
		+ 2^(2/3)3^(1/3)a^(2/3)(-1875-969l+476l^2-102l^3+17l^4)
		- 2a(1803+6655l-7857l^2+1746l^3))/2304),$MinPrecision];
	{\[Omega],Alm}
]


(* ::Subsection::Closed:: *)
(*Set SpinWeight, SelectMode, and Data-Variable Names*)


SetSpinWeight[s_Integer]:=
Module[{},
	SetSpinWeight::spinweight="Invalid TTMR Spin Weight : `1`";
	SetSpinWeight::confirm="All KerrMode routines (TTMR) set for Spin-Weight s = `1`";
	If[s!=0 && s!=1 && s!=2,Message[SetSpinWeight::spinweight,s]];
	modeName:=GetKerrName[TTMR,s];
	SchTable:=GetSchName[TTMR,s];
	SetOptions[KerrTTMR`SchwarzschildTTMR,SpinWeight->s];
	SetOptions[KerrTTMR`KerrTTMRSequence,SpinWeight->s];
	SetOptions[KerrTTMR`KerrTTMRRefineSequence,SpinWeight->s];
	SetOptions[KerrModes`SchwarzschildOmega,SpinWeight->s];
	SetOptions[KerrModes`KerrOmegaList,SpinWeight->s];
	SetOptions[KerrModes`KerrOmegaListS,SpinWeight->s];
	SetOptions[KerrModes`KerraOmegaList,SpinWeight->s];
	SetOptions[KerrModes`KerraOmegaListS,SpinWeight->s];
	SetOptions[KerrModes`KerrAList,SpinWeight->s];
	SetOptions[KerrModes`KerrAListS,SpinWeight->s];
	SetOptions[KerrModes`ModePlotOmega,SpinWeight->s];
	SetOptions[KerrModes`ModePlotA,SpinWeight->s];
	SetOptions[KerrModes`ModePlotOmegaTones,SpinWeight->s];
	SetOptions[KerrModes`SWSFLists,SpinWeight->s];
	SetOptions[KerrModes`AngularModeRealPath,SpinWeight->s];
	Print[Style[StringForm["All KerrMode routines (TTMR) set for Spin-Weight s = `1`",s],{Medium,Darker[Green]}]];
]


SelectMode[funcmodechoice_]:=
Module[{},
	SelectMode::polynomial="Mode set to find Polynomial solutions";
	SelectMode::continuedfraction="Mode set to find Continued Fraction solutions";
	SelectMode::choose="Must Choose PolynomialMode or ContinuedFractionMode for respective solutions";
	Switch[funcmodechoice,
			PolynomialMode,
				Print[Style["Mode set to find Polynomial solutions",{Medium,Darker[Green]}]];
				Unprotect[RunCFConvergence];RunCFConvergence=False;Protect[RunCFConvergence];
				ModeFunction[n_,s_,m_,a_,Alm_,\[Omega]_,Nrcf_]=Starobinsky[s,m,a,Alm,\[Omega]],
			ContinuedFractionMode,
				Print[Style["Mode set to find Continued Fraction solutions",{Medium,Darker[Green]}]];
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


KerrTTMRSequence::usage=
"KerrTTMRSequence[l,m,n,\[Epsilon]] "<>
"extend a sequence of Kerr right total-transmission mode solutions with harmonic "<>
"index l, azimuthal index m, and overtone n.  The solution accuracy is refined until "<>
"it is below an absolute accuracy of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\)."


KerrTTMRRefineSequence::usage=
"KerrTTMRRefineSequence[l,m,n,\[Epsilon]] "<>
"refine a sequence of Kerr right total-transmission mode solutions with harmonic "<>
"index l, azimuthal index m, and overtone n.  The solution accuracy is refined until "<>
"it is below an absolute accuracy of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\)."


SchwarzschildTTMR::usage=
"SchwarzschildTTMR[l,n] "<>
"computes the Schwarzschild right total-transmission mode solution for overtone n "<>
"and harmonic index l.  The mode is computed to an accuracy of \!\(\*SuperscriptBox[\(10\), \(-14\)]\) and stored in "<>
"the global variable SchTable[l,n]."


PlotSchTTMR::usage=
"PlotSchTTMR[l] "<>
"plots both the \"positive\" and \"negative\" frequency Schwarzschild right "<>
"total-transmission modes with harmonic mode l."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Kerr TTMR methods*)


(* ::Subsection::Closed:: *)
(*Adaptive Bisection sequencer*)


Options[KerrTTMRSequence]=Options[KerrModes`Private`KerrModeSequence];


KerrTTMRSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{ModeSavePrecision=$MinPrecision,saneopts},
	KerrTTMRSequence::spinweight="SpinWeight has not been set.  You must call SetSpinWeight[]";
	KerrTTMRSequence::argl="The order l is set to `1`, but must be \[GreaterEqual] |`2`|";
	KerrTTMRSequence::argm="The index m is set to `1`, but must be between -`2` and `2`";
	KerrTTMRSequence::argn="The overtone n is set to `1`, but must be \[GreaterEqual] 0";
	If[OptionValue[SpinWeight]==Null[],Message[SchwarzschildTTMR::spinweight];Abort[]];
	If[l<Abs[OptionValue[SpinWeight]],
			Message[KerrTTMRSequence::argl,l,OptionValue[SpinWeight]];Abort[]];
	If[Abs[m]>l,
			Message[KerrTTMRSequence::argm,m,l];Abort[]];
	If[n<0,Message[KerrTTMRSequence::argn,n];Abort[]];
	(* saneopts ensures options set via SetOptions[KerQNMSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrTTMRSequence],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`KerrModeSequence[l,m,n,\[Epsilon],FilterRules[saneopts,Options[KerrTTMRSequence]]],
				$MinPrecision=ModeSavePrecision;Abort[]];
	$MinPrecision=ModeSavePrecision;
]


Options[KerrTTMRRefineSequence]=Options[KerrModes`Private`KerrModeRefineSequence];


KerrTTMRRefineSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{SavePrecision=$MinPrecision,saneopts},
	(* saneopts ensures options set via SetOptions[KerTTMRRefineSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrTTMRRefineSequence],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`KerrModeRefineSequence[l,m,n,\[Epsilon],saneopts],
				$MinPrecision=SavePrecision;Abort[]];
	$MinPrecision=SavePrecision;
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


Options[PlotSchTTMR]=Union[{SpinWeight->2},Options[KerrModes`Private`PlotSchModes]];


PlotSchTTMR[l_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],ptable},
	PlotSchTTMR::spinweight="Invalid TTMR Spin Weight : `1`";
	ptable:=GetSchName[TTMR,s];
	KerrModes`Private`PlotSchModes[l,PlotTable->ptable,
									FilterRules[{opts},Options[KerrModes`Private`PlotSchModes]]]
]


(* ::Section::Closed:: *)
(*End of KerrTTMR Package*)


End[] (* `Private` *)


EndPackage[]
