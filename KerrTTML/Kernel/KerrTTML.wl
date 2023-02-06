(* ::Package:: *)

(* ::Title:: *)
(*Left Total Transmission Modes of Kerr*)


(* ::Section::Closed:: *)
(*Begin KerrTTML Package*)


BeginPackage["KerrTTML`",{"KerrModes`"}]


Unprotect[KerrTTMLDebug];
KerrTTMLDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[KerrTTMLDebug,Unprotect["KerrTTML`*"];Unprotect["KerrTTML`Private`*"]];
Protect[KerrTTMLDebug];


(* ::Section::Closed:: *)
(*Documentation of External Functions in KerrModes Namespace*)


SelectMode::usage=
"SelectMode[type] "<>
"specifies the type of mode function when solving for TTMs.  This function must be "<>
"called to set the type of the mode function before any calls to functions in the TTML` "<>
"package."


(* ::Section::Closed:: *)
(*Definitions for KerrModes Namespace*)


Begin["KerrModes`Private`"]


(* ::Subsection::Closed:: *)
(*Define the TTML Radial Equation recurrence relation coefficients*)


rp=1+Sqrt[1-a^2]; rm=1-Sqrt[1-a^2]; \[Sigma]p=(2\[Omega]*rp-m a)/(rp-rm); \[Sigma]m=(2\[Omega]*rm-m a)/(rp-rm);


\[Zeta]=I*\[Omega]; \[Xi]=I*\[Sigma]p; \[Eta]=-I*\[Sigma]m;


p=(rp-rm)\[Zeta]/2; \[Alpha]=1+s+\[Xi]+\[Eta]-2\[Zeta]+s*I*\[Omega]/\[Zeta]; \[Gamma]=1+s+2\[Eta]; \[Delta]=1+s+2\[Xi]; \[Sigma]=Alm+a^2*\[Omega]^2-8\[Omega]^2+p(2\[Alpha]+\[Gamma]-\[Delta])+(1+s-(\[Gamma]+\[Delta])/2)(s+(\[Gamma]+\[Delta])/2);


D0=Simplify[\[Delta]]; D1=Simplify[-2+4p-2\[Alpha]+\[Gamma]-\[Delta]]; D2=Simplify[2+2\[Alpha]-\[Gamma]];
D3=Simplify[4p*\[Alpha]-\[Alpha]*\[Delta]-\[Sigma]]; D4=Simplify[\[Alpha](1+\[Alpha]-\[Gamma])];


\[Alpha]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D0+1)n+D0;
\[Beta]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=-2n^2+(D1+2)n+D3;
\[Gamma]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D2-3)n+D4-D2+2;


Clear[rp,rm,\[Sigma]p,\[Sigma]m,\[Zeta],\[Xi],\[Eta],p,\[Alpha],\[Gamma],\[Delta],\[Sigma],D0,D1,D2,D3,D4];
If[!KerrTTMLDebug,Protect[\[Alpha]r,\[Beta]r,\[Gamma]r]];


(* ::Subsection::Closed:: *)
(*Approximation of remainder at the Nmax element of the continued fraction*)


RadialCFRemainder[s_Integer,m_Integer,a_Rational|a_Integer,
				  Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:=
Module[{C12tmp,C1,C2,C3,C4,C5,Rem,Err},
	C12tmp=-4I Sqrt[1-a^2]\[Omega];
	C1 = Sqrt[C12tmp];
	If[Re[C1]>0,C1=-C1];
	C2=-(3-4s+8I(\[Omega]+Sqrt[1-a^2]\[Omega]))/4;
	C3=(3 + 16*Alm + 16*s + 16*s^2 + 32*a*m*\[Omega] - 256*\[Omega]^2 + 80*a^2*\[Omega]^2 - 256*Sqrt[1 - a^2]*\[Omega]^2 - (32*I)*(-5*Sqrt[1 - a^2]*\[Omega] + 2*Sqrt[1 - a^2]*s*\[Omega]))/(32*C1);
	C4=((-3*I)*Sqrt[1 - a^2] - (16*I)*Sqrt[1 - a^2]*s^2 - 288*\[Omega] + 288*a^2*\[Omega] +(96*I)*a*Sqrt[1 - a^2]*m*\[Omega] - (896*I)*\[Omega]^2 + (896*I)*a^2*\[Omega]^2 - (896*I)*Sqrt[1 - a^2]*\[Omega]^2 + (816*I)*a^2*Sqrt[1 - a^2]*\[Omega]^2 - 256*a*m*\[Omega]^2 +256*a^3*m*\[Omega]^2 - 512*a*Sqrt[1 - a^2]*m*\[Omega]^2 + 2048*\[Omega]^3 - 2176*a^2*\[Omega]^3 + 128*a^4*\[Omega]^3 + 2048*Sqrt[1 - a^2]*\[Omega]^3 - 1024*a^2*Sqrt[1 - a^2]*\[Omega]^3 +16*s*((-I)*Sqrt[1 - a^2] + 8*\[Omega] - 8*a^2*\[Omega]) +16*Alm*((-I)*Sqrt[1 - a^2] + 8*(-1 + a^2)*\[Omega]))/(256*(-1 + a^2)*\[Omega]);
	C5=(63*I - (256*I)*Alm^2 - (256*I)*s^4 - 1344*Sqrt[1 - a^2]*\[Omega] + (576*I)*a*m*\[Omega] -(59904*I)*\[Omega]^2 + (56736*I)*a^2*\[Omega]^2 - (1536*I)*Sqrt[1 - a^2]*\[Omega]^2 - 30720*a*Sqrt[1 - a^2]*m*\[Omega]^2 + (7168*I)*a^2*m^2*\[Omega]^2 + 147456*\[Omega]^3 - 147456*a^2*\[Omega]^3 + 147456*Sqrt[1 - a^2]*\[Omega]^3 - 35840*a^2*Sqrt[1 - a^2]*\[Omega]^3 - (81920*I)*a*m*\[Omega]^3 + (15360*I)*a^3*m*\[Omega]^3 - (81920*I)*a*Sqrt[1 - a^2]*m*\[Omega]^3 + (262144*I)*\[Omega]^4 - (172032*I)*a^2*\[Omega]^4 + (1792*I)*a^4*\[Omega]^4 + (262144*I)*Sqrt[1 - a^2]*\[Omega]^4 - (40960*I)*a^2*Sqrt[1 - a^2]*\[Omega]^4 - 512*s^3*(I + 4*Sqrt[1 - a^2]*\[Omega]) - (32*I)*s^2*(-1 + 32*((-9*I)*Sqrt[1 - a^2] + a*m)*\[Omega] +  16*(-24 + 13*a^2 + 16*Sqrt[1 - a^2])*\[Omega]^2) - 32*s*(-9*I + 4*(59*Sqrt[1 - a^2] - (24*I)*a*m)*\[Omega] + 16*((13*I)*a^2 + (8*I)*(-1 + 2*Sqrt[1 - a^2]) + 8*a*Sqrt[1 - a^2]*m)*\[Omega]^2 + 64*(-16*(1 + Sqrt[1 - a^2]) + a^2*(16 + 5*Sqrt[1 - a^2]))*\[Omega]^3) - (32*I)*Alm*(-9 + 16*s^2 + ((-224*I)*Sqrt[1 - a^2] - 96*a*m)*\[Omega] + 16*(16 - 11*a^2 + 16*Sqrt[1 - a^2])*\[Omega]^2 + 16*s*(1 - (4*I)*Sqrt[1 - a^2]*\[Omega])))/(8192*Sqrt[1 - a^2]*C1*\[Omega]);
	Rem=1+C1/Sqrt[Nmax]+C2/Nmax+C3/(Nmax^(3/2))+C4/Nmax^2+C5/(Nmax^(5/2));
	Err=Abs[(C5/(Nmax^(5/2)))/Rem];
	{Rem,Err}
]


If[!KerrTTMLDebug,Protect[RadialCFRemainder]];


(* ::Subsection::Closed:: *)
(*AsymptoteFunction*)


AsymptoteFunction[s_Integer,l_Integer,m_Integer,a_Rational|a_Integer]:=
Module[{\[Omega],Alm},
	\[Omega]=SetPrecision[-(2^(2/3)3^(1/3)a^(-4/3)+(4/3l-2)/a
		-(39+6l-2l^2)/(18 2^(2/3) 3^(1/3))a^(-2/3)
		+(6561-4545l+171l^2-38l^3)/(1296 2^(1/3) 3^(2/3))a^(-1/3)
		+(-45+27l-6l^3+l^4)/384)I,$MinPrecision];
	Alm=SetPrecision[((2l-3)(I \[Omega] a)
		+ 17/4+3/2*l-1/2l^2
		+ (-24 2^(1/3) 3^(2/3) a^(1/3)(177-109l-9l^2+2l^3)
		+ 2^(2/3)3^(1/3)a^(2/3)(-1875-969l+476l^2-102l^3+17l^4)
		- 2a(1803+6655l-7857l^2+1746l^3))/2304),$MinPrecision];
	{\[Omega],Alm}
]


(* ::Subsection::Closed:: *)
(*Set SpinWeight, SelectMode, and Data-Variable Names*)


SetSpinWeight[s_Integer]:=
Module[{},
	SetSpinWeight::spinweight="Invalid TTML Spin Weight : `1`";
	SetSpinWeight::confirm="All KerrMode routines (TTML) set for Spin-Weight s = `1`";
	If[s!=0 && s!=-1 && s!=-2,Message[SetSpinWeight::spinweight,s]];
	modeName:=GetKerrName[TTML,s];
	SchTable:=GetSchName[TTML,s];
	SetOptions[KerrTTML`SchwarzschildTTML,SpinWeight->s];
	SetOptions[KerrTTML`KerrTTMLSequence,SpinWeight->s];
	SetOptions[KerrTTML`KerrTTMLRefineSequence,SpinWeight->s];
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
	Print[Style[StringForm["All KerrMode routines (TTML) set for Spin-Weight s = `1`",s],{Medium,Darker[Green]}]];
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
(*Documentation of External Functions in KerrTTML Namespace*)


KerrTTMLSequence::usage=
"KerrTTMLSequence[l,m,n,\[Epsilon]] "<>
"extend a sequence of Kerr left total-transmission mode solutions with harmonic "<>
"index l, azimuthal index m, and overtone n.  The solution accuracy is refined until "<>
"it is below an absolute accuracy of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\)."


KerrTTMLRefineSequence::usage=
"KerrTTMLRefineSequence[l,m,n,\[Epsilon]] "<>
"refine a sequence of Kerr left total-transmission mode solutions with harmonic "<>
"index l, azimuthal index m, and overtone n.  The solution accuracy is refined until "<>
"it is below an absolute accuracy of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\)."


SchwarzschildTTML::usage=
"SchwarzschildTTML[l,n] "<>
"computes the Schwarzschild left total-transmission mode solution for overtone n "<>
"and harmonic index l.  The mode is computed to an accuracy of \!\(\*SuperscriptBox[\(10\), \(-14\)]\) and stored in "<>
"the global variable SchTable[l,n]."


PlotSchTTML::usage=
"PlotSchTTML[l] "<>
"plots both the \"positive\" and \"negative\" frequency Schwarzschild left "<>
"total-transmission modes with harmonic mode l."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Kerr TTML methods*)


(* ::Subsection::Closed:: *)
(*Adaptive Bisection sequencer*)


Options[KerrTTMLSequence]=Options[KerrModes`Private`KerrModeSequence];


KerrTTMLSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{ModeSavePrecision=$MinPrecision,saneopts},
	KerrTTMLSequence::spinweight="SpinWeight has not been set.  You must call SetSpinWeight[]";
	KerrTTMLSequence::argl="The order l is set to `1`, but must be \[GreaterEqual] |`2`|";
	KerrTTMLSequence::argm="The index m is set to `1`, but must be between -`2` and `2`";
	KerrTTMLSequence::argn="The overtone n is set to `1`, but must be \[GreaterEqual] 0";
	If[OptionValue[SpinWeight]==Null[],Message[SchwarzschildTTML::spinweight];Abort[]];
	If[l<Abs[OptionValue[SpinWeight]],
			Message[KerrTTMLSequence::argl,l,OptionValue[SpinWeight]];Abort[]];
	If[Abs[m]>l,
			Message[KerrTTMLSequence::argm,m,l];Abort[]];
	If[n<0,Message[KerrTTMLSequence::argn,n];Abort[]];
	(* saneopts ensures options set via SetOptions[KerQNMSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrTTMLSequence],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`KerrModeSequence[l,m,n,\[Epsilon],FilterRules[saneopts,Options[KerrTTMLSequence]]],
				$MinPrecision=ModeSavePrecision;Abort[]];
	$MinPrecision=ModeSavePrecision;
]


Options[KerrTTMLRefineSequence]=Options[KerrModes`Private`KerrModeRefineSequence];


KerrTTMLRefineSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{SavePrecision=$MinPrecision,saneopts},
	(* saneopts ensures options set via SetOptions[KerTTMLRefineSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrTTMLRefineSequence],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`KerrModeRefineSequence[l,m,n,\[Epsilon],saneopts],
				$MinPrecision=SavePrecision;Abort[]];
	$MinPrecision=SavePrecision;
]


(* ::Section::Closed:: *)
(*Initial Guesses*)


Options[SchwarzschildTTML]=Options[KerrModes`Private`SchwarzschildMode];


SchwarzschildTTML[l_Integer,n_Integer,
				opts:OptionsPattern[]]:=
Module[{SavePrecision=$MinPrecision,saneopts},
	SchwarzschildTTML::spinweight="SpinWeight has not been set.  You must call SetSpinWeight[]";
	SchwarzschildTTML::argl="The order l is set to `1`, but must be \[GreaterEqual] |`2`|";
	SchwarzschildTTML::argn="The overtone n is set to `1`, but must be \[GreaterEqual] 0";
	If[OptionValue[SpinWeight]==Null[],Message[SchwarzschildTTML::spinweight];Abort[]];
	If[l<Abs[OptionValue[SpinWeight]],
			Message[SchwarzschildTTML::argl,l,OptionValue[SpinWeight]];Abort[]];
	If[n<0,Message[SchwarzschildTTML::argn,n];Abort[]];
	(* saneopts ensures options set via SetOptions[SchwarzschildTTML,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[SchwarzschildTTML],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`SchwarzschildMode[l,n,FilterRules[saneopts,Options[SchwarzschildTTML]]],
				$MinPrecision=SavePrecision;Abort[]];
	$MinPrecision=SavePrecision;
]


(* ::Section::Closed:: *)
(*Graphics*)


Options[PlotSchTTML]=Union[{SpinWeight->-2},Options[KerrModes`Private`PlotSchModes]];


PlotSchTTML[l_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],ptable},
	PlotSchTTML::spinweight="Invalid TTML Spin Weight : `1`";
	ptable:=GetSchName[TTML,s];
	KerrModes`Private`PlotSchModes[l,PlotTable->ptable,
									FilterRules[{opts},Options[KerrModes`Private`PlotSchModes]]]
]


(* ::Section::Closed:: *)
(*End of KerrTTML Package*)


End[] (* `Private` *)


EndPackage[]
