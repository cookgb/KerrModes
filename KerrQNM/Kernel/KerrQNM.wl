(* ::Package:: *)

(* ::Title:: *)
(*QuasiNormal Modes of Kerr*)


(* ::Section::Closed:: *)
(*Begin KerrQNM Package*)


BeginPackage["KerrQNM`",{"KerrModes`"}]


Unprotect[KerrQNMDebug];
KerrQNMDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[KerrQNMDebug,Unprotect["KerrQNM`*"];Unprotect["KerrQNM`Private`*"]];
Protect[KerrQNMDebug];


(* ::Section:: *)
(*Documentation of External Functions in KerrModes Namespace*)


(* ::Section::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Set SpinWeight, and Data-Variable Names*)


SetSpinWeight[s_Integer]:=
Module[{},
(* Mode selection preset for QNM *)
	Unprotect[RunCFConvergence];RunCFConvergence=True;Protect[RunCFConvergence];
	ModeFunction[n_,s1_,m_,a_,Alm_,\[Omega]_,Nrcf_]=RadialCFRem[n,s1,m,a,Alm,\[Omega],Nrcf];
	SetSpinWeight::spinweight="Invalid QNM Spin Weight : `1`";
	SetSpinWeight::confirm="All KerrMode routines (QNM) set for Spin-Weight s = `1`";
	If[s!=0 && s!=-1 && s!=-2,Message[SetSpinWeight::spinweight,s]];
	modeName:=GetKerrName[QNM,s];
	SchTable:=GetSchName[QNM,s];
	SetOptions[KerrQNM`SchwarzschildQNM,SpinWeight->s];
	SetOptions[KerrQNM`KerrQNMSequence,SpinWeight->s];
	SetOptions[KerrQNM`KerrQNMRefineSequence,SpinWeight->s];
	SetOptions[KerrModes`KerrOmegaList,SpinWeight->s];
	SetOptions[KerrModes`KerrOmegaListS,SpinWeight->s];
	SetOptions[KerrModes`SWSFLists,SpinWeight->s];
	SetOptions[KerrModes`AngularModeRealPath,SpinWeight->s];
	Print[Style[StringForm["All KerrMode routines (QNM) set for Spin-Weight s = `1`",s],{Medium,Darker[Green]}]];
]


If[!KerrModeDebug,Protect[SetSpinWeight]];


End[] (* KerrModes`Private` *)


(* ::Section::Closed:: *)
(*Documentation of External Functions in KerrQNM Namespace*)


KerrQNMSequence::usage=
"KerrQNMSequence[l,m,n,\[Epsilon]]  "<>
"extend a sequence of Kerr QNM solutions with harmonic index l, azimuthal index m, "<>
"and overtone n.  The solution accuracy is refined until it is below an absolute "<>
"accuracy of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\)."


KerrQNMRefineSequence::usage=
"KerrQNMRefineSequence[l,m,n,\[Epsilon]]  "<>
"refine a sequence of Kerr QNM solutions with harmonic index l, azimuthal index m, "<>
"and overtone n.  The solution accuracy is refined until it is below an absolute "<>
"accuracy of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\)."


SchwarzschildQNM::usage=
"SchwarzschildQNM[l,n]  "<>
"computes the Schwarzschild quasi-normal mode solution for overtone n and harmonic "<>
"index l.  The mode is computed to an accuracy of \!\(\*SuperscriptBox[\(10\), \(-14\)]\) and stored in the global "<>
"variable SchTable[l,n]."


PlotSchQNM::usage=
"PlotSchQNM[l] "<>
"plots both the \"positive\" and \"negative\" frequency Schwarzschild QNMs "<>
"with harmonic mode l."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Kerr QNM methods*)


(* ::Subsection::Closed:: *)
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


Options[KerrQNMRefineSequence]=Options[KerrModes`Private`KerrModeRefineSequence];


KerrQNMRefineSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{SavePrecision=$MinPrecision,saneopts},
	(* saneopts ensures options set via SetOptions[KerQNMRefineSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrQNMRefineSequence],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrModes`Private`KerrModeRefineSequence[l,m,n,\[Epsilon],saneopts],
				$MinPrecision=SavePrecision;Abort[]];
	$MinPrecision=SavePrecision;
]


(* ::Section::Closed:: *)
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


(* ::Section::Closed:: *)
(*Graphics*)


Options[PlotSchQNM]=Union[{SpinWeight->-2},Options[KerrModes`Private`PlotSchModes]];


PlotSchQNM[l_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[SpinWeight],ptable},
	PlotSchQNM::spinweight="Invalid QNM Spin Weight : `1`";
	ptable:=KerrModes`Private`GetSchName[QNM,s];
	KerrModes`Private`PlotSchModes[l,PlotTable->ptable,
									FilterRules[{opts},Options[KerrModes`Private`PlotSchModes]]]
]


(* ::Section::Closed:: *)
(*End of KerrQNM Package*)


End[] (* `Private` *)


EndPackage[]
