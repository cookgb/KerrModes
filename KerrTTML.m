(* ::Package:: *)

(* ::Title:: *)
(*Total Transmission Modes of Kerr*)


(* ::Section:: *)
(*Begin KerrTTML Package*)


BeginPackage["KerrTTML`",{"KerrModes`"}]


Unprotect[KerrTTMLDebug];
KerrTTMLDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[KerrTTMLDebug,Unprotect["KerrTTML`*"];Unprotect["KerrTTML`Private`*"]];
Protect[KerrTTMLDebug];


(* ::Section:: *)
(*Documentation of External Functions in KerrModes Namespace*)


SetSpinWeight::usage=
	"SetSpinWeight[s] sets the value of the spin-weight used in all subsequent "<>
	"TTML computations:\n"<>
	"\t s=-2 : Gravitational perturbations\n"<>
	"\t s=-1 : Electro-Magnetic perturbations\n"<>
	"\t s= 0 : Scalar perturbations."


(* ::Section:: *)
(*Definitions for KerrModes Namespace*)


Begin["KerrModes`Private`"]


(* ::Subsection:: *)
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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Set SpinWeight and Data-Variable Names*)


SetSpinWeight[s_Integer]:=
Module[{},
	SetSpinWeight::spinweight="Invalid TTML Spin Weight : `1`";
	Switch[s,
		   -2,modeName:=Global`KerrTTML; SchTable:=Global`SchTTMLTable,
		   -1,modeName:=Global`KerrTTMLe; SchTable:=Global`SchTTMLeTable,
		    0,modeName:=Global`KerrTTMLs; SchTable:=Global`SchTTMLsTable,
			_,Message[SetSpinWeight::spinweight,s];Abort[]
		  ];
	SetOptions[KerrTTML`SchwarzschildTTML,SpinWeight->s];
	Print["All KerrMode routines (TTML) set for Spin-Weight s = ",s];
]


If[!KerrModeDebug,Protect[SetSpinWeight]];


End[] (* KerrModes`Private` *)


(* ::Section:: *)
(*Documentation of External Functions in KerrTTML Namespace*)


SchwarzschildTTML::usage=
	"SchwarzschildTTML[l,n] computes the Quasi-Normal Mode solutions for overtone n of mode l.  "<>
	"The mode is computed to an accuracy of \!\(\*SuperscriptBox[\(10\), \(-14\)]\).  For given 'l', if solutions with overtones "<>
	"(n-1) and (n-2) have not been computed, then the routine is recersively called for overtone "<>
	"'(n-1).  If no solutions exist for mode l, then the first two overtones of (l-1) and (l-2) "<>
	"are used to extrapolate initial guesses for these modes.\n\n"<>
	"Options:\n"<>
	"\t SpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set via a call to SetSpinWeight before any KerrTTML\n"<>
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


PlotSchTTML::usage=
	"SchTTMLTable[l] plots both the \"positive\" and \"negative\" frequency TTMLs.  "<>
	"By default, the gravitational modes are plotted, but the PlotSpinWeight option "<>
	"can be set to change this.\n\n"<>
	"Options:\n"<>
	"\t PlotSpinWeight->-2 : -2,-1,0\n\n"<>
	"SchTTMLTable also take all of the options available to ListPlot.\n"


(* ::Subsection:: *)
(*Reserved Globals*)


Protect[PlotSpinWeight];


Begin["`Private`"]


(* ::Section:: *)
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


(* ::Section:: *)
(*Graphics*)


Options[PlotSchTTML]=Union[{PlotSpinWeight->-2},Options[KerrModes`Private`PlotSchModes]];


PlotSchTTML[l_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[PlotSpinWeight],ptable},
	PlotSchTTML::spinweight="Invalid TTML Spin Weight : `1`";
	Switch[s,
		   -2,ptable:=Global`SchTTMLTable,
		   -1,ptable:=Global`SchTTMLeTable,
		    0,ptable:=Global`SchTTMLsTable,
			_,Message[PlotSchTTML::spinweight,s];Abort[]
		  ];
	KerrModes`Private`PlotSchModes[l,PlotTable->ptable,
									FilterRules[{opts},Options[KerrModes`Private`PlotSchModes]]]
]


(* ::Section::Closed:: *)
(*End of KerrTTML Package*)


End[] (* `Private` *)


EndPackage[]