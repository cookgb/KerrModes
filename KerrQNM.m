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
(*Documentation of External Functions*)


SetSpinWeight::usage=
	"SetSpinWeight[s] sets the value of the spin-weight used in all subsequent "<>
	"QNM computations:\n"<>
	"\t s=-2 : Gravitational perturbations\n"<>
	"\t s=-1 : Electro-Magnetic perturbations\n"<>
	"\t s= 0 : Scalar perturbations."


(* ::Section:: *)
(*Definitions for KerrModes Namespace*)


Begin["KerrModes`Private`"]


(* ::Subsection:: *)
(*Define the QNM Radial Equation recurrence relation coefficients*)


rp=1+Sqrt[1-a^2];rm=1-Sqrt[1-a^2];\[Sigma]p=(2\[Omega] rp-m a)/(rp-rm);


C0=1-s-2I \[Sigma]p;C1=2(1+s-2I \[Omega]);C2=2I(rp-rm)\[Omega];C3=(1-4I \[Omega])(1+s-2I \[Sigma]p);C4=(4rp-a^2)\[Omega]^2 + 2((rp-rm)\[Sigma]p+I rp-2I(1+s))\[Omega]-Alm;


D0=Simplify[C0];D1=Simplify[-2C0-C1+C2];D2=Simplify[C0+C1];D3=Simplify[-C3+C4];D4=C3;


\[Alpha]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D0+1)n+D0;
\[Beta]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=-2n^2+(D1+2)n+D3;
\[Gamma]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D2-3)n+D4-D2+2;


Clear[rp,rm,\[Sigma]p,C0,C1,C2,C3,C4,D0,D1,D2,D3,D4];
If[!KerrQNMDebug,Protect[\[Alpha]r,\[Beta]r,\[Gamma]r]];


(* ::Subsection:: *)
(*Approximation of remainder at the Nmax element of the continued fraction*)


RadialCFRemainder[s_Integer,m_Integer,a_Rational|a_Integer,
				  Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:=
Module[{C12tmp,C1,C2,C3,C4,C5,Rem,Err},
	C12tmp=-4I Sqrt[1-a^2]\[Omega];
	C1 = Sqrt[C12tmp];
	If[Re[C1]>0,C1=-C1];
	C2=-(3-4s+8I(\[Omega]+Sqrt[1-a^2]\[Omega]))/4;
	C3=-(3+16Alm+16s+16s^2-32a m \[Omega]-128\[Omega]^2+80a^2\[Omega]^2-128Sqrt[1-a^2]\[Omega]^2-32I(-5Sqrt[1-a^2]\[Omega]+4Sqrt[1-a^2]s \[Omega]))/(32C1);
	C4=(32\[Omega](9+4Alm+2s(-5+4s)+4a \[Omega](-2m+a \[Omega]))+(I(3+16Alm+96a m \[Omega]+16(s^2+(32-51a^2+32Sqrt[1-a^2])\[Omega]^2+s(1+16a \[Omega](-m+2a \[Omega])))))/Sqrt[1-a^2])/(256\[Omega]);
	C5=(1/(8192(-1+a^2)C1 \[Omega]))I(-256Sqrt[1-a^2]Alm^2+Sqrt[1-a^2](63+288s+32s^2-512s^3-256s^4)+64I(21+100s+48s^2-64s^3+I a Sqrt[1-a^2]m(9+48s+112s^2)+a^2(-21-100s-48s^2+64s^3))\[Omega]+32(320I a m(-3+4s)-320I a^3m(-3+4s)-8(-3+225Sqrt[1-a^2]-16(1+19Sqrt[1-a^2])s+16(-1+3Sqrt[1-a^2])s^2)+a^2(3(-8+591Sqrt[1-a^2])+224Sqrt[1-a^2]m^2-16(8+133Sqrt[1-a^2])s+16(-8+59Sqrt[1-a^2])s^2))\[Omega]^2-1024I(-8I a(1+Sqrt[1-a^2])m-I a^3(-8+15Sqrt[1-a^2])m-5a^4(-7+4s)+8(1+Sqrt[1-a^2])(3+4s)-a^2(59+24Sqrt[1-a^2]+4(3+8Sqrt[1-a^2])s))\[Omega]^3+256(-128(1+Sqrt[1-a^2])+a^4(-80+7Sqrt[1-a^2])+16a^2(13+9Sqrt[1-a^2]))\[Omega]^4-32Alm(Sqrt[1-a^2](-9+16s+16s^2)-32I(7+3I a Sqrt[1-a^2]m-4s+a^2(-7+4s))\[Omega]-16(8(1+Sqrt[1-a^2])+a^2(-8+11Sqrt[1-a^2]))\[Omega]^2));
	Rem=1+C1/Sqrt[Nmax]+C2/Nmax+C3/(Nmax^(3/2))+C4/Nmax^2+C5/(Nmax^(5/2));
	Err=Abs[(C5/(Nmax^(5/2)))/Rem];
	{Rem,Err}
]


If[!KerrQNMDebug,Protect[RadialCFRemainder]];


(* ::Subsection:: *)
(*Set SpinWeight and Data-Variable Names*)


SetSpinWeight[s_Integer]:=
Module[{},
	Switch[s,
		   -2,modeName:=Global`KerrQNM; SchTable:=Global`SchQNMTable,
		   -1,modeName:=Global`KerrQNMe; SchTable:=Global`SchQNMeTable,
		    0,modeName:=Global`KerrQNMs; SchTable:=Global`SchQNMsTable,
			_,Print["Invalid QNM Spin Weight : ",s];Abort[]
		  ];
	SetOptions[SchwarzschildMode,ModeSpinWeight->s];
	Print["All KerrMode routines (QNM) set for Spin-Weight s = ",s];
]


If[!KerrQNMDebug,Protect[SetSpinWeight]];


End[] (* KerrModes`Private` *)


(* ::Section:: *)
(*Documentation of External Functions*)


(* ::Subsection:: *)
(*Reserved Globals*)


Protect[];


Begin["`Private`"]


(* ::Section::Closed:: *)
(*End of KerrQNM Package*)


End[] (* `Private` *)


EndPackage[]
