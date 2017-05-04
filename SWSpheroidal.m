(* ::Package:: *)

(* ::Title:: *)
(*Spin-Weighted Spheroidal Functions*)


(* ::Section:: *)
(*Documentation *)


(* ::Section::Closed:: *)
(*Begin SWSpheroidal Package*)


BeginPackage["SWSpheroidal`"]


Unprotect[SWSphDebug];
SWSphDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[SWSphDebug,Unprotect["SWSpheroidal`*"];Unprotect["SWSpheroidal`Private`*"]];
Protect[SWSphDebug];


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


SpinWeightedSpheroidal::usage=
	"SpinWeightedSpheroidal[m,s,c,N] solves the N-dimensional discrete "<>
	"approximation for the spin-weighted functions of spin-weight s, "<>
	"'magnetic' index m, and oblateness parameter c.  The solution returns "<>
	"a list containing {values, vectors}, where values is a list of all "<>
	"the eigenvalues and vectors is a list of all the corresponding spectral "<>
	"coefficients.  These lists are sorted in order of the real part of the "<>
	"eigenvalues."


AngularSpectralRoot::usage=
	"AngularSpectralRoot[s,m,c,Alm,N] solves the N-dimensional discrete "<>
	"approximation for the spin-weighted functions of spin-weight s, "<>
	"'magnetic' index m, and oblateness parameter c.  The solution returns "<>
	"the eigenvalue closest to Alm, and the corresponding spectral coefficients."


AngularSpectralRootIndex::usage=
	"AngularSpectralRoot[s,m,c,index,N] solves the N-dimensional discrete "<>
	"approximation for the spin-weighted functions of spin-weight s, "<>
	"'magnetic' index m, and oblateness parameter c.  The solution returns "<>
	"the eigenvalue and spectral coefficients for the index-th eigenvalue."


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Spin-Weighted Spheroidal Functions : Spectral Method*)


(* ::Text:: *)
(*Evaluate the Spin-weighted Oblate Spheroidal Functions using a spectral method.*)


(* ::Subsection::Closed:: *)
(*Matrix Coefficients :*)


Flms[l_Integer,m_Integer,s_Integer]:=Flms[l,m,s]=
									Piecewise[{{Sqrt[(l+m+1)(l-m+1)/((2l+3)(2l+1))],s==0}},
												Sqrt[(l+m+1)(l-m+1)/((2l+3)(2l+1)) (l+s+1)(l-s+1)/((l+1)(l+1))]];


Glms[l_Integer,m_Integer,s_Integer]:=Glms[l,m,s]=
									Piecewise[{{0,l==0}},
												Sqrt[(l+m)(l-m)/((2l+1)(2l-1)) (l+s)(l-s)/l^2]];


Hlms[l_Integer,m_Integer,s_Integer]:=Hlms[l,m,s]=
									Piecewise[{{0,l==0||s==0}},
												-m s/(l(l+1))];


Alms[l_Integer,m_Integer,s_Integer]:=Alms[l,m,s]=
									Flms[l+1,m,s]Flms[l,m,s];


Dlms[l_Integer,m_Integer,s_Integer]:=Dlms[l,m,s]=
									Hlms[l+1,m,s]Flms[l,m,s]+Hlms[l,m,s]Flms[l,m,s];


Blms[l_Integer,m_Integer,s_Integer]:=Blms[l,m,s]=
									Glms[l+1,m,s]Flms[l,m,s]+Glms[l,m,s]Flms[l-1,m,s] +Hlms[l,m,s]^2;


Elms[l_Integer,m_Integer,s_Integer]:=Elms[l,m,s]=
									Glms[l,m,s]Hlms[l-1,m,s]+Glms[l,m,s]Hlms[l,m,s];


Clms[l_Integer,m_Integer,s_Integer]:=Clms[l,m,s]=
									Glms[l,m,s]Glms[l-1,m,s];


If[!SWSphDebug,Protect[Alms,Blms,Clms,Dlms,Elms,Flms,Glms,Hlms]];


(* ::Subsection::Closed:: *)
(*Construct and solve the spectral matrix for Spin-weighted Oblate Spheroidal functions with azimuthal index m*)


Mat[i_Integer,j_Integer,m_Integer,s_Integer,c_?NumberQ]:=
Module[{lmin,l},
	lmin=Max[Abs[m],Abs[s]];
	l=lmin+(j-1);
	Return[Piecewise[
				{{l(l+1)-s(s+1)- c^2Blms[l,m,s]+2c s Hlms[l,m,s],i==j},
				{-c^2 Dlms[l,m,s]+2c s Flms[l,m,s],i-1==j},
				{-c^2Alms[l,m,s],i-2==j},
				{-c^2 Elms[l,m,s]+2c s Glms[l,m,s],i+1==j},
				{-c^2Clms[l,m,s],i+2==j}},
				0]];
]


SpinWeightedSpheroidal[m_Integer,s_Integer,c_?NumberQ,N_Integer] := 
Module[{unsorted,sortorder,i,j},
	unsorted=Eigensystem[Table[Mat[i,j,m,s,c],{i,1,N},{j,1,N}]];
	sortorder= Ordering[Re[unsorted[[1]]]];
	{unsorted[[1]][[sortorder]],unsorted[[2]][[sortorder]]}
]


AngularSpectralRoot[s_Integer,m_Integer,c_?NumberQ,Alm_?NumberQ,N_Integer]:=
Module[{sol,diff,index},
	sol=SpinWeightedSpheroidal[m,s,c,N];
	diff = Abs[Alm - sol[[1]]];
	index=Position[diff,Min[diff]][[1,1]];
	{sol[[1,index]],N,sol[[2,index]]}
]


AngularSpectralRootIndex[s_Integer,m_Integer,c_?NumberQ,index_Integer,N_Integer]:=
Module[{sol},
	sol=SpinWeightedSpheroidal[m,s,c,N];
	{sol[[1,index]],N,sol[[2,index]]}
]


If[!QNMDebug,Protect[Mat,SpinWeightedSpheroidal,AngularSpectralRoot,AngularSpectralRootForl]];


(* ::Section::Closed:: *)
(*End of SWSpheroidal Package*)


End[] (* `Private` *)


EndPackage[]
