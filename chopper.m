(* ::Package:: *)

BeginPackage["Chopper`"]


Chopper::usage=
	""


FatCamp::usage=
	""


Begin["A`Private`"]


FatCamp[comp_, \[Epsilon]_]:=Module[{num=comp, eps=\[Epsilon]},
	rea = Re[num];
	imag = Im[num];
	If[rea != 0 && imag != 0,
		If[Abs[rea] < Abs[imag] && Abs[(rea / imag)]<(10^(eps)), num=imag*I];
		If[Abs[imag] < Abs[rea] && Abs[(imag / rea)]<(10^(eps)), num=rea];
	];
	Return[num]
]


Options[Chopper]={\[Epsilon]-> -20}


Chopper[packet_,OptionsPattern[]]:=Module[{dat=packet,eps=OptionValue[\[Epsilon]]},
	Do[
		dat=ReplacePart[dat,{i,2,1}-> FatCamp[dat[[i,2,1]],eps]];
		dat=ReplacePart[dat,{i,3,1}-> FatCamp[dat[[i,3,1]],eps]];
		dat=ReplacePart[dat,{i,3,3}-> (FatCamp[#,eps] & /@ dat[[i,3,3]])];,
		{i, 1, Length[dat]}
	];
	Return[dat]
]


End[]


EndPackage[]
