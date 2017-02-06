(* ::Package:: *)

(* ::Title:: *)
(*Mode Plot routines for HDF5 data sets*)


(* ::Section::Closed:: *)
(*Begin ModePlotHDF5 Package*)


BeginPackage["ModePlotHDF5`"]


Unprotect[MPHDF5Debug];
MPHDF5Debug=True; (* Set this to True to allow reloading of Package with changes *)
If[MPHDF5Debug,Unprotect["MPHDF5Debug`*"];Unprotect["MPHDF5Debug`Private`*"]];
Protect[MPHDF5Debug];


(* ::Section:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*Plotting Routines*)


ModePlotHDF5Dir::usage=""


Read\[Omega]KerrQNM::usage=""


ReadAlmKerrQNM::usage=""


ReadSpectralCoefs::usage=""


QNMPlot\[Omega]::usage=""


QNMPlotAlm::usage=""


QNMPlotAll\[Omega]Tones::usage=""


QNMPlotAllAlmTones::usage=""


QNMPlotmTones::usage=""


PlotAccumulation\[Omega]::usage=""


PlotAccumulationAlm::usage=""


SpectralCoef::usage=""


SpectralCoefAll::usage=""


LoopZeros::usage=""


\[Sigma]ptestKerrQNM::usage=""


NonGenericMode::usage=""


(* ::Subsection:: *)
(*Reserved Globals*)


Protect[PlotSchwarzschild,SchwarzschildStyle,SchwarzschildOnly,ModePlotRange,OvertoneList,ModeMaxl,
		OTskip,OTmultiple,MarkerSize,LineThickness,TruncateLongSequences,ZeroTol,MaxFail,
        ModReal];


Begin["`Private`"]


(* ::Section:: *)
(*Plotting Routines*)


(* ::Subsection::Closed:: *)
(*Set HDF5 data directory*)


ModePlotHDF5DIR="";
ModePlotHDF5Dir[dir_]:=ModePlotHDF5dir=dir;


(* ::Subsection::Closed:: *)
(*Low Level reoutines to read data from HDF file*)


ReadKerrQNM[qnm_List]:=Module[{nname,mnamea,mname,qnmname},
If[Head[qnm[[3]]]==List,
If[Length[qnm[[3]]]!=2 &&Not[IntegerQ[qnm[[3]]]],Print["Improper overtone index ",qnm];Abort[]];
nname=If[qnm[[3,1]]<10,"0"<>ToString[qnm[[3,1]]],ToString[qnm[[3,1]]]],
Null[],
nname=If[qnm[[3]]<10,"0"<>ToString[qnm[[3]]],ToString[qnm[[3]]]]
];
If[Not[IntegerQ[qnm[[2]]]],Print["Improper m index"];Abort[]];
mnamea=If[Abs[qnm[[2]]]<10,"0"<>ToString[Abs[qnm[[2]]]],ToString[Abs[qnm[[2]]]]];
mname=If[qnm[[2]]<0,"-"<>mnamea,"+"<>mnamea];
If[Not[IntegerQ[qnm[[1]]]],Print["Improper l index"];Abort[]];
qnmname = StringReplace[ToString[qnm]," "->""];
Import[ModePlotHDF5dir<>"KerrQNM_"<>nname<>".h5",{"HDF5","Datasets",{"/n"<>nname<>"/m"<>mname<>"/"<>qnmname}}]
]
ReadKerrQNM[l_Integer,m_Integer,n_Integer|n_List]:=ReadKerrQNM[{l,m,n}]


(* ::Subsection::Closed:: *)
(*Routines to read \[Omega], Alm, Spectral Coefficients*)


Options[Read\[Omega]KerrQNM]={TruncateLongSequences->True};
Read\[Omega]KerrQNM[qnm_List,range_List:{1,-1},opts:OptionsPattern[]]:=
Module[{rawdat,Nelems,full,amod,short,a,lr=range,trunc=OptionValue[TruncateLongSequences]},
Off[Import::dataset];
rawdat=ReadKerrQNM[qnm];
On[Import::dataset];
If[rawdat==$Failed,Return[$Failed]];
If[lr=={1,-1} && trunc,
If[qnm=={2,0,{11,1}},lr={1,20000}];If[qnm=={2,0,{12,1}},lr={1,22000}];If[qnm=={2,0,{13,1}},lr={1,35000}];If[qnm=={2,0,{14,1}},lr={1,35000}];
If[qnm=={2,0,{15,1}},lr={1,35000}];If[qnm=={2,0,{16,1}},lr={1,25000}];If[qnm=={2,0,{17,1}},lr={1,25000}];If[qnm=={2,0,{18,1}},lr={1,25000}];
If[qnm=={2,0,{19,1}},lr={1,25000}];If[qnm=={2,0,{20,1}},lr={1,25000}];If[qnm=={2,0,{21,1}},lr={1,25000}];If[qnm=={2,0,{22,1}},lr={1,25000}];
If[qnm=={2,0,{23,1}},lr={1,25000}];If[qnm=={2,0,{24,1}},lr={1,25000}];If[qnm=={2,0,{25,1}},lr={1,25000}];If[qnm=={2,0,{26,1}},lr={1,25000}];
If[qnm=={3,0,19},lr={1,20000}];If[qnm=={3,0,20},lr={1,20000}];If[qnm=={3,0,21},lr={1,22000}];If[qnm=={3,0,22},lr={1,35000}];
If[qnm=={3,0,23},lr={1,35000}];If[qnm=={3,0,24},lr={1,25000}];If[qnm=={3,0,25},lr={1,25000}];If[qnm=={3,0,26},lr={1,25000}];
If[qnm=={3,0,27},lr={1,25000}];If[qnm=={3,0,28},lr={1,25000}];If[qnm=={3,0,29},lr={1,25000}];If[qnm=={3,0,30},lr={1,25000}];
If[qnm=={3,0,31},lr={1,25000}];
If[qnm=={4,0,26},lr={1,22000}];If[qnm=={4,0,27},lr={1,25000}];If[qnm=={4,0,28},lr={1,35000}];If[qnm=={4,0,29},lr={1,35000}];
If[qnm=={4,0,30},lr={1,25000}];If[qnm=={4,0,31},lr={1,25000}];
];
full=Function[x,{x[[1]],-x[[2]]}]/@Take[rawdat,lr,{2,3}];
a=Flatten[Take[rawdat,lr,1]];
amod = Position[Mod[#,1/20]&/@Function[x,x/10^16]/@ Function[x,IntegerPart[10^16 x]]/@a,0];
If[1-a[[-1]]<10^(-6),AppendTo[amod,Length[full]]];
short=full[[Flatten[amod]]];
{full,short}
]
Read\[Omega]KerrQNM[l_Integer,m_Integer,n_Integer|n_List,range_List:{1,-1},opts:OptionsPattern[]]:=Read\[Omega]KerrQNM[{l,m,n},range,opts]


Options[ReadAlmKerrQNM]={TruncateLongSequences->True};
ReadAlmKerrQNM[qnm_List,range_List:{1,-1},opts:OptionsPattern[]]:=
Module[{rawdat,Nelems,full,amod,short,a,lr=range,trunc=OptionValue[TruncateLongSequences]},
Off[Import::dataset];
rawdat=ReadKerrQNM[qnm];
On[Import::dataset];
If[rawdat==$Failed,Return[$Failed]];
If[lr=={1,-1} && trunc,
If[qnm=={2,0,{11,1}},lr={1,20000}];If[qnm=={2,0,{12,1}},lr={1,22000}];If[qnm=={2,0,{13,1}},lr={1,35000}];If[qnm=={2,0,{14,1}},lr={1,35000}];
If[qnm=={2,0,{15,1}},lr={1,35000}];If[qnm=={2,0,{16,1}},lr={1,25000}];If[qnm=={2,0,{17,1}},lr={1,25000}];If[qnm=={2,0,{18,1}},lr={1,25000}];
If[qnm=={2,0,{19,1}},lr={1,25000}];If[qnm=={2,0,{20,1}},lr={1,25000}];If[qnm=={2,0,{21,1}},lr={1,25000}];If[qnm=={2,0,{22,1}},lr={1,25000}];
If[qnm=={2,0,{23,1}},lr={1,25000}];If[qnm=={2,0,{24,1}},lr={1,25000}];If[qnm=={2,0,{25,1}},lr={1,25000}];If[qnm=={2,0,{26,1}},lr={1,25000}];
If[qnm=={3,0,19},lr={1,20000}];If[qnm=={3,0,20},lr={1,20000}];If[qnm=={3,0,21},lr={1,22000}];If[qnm=={3,0,22},lr={1,35000}];
If[qnm=={3,0,23},lr={1,35000}];If[qnm=={3,0,24},lr={1,25000}];If[qnm=={3,0,25},lr={1,25000}];If[qnm=={3,0,26},lr={1,25000}];
If[qnm=={3,0,27},lr={1,25000}];If[qnm=={3,0,28},lr={1,25000}];If[qnm=={3,0,29},lr={1,25000}];If[qnm=={3,0,30},lr={1,25000}];
If[qnm=={3,0,31},lr={1,25000}];
If[qnm=={4,0,26},lr={1,22000}];If[qnm=={4,0,27},lr={1,25000}];If[qnm=={4,0,28},lr={1,35000}];If[qnm=={4,0,29},lr={1,35000}];
If[qnm=={4,0,30},lr={1,25000}];If[qnm=={4,0,31},lr={1,25000}];
];
full=Take[rawdat,lr,{4,5}];
a=Flatten[Take[rawdat,lr,1]];
amod = Position[Mod[#,1/20]&/@Function[x,x/10^16]/@ Function[x,IntegerPart[10^16 x]]/@a,0];
If[1-a[[-1]]<10^(-6),AppendTo[amod,Length[full]]];
short=full[[Flatten[amod]]];
{full,short}
]
ReadAlmKerrQNM[l_Integer,m_Integer,n_Integer|n_List,range_List:{1,-1},opts:OptionsPattern[]]:=ReadAlmKerrQNM[{l,m,n},range,opts]


ReadSpectralCoefs[qnm_List,range_List:{1,-1}]:=Module[{rawdat},
Off[Import::dataset];
rawdat=ReadKerrQNM[qnm];
On[Import::dataset];
If[rawdat==$Failed,Return[$Failed]];
Drop[Take[rawdat,range],0,{2,5}]
]
ReadSpectralCoefs[l_Integer,m_Integer,n_Integer|n_List,range_List:{1,-1}]:=ReadSpectralCoefs[{l,m,n},range]


(* ::Subsection:: *)
(*Routines to Plot \[Omega] and Alm for specific lists of modes {l, m, n}*)


Options[QNMPlot\[Omega]]=Union[{ModePlotRange->{1,-1},PlotSchwarzschild->True,SchwarzschildStyle->{{Gray,Dashed}},SchwarzschildOnly->False},
						Options[Read\[Omega]KerrQNM],Options[ListLinePlot],Options[ListPlot]];
SetOptions[QNMPlot\[Omega],AxesLabel->Automatic];
QNMPlot\[Omega][qnmlist_List,opts:OptionsPattern[]]:=
Module[{alllists,onlyovertones,plotlist,
		range=OptionValue[ModePlotRange],alabel=OptionValue[AxesLabel],
		pos,qnml=qnmlist,schplot,plotsch=OptionValue[PlotSchwarzschild],onlysch=OptionValue[SchwarzschildOnly]},
SetOptions[ListLinePlot,ImageSize->800,TicksStyle->Directive[14],PlotRange->All];
SetOptions[ListPlot,ImageSize->800,TicksStyle->Directive[14],PlotRange->All,PlotMarkers->{Automatic,Small}];
alllists=Read\[Omega]KerrQNM[#,range,FilterRules[{opts},Options[Read\[Omega]KerrQNM]]]& /@ qnmlist;
(* onlyovertones makes sure all qnm in alllists have same value of l and m *)
onlyovertones=Length[DeleteCases[qnmlist,{l_,m_,n_}/;l==qnmlist[[1,1]]&&m==qnmlist[[1,2]]]]==0;
If[MemberQ[alllists,$Failed],
Print["No data for QNMs :",qnmlist[[Flatten[Position[alllists,$Failed],1]]]];Return[$Failed];
];
plotlist= {ListLinePlot[Part[#,1] & /@ alllists,FilterRules[{opts},FilterRules[Options[ListLinePlot],Except[{PlotMarkers,Joined}]]]]};
AppendTo[plotlist,
ListPlot[Part[#,2] & /@ alllists,FilterRules[FilterRules[{opts},Options[ListPlot]],Except[{Joined,PlotLegends}]]]];
If[onlysch,
schplot=Part[#,2,1] & /@ alllists,
If[plotsch&&onlyovertones && range[[1]]==1,
axismodes=Transpose[{{{2,_,8},{0,2}},{{3,_,40},{0,10}}}];
schplot=Part[#,2,1] & /@ alllists;
pos=Position[qnml,{_,_,{_,1}}];
For[i=Length[pos],i>0,--i,qnml=Drop[qnml,pos[[i]]];schplot=Drop[schplot,pos[[i]]]];
qnml = qnml/.{l_,m_,{n_,0}}->{l,m,n};
For[i=1,i<=Length[axismodes[[1]]],++i,
pos=Position[qnml,axismodes[[1,i]]];
If[Length[pos]>=1,
schplot=ReplacePart[schplot,#->axismodes[[2,i]]]&/@Position[qnml,axismodes[[1,i]]]
];
];
AppendTo[plotlist,ListLinePlot[schplot,PlotStyle->OptionValue[SchwarzschildStyle]]];
];
];
Clear[alllists];
If[alabel==Automatic,
alabel = {Style["Re(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))",16],Style["-Im(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))",16]};
];
If[onlysch,
ListPlot[schplot,FilterRules[FilterRules[{opts},Options[ListPlot]],Except[{Joined,PlotLegends}]]],
Show[plotlist,AxesLabel->alabel,PlotRange->OptionValue[PlotRange]]
]
]


Options[QNMPlotAlm]=Union[{ModePlotRange->{1,-1}},Options[Read\[Omega]KerrQNM],Options[ListLinePlot],Options[ListPlot]];
SetOptions[QNMPlotAlm,AxesLabel->Automatic];
QNMPlotAlm[qnmlist_List,opts:OptionsPattern[]]:=
Module[{alllists,plotlist,
		range=OptionValue[ModePlotRange],alabel=OptionValue[AxesLabel]},
SetOptions[ListLinePlot,ImageSize->800,TicksStyle->Directive[14],PlotRange->All];
SetOptions[ListPlot,ImageSize->800,TicksStyle->Directive[14],PlotRange->All,PlotMarkers->{Automatic,Small}];
alllists=ReadAlmKerrQNM[#,range,FilterRules[{opts},Options[ReadAlmKerrQNM]]]& /@ qnmlist;
If[MemberQ[alllists,$Failed],
Print["No data for QNMs :",qnmlist[[Flatten[Position[alllists,$Failed],1]]]];Return[$Failed];
];
plotlist= {ListLinePlot[Part[#,1] & /@ alllists,FilterRules[{opts},FilterRules[Options[ListLinePlot],Except[{PlotMarkers,Joined}]]]]};
AppendTo[plotlist,
ListPlot[Part[#,2] & /@ alllists,FilterRules[FilterRules[{opts},Options[ListPlot]],Except[{Joined,PlotLegends}]]]];
Clear[alllists];
If[alabel==Automatic,
alabel={Style[Re["\!\(\*SubscriptBox[\(A\), \(lm\)]\)"],16],Style[Im["\!\(\*SubscriptBox[\(A\), \(lm\)]\)"],16]};
];
Show[plotlist,AxesLabel->alabel,PlotRange->OptionValue[PlotRange]]
]


(* ::Subsection::Closed:: *)
(*Routines to Plot \[Omega] and Alm for sets of modes for a specific "m"*)


(* ::Subsubsection::Closed:: *)
(*Lists of known Overtone Multiplets*)


QNMOTlist=Join[Table[{2,-2,n},{n,13,16}],Table[{3,-2,n},{n,26,29}],{{2,1,8},{2,2,8}},Table[{2,0,n},{n,8,26}],Table[{2,1,n},{n,34,39}]];
QNMOTSlist={{2,-2,15},{3,-2,28}}; (* Special cases of "single" multiplets *)


(* ::Subsubsection::Closed:: *)
(*Ploting routines*)


Options[QNMPlotAll\[Omega]Tones]=Union[{OvertoneList->Range[0,15],MarkerSize->10,LineThickness->Automatic},Options[QNMPlot\[Omega]]];
SetOptions[QNMPlotAll\[Omega]Tones,PlotMarkers->Automatic];
QNMPlotAll\[Omega]Tones[m_Integer,llist_List,opts:OptionsPattern[]]:=Module[{lmin,lvals=llist,nlist=OptionValue[OvertoneList],qnms,i,j,qnmot,label,plots={},colors,marks,pstyle=OptionValue[PlotStyle],pmarks=OptionValue[PlotMarkers],msize=OptionValue[MarkerSize]},
lmin = Max[2,Abs[m]];
lvals=DeleteCases[lvals,l_/;l<lmin];
label=OptionValue[PlotLabel];
If[label==Automatic,
label = Style["m="<>ToString[m]<>" : (l="<>ToString[lvals]<>", n="<>ToString[nlist]<>")",20];
];
For[i=1,i<=Length[lvals],++i,
qnms=Table[{lvals[[i]],m,nlist[[j]]},{j,1,Length[nlist]}];
colors=Range[1,Length[qnms]];
marks=Range[1,Length[qnms]];
For[j=Length[qnms],j>0,--j,
qnmot=qnms[[j]];
If[MemberQ[QNMOTlist,qnmot],
If[MemberQ[QNMOTSlist,qnmot],
qnms=ReplacePart[qnms,j->{qnmot[[1]],qnmot[[2]],{qnmot[[3]],0}}],
qnms=ReplacePart[qnms,j->{qnmot[[1]],qnmot[[2]],{qnmot[[3]],1}}];
qnms=Insert[qnms,{qnmot[[1]],qnmot[[2]],{qnmot[[3]],0}},j];
colors=Insert[colors,colors[[j]],j];
marks=Insert[marks,marks[[j]],j];
];
];
];
If[OptionValue[LineThickness]==Automatic,
If[OptionValue[PlotStyle]==Automatic,pstyle=ColorData[1,#]& /@colors],Null,
If[OptionValue[PlotStyle]==Automatic,pstyle={ColorData[1,#],Thickness[OptionValue[LineThickness]]}& /@colors]
];
If[OptionValue[PlotMarkers]==Automatic,pmarks={{Global`\[FilledCircle],Global`\[FilledSmallSquare],Global`\[FilledDiamond],Global`\[FilledUpTriangle],Global`\[FilledDownTriangle],
					Global`\[EmptyCircle],Global`\[EmptySquare],Global`\[EmptyDiamond],Global`\[EmptyUpTriangle],Global`\[EmptyDownTriangle]}[[Mod[#-1,10]+1]],msize}& /@marks];
AppendTo[plots,QNMPlot\[Omega][qnms,PlotStyle->pstyle,PlotMarkers->pmarks,FilterRules[FilterRules[{opts},Options[QNMPlot\[Omega]]],Except[{PlotStyle,PlotMarkers}]]]];
];
Show[plots,PlotLabel->label,PlotRange->OptionValue[PlotRange]]
]


Options[QNMPlotAllAlmTones]=Union[{OvertoneList->Range[0,15],MarkerSize->10,LineThickness->Automatic},Options[QNMPlotAlm]];
SetOptions[QNMPlotAllAlmTones,PlotMarkers->Automatic];
QNMPlotAllAlmTones[m_Integer,llist_List,opts:OptionsPattern[]]:=Module[{lmin,lvals=llist,nlist=OptionValue[OvertoneList],qnms,i,j,qnmot,label,plots={},colors,marks,pstyle=OptionValue[PlotStyle],pmarks=OptionValue[PlotMarkers],msize=OptionValue[MarkerSize]},
lmin = Max[2,Abs[m]];
lvals=DeleteCases[lvals,l_/;l<lmin];
label=OptionValue[PlotLabel];
If[label==Automatic,
label = Style["m="<>ToString[m]<>" : (l="<>ToString[lvals]<>", n="<>ToString[nlist]<>")",20];
];
For[i=1,i<=Length[lvals],++i,
qnms=Table[{lvals[[i]],m,nlist[[j]]},{j,1,Length[nlist]}];
colors=Range[1,Length[qnms]];
marks=Range[1,Length[qnms]];
For[j=Length[qnms],j>0,--j,
qnmot=qnms[[j]];
If[MemberQ[QNMOTlist,qnmot],
If[MemberQ[QNMOTSlist,qnmot],
qnms=ReplacePart[qnms,j->{qnmot[[1]],qnmot[[2]],{qnmot[[3]],0}}],
qnms=ReplacePart[qnms,j->{qnmot[[1]],qnmot[[2]],{qnmot[[3]],1}}];
qnms=Insert[qnms,{qnmot[[1]],qnmot[[2]],{qnmot[[3]],0}},j];
colors=Insert[colors,colors[[j]],j];
marks=Insert[marks,marks[[j]],j];
];
];
];
If[OptionValue[LineThickness]==Automatic,
If[OptionValue[PlotStyle]==Automatic,pstyle=ColorData[1,#]& /@colors],Null,
If[OptionValue[PlotStyle]==Automatic,pstyle={ColorData[1,#],Thickness[OptionValue[LineThickness]]}& /@colors]
];
If[OptionValue[PlotMarkers]==Automatic,pmarks={{Global`\[FilledCircle],Global`\[FilledSmallSquare],Global`\[FilledDiamond],Global`\[FilledUpTriangle],Global`\[FilledDownTriangle],
					Global`\[EmptyCircle],Global`\[EmptySquare],Global`\[EmptyDiamond],Global`\[EmptyUpTriangle],Global`\[EmptyDownTriangle]}[[Mod[#-1,10]+1]],msize}& /@marks];
AppendTo[plots,QNMPlotAlm[qnms,PlotStyle->pstyle,PlotMarkers->pmarks,FilterRules[FilterRules[{opts},Options[QNMPlot\[Omega]]],Except[{PlotStyle,PlotMarkers}]]]];
];
Show[plots,PlotLabel->label,PlotRange->OptionValue[PlotRange]]
]


(* ::Subsubsection::Closed:: *)
(*Routine to plot all modes for specified "m"*)


Options[QNMPlotmTones]=Union[{ModeMaxl->16},Options[QNMPlotAll\[Omega]Tones]];
QNMPlotmTones[m_,opts:OptionsPattern[]]:=Module[{lmax=OptionValue[ModeMaxl]},
Flatten[Table[{QNMPlotAll\[Omega]Tones[m,{l},FilterRules[{opts},Options[QNMPlotAll\[Omega]Tones]]],QNMPlotAllAlmTones[m,{l},FilterRules[{opts},Options[QNMPlotAll\[Omega]Tones]]]},{l,Max[2,Abs[m]],lmax}]]
]


(* ::Subsection::Closed:: *)
(*Routines for Accumulation Plots*)


Options[OTLists]={OTskip->{},OTmultiple->{}};
OTLists[n_Integer|n_List,overtones_List,OptionsPattern[]]:=Module[{skip=DeleteDuplicates[OptionValue[OTskip]],
multiple=DeleteDuplicates[OptionValue[OTmultiple]],otsort,fitot=n,fitotsave=n,
fitind,nf=DeleteDuplicates[overtones],nfus,intnf={},i,multints,posnf,possk,pos,
totlist,nprime,nrm,intlist},
otsort[a_,b_]:=If[Head[a]==Head[b],OrderedQ[{a,b}],
					Null[],
					If[Head[a]==List,a[[1]]<=b,Null[],a<=b[[1]]]];
nfus=nf;
multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
skip=Sort[skip,otsort];
multiple=Sort[multiple,otsort];
nf=Sort[nf,otsort];
(* Abort if any multiple overtone has length one in OTmultiple *)
For[i=1,i<=Length[multiple],++i,If[multiple[[i,2]]<2,Print["Otmultiple ",multiple[[i]]," is not a multiple."];Abort[]]
];
(* Abort if any multiple overtone of List type are not in OTmultiple *)
posnf=Flatten[Drop[Position[Head /@ nf,List],1]];
possk=Flatten[Drop[Position[Head /@ skip,List],1]];
For[i=1,i<=Length[posnf],++i,
pos=Flatten[Position[multints,nf[[posnf[[i]],1]]]];
If[Length[pos]==0,Print["overtone ",nf[[posnf[[i]]]]," is not in OTmultiple"];Abort[]];
];
For[i=1,i<=Length[possk],++i,
pos=Flatten[Position[multints,skip[[possk[[i]],1]]]];
If[Length[pos]==0,Print["OTskip ",skip[[possk[[i]]]]," is not in OTmultiple"];Abort[]];
];
If[Head[fitot]==List,
pos=Flatten[Position[multints,fitot[[1]]]];
If[Length[pos]==0,Print["Fit overtone ",fitot," is not in OTmultiple"];Abort[]];
];
(* Convert any multiple overtone of Integer type to List type. eg 8\[Rule]{8,0} *)
(* Save converted Integer elements of nf so they can be returned to original form *)
posnf=Flatten[Position[Head /@ nf,Integer]];
possk=Flatten[Position[Head /@ skip,Integer]];
For[i=1,i<=Length[multints],++i,
pos=posnf[[Flatten[Position[nf[[posnf]],multints[[i]]]]]];
If[Length[pos]==1,nf[[posnf[[pos[[1]]]]]]={multints[[i]],0};AppendTo[intnf,multints[[i]]]];
pos=possk[[Flatten[Position[skip[[possk]],multints[[i]]]]]];
If[Length[pos]==1,skip[[possk[[pos[[1]]]]]]={multints[[i]],0}];
];
If[Head[fitot]==Integer && fitot>=0,
pos=Flatten[Position[multints,fitot]];
If[Length[pos]==1,fitot={fitot,0}];
];
(* Get list of all overtones used in fit, removing specified overtones *)
nrm=Flatten[Position[nf,#] & /@ skip,1];
nf=Delete[nf,nrm];
nrm=Flatten[Position[nfus,#] & /@ skip,1];
nfus=Delete[nfus,nrm];
(* Get list of all remaining non-multiple overtones in list *)
nrm=Drop[Position[Head /@ nf,List],1];
intlist=Delete[nf,nrm];
(* Get list of all possible integer overtones *)
totlist=Table[i,{i,Min[intlist,multints],Max[intlist,multints]}];
nrm=Flatten[Position[totlist,#] & /@ multints,1];
totlist=Delete[totlist,nrm];
(* Add list of all possible multiple overtones *)
For[i=1,i<=Length[multiple],++i,
totlist=Join[totlist,Table[{multiple[[i,1]],j},{j,0,multiple[[i,2]]-1}]];
];
(* remove specified overtones *)
nrm=Flatten[Position[totlist,#] & /@ skip,1];
totlist=Sort[Delete[totlist,nrm],otsort];
nprime=Table[Position[totlist,nf[[i]]][[1,1]]-1,{i,1,Length[nf]}];
(* return integer overtones in nf *)
For[i=1,i<=Length[intnf],++i,
pos=Flatten[Position[nf,{intnf[[i]],0}]];
If[Length[pos]>0,nf[[pos[[1]]]]=intnf[[i]]];
];
pos=Flatten[Position[totlist,fitot]];
If[Length[pos]>0,fitind=pos[[1]]-1,fitind=-1];
{fitotsave,fitind,nfus,nprime}
]


Options[PlotAccumulation\[Omega]]=Union[{OTskip->{},OTmultiple->{}},Options[Plot],Options[ListPlot]];
SetOptions[PlotAccumulation\[Omega],AxesLabel->Automatic,AxesOrigin->Automatic,PlotLabel->Automatic,PlotRange->Automatic];
PlotAccumulation\[Omega][l_Integer,m_Integer,n_List,Nv_Integer,
					a0_Real|a0_Rational|a0_Integer,
					opts:OptionsPattern[]]:=
Module[{s=-2,KerrSEQ,nf,nprime,fitlists,i,j,qnmend,qnmRe,qnmIm,ReList={},ImList={},Refit,ReLOF,ReParams,ReFitData,Imfit,ImLOF,ImParams,ImFitData,rplot,iplot,label,prange=OptionValue[PlotRange],aorig=OptionValue[AxesOrigin],plabel=OptionValue[PlotLabel],realabel=OptionValue[AxesLabel],imalabel=OptionValue[AxesLabel]},
SetOptions[ListPlot,ImageSize->800,TicksStyle->Directive[14],PlotRange->All,PlotMarkers->Automatic];
fitlists=OTLists[-1,n,FilterRules[{opts},Options[OTLists]]];
nf=fitlists[[3]];
nprime=fitlists[[4]];
For[i=1,i<=Length[nf],++i,
qnmend=Function[x,{1-x[[1]],x[[2]]-m/2,x[[3]]}]/@Take[ReadKerrQNM[l,m,nf[[i]]],{-Nv,-1},{1,3}];
qnmRe=Function[x,{x[[1]],x[[2]]}]/@qnmend;
qnmIm=Function[x,{x[[1]],x[[3]]}]/@qnmend;
qnmRe=DeleteCases[qnmRe,{eps_,_}/;eps>a0];
qnmIm=DeleteCases[qnmIm,{eps_,_}/;eps>a0];
AppendTo[ReList,qnmRe];AppendTo[ImList,qnmIm];
];
ReFitData=Flatten[Table[{nprime[[i]],ReList[[i,j,1]],ReList[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ReList[[i]]]}],1];
ImFitData=Flatten[Table[{nprime[[i]],ImList[[i,j,1]],ImList[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ImList[[i]]]}],1];
If[nf[[1]]==0,
Refit=NonlinearModelFit[ReFitData,-\[Alpha]1 Sqrt[\[Epsilon]1/2]+(\[Alpha]2+\[Alpha]3 n1)\[Epsilon]1,{\[Alpha]1,\[Alpha]2,\[Alpha]3},{n1,\[Epsilon]1}];
Imfit=NonlinearModelFit[ImFitData,-(n1+1/2)(Sqrt[\[Epsilon]1/2]-\[Alpha]4 \[Epsilon]1),{\[Alpha]4},{n1,\[Epsilon]1}],
Refit=NonlinearModelFit[ReFitData,+(\[Alpha]1+\[Alpha]2 n1)\[Epsilon]1,{\[Alpha]1,\[Alpha]2},{n1,\[Epsilon]1}];
Imfit=NonlinearModelFit[ImFitData,-1(\[Alpha]3+n1+1/2)Sqrt[\[Epsilon]1/2]+(\[Alpha]4+\[Alpha]5 n1)\[Epsilon]1,{\[Alpha]3,\[Alpha]4,\[Alpha]5},{n1,\[Epsilon]1}]
];

Print["Re(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))= ",Normal[Refit]];
Print["Im(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))= ",Normal[Imfit]];
Off[SetPrecision::"precsm"];Off[N::"precsm"];
Print[MatrixForm[{{"Re(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))",Refit["ParameterConfidenceIntervalTable"]},{"Im(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))",Imfit["ParameterConfidenceIntervalTable"]}}]];
On[SetPrecision::"precsm"];On[N::"precsm"];

(* Modify nf for multiplet labels *)
For[i=1,i<=Length[nf],++i,
If[Length[nf[[i]]]==2,nf[[i]]=Subscript[nf[[i,1]], nf[[i,2]]]];
];
ReParams=Refit["BestFitParameters"];
ImParams=Imfit["BestFitParameters"];
If[nf[[1]]==0,
ReLOF[nf_,\[Epsilon]f_]:=-\[Alpha]1 Sqrt[\[Epsilon]f/2]/.ReParams;
ImLOF[nf_,\[Epsilon]f_]:=-(nf+1/2)Sqrt[\[Epsilon]f/2]/.ImParams,
ReLOF[nf_,\[Epsilon]f_]:=+(\[Alpha]1+\[Alpha]2 nf)\[Epsilon]f/.ReParams;
ImLOF[nf_,\[Epsilon]f_]:=-1(\[Alpha]3+nf+1/2)Sqrt[\[Epsilon]f/2]/.ImParams
];
If[prange==Automatic,prange={{0,a0},All}];
If[aorig==Automatic,aorig={0,0}];
If[plabel==Automatic,plabel = Style[DisplayForm[RowBox[{"l=",l," m=",m," n=",nf}]],16];];
If[realabel==Automatic,realabel = {Style["1-\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],Style["Re(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))-m/2",16]};];
rplot=Show[ListPlot[ReList,PlotRange->prange,AxesOrigin->aorig,AxesLabel->realabel,PlotLabel->plabel,FilterRules[FilterRules[{opts},Options[ListPlot]],Except[{Joined,PlotRange,AxesOrigin,AxesLabel,PlotLabel}]]],Plot[ReLOF[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0},PlotStyle->{{Red,Dashed}}],Plot[Refit[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0}]];
If[imalabel==Automatic,imalabel = {Style["1-\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],Style["Im(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))",16]};];
iplot=Show[ListPlot[ImList,PlotRange->prange,AxesOrigin->aorig,AxesLabel->imalabel,PlotLabel->plabel,FilterRules[FilterRules[{opts},Options[ListPlot]],Except[{Joined,PlotRange,AxesOrigin,AxesLabel,PlotLabel}]]],Plot[ImLOF[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0},PlotStyle->{{Red,Dashed}}],Plot[Imfit[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0}]];
{rplot,iplot}
]


Options[PlotAccumulationAlm]=Union[{OTskip->{},OTmultiple->{}},Options[Plot],Options[ListPlot]];
SetOptions[PlotAccumulationAlm,AxesLabel->Automatic,AxesOrigin->Automatic,PlotLabel->Automatic,PlotRange->Automatic];
PlotAccumulationAlm[l_Integer,m_Integer,n_List,Nv_Integer,
					a0_Real|a0_Rational|a0_Integer,
					opts:OptionsPattern[]]:=
Module[{s=-2,KerrSEQ,nf,nprime,fitlists,i,j,qnmend,qnmRe,qnmIm,ReList={},ImList={},Refit,ReLOF,ReParams,ReFitData,Imfit,ImLOF,ImParams,ImFitData,rplot,iplot,label,A0,intercept,\[Delta],prange=OptionValue[PlotRange],raorig=OptionValue[AxesOrigin],iaorig=OptionValue[AxesOrigin],plabel=OptionValue[PlotLabel],realabel=OptionValue[AxesLabel],imalabel=OptionValue[AxesLabel]},
SetOptions[ListPlot,ImageSize->800,TicksStyle->Directive[14],PlotRange->All,PlotMarkers->Automatic];
fitlists=OTLists[-1,n,FilterRules[{opts},Options[OTLists]]];
nf=fitlists[[3]];
nprime=fitlists[[4]];
For[i=1,i<=Length[nf],++i,
qnmend=Function[x,{1-x[[1]],x[[4]],x[[5]]}]/@Take[ReadKerrQNM[l,m,nf[[i]]],{-Nv,-1},{1,5}];
qnmRe=Function[x,{x[[1]],x[[2]]}]/@qnmend;
qnmIm=Function[x,{x[[1]],x[[3]]}]/@qnmend;
qnmRe=DeleteCases[qnmRe,{eps_,_}/;eps>a0];
qnmIm=DeleteCases[qnmIm,{eps_,_}/;eps>a0];
AppendTo[ReList,qnmRe];AppendTo[ImList,qnmIm];
];
ReFitData=Flatten[Table[{nprime[[i]],ReList[[i,j,1]],ReList[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ReList[[i]]]}],1];
ImFitData=Flatten[Table[{nprime[[i]],ImList[[i,j,1]],ImList[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ImList[[i]]]}],1];
If[nf[[1]]==0,
Refit=NonlinearModelFit[ReFitData,l(l+1)-2+\[Beta]1+\[Beta]2 Sqrt[\[Epsilon]1/2]+(\[Beta]3+\[Beta]4 n1+\[Beta]5 n1^2)\[Epsilon]1,{\[Beta]1,\[Beta]2,\[Beta]3,\[Beta]4,\[Beta]5},{n1,\[Epsilon]1}];
Imfit=NonlinearModelFit[ImFitData,(n1+1/2)(\[Beta]6 Sqrt[\[Epsilon]1/2]+\[Beta]7 \[Epsilon]1),{\[Beta]6,\[Beta]7},{n1,\[Epsilon]1}],
Refit=NonlinearModelFit[ReFitData,l(l+1)-2+\[Beta]1+(\[Beta]2+\[Beta]3 n1+\[Beta]4 n1^2)\[Epsilon]1,{\[Beta]1,\[Beta]2,\[Beta]3,\[Beta]4},{n1,\[Epsilon]1}];
Imfit=NonlinearModelFit[ImFitData,(\[Beta]5+n1+1/2)\[Beta]6 Sqrt[\[Epsilon]1/2]+(\[Beta]7+\[Beta]8 n1)\[Epsilon]1,{\[Beta]5,\[Beta]6,\[Beta]7,\[Beta]8},{n1,\[Epsilon]1}]
];
A0=Refit["BestFitParameters"];
intercept=l(l+1)-2+\[Beta]1/.A0;
\[Delta]=Sqrt[7/4 m^2 - (-2+1/2)^2-(l(l+1)-2+\[Beta]1)/.A0];
Print["\[Delta] = ",\[Delta]];
Print["Re(\!\(\*SubscriptBox[\(A\), \(lm\)]\))= ",Normal[Refit]];
Print["Im(\!\(\*SubscriptBox[\(A\), \(lm\)]\))= ",Normal[Imfit]];
Off[SetPrecision::"precsm"];Off[N::"precsm"];
Print[MatrixForm[{{"Re(\!\(\*SubscriptBox[\(A\), \(lm\)]\))",Refit["ParameterConfidenceIntervalTable"]},{"Im(\!\(\*SubscriptBox[\(A\), \(lm\)]\))=",Imfit["ParameterConfidenceIntervalTable"]}}]];
On[SetPrecision::"precsm"];On[N::"precsm"];

(* Modify nf for multiplet labels *)
For[i=1,i<=Length[nf],++i,
If[Length[nf[[i]]]==2,nf[[i]]=Subscript[nf[[i,1]], nf[[i,2]]]];
];
ReParams=Refit["BestFitParameters"];
ImParams=Imfit["BestFitParameters"];
If[nf[[1]]==0,
ReLOF[nf_,\[Epsilon]f_]:=l(l+1)-2+\[Beta]1+\[Beta]2 Sqrt[\[Epsilon]f/2]/.ReParams;
ImLOF[nf_,\[Epsilon]f_]:=(nf+1/2)\[Beta]6 Sqrt[\[Epsilon]f/2]/.ImParams,
ReLOF[nf_,\[Epsilon]f_]:=l(l+1)-2+\[Beta]1+(\[Beta]2+\[Beta]3 nf+\[Beta]4 nf^2)\[Epsilon]f/.ReParams;
ImLOF[nf_,\[Epsilon]f_]:=(\[Beta]5+nf+1/2)\[Beta]6 Sqrt[\[Epsilon]f/2]/.ImParams
];
If[prange==Automatic,prange={{0,a0},All}];
If[raorig==Automatic,raorig={0,intercept}];
If[plabel==Automatic,plabel = Style[DisplayForm[RowBox[{"l=",l," m=",m," n=",nf}]],16];];
If[realabel==Automatic,realabel = {Style["1-\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],Style["Re(\!\(\*SubscriptBox[\(A\), \(lm\)]\))",16]};];
rplot=Show[ListPlot[ReList,PlotRange->prange,AxesOrigin->raorig,AxesLabel->realabel,PlotLabel->plabel,FilterRules[FilterRules[{opts},Options[ListPlot]],Except[{Joined,PlotRange,AxesOrigin,AxesLabel,PlotLabel}]]],Plot[ReLOF[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0},PlotStyle->{{Red,Dashed}}],Plot[Refit[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0}]];
If[iaorig==Automatic,iaorig={0,0}];
If[imalabel==Automatic,imalabel = {Style["1-\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],Style["Im(\!\(\*SubscriptBox[\(A\), \(lm\)]\))",16]};];
iplot=Show[ListPlot[ImList,PlotRange->prange,AxesOrigin->iaorig,AxesLabel->imalabel,PlotLabel->plabel,FilterRules[FilterRules[{opts},Options[ListPlot]],Except[{Joined,PlotRange,AxesOrigin,AxesLabel,PlotLabel}]]],Plot[ImLOF[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0},PlotStyle->{{Red,Dashed}}],Plot[Imfit[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0}]];
{rplot,iplot}
]


(* ::Subsection::Closed:: *)
(*Routines for Spectral Coefficients*)


Options[SpectralCoef]=Options[ListPointPlot3D];
SetOptions[SpectralCoef,PlotLabel->Automatic,AxesLabel->Automatic,PlotRange->Automatic];
SpectralCoef[l_,m_,n_,opts:OptionsPattern[]]:=Module[{rawdat,Na,Nl,lmin,lmax,lind,amin,aind,alist,llist,plotline,plotlist={},label=OptionValue[PlotLabel],alabel=OptionValue[AxesLabel],range=OptionValue[PlotRange]},
SetOptions[ListPointPlot3D,TicksStyle->Directive[14]];
rawdat=ReadSpectralCoefs[l,m,n];
alist=Flatten[Take[rawdat,{1,-1},1]];
llist=Complex[#[[1]],#[[2]]]&/@Partition[#,2]&/@Drop[rawdat,0,1];
Na=Length[alist];
Nl=Length[llist[[1]]];
lmin=Max[Abs[m],2];
lmax=Nl+lmin-1;
Print["lmin = ",lmin," : lmax = ",lmax];
For[lind=lmin,lind<=lmax,++lind,
amin=1;
While[True,
If[Abs[llist[[amin,lind-lmin+1]]]>0,Break[],If[amin==Na,Break[],++amin]];
];
(* Print["**lind = ",lind," amin = ",amin]; *)
plotline=Table[{lind,alist[[i]],Log10[Abs[llist[[i,lind-lmin+1]]]]},{i,amin,Na}];
AppendTo[plotlist,plotline];
];
If[label==Automatic,label=Style["l="<>ToString[l]<>" m="<>ToString[m]<> " n="<>ToString[n],16]];
If[alabel==Automatic,alabel={Style["\!\(\*SuperscriptBox[\(\[ScriptL]\), \(\[Prime]\)]\)",20],Style["\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],""}];
If[range==Automatic,range={-15,0}];
ListPointPlot3D[plotlist,PlotRange->range,AxesLabel->alabel,PlotLabel->label,FilterRules[FilterRules[{opts},Options[ListPointPlot3D]],Except[{PlotRange,AxesLabel,PlotLabel}]]]
]


SpectralCoefAll[n_]:=Module[{l,m,rawdat,Nl,lmin,lmax},
For[l=2,l<=12,++l,
For[m=-l,m<=l,++m,
rawdat=ReadSpectralCoefs[l,m,n];
Nl=(Length[rawdat[[1]]]-1)/2;
lmin=Max[Abs[m],2];
lmax=Nl+lmin-1;
Print["(",l,",",m,",",n,") lmin = ",lmin," : lmax = ",lmax,"  [",Nl,"]"];
];
];
]


(* ::Subsection::Closed:: *)
(*Routine to find "Loop Zeros"*)


Options[LoopZeros]={ZeroTol->10^(-6),MaxFail->2};
LoopZeros[l_Integer,m_Integer,n_Integer,OptionsPattern[]]:=Module[{dat0,dat1,term0,term1,redat\[Omega],imdat\[Omega],redatAlm,imdatAlm,
refit\[Omega],imfit\[Omega],refitAlm,imfitAlm,a,azero,zeros,
i,lmindex,retst,decflag=False,failsmall=0,
tol=OptionValue[ZeroTol],maxfail=OptionValue[MaxFail]},
If[l==2 && 9<=n<=26,
(* Find the interpolation to the termination of the 0 sequence *)
dat0=Take[ReadKerrQNM[l,m,{n,0}],{1,-1},{1,5}];
redat\[Omega]=Function[x,{x[[1]],x[[2]]}]/@Take[dat0,{-5,-1}];
imdat\[Omega]=Function[x,{x[[1]],x[[3]]}]/@Take[dat0,{-5,-1}];
refit\[Omega]=Fit[redat\[Omega],{1,a},a];
(*Print[Show[ListLinePlot[redat\[Omega]],Plot[refit\[Omega],{a,redat\[Omega][[1,1]],redat\[Omega][[5,1]]}]]];*)
azero=Solve[refit\[Omega]==0,a][[1]];
imfit\[Omega]=Fit[imdat\[Omega],{1,a},a];
(*Print[Show[ListLinePlot[imdat\[Omega]],Plot[imfit\[Omega],{a,imdat\[Omega][[1,1]],imdat\[Omega][[5,1]]}]]];*)
redatAlm=Function[x,{x[[1]],x[[4]]}]/@Take[dat0,{-5,-1}];
imdatAlm=Function[x,{x[[1]],x[[5]]}]/@Take[dat0,{-5,-1}];
refitAlm=Fit[redatAlm,{1,a},a];
imfitAlm=Fit[imdatAlm,{1,a},a];
(*Print[Show[ListLinePlot[redat\[Omega]],Plot[refit\[Omega],{a,redat\[Omega][[1,1]],redat\[Omega][[5,1]]}]]];*)
(*Print[Show[ListLinePlot[imdat\[Omega]],Plot[imfit\[Omega],{a,imdat\[Omega][[1,1]],imdat\[Omega][[5,1]]}]]];*)
(*Print[Show[ListLinePlot[redatAlm],Plot[refitAlm,{a,redatAlm[[1,1]],redatAlm[[5,1]]}]]];*)
(*Print[Show[ListLinePlot[imdatAlm],Plot[imfitAlm,{a,imdatAlm[[1,1]],imdatAlm[[5,1]]}]]];*)
zeros={{a,imfit\[Omega],refitAlm,Length[dat0],refit\[Omega],imfitAlm,Abs[a-redat\[Omega][[-1,1]]],
	Sqrt[(refit\[Omega]-redat\[Omega][[-1,2]])^2+(imfit\[Omega]-imdat\[Omega][[-1,2]])^2],
	Sqrt[(refitAlm-redatAlm[[-1,2]])^2+(imfitAlm-imdatAlm[[-1,2]])^2]}}/.azero;
Clear[dat0]
,
zeros={};
];

(* Find the interpolation to the beginning of the 1 sequence *)
If[l==2 && 8<=n<=26,
dat1=Take[ReadKerrQNM[l,m,{n,1}],{1,-1},{1,5}];
redat\[Omega]=Function[x,{x[[1]],x[[2]]}]/@Take[dat1,{1,5}];
imdat\[Omega]=Function[x,{x[[1]],x[[3]]}]/@Take[dat1,{1,5}];
refit\[Omega]=Fit[redat\[Omega],{1,a},a];
azero=Solve[refit\[Omega]==0,a][[1]];
imfit\[Omega]=Fit[imdat\[Omega],{1,a},a];
redatAlm=Function[x,{x[[1]],x[[4]]}]/@Take[dat1,{1,5}];
imdatAlm=Function[x,{x[[1]],x[[5]]}]/@Take[dat1,{1,5}];
refitAlm=Fit[redatAlm,{1,a},a];
imfitAlm=Fit[imdatAlm,{1,a},a];
(*Print[Show[ListLinePlot[redat\[Omega]],Plot[refit\[Omega],{a,redat\[Omega][[1,1]],redat\[Omega][[5,1]]}]]];*)
(*Print[Show[ListLinePlot[imdat\[Omega]],Plot[imfit\[Omega],{a,imdat\[Omega][[1,1]],imdat\[Omega][[5,1]]}]]];*)
(*Print[Show[ListPlot[redatAlm],Plot[refitAlm,{a,redatAlm[[1,1]],redatAlm[[5,1]]}]]];*)
(*Print[Show[ListPlot[imdatAlm],Plot[imfitAlm,{a,imdatAlm[[1,1]],imdatAlm[[5,1]]}]]];*)
AppendTo[zeros,{a,imfit\[Omega],refitAlm,1,refit\[Omega],imfitAlm,Abs[a-redat\[Omega][[1,1]]],
	Sqrt[(refit\[Omega]-redat\[Omega][[1,2]])^2+(imfit\[Omega]-imdat\[Omega][[1,2]])^2],
	Sqrt[(refitAlm-redatAlm[[1,2]])^2+(imfitAlm-imdatAlm[[1,2]])^2]}/.azero],
dat1=Take[ReadKerrQNM[l,m,n],{1,-1},{1,5}];
];

(* Find the interpolations to the loops of the 1 sequence *)
redat\[Omega]=Function[x,{x[[1]],x[[2]]}]/@dat1;
imdat\[Omega]=Function[x,{x[[1]],x[[3]]}]/@dat1;
retst=redat\[Omega][[1,2]];
lmindex={};
For[i=2,i<=Length[redat\[Omega]],++i,
If[redat\[Omega][[i,2]]>retst && decflag,
If[Abs[retst]<tol,
AppendTo[lmindex,i-1];
If[failsmall>0,Print["Warning false pullaway at ",lmindex[[-1]]]]
,
If[Length[lmindex]>0,++failsmall]
];
decflag=False
,
If[redat\[Omega][[i,2]]<retst && Not[decflag],decflag=True];
];
If[failsmall>maxfail,Break[]];
retst=redat\[Omega][[i,2]];
];
For[i=1,i<=Length[lmindex],++i,
redat\[Omega]=Function[x,{x[[1]],x[[2]]}]/@Take[dat1,{lmindex[[i]]-2,lmindex[[i]]+2}];
imdat\[Omega]=Function[x,{x[[1]],x[[3]]}]/@Take[dat1,{lmindex[[i]]-2,lmindex[[i]]+2}];
refit\[Omega]=Interpolation[redat\[Omega],InterpolationOrder->3,Method->"Spline"];
(*Print[Show[ListLinePlot[redat],Plot[refit[a],{a,redat[[1,1]],redat[[5,1]]}]]];*)
azero=FindRoot[refit\[Omega]'[a]==0,{a,redat\[Omega][[3,1]]}][[1]];
imfit\[Omega]=Fit[imdat\[Omega],{1,a},a];
redatAlm=Function[x,{x[[1]],x[[4]]}]/@Take[dat1,{lmindex[[i]]-2,lmindex[[i]]+2}];
imdatAlm=Function[x,{x[[1]],x[[5]]}]/@Take[dat1,{lmindex[[i]]-2,lmindex[[i]]+2}];
refitAlm=Fit[redatAlm,{1,a},a];
imfitAlm=Interpolation[imdatAlm,InterpolationOrder->3,Method->"Spline"];
AppendTo[zeros,{a,imfit\[Omega],refitAlm,lmindex[[i]],refit\[Omega][a],imfitAlm[a],Abs[a-redat\[Omega][[3,1]]],
	Sqrt[(refit\[Omega][a]-redat\[Omega][[3,2]])^2+(imfit\[Omega]-imdat\[Omega][[3,2]])^2],
	Sqrt[(refitAlm-redatAlm[[3,2]])^2+(imfitAlm[a]-imdatAlm[[3,2]])^2]}/.azero];
];
zeros
]


(* ::Subsection::Closed:: *)
(*Routines to test for Non-generic modes*)


Options[\[Sigma]ptestKerrQNM]={ModReal->True,TruncateLongSequences->True};
\[Sigma]ptestKerrQNM[qnm_List,range_List:{1,-1},opts:OptionsPattern[]]:=
Module[{rawdat,Nelems,\[Omega],a,m,\[Sigma]ptest,lr=range,modreal=OptionValue[ModReal],trunc=OptionValue[TruncateLongSequences]},
Off[Import::dataset];
rawdat=ReadKerrQNM[qnm];
On[Import::dataset];
If[rawdat==$Failed,Return[$Failed]];
If[lr=={1,-1} && trunc,
If[qnm=={2,0,{11,1}},lr={1,20000}];If[qnm=={2,0,{12,1}},lr={1,22000}];If[qnm=={2,0,{13,1}},lr={1,35000}];If[qnm=={2,0,{14,1}},lr={1,35000}];
If[qnm=={2,0,{15,1}},lr={1,35000}];If[qnm=={2,0,{16,1}},lr={1,25000}];If[qnm=={2,0,{17,1}},lr={1,25000}];If[qnm=={2,0,{18,1}},lr={1,25000}];
If[qnm=={2,0,{19,1}},lr={1,25000}];If[qnm=={2,0,{20,1}},lr={1,25000}];If[qnm=={2,0,{21,1}},lr={1,25000}];If[qnm=={2,0,{22,1}},lr={1,25000}];
If[qnm=={2,0,{23,1}},lr={1,25000}];If[qnm=={2,0,{24,1}},lr={1,25000}];If[qnm=={2,0,{25,1}},lr={1,25000}];If[qnm=={2,0,{26,1}},lr={1,25000}];
If[qnm=={3,0,19},lr={1,20000}];If[qnm=={3,0,20},lr={1,20000}];If[qnm=={3,0,21},lr={1,22000}];If[qnm=={3,0,22},lr={1,35000}];
If[qnm=={3,0,23},lr={1,35000}];If[qnm=={3,0,24},lr={1,25000}];If[qnm=={3,0,25},lr={1,25000}];If[qnm=={3,0,26},lr={1,25000}];
If[qnm=={3,0,27},lr={1,25000}];If[qnm=={3,0,28},lr={1,25000}];If[qnm=={3,0,29},lr={1,25000}];If[qnm=={3,0,30},lr={1,25000}];
If[qnm=={3,0,31},lr={1,25000}];
If[qnm=={4,0,26},lr={1,22000}];If[qnm=={4,0,27},lr={1,25000}];If[qnm=={4,0,28},lr={1,35000}];If[qnm=={4,0,29},lr={1,35000}];
If[qnm=={4,0,30},lr={1,25000}];If[qnm=={4,0,31},lr={1,25000}];
];
\[Omega]=Function[x,x[[1]]+I x[[2]]]/@Take[rawdat,lr,{2,3}];
a=Flatten[Take[rawdat,lr,1]];
m=qnm[[2]];
\[Sigma]ptest=I (2\[Omega](1+Sqrt[1-a^2])-m a)/Sqrt[1-a^2];
If[modreal,{Transpose[{a,Mod[Re[\[Sigma]ptest],1]}],Transpose[{a,Im[\[Sigma]ptest]}]},
{Transpose[{a,Re[\[Sigma]ptest]}],Transpose[{a,Im[\[Sigma]ptest]}]}]
]
\[Sigma]ptestKerrQNM[l_Integer,m_Integer,n_Integer|n_List,range_List:{1,-1},opts:OptionsPattern[]]:=\[Sigma]ptestKerrQNM[{l,m,n},range,opts]


NonGenericMode[qnm_List,\[Epsilon]_]:=
Module[{\[Sigma]list,len,i,sign,lastsign,dec,lastdec,imin,imax,a,fit,coefs,amin},
\[Sigma]list=\[Sigma]ptestKerrQNM[qnm,{1,1000},ModReal->False,TruncateLongSequences->False];
sign=Sign[\[Sigma]list[[2,1,2]]];
dec=\[Sigma]list[[2,2,2]]<\[Sigma]list[[2,1,2]];
len=Length[\[Sigma]list[[1]]];
For[i=2,i<=len,++i,
(* look for sign flip *)
lastsign=sign;
lastdec=dec;
sign=Sign[\[Sigma]list[[2,i,2]]];
dec=\[Sigma]list[[2,i,2]]<\[Sigma]list[[2,i-1,2]];
(*Print[i" (",sign,",",lastsign,") (",dec,",",lastdec,")"];*)
If[lastsign sign<=0,
sign=Sign[\[Sigma]list[[1,i,2]]];
Print["Sign change at ",i];
If[\[Sigma]list[[1,i,2]]<\[Epsilon]||1-\[Sigma]list[[1,i,2]]<1\[Epsilon],
Print["Possible non-generic mode of ",qnm," at a = ",\[Sigma]list[[1,i,1]]," : (",i,")"];
],
(* look for minimum of zero *)
If[lastdec && !dec,
imin=Max[1,i-2];imax=Min[len,imin+2];imin=Max[1,imax-2];
fit=Fit[Take[\[Sigma]list[[2]],{imin,imax}],{1,a,a^2},a];
Print[Show[ListPlot[Take[\[Sigma]list[[2]],{imin,imax}],PlotRange->All],Plot[fit,{a,\[Sigma]list[[2,imin,1]],\[Sigma]list[[2,imax,1]]}]]];
coefs=CoefficientList[fit,a];
amin=-coefs[[2]]/(2coefs[[3]]);
Print[" a = ",amin," : coefs = ",coefs];
If[\[Sigma]list[[1,i,2]]<\[Epsilon]||1-\[Sigma]list[[1,i,2]]<1\[Epsilon],
Print["Possible non-generic mode of ",qnm," at a = ",\[Sigma]list[[1,i,1]]," : (",i,")"];
];
],
Print["Bad Search"];
];
];
]


(* ::Section::Closed:: *)
(*End of ModePlotHDF5 Package*)


End[] (* `Private` *)


EndPackage[]
