(* ::Package:: *)

(* ::Title:: *)
(* Quasi-Normal Modes of Kerr*)


(* ::Section::Closed:: *)
(*Documentation *)


(* ::Item:: *)
(*Overview*)


(* ::ItemParagraph:: *)
(*The KerrQNM Mathematica package provides functions to compute and plot the Quasi-Normal Modes of gravitational (s=-2), electromagnetic (s=-1), and scalar (s=0) perturnbations of the Kerr geometry.*)


(* ::ItemParagraph:: *)
(*QNM sequences for a specific (l,m) mode and overtone n are parameterized by the dimensionless spin parameter a \[Element][0,1).  Sequences usually start with a=0 using initial guesses stored in tables of the Schwarzschild QNMs.*)


(* ::Item:: *)
(*Spin Weight*)


(* ::ItemParagraph:: *)
(*All Package functions take the option QNMSpinWeight.  By default, QNMSpinWeight\[RightArrow]Null and all functions require the spin weight to be set before they will work.  The Package function:*)


(* ::Subsubitem:: *)
(*SetSpinWeight[s]*)


(* ::ItemParagraph:: *)
(*will set the default value for QNMSpinWeight to s for all Package functions.  Normal workflow is to load the package and then call SetSpinWeight.  For example:*)


(* ::Program:: *)
(*<< "full_path/KerrQNM.m";*)
(*SetSpinWeight[-2]*)


(* ::Item:: *)
(*Data Storage and Indexing:*)


(* ::ItemParagraph:: *)
(*Data is stored in Global variables and it is the responsibility of the user to save and retreive data as needed. *)


(* ::Subitem:: *)
(*Schwarzschild Data:*)


(* ::SubitemParagraph:: *)
(*Data for the Schwarzschild QNMs are stored in one of three global variables depending on the spin-weight s being considered.*)


(* ::Subsubitem:: *)
(*s = -2 : SchQNMTable*)


(* ::Subsubitem:: *)
(*s = -1 : SchQNMeTable*)


(* ::Subsubitem:: *)
(*s = 0 : SchQNMsTable*)


(* ::Subitem:: *)
(*Schwarzschild Data Indexing:*)


(* ::SubitemParagraph:: *)
(*The Schwarzschild data tables are indexed by both single and double integer indices:*)


(* ::SubsubitemParagraph:: *)
(*SchQNMTable[l] : stores the total number of overtones stored for the specified mode l.*)


(* ::SubsubitemParagraph:: *)
(*SchQNMTable[l,n] : stores the QNM for mode l, overtone n.  l \[Element] Integers >= |s| ; n \[Element] {0,1,...,SchQNMTable[l]-1}*)


(* ::Subitem:: *)
(*Kerr Data:*)


(* ::SubitemParagraph:: *)
(*Data for the Kerr QNMs sequences are stored in one of three global variables depending on the spin-weight s being considered.*)


(* ::Subsubitem:: *)
(*s = -2 : KerrQNM*)


(* ::Subsubitem:: *)
(*s = -1 : KerrQNMe*)


(* ::Subsubitem:: *)
(*s = 0 : KerrQNMs*)


(* ::Subitem:: *)
(*Kerr Data Indexing:*)


(* ::SubitemParagraph:: *)
(*The Kerr data tables are indexed by an index triplet.*)


(* ::SubitemParagraph:: *)
(*KerrQNMTable[l,m,n] : stores a List of QNMs for mode (l,m) and overtone n.  l \[Element]Integers >= Max(|s|,|m|) ; n is usually an integer >= 0, but can be a two element List {n,i} denoting an element of an overtone multiplet.  Each element in the QNM List specifies the QNM for a specific value of the dimensionless spin a.*)


(* ::ItemParagraph:: *)
(*All data currently in memory can be written to a file using the Save command and retreived from a file using <<.  For example:*)


(* ::Program:: *)
(* Save["full_path_and_filename",SchQNMTable];*)


(* ::ItemParagraph:: *)
(*and*)


(* ::Program:: *)
(*<< "full_path_and_filename";*)


(* ::ItemParagraph:: *)
(*The KerrQNMTable files can become quite large, so it is recommended that users split the storage among multiple files.*)


(* ::Item:: *)
(*Schwarzschild Tables*)


(* ::ItemParagraph:: *)
(*The construction of Kerr QNM sequences requires initial quesses which are usually obtained from tables of Schwarzschild QNMs as described above.  These tables can be constructed using the Package function:*)


(* ::Subsubitem:: *)
(*SchwarzschildQNM[l,n]*)


(* ::SubitemParagraph:: *)
(*If a table entry exists for mode l and overtone n, then calling this function uses the existing entry as a starting guess and recomputes the mode to an accuracy of (10^-14).  If a table entry does not exist, then extrapolation is used to obtain an starting guess which is then refined.*)


(* ::Subitem:: *)
(*Initial Guesses:*)


(* ::SubitemParagraph:: *)
(*At least 4 initial guesses must be manually entered in the Schwarzschild table before SchwarzschildQNM[l,n] will work properly.  Entries for l=|s|, n=0, 1 and l=|s|+1, n=0, 1 are required.  Additional initial guesses may be required if extrapolation does not produce guesses sufficiently close to the given mode.*)


(* ::Subitem:: *)
(*Algebraically Special Modes:*)


(* ::SubitemParagraph:: *)
(*SchwarzschildQNM[l,n] cannot obtain solutions at the algebraically speical frequencies.  At these frequencies, which are characterized by a purely imaginary frequency, the Schwarzschild table entry must be entered manually.*)


(* ::Subitem:: *)
(*Plotting*)


(* ::SubitemParagraph:: *)
(*A plot of all overtones for mode l can be obtained by calling the Package function:*)


(* ::Subsubitem:: *)
(*QNMPlotSch[l]*)


(* ::Item:: *)
(*Kerr QNM Sequences*)


(* ::ItemParagraph:: *)
(*Sequences usually start with a=0 using initial guesses stored in tables of the Schwarzschild QNMs.  Such sequences can be constructed using the Package function:*)


(* ::Subsubitem:: *)
(*KerrQNMSequence[l,m,n,\[Epsilon]]*)


(* ::SubitemParagraph:: *)
(*If a table entry exists for mode (l,m) and overtone n, then the sequence is extended starting from the largest value of the dimensionless spin parameter a.  If no table entry exists, a new sequence is started from a=0.  Solutions are found to an accuracy of (10^\[Epsilon]).  A step size of \[CapitalDelta]a=10^-3 is used unless finer steps are required.  Adaptive step sizing will reduce the step size to as small as  \[CapitalDelta]a=10^-6 as needed.*)


(* ::Subitem:: *)
(*Special Sequences*)


(* ::SubitemParagraph:: *)
(*Under certain circumstances, a sequence cannot start at a=0.  If a table entry does not exist for mode (l,m) and overtone n, then the option*)


(* ::Subsubitem:: *)
(*QNMaStart \[RightArrow] {a,\[Omega],Alm}*)


(* ::SubitemParagraph:: *)
(*can be given to KerrQNMSequence[l,m,n,\[Epsilon]].  This specifies that the sequence should be started with the specified value of the dimensionless spin parameter, and the given initial guesses for the complex frequency, and separation constant.  An optional 4th arguement is allowed, specifying the initial Integer size of the spectral matrix used to solve the angular Teukolsky Equation(see below).*)


(* ::SubitemParagraph:: *)
(*Any sequence that does not begin with a=0 can be extended to smaller values of a using the Package function:*)


(* ::Subsubitem:: *)
(*KerrQNMSequenceReverse[l,m,n,\[Epsilon]]*)


(* ::SubitemParagraph:: *)
(*KerrQNMSequenceReverse will refine its stepsize as it approaches a=0, but will not attempt to evaluate a solutions at a=0. Note that currently, KerrQNMSequenceReverse will use as the initial stepsize, the value of \[CapitalDelta]a found between the first two entries in the sequence List, requires that there be at least 3 entries in the List, and assumes that the stepsize between the first three entries is the same.   A factor of 10 refinement in the stepsize can be forced by setting the optional last argument to True:*)


(* ::Subsubitem:: *)
(*KerrQNMSequenceReverse[l,m,n,\[Epsilon],True]*)


(* ::SubitemParagraph:: *)
(*Some mode sequences asymptotically approach the real limiting frequency of \[Omega]=m/2 as a\[RightArrow]1.  In these cases, a large number of overtones are in close proximity requiring very accurate initial guesses in order to assure solutions are found with the desired overtone.  An asymptotic fitting function can be used to obtain accurate initial guesses.  The Package function:*)


(* ::Subsubitem:: *)
(*KerrQNMAccumulation[l,m,n,\[Epsilon],overtones,Nv,a0]*)


(* ::SubitemParagraph:: *)
(*will extend the sequence for mode (l,m) and overtone n in the direction of increasing a.  Solutions are found to an accuracy of (10^\[Epsilon]).  overtones is a List of overtone sequences of the mode (l,m) that are used to set the parameters of the fitting functions.  The mode being extended can be an element of overtones.  Only elements of the sequences listed in overtones with a>=1-a0 will be used to define the fit, and at most the last Nv elements will be used.*)


(* ::SubitemParagraph:: *)
(*The fitting function used in KerrQNMAccumulation can be of two forms depending on whether or not the fundamental n=0 mode asymptotes to \[Omega]=m/2.  In order to obtain the correct fitting function, it is essential that the lowest overtone that asymptotes to \[Omega]=m/2 for mode (l,m) must be the first element of the List overtones.*)


(* ::Item:: *)
(*Eigenvalue Solvers*)


(* ::ItemParagraph:: *)
(*Kerr QNM sequences are constructed from individual solutions at specific values of the dimensionless angular momentum 'a'.  These individual solutions require solving a coupled pair of eigenvalue problems.  These are solved iteratively by the Package function*)


(* ::Subitem:: *)
(*QNMSolution[n,s,l,m,a,\[Omega]g,Almg,\[Epsilon],relax,Nrcf,Nm,\[Omega]0,Alm0,rl,rt]*)


(* ::ItemParagraph:: *)
(*The output of QNMSolution is a 3-element List*)


(* ::Program:: *)
(*{Valid,\[Alpha],solution}*)


(* ::Subsubitem:: *)
(*Valid is a Boolean.  If False, the iterative solver failed to converge to the desired accuracy.*)


(* ::Subsubitem:: *)
(*\[Alpha] is the final value of the under-relaxation parameter used in the iterative solver.*)


(* ::Subsubitem:: *)
(*solution is a 3-element list {a,radial_solution,angular_solution}*)


(* ::Subsubitem:: *)
(*a is the value of the dimensionless angular momentum used in computing this QNM.*)


(* ::Subsubitem:: *)
(*radial_solution is the second element of the output of RadialLentzRoot described below.*)


(* ::Subsubitem:: *)
(*angular_solution is the output of AngularSpectralRoot described below.*)


(* ::ItemParagraph:: *)
(*QNMSolution calls the individual eigenvalue solvers for the angular and radial Teukolsky equations.  The angular equation is solve directly using a spectral eigenvalue method implemented by the Package function*)


(* ::Subitem:: *)
(*AngularSpectralRoot[s,m,c,Alm,N]*)


(* ::ItemParagraph:: *)
(*The output of AngularSpectralRoot is the 3-element List*)


(* ::Program:: *)
(*{Alm, N, Coefficients}*)


(* ::Subsubitem:: *)
(*Alm is the complex separation constant*)


(* ::Subsubitem:: *)
(*N is the dimension of the spectral eigenvalue matrix*)


(* ::Subsubitem:: *)
(*Coefficients is an N-element list of expansion coefficients*)


(* ::SubsubitemParagraph:: *)
(*The expansion of the solution is in terms of spin-weighted spherical harmonics.  For any overtone of mode (l,m), the expansion includes harmonics beginning at l'=Max(|m|,|s|) and are truncated at l'=N+Max(|m|,|s|)-1.  The solution is normalized so that the coefficient for l'=l is real, and the sum of the magnitudes of the coefficients is one.*)
(**)


(* ::ItemParagraph:: *)
(*The radial equation is solve by finding the roots of a continued fraction.  The roots are obtained via a 2-dimensional Newton's method implemented by the Package function*)


(* ::Subitem:: *)
(*RadialLentzRoot[n,s,m,a,Alm,\[Omega],N,\[Epsilon],Radius]*)


(* ::ItemParagraph:: *)
(*The output of RadialLentzRoot is the 2-element List containg diagnostic information from the Newton solver and the solution*)


(* ::Program:: *)
(*{diagnostics, radial_solution}*)


(* ::ItemParagraph:: *)
(*diagnostics is a 2-element List of Booleans {converged,slow}.*)


(* ::Subsubitem:: *)
(*converged is True when the solution has converged to the desired accuracy.*)


(* ::Subsubitem:: *)
(*slow is False for normal convergence.  If True, then the convergence is slow.*)


(* ::ItemParagraph:: *)
(*radial_solution is a 5-element list {\[Omega],Ninv,Nrcf,\[Epsilon],residual}.*)


(* ::Subsubitem:: *)
(*\[Omega] is the complex frequency of the QNM*)


(* ::Subsubitem:: *)
(*Ninv is the 'inversion' of the continued fraction used*)


(* ::Subsubitem:: *)
(*Nrcf is the depth at which the radial continued fraction is truncated.*)


(* ::Subsubitem:: *)
(*\[Epsilon] is the requested accuracy of the solution \[Omega].  A solution will have converged to the point where the magnitude of successive changes to \[Omega] are less than 10^\[Epsilon]. *)


(* ::Subsubitem:: *)
(*residual is the magnitude of the last change made to \[Omega].*)


(* ::Item:: *)
(*Plotting*)


(* ::ItemParagraph:: *)
(*Several Package functions are available to plot the the various mode sequences.  *)
(**)
(*The mode frequencies \[Omega] are usually plotted with the imaginary axis inverted.  That is, the complex conjugate of \[Omega] is plotted.  The functions used to plot the frequency are:*)


(* ::Subsubitem:: *)
(*QNMPlotOmega[l,m,n]*)


(* ::SubsubitemParagraph:: *)
(*Plots the single sequence for overtone n (can be a multiplet index) of mode (l,m).*)


(* ::Subsubitem:: *)
(*QNMPlotOmega[l,n]*)


(* ::SubsubitemParagraph:: *)
(*Plots overtone n of all modes (l,m) with |m|<=l.  Overtone multiplets must be designated.*)


(* ::Subsubitem:: *)
(*QNMPlotOmegaTones[l,m]*)


(* ::SubsubitemParagraph:: *)
(*Plots all overtones of modes (l,m).  Overtone multiplets must be designated.*)


(* ::Subsubitem:: *)
(*QNMPlotAccumulation\[Omega][l,m,overtones,Nv,a0]*)


(* ::SubsubitemParagraph:: *)
(*Computes and plots fitting functions for the real and imaginary parts of the specified overtones of modes (l,m) versus 1-a over the range (0,a0).  Overtone multiplets must be designated.  The fitting functions are*)


(* ::SubsubitemParagraph:: *)
(*	If 0\[Element]{overtones}*)


(* ::Program:: *)
(*			Re[\[Omega]] : m/2 - \[Alpha]1*Sqrt[\[Epsilon]/2] + (\[Alpha]2 + \[Alpha]3*n)*\[Epsilon]*)
(*			Im[\[Omega]] : -(n+1/2)*(Sqrt[\[Epsilon]/2] - \[Alpha]4*\[Epsilon])*)


(* ::SubsubitemParagraph:: *)
(*	otherwise*)


(* ::Program:: *)
(*			Re[\[Omega]] : m/2 + (\[Alpha]1 + \[Alpha]2*n1)*\[Epsilon]*)
(*			Im[\[Omega]] : -(\[Alpha]3 + n + 1/2)*Sqrt[\[Epsilon]/2] + (\[Alpha]4 + \[Alpha]5*n)*\[Epsilon]*)


(* ::SubsubitemParagraph:: *)
(*Symbols denote the data points used in the fit.  Solid lines denote the full fit function, and dashed red lines indicate the fit function including only the terms through leading order in \[Epsilon]=1-a.*)


(* ::Subsubitem:: *)
(*KerrOmegaList[l,m,n]*)


(* ::SubsubitemParagraph:: *)
(*Utility routine that returns the List of points corresponding to the complex conjugate of \[Omega] for overtone n of mode (l,m).  All data points in the sequence are returned.*)


(* ::Subsubitem:: *)
(*KerrOmegaListS[l,m,n]*)


(* ::SubsubitemParagraph:: *)
(*Utility routine that returns the List of points corresponding to the complex conjugate of \[Omega] for overtone n of mode (l,m). Only data points where 'a' is a multiple of 0.05 are returned (and the last point near a=1).*)


(* ::Subsubitem:: *)
(*SchwarzschildOmega[l,m,n]*)


(* ::SubsubitemParagraph:: *)
(*Utility routine that returns the single point corresponding to the complex conjugate of \[Omega] in the Schwarzschild limit for overtone n of mode (l,m).*)


(* ::ItemParagraph:: *)
(*The corresponding functions used to plot the separation constant Alm are:*)


(* ::Subsubitem:: *)
(*QNMPlotA[l,m,n]*)


(* ::Subsubitem:: *)
(*QNMPlotA[l,n]*)


(* ::Subsubitem:: *)
(*QNMPlotATones[l,m]*)


(* ::Subsubitem:: *)
(*QNMPlotAccumulationAlm[l,m,overtones,Nv,a0]*)


(* ::SubsubitemParagraph:: *)
(*The fitting functions are*)


(* ::SubsubitemParagraph:: *)
(*	If 0\[Element]{overtones}*)


(* ::Program:: *)
(*			Re[Alm] : l(l+1)-s(s+1) + \[Beta]1 + \[Beta]2*Sqrt[\[Epsilon]/2]+(\[Beta]3 + \[Beta]4*n + \[Beta]5*n^2)*\[Epsilon]*)
(*			Re[Alm] : (n+1/2)*(\[Beta]6*Sqrt[\[Epsilon]/2] + \[Beta]7*\[Epsilon])*)


(* ::SubsubitemParagraph:: *)
(*	otherwise*)


(* ::Program:: *)
(*			Re[Alm] : l(l+1)-s(s+1) + \[Beta]1 + (\[Beta]2 + \[Beta]3*n + \[Beta]4*n^2)*\[Epsilon]*)
(*			Im[Alm] : (\[Beta]5 + n + 1/2)*\[Beta]6*Sqrt[\[Epsilon]/2] + (\[Beta]7 + \[Beta]8*n)*\[Epsilon]*)


(* ::Subsubitem:: *)
(*KerrAList[l,m,n]*)


(* ::Subsubitem:: *)
(*KerrAListS[l,m,n]*)


(* ::Section::Closed:: *)
(*Begin KerrQNM Package*)


BeginPackage["KerrQNM`"]


Unprotect[QNMDebug];
QNMDebug=True; (* Set this to True to allow reloading of Package with changes *)
If[QNMDebug,Unprotect["KerrQNM`*"];Unprotect["KerrQNM`Private`*"]];
Protect[QNMDebug];


(* ::Section:: *)
(*Documentation of External Functions*)


(* ::Subsection:: *)
(*Sequencers*)


KerrQNMSequenceB::usage=""


KerrQNMRefineSequenceB::usage=""KerrQNMSequence[l,m,n,\[Epsilon]] computes a sequence of Quasi-Normal Mode solutions "<>
	"for overtone n of mode (l,m).  The solutions are computed to an absolute "<>
	"accuracty of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  The sequence is "<>
	"parameterized by increasing values of the dimensionless angular momentum "<>
	"'a' starting at a=0 (or the largest value of 'a' already computed) up to "<>
	"(but not including) a=1.\n \n "<>


KerrQNMSequence::usage=
	"KerrQNMSequence[l,m,n,\[Epsilon]] computes a sequence of Quasi-Normal Mode solutions "<>
	"for overtone n of mode (l,m).  The solutions are computed to an absolute "<>
	"accuracty of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  The sequence is "<>
	"parameterized by increasing values of the dimensionless angular momentum "<>
	"'a' starting at a=0 (or the largest value of 'a' already computed) up to "<>
	"(but not including) a=1.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  'n' "<>
	"can be either an Integer or an overtone multiplet index.  An overtone multilplet index "<>
	"is a 2 element list {n,mult}, where 'n' is the Integer overtone number, and 'mult' is "<>
	"in the range 0,1,...,(Nmult-1), with 'Nmult' the number of sequences with the same overtone "<>
	"index.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t Min\[CapitalDelta]alevel\[Rule]1 : 1,2,3,4,5,6\n"<>
	"\t Max\[CapitalDelta]alevel\[Rule]4 : 1,2,3,4,5,6\n"<>
	"\t\t The step size \[CapitalDelta]a is set to \!\(\*SuperscriptBox[\(10\), \(-\((2 + \[CapitalDelta]alevel)\)\)]\)\n"<>
	"\t QNMaStart\[Rule]0 : {a,\[Omega],\!\(\*SubscriptBox[\(A\), \(lm\)]\)}\n"<>
	"\t\t This option is only used when starting a new sequence.\n"<>
	"\t\t The List contains the initial value of 'a', and initial guesses for \[Omega] and \!\(\*SubscriptBox[\(A\), \(lm\)]\).\n"<>
	"\t\t An option 4th argument, specifying the initial Integer size of the spectral\n"<>
	"\t\t matrix used to solve the angular Teukolsky equation.\n"<>
	"\t QNMPrecision\[Rule]24\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t The radial solver make use of numerical derivatives.  The relative step size\n"<>
	"\t\t is set to \!\(\*SuperscriptBox[\(10\), \(JacobianStep\)]\).\n"<>
	"\t RadialCFDepth\[Rule]1\n"<>
	"\t\t The initial Radial Continued Fraction Depth is usually taken from the prior \n"<>
	"\t\t solution in the sequence.  Fractional values reduce this initial vaue by that\n"<>
	"\t\t fraction.  Integer values larger than 300 are used as the new initial value.\n"<>
	"\t RadialCFDigits\[Rule]8\n"<>
	"\t\t Approximate number of digits of agreement when testing the continued fraction\n"<>
	"\t\t pproximation.  The continued fraction is evaluated in the proximity of the \n"<>
	"\t\t guessed root, then evaluated again with a depth half-again deeper.  The first\n"<>
	"\t\t RadialCFDigits non-vanishing digits must agree for the solution to be accepted.\n"<>
	"\t SolutionRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial/Angular iterations.\n"<>
	"\t RadialRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial Newton iterations.\n"<>
	"\t SolutionDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial/Angular iterations.\n"<>
	"\t RadialDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial Newton iterations.\n"<>
	"\t SolutionWindowl\[Rule]1/2\n"<>
	"\t SolutionWindowt\[Rule]1/3\n"<>
	"\t\t Set the size of the solution window.  Solutions must fall within this window\n"<>
	"\t\t to be accepted.  Setting either to 0 causes all solutions to be accepted."


KerrQNMSequenceReverse::usage=
	"KerrQNMSequenceReverse[l,m,n,\[Epsilon]] or\n"<>
	"KerrQNMSequenceReverse[l,m,n,\[Epsilon],refine] computes "<>
	"a sequence of Quasi-Normal Mode solutions for overtone n of mode (l,m).  The solutions "<>
	"are computed to an absolute accuracty of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  The "<>
	"sequence is parameterized by decreasing values of the dimensionless angular momentum 'a' "<>
	"starting at the smallest value of 'a' already computed and continuing down to a=\!\(\*SuperscriptBox[\(10\), \(-8\)]\).  "<>
	"If 'refine' is present and set to True, the step size, \[CapitalDelta]a, is reduced by 1/10.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  'n' "<>
	"can be either an Integer or an overtone multiplet index.  An overtone multilplet index "<>
	"is a 2 element list {n,mult}, where 'n' is the Integer overtone number, and 'mult' is "<>
	"in the range 0,1,...,(Nmult-1), with 'Nmult' the number of sequences with the same overtone "<>
	"index.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t QNMPrecision\[Rule]24\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t The radial solver make use of numerical derivatives.  The relative step size\n"<>
	"\t\t is set to \!\(\*SuperscriptBox[\(10\), \(JacobianStep\)]\).\n"<>
	"\t RadialCFDepth\[Rule]1\n"<>
	"\t\t The initial Radial Continued Fraction Depth is usually taken from the prior \n"<>
	"\t\t solution in the sequence.  Fractional values reduce this initial vaue by that\n"<>
	"\t\t fraction.  Integer values larger than 300 are used as the new initial value.\n"<>
	"\t SolutionRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial/Angular iterations.\n"<>
	"\t RadialRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial Newton iterations.\n"<>
	"\t SolutionDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial/Angular iterations.\n"<>
	"\t RadialDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial Newton iterations."


KerrQNMInsert::usage=
	"KerrQNMInsert[l,m,n,\[Epsilon],a0] computes additional Quasi-Normal Mode solutions "<>
	"for overtone n of mode (l,m).  The solutions are computed to an absolute "<>
	"accuracty of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  The sequence is "<>
	"parameterized by increasing values of the dimensionless angular momentum "<>
	"'a'.  Nine new solutions are computed between a=a0 and the next computed value of a.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  'n' "<>
	"can be either an Integer or an overtone multiplet index.  An overtone multilplet index "<>
	"is a 2 element list {n,mult}, where 'n' is the Integer overtone number, and 'mult' is "<>
	"in the range 0,1,...,(Nmult-1), with 'Nmult' the number of sequences with the same overtone "<>
	"index.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t QNMPrecision\[Rule]24\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t The radial solver make use of numerical derivatives.  The relative step size\n"<>
	"\t\t is set to \!\(\*SuperscriptBox[\(10\), \(JacobianStep\)]\).\n"<>
	"\t RadialCFDepth\[Rule]1\n"<>
	"\t\t The initial Radial Continued Fraction Depth is usually taken from the prior \n"<>
	"\t\t solution in the sequence.  Fractional values reduce this initial vaue by that\n"<>
	"\t\t fraction.  Integer values larger than 300 are used as the new initial value.\n"<>
	"\t RadialCFDigits\[Rule]8\n"<>
	"\t\t Approximate number of digits of agreement when testing the continued fraction\n"<>
	"\t\t pproximation.  The continued fraction is evaluated in the proximity of the \n"<>
	"\t\t guessed root, then evaluated again with a depth half-again deeper.  The first\n"<>
	"\t\t RadialCFDigits non-vanishing digits must agree for the solution to be accepted.\n"<>
	"\t SolutionRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial/Angular iterations.\n"<>
	"\t RadialRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial Newton iterations.\n"<>
	"\t SolutionDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial/Angular iterations.\n"<>
	"\t RadialDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial Newton iterations."


KerrQNMAccumulation::usage=
	"KerrQNMAccumulation[l,m,n,\[Epsilon],overtones,Nv,a0] computes a sequence of Quasi-Normal "<>
	"Mode solutions for overtone n of mode (l,m).  The solutions are computed to an "<>
	"absolute accuracty of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  The sequence is "<>
	"parameterized by increasing values of the dimensionless angular momentum 'a', "<>
	"starting at the largest value of 'a' already computed up to a=1-\!\(\*SuperscriptBox[\(10\), \(\(-\)\(8.\)\(\\\ \)\)]\) "<>
	"Initial guesses are based on expansions for \[Omega] and \!\(\*SubscriptBox[\(A\), "<>
	"\(lm\)]\) valid near a=1.  The fitting function is based on existing data for the "<>
	"(l,m) modes with overtone values given in the List 'overtones'.  Nv and a0 are used "<>
	"to restrict the number of points in each (l,m,n) sequence that are used to "<>
	"determine the fit.  Only the last Nv points in each sequence are used, and only if "<>
	"the value of 'a' is greater than 1-a0.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  'n' "<>
	"and the elements of 'overtones' can be either an Integer or an overtone multiplet "<>
	"index.  An overtone multiplet index is a 2 element list {n,mult}, where 'n' is the "<>
	"Integer overtone number, and 'mult' is in the range 0,1,...,(Nmult-1), with 'Nmult' "<>
	"the number of sequences with the same overtone index.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t OTskip\[Rule]{}\n"<>
	"\t\t List of overtone values or multiplet indices that do not approach \n"<>
	"\t\t the limit point.\n"<>
	"\t OTmultiple\[Rule]{}\n"<>
	"\t\t List of overtone multiplets.  An overtone multiplet is a List {n,Nmult},\n"<>
	"\t\t where 'n' is the overtone index,and 'Nmult' is the number of sequences\n"<>
	"\t\t with the same overtone index.\n"<>
	"\t QNMPrecision\[Rule]24\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t The radial solver make use of numerical derivatives.  The relative step size\n"<>
	"\t\t is set to \!\(\*SuperscriptBox[\(10\), \(JacobianStep\)]\).\n"<>
	"\t RadialCFDepth\[Rule]1\n"<>
	"\t\t The initial Radial Continued Fraction Depth is usually taken from the prior \n"<>
	"\t\t solution in the sequence.  Fractional values reduce this initial vaue by that\n"<>
	"\t\t fraction.  Integer values larger than 300 are used as the new initial value.\n"<>
	"\t SolutionRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial/Angular iterations.\n"<>
	"\t RadialRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial Newton iterations.\n"<>
	"\t SolutionDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial/Angular iterations.\n"<>
	"\t RadialDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial Newton iterations."


SchwarzschildQNM::usage=
	"SchwarzschildQNM[l,n] computes the Quasi-Normal Mode solutions for overtone n of mode l.  "<>
	"The mode is computed to an accuracy of \!\(\*SuperscriptBox[\(10\), \(-14\)]\).  For given 'l', if solutions with overtones "<>
	"(n-1) and (n-2) have not been computed, then the routine is recersively called for overtone "<>
	"'(n-1).  If no solutions exist for mode l, then the first two overtones of (l-1) and (l-2) "<>
	"are used to extrapolate initial guesses for these modes.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t QNMPrecision\[Rule]24\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t The radial solver make use of numerical derivatives.  The relative step size\n"<>
	"\t\t is set to \!\(\*SuperscriptBox[\(10\), \(JacobianStep\)]\).\n"<>
	"\t RadialCFDepth\[Rule]1\n"<>
	"\t\t The initial Radial Continued Fraction Depth is usually taken from the prior \n"<>
	"\t\t solution in the sequence.  Fractional values reduce this initial vaue by that\n"<>
	"\t\t fraction.  Integer values larger than 300 are used as the new initial value.\n"<>
	"\t RadialRelax\[Rule]1\n"<>
	"\t\t Initial under-relaxation parameter used in Radial Newton iterations.\n"<>
	"\t SchDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Schwarzschile iterations.\n"<>
	"\t RadialDebug\[Rule]0 : 0,1,2,3,...\n"<>
	"\t\t 'Verbosity' level during the Radial Newton iterations.\n"<>
	"\t SchAnSol\[Rule]False\n"<>
	"\t\t Call SchwarzschildQNM[l,n,Nmax] to us Mathematica routines to obtain a \n"<>
	"\t\t potientially superior initial guess for the solution.\n\n"<>
	"SchwarzschildQNM[l,n,Nmax] computes the Quasi-Normal Mode solutions for overtone n of mode l.  "<>
	"The mode is computed with the radial continued fraction is truncated at a depth of Nmax.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t QNMPrecision\[Rule]24"


(* ::Subsection::Closed:: *)
(*Eigenvalue Solvers*)


QNMSolution::usage=
	"QNMSolution[n,s,l,m,a,\[Omega]g,Almg,\[Epsilon],relax,Nrcf,Nm,\[Omega]0,Alm0,rl,rt] finds a "<>
	"solution of the coupled radial and angular Teukolsky equations with "<>
	"spin-weight s, 'magnetic' index m, and dimensionless angular momentum a.  "<>
	"\[Omega]g and Almg are initial guesses for the frequency and separation constant, "<>
	"which are assumed to be associated with azimuthal index l and overtone "<>
	"index n.\n\n"<>
	"\tThe solution is found by iteration, calling AngularSpectralRoot and "<>
	"RadialLentzRoot until the magnitude of the change produced by either "<>
	"is less than \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  'relax' is an initial "<>
	"under-relaxation parameter used in updating the current guesses. For the "<>
	"radial equation, the continued fraction is truncated at the "<>
	"\!\(\*SuperscriptBox[\(Nrcf\),\(th\)]\) term.  For the angular equation, "<>
	"Nm sets the size of the spectral approximation matrix.\n\n"<>
	"\t\[Omega]0, Alm0, rl, and rt are used to specify solution windows centered "<>
	"around the initial guesses \[Omega]g and Almg.  \[Omega]0 and Alm0 represent the prior "<>
	"solutions in a sequence of solutions.  |\[Omega]0-\[Omega]g| and |Alm0-Almg| set a "<>
	"length scale d for each solution.  The window for each solution is a "<>
	"portion of an anulus centered on the prior solution.  The difference "<>
	"between the inner and outer rings of the anulus is 2*rl*d with the "<>
	"guess solution centered between.  The width of the wedge is 2*rt*d "<>
	"along the outer ring.  Solutions falling outside the solution window "<>
	"for either quantity are rejected.  If rl=0 or rt=0, the solution is "<>
	"always accepted.\n\n"<>
	"Options:\n"<>
	"\t SolutionDebug\[Rule]0 : Integer\n"<>
	"\t\t Verbosity of debugging output for solutions iteration.  Increasing\n"<>
	"\t\t value increases verbosity.\n"<>
	"\t RadialDebug\[Rule]0 : Integer\n"<>
	"\t\t Verbosity of debugging output.  Increasing value increases verbosity.\n"<>
	"\t RadialRelax\[Rule]1: Rational(or Integer)\n"<>
	"\t\t Under-relaxation parameter.\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t Log10 of relative step size used in numerical evaluation of the \n"<>
	"\t\t Jacobian in Newton's method."


SpinWeightedSpheroidal::usage=""


AngularSpectralRoot::usage=
	"AngularSpectralRoot[s,m,c,Alm,N] solves the N-dimensional discrete "<>
	"approximation for the spin-weighted functions of spin-weight s, "<>
	"'magnetic' index m, and oblateness parameter c.  The solution returns "<>
	"the eigenvalue closest to Alm, and the corresponding spectral coefficients."


RadialCFRem::usage=""


AngularSpectralRootForl::usage=
	"AngularSpectralRoot[s,l,m,c,N] solves the N-dimensional discrete "<>
	"approximation for the spin-weighted functions of spin-weight s, "<>
	"'magnetic' index m, and oblateness parameter c.  The solution returns "<>
	"the eigenvalue and spectral coefficients for given multipole index l."


RadialLentzRoot::usage=
	"RadialLentzRoot[n,s,m,a,Alm,\[Omega],N,\[Epsilon],Radius] solves the radial Teukolsky "<>
	"equation with spin-weight s, 'magnetic' index m, dimensionless angular "<>
	"momentum a, and separation constant Alm.  The solution is obtained by "<>
	"finding the root of the \!\(\*SuperscriptBox[\(n\), \(th\)]\) inversion "<>
	"of the associated continued fraction equation.  Newton's method is used "<>
	"and \[Omega] is taken as the initial guess.  The continued fraction is evalueated "<>
	"'bottom up' starting with the approximate remainder for the "<>
	"\!\(\*SuperscriptBox[\(N\), \(th\)]\) term.  Newton's method terminates "<>
	"when the corrections are smaller than \!\(\*SuperscriptBox[\(10\),\(\[Epsilon]\)]\).  "<>
	"Radius sets an upper limit to the magnitude of the Newton correction "<>
	"during each iteration.\n\n"<>
	"Options:\n"<>
	"\t RadialDebug\[Rule]0 : Integer\n"<>
	"\t\t Verbosity of debugging output.  Increasing value increases verbosity.\n"<>
	"\t RadialRelax\[Rule]1: Rational(or Integer)\n"<>
	"\t\t Under-relaxation parameter.\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t Log10 of relative step size used in numerical evaluation of the \n"<>
	"\t\t Jacobian in Newton's method."


RadialLentzRoot2::usage=
	"RadialLentzRoot[n,s,m,a,Alm,\[Omega],Nrcf,Nm,\[Epsilon],Radius] solves the radial Teukolsky "<>
	"equation with spin-weight s, 'magnetic' index m, dimensionless angular "<>
	"momentum a, and separation constant Alm.  The solution is obtained by "<>
	"finding the root of the \!\(\*SuperscriptBox[\(n\), \(th\)]\) inversion "<>
	"of the associated continued fraction equation.  Newton's method is used "<>
	"and \[Omega] is taken as the initial guess.  For each value of \[Omega] considered, Alm is "<>
	"recomputed using an Nm dimensional spectral representation."<>
	"The continued fraction is evalueated "<>
	"'bottom up' starting with the approximate remainder for the "<>
	"\!\(\*SuperscriptBox[\(Nrcf\), \(th\)]\) term.  Newton's method terminates "<>
	"when the corrections are smaller than \!\(\*SuperscriptBox[\(10\),\(\[Epsilon]\)]\).  "<>
	"Radius sets an upper limit to the magnitude of the Newton correction "<>
	"during each iteration.\n\n"<>
	"Options:\n"<>
	"\t RadialDebug\[Rule]0 : Integer\n"<>
	"\t\t Verbosity of debugging output.  Increasing value increases verbosity.\n"<>
	"\t RadialRelax\[Rule]1: Rational(or Integer)\n"<>
	"\t\t Under-relaxation parameter.\n"<>
	"\t JacobianStep\[Rule]-10\n"<>
	"\t\t Log10 of relative step size used in numerical evaluation of the \n"<>
	"\t\t Jacobian in Newton's method."


(* ::Subsection::Closed:: *)
(*Utility Routines*)


VerifyExpansion::usage=
	"VerifyExpansion[l,m,n] verify the expansion coefficients for the spin-weighted "<>
	"spheroidal function associated with each QNM solution in the Kerr QNM sequence with "<>
	"mode (l,m) and overtone n.  Each solution is checked to ensure that the function is "<>
	"properly normalized, and that the \"l-th\" coefficient is real.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  'n' "<>
	"can be either an Integer or an overtone multiplet index.  An overtone multiplet index "<>
	"is a 2 element list {n,mult}, where 'n' is the Integer overtone number, and 'mult' is "<>
	"in the range 0,1,...,(Nmult-1), with 'Nmult' the number of sequences with the same "<>
	"overtone index.\n\n"<>
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call."


RadialCFErrorEst::usage=""


KerrQNMMakeMultiplet::usage=
	"KerrQNMMakeMultiplet[l,m,n] converts an Integer overtone index 'n' for the "<>
	"existing Kerr QNM sequence with mode (l,m) into an overtone multiplet index.  "<>
	"By default, the sequence number is set to 0.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  An "<>
	"overtone multiplet index is a 2 element list {n,mult}, where 'n' is the Integer "<>
	"overtone number, and the sequence number 'mult' is in the range 0,1,...,(Nmult-1), "<>
	"with 'Nmult' the number of sequences with the same overtone index.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t OTmultiple\[Rule]0\n"<>
	"\t\t Sequence number to be used in creating the overtone multiplet index."


ShortenQNMSequence::usage=
	"ShortenQNMSequence[l,m,n,N] removes the last N elements of the Quasi-Normal "<>
	"Mode sequence (l,m,n)."


MergeQNMSequence::usage=""


SetSpinWeight::usage=
	"SetSpinWeight[s] sets the value of the spin-weight used in all subsequent "<>
	"computations:\n"<>
	"\t s=-2 : Gravitational perturbations\n"<>
	"\t s=-1 : Electro-Magnetic perturbations\n"<>
	"\t s= 0 : Scalar perturbations."


(* ::Subsection::Closed:: *)
(*Plotting Routines*)


QNMPlotSch::usage=
	"QNMPlotSch[l] plots the Schwarzschild QNM frequency for all computed "<>
	"overtones for mode l.  The imaginary axis is inverted.\n\n"<>
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call."


SchwarzschildOmega::usage=
	"SchwarzschildOmega[l,m,n] returns the conjucate of the Schwarzschild QNM "<>
	"frequency for mode l and overtone n."


KerrOmegaList::usage=
	"KerrOmegaList[l,m,n] returns the list of the conjugate of the Kerr QNM "<>
	"frequencies for mode (l,m) and overtone n."


KerrOmegaListS::usage=
	"KerrOmegaListS[l,m,n] returns the short list of the conjugate of the Kerr QNM "<>
	"frequencies for mode (l,m) and overtone n.  Only frequencies where the dimensionless "<>
	"angular momentum a is a multiple of 0.05 are included, along with the frequency "<>
	"for the largest a if 1-a \[LessEqual] \!\(\*SuperscriptBox[\(10\), \(-6\)]\)."


KerrAList::usage=
	"KerrAList[l,m,n] returns the list of the Kerr QNM separation constants for "<>
	"mode (l,m) and overtone n."


KerrAListS::usage=
	"KerrAListS[l,m,n] returns the short list of the Kerr QNM separation constants for "<>
	"mode (l,m) and overtone n.  Only separation constants where the dimensionless "<>
	"angular momentum a is a multiple of 0.05 are included, along with the constant "<>
	"for the largest a if 1-a \[LessEqual] \!\(\*SuperscriptBox[\(10\), \(-6\)]\)."


QNMPlotOmega::usage=
	"QNMPlotOmega[l,m,n] plots the Kerr QNM frequency for mode (l,m)"<>
	"and overtone n.  The imaginary axis is inverted.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  'n' "<>
	"can be either an Integer or an overtone multiplet index.  An overtone multiplet index "<>
	"is a 2 element list {n,mult}, where 'n' is the Integer overtone number, and 'mult' is "<>
	"in the range 0,1,...,(Nmult-1), with 'Nmult' the number of sequences with the same "<>
	"overtone index.\n\n"<>
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n\n\n"<>
	"QNMPlotOmega[l,n] plots Kerr QNM frequency for all modes (l,m) with |m|\[LessEqual]l "<>
	"and overtone n.  The imaginary axis is inverted.\n\n"<>
	"Overtone Multiplets: See above.  In this case n must be an Integer.\n\n"
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t OTmultiple\[Rule]{}\n"<>
	"\t\t List of overtone multiplets.  An overtone multiplet is a List {m,Nmult},\n"<>
	"\t\t where -l\[LessEqual]m\[LessEqual]l is an Integer and mode (l,m) with overtone n is an Overtone \n"<>
	"\t\t Multiplet."


QNMPlotA::usage=
	"QNMPlotA[l,m,n] plots the Kerr QNM separation constant for mode (l,m) "<>
	"and overtone n.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  'n' "<>
	"can be either an Integer or an overtone multiplet index.  An overtone multiplet index "<>
	"is a 2 element list {n,mult}, where 'n' is the Integer overtone number, and 'mult' is "<>
	"in the range 0,1,...,(Nmult-1), with 'Nmult' the number of sequences with the same "<>
	"overtone index.\n\n"<>
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n\n\n"<>
	"QNMPlotOmega[l,n] plots Kerr QNM separation constant for all modes (l,m) with |m|\[LessEqual]l "<>
	"and overtone n.\n\n"<>
	"Overtone Multiplets: See above.  In this case n must be an Integer.\n\n"
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t OTmultiple\[Rule]{}\n"<>
	"\t\t List of overtone multiplets.  An overtone multiplet is a List {m,Nmult},\n"<>
	"\t\t where -l\[LessEqual]m\[LessEqual]l is an Integer and mode (l,m) with overtone n is an Overtone \n"<>
	"\t\t Multiplet."


QNMPlotOmegaTones::usage=
	"QNMPlotOmegaTones[l,m] plots the Kerr QNM frequency for all modes (l,m) "<>
	"and overtone n.  The imaginary axis is inverted.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.\n\n"<>
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t OTmultiple\[Rule]{}\n"<>
	"\t\t List of overtone multiplets.  An overtone multiplet is a List {n,Nmult},\n"<>
	"\t\t where 'n' is the overtone index,and 'Nmult' is the number of sequences\n"<>
	"\t\t with the same overtone index."


QNMPlotATones::usage=
	"QNMPlotATones[l,m] plots the Kerr QNM separation constant for all modes (l,m) "<>
	"and overtone n.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.\n\n"<>
	"Options:  All options for ListPlot are allowed.\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t OTmultiple\[Rule]{}\n"<>
	"\t\t List of overtone multiplets.  An overtone multiplet is a List {n,Nmult},\n"<>
	"\t\t where 'n' is the overtone index,and 'Nmult' is the number of sequences\n"<>
	"\t\t with the same overtone index."


QNMPlotAccumulation\[Omega]::usage=
	"QNMPlotAccumulation\[Omega][l,m,overtones,Nv,a0] separately fits and plots the real and "<>
	"imaginary parts of the Kerr QNM frequency versus 1-a, where 'a' is the dimensionless "<>
	"angular momentum.  The fitting functions are based on data from the (l,m) modes "<>
	"with overtone values given in the List 'overtones'.  Nv and a0 are used to restrict the "<>
	"number of points in each (l,m,n) sequence that are used to determine the fit.  Only "<>
	"the last Nv points in each sequence are used, and only if the value of 'a' is greater "<>
	"than 1-a0.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  The "<>
	"elements of 'overtones' can be either an Integer or an overtone multiplet index.  An "<>
	"overtone multiplet index is a 2 element list {n,mult}, where 'n' is the Integer "<>
	"overtone number, and 'mult' is in the range 0,1,...,(Nmult-1), with 'Nmult' the number "<>
	"of sequences with the same overtone index.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t OTskip\[Rule]{}\n"<>
	"\t\t List of overtone values or multiplet indices that do not approach the \n"<>
	"\t\t limit point.\n"<>
	"\t OTmultiple\[Rule]{}\n"<>
	"\t\t List of overtone multiplets.  An overtone multiplet is a List {n,Nmult},\n"<>
	"\t\t where 'n' is the overtone index,and 'Nmult' is the number of sequences\n"<>
	"\t\t with the same overtone index."


QNMPlotAccumulationAlm::usage=
	"QNMPlotAccumulationAlm[l,m,overtones,Nv,a0] separately fits and plots the real and "<>
	"imaginary parts of the Kerr QNM separation constant versus 1-a, where 'a' is the "<>
	"dimensionless angular momentum.  The fitting functions are based on data from the "<>
	"(l,m) modes with overtone values given in the List 'overtones'.  Nv and a0 are used "<>
	"to restrict the number of points in each (l,m,n) sequence that are used to determine "<>
	"the fit.  Only the last Nv points in each sequence are used, and only if the value "<>
	"of 'a' is greater than 1-a0.\n\n"<>
	"Overtone Multiplets: There are cases where more than one sequence is associated with "<>
	"the same overtone n of mode (l,m).  Such sets are called overtone multiplets.  The "<>
	"elements of 'overtones' can be either an Integer or an overtone multiplet index.  An "<>
	"overtone multiplet index is a 2 element list {n,mult}, where 'n' is the Integer "<>
	"overtone number, and 'mult' is in the range 0,1,...,(Nmult-1), with 'Nmult' the number "<>
	"of sequences with the same overtone index.\n\n"<>
	"Options:\n"<>
	"\t QNMSpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set before any KerrQNM function call.\n"<>
	"\t OTskip\[Rule]{}\n"<>
	"\t\t List of overtone values or multiplet indices that do not approach the \n"<>
	"\t\t limit point.\n"<>
	"\t OTmultiple\[Rule]{}\n"<>
	"\t\t List of overtone multiplets.  An overtone multiplet is a List {n,Nmult},\n"<>
	"\t\t where 'n' is the overtone index,and 'Nmult' is the number of sequences\n"<>
	"\t\t with the same overtone index."


QNMColor::usage=
	"QNMColor[i] returns the entry i in the default color scheme used in Mathematica \n"<>
	"Plotting.  The scheme has 4 entries which are repeated indefinitely."


QNMMark::usage=
	"QNMMark[i] returns the entry i in the default PlotMarker scheme used in Mathematica \n"<>
	"Plotting.  The scheme has 8 entries which are repeated indefinitely."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[QNMSpinWeight,Min\[CapitalDelta]alevel,Max\[CapitalDelta]alevel,Max\[CapitalDelta]\[Phi],Max\[CapitalDelta]\[Omega],\[Omega]poly,QNMaStart,QNMGuess,QNMPrecision,
		NoNeg\[Omega],SolutionWindowl,SolutionWindowt,SolutionDebug,SolutionRelax,SolutionAsym,
		SolutionSlow,SolutionOscillate,SolutionIter,SolutionSpiral,SolutionSpiralAuto,
		RadialCFMinDepth,RadialCFMaxGuess,RadialCFDepth,RadialCFDigits,RCFPower,
		RadialDebug,RadialRelax,JacobianStep,Root\[Epsilon],ImaginaryAxis,
		OTskip,OTmultiple,OTfit,SchAnSol,SchDebug,InterpOrder,NInsert,ExtrapolationOrder,
		SeqDirection,Maximala\[Epsilon],Minblevel,Maxblevel,CurvatureRatio,
		Refinement,Index,RefinementPlot,SeqLevel,RadialCFLevel,AccuracyLevel,PrecisionLevel,StepRatio,CurveRatio,
		RefinementAction,ForceRefinement,RefineAccuracy,RefinePrecision,RefineAdapt,FixAdapt,RemoveLevels,
		LimitRefinement,Minima];


Protect[\[Delta]\[Omega]r,\[Delta]\[Omega]i,\[Epsilon]1,n1,\[Alpha]1,\[Alpha]2,\[Alpha]3,\[Alpha]4,\[Alpha]5,\[Beta]1,\[Beta]2,\[Beta]3,\[Beta]4,\[Beta]5,\[Beta]6,\[Beta]7,\[Beta]8];


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Angular Equation : Spectral Method*)


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


(*If[!QNMDebug,Protect[Alms,Blms,Clms,Dlms,Elms,Flms,Glms,Hlms]];*)


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


AngularSpectralRootForl[s_Integer,index_Integer,m_Integer,c_?NumberQ,N_Integer]:=
Module[{sol},
	sol=SpinWeightedSpheroidal[m,s,c,N];
	{sol[[1,index]],N,sol[[2,index]]}
]


If[!QNMDebug,Protect[Mat,SpinWeightedSpheroidal,AngularSpectralRoot,AngularSpectralRootForl]];


(* ::Section::Closed:: *)
(*Radial Equation : Modified Leaver' s Method*)


(* ::Subsection::Closed:: *)
(*Define the Radial Equation recurrence relation coefficients*)


rp=1+Sqrt[1-a^2];rm=1-Sqrt[1-a^2];\[Sigma]p=(2\[Omega] rp-m a)/(rp-rm);


C0=1-s-2I \[Sigma]p;C1=2(1+s-2I \[Omega]);C2=2I(rp-rm)\[Omega];C3=(1-4I \[Omega])(1+s-2I \[Sigma]p);C4=(4rp-a^2)\[Omega]^2 + 2((rp-rm)\[Sigma]p+I rp-2I(1+s))\[Omega]-Alm;


D0=Simplify[C0];D1=Simplify[-2C0-C1+C2];D2=Simplify[C0+C1];D3=Simplify[-C3+C4];D4=C3;


\[Alpha]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D0+1)n+D0;
\[Beta]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=-2n^2+(D1+2)n+D3;
\[Gamma]r[n_Integer|n_Plus|n_Symbol,s_Integer,m_Integer,a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ]=n^2+(D2-3)n+D4-D2+2;


Clear[rp,rm,\[Sigma]p,C0,C1,C2,C3,C4,D0,D1,D2,D3,D4];
If[!QNMDebug,Protect[\[Alpha]r,\[Beta]r,\[Gamma]r]];


(* ::Subsection::Closed:: *)
(*Newton' s Method for finding roots of Radial Equation*)


Options[RadialLentzStep]={RadialDebug->0}


RadialLentzStep[n_Integer,s_Integer,m_Integer,
				a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ,
				\[Omega]step_Real|\[Omega]step_Rational|\[Omega]step_Integer,
				Nrcf_Integer,\[Epsilon]_Integer,OptionsPattern[]]:= 
Module[{sol,solpr,solpi,\[CapitalDelta]\[Omega]r,\[CapitalDelta]\[Omega]i,RedRd\[Omega]r,RedRd\[Omega]i,ImdRd\[Omega]r,ImdRd\[Omega]i,\[Delta]\[Omega],WorkPrec,
		radialdebug=OptionValue[RadialDebug]},
	(* Alternative for Lentz's method
		sol=RadialCFLentz[n,s,m,a,Alm,\[Omega],10^\[Epsilon]];
	*)
	sol=RadialCFRem[n,s,m,a,Alm,\[Omega],Nrcf];
	If[Abs[sol[[1]]]==0,
		Return[{Solve[{\[Delta]\[Omega]r==0,\[Delta]\[Omega]i==0},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol}];
	]; (* Avoid error case *)
	\[CapitalDelta]\[Omega]r=Re[\[Omega]]\[Omega]step;\[CapitalDelta]\[Omega]i=Im[\[Omega]]\[Omega]step;
	(* Alternative for Lentz's method
		solpr=RadialCFLentz[n,s,m,a,Alm,Re[\[Omega]](1+\[Omega]step)+I Im[\[Omega]],10^\[Epsilon]];
		solpi=RadialCFLentz[n,s,m,a,Alm,Re[\[Omega]]+I Im[\[Omega]](1+\[Omega]step),10^\[Epsilon]];
	*)
	solpr=RadialCFRem[n,s,m,a,Alm,Re[\[Omega]](1+\[Omega]step)+I Im[\[Omega]],Nrcf];
	solpi=RadialCFRem[n,s,m,a,Alm,Re[\[Omega]]+I Im[\[Omega]](1+\[Omega]step),Nrcf];
	RedRd\[Omega]r=Re[solpr[[1]]-sol[[1]]]/\[CapitalDelta]\[Omega]r;
	RedRd\[Omega]i=Re[solpi[[1]]-sol[[1]]]/\[CapitalDelta]\[Omega]i;
	ImdRd\[Omega]r=Im[solpr[[1]]-sol[[1]]]/\[CapitalDelta]\[Omega]r;
	ImdRd\[Omega]i=Im[solpi[[1]]-sol[[1]]]/\[CapitalDelta]\[Omega]i;
	If[radialdebug>4,
		Print[RedRd\[Omega]r,"\[Delta]\[Omega]r + ",RedRd\[Omega]i,"\[Delta]\[Omega]i==",-Re[sol[[1]]]];
		Print[ImdRd\[Omega]r,"\[Delta]\[Omega]r + ",ImdRd\[Omega]i,"\[Delta]\[Omega]i==",-Im[sol[[1]]]]
	];
	If[Accuracy[RedRd\[Omega]r]<=0,RedRd\[Omega]r=0];
	If[Accuracy[RedRd\[Omega]i]<=0,RedRd\[Omega]i=0];
	If[Accuracy[ImdRd\[Omega]r]<=0,ImdRd\[Omega]r=0];
	If[Accuracy[ImdRd\[Omega]i]<=0,ImdRd\[Omega]i=0];
	If[RedRd\[Omega]r==0 && RedRd\[Omega]i==0,
		Return[{Solve[{\[Delta]\[Omega]r==2 10^(\[Epsilon]+2),\[Delta]\[Omega]i==-3 10^(\[Epsilon]+2)},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol}];
	]; (* Avoid error case *)
	If[ImdRd\[Omega]r==0 && ImdRd\[Omega]i==0,
		Return[{Solve[{\[Delta]\[Omega]r==-2 10^(\[Epsilon]+2),\[Delta]\[Omega]i==3 10^(\[Epsilon]+2)},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol}];
	]; (* Avoid error case *)
	WorkPrec = 2Max[$MinPrecision,Quiet[Precision[sol[[1]]]]];
	\[Delta]\[Omega]=Solve[{RedRd\[Omega]r \[Delta]\[Omega]r + RedRd\[Omega]i \[Delta]\[Omega]i==-Re[sol[[1]]], ImdRd\[Omega]r \[Delta]\[Omega]r + ImdRd\[Omega]i \[Delta]\[Omega]i==-Im[sol[[1]]]},{\[Delta]\[Omega]r,\[Delta]\[Omega]i},WorkingPrecision->WorkPrec][[1]];
	If[radialdebug>3,Print["RadialLentzStep : ",\[Delta]\[Omega]," : ",sol]];
	If[radialdebug>5,Print["Precision - \[Delta]\[Omega] : ",Precision[\[Delta]\[Omega][[1]]]," " ,Precision[\[Delta]\[Omega][[2]]]]];
	{\[Delta]\[Omega],sol}
]


Options[RadialLentzStep2D]={RadialDebug->0}


RadialLentzStep2D[n_Integer,s_Integer,m_Integer,
				a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ,
				\[Omega]step_Real|\[Omega]step_Rational|\[Omega]step_Integer,
				Nrcf_Integer,Nm_Integer,\[Epsilon]_Integer,OptionsPattern[]]:= 
Module[{\[Omega]r,\[Omega]i,Almr,sol,solpr,\[CapitalDelta]sol,pcount,pflag,\[CapitalDelta]\[Omega]r,\[CapitalDelta]\[Omega]i,RedRd\[Omega]r,RedRd\[Omega]i,ImdRd\[Omega]r,ImdRd\[Omega]i,\[Delta]\[Omega],WorkPrec,
		radialdebug=OptionValue[RadialDebug]},
	sol=RadialCFRem[n,s,m,a,Alm,\[Omega],Nrcf];
	\[CapitalDelta]\[Omega]r=Re[\[Omega]]\[Omega]step;\[CapitalDelta]\[Omega]i=Im[\[Omega]]\[Omega]step;
	\[Omega]r=Re[\[Omega]](1+\[Omega]step)+I Im[\[Omega]];
	\[Omega]i=Re[\[Omega]]+I Im[\[Omega]](1+\[Omega]step);
	Almr=AngularSpectralRoot[s,m,a*\[Omega]r,Alm,Nm][[1]];
	pcount=0;pflag=False;
	While[True,
		solpr=RadialCFRem[n,s,m,a,Almr,\[Omega]r,Nrcf];
		\[CapitalDelta]sol=solpr[[1]]-sol[[1]];
		Off[Precision::mnprec,Accuracy::mnprec];
		If[MyPrecision[\[CapitalDelta]sol]==Precision[\[CapitalDelta]sol],
			On[Precision::mnprec,Accuracy::mnprec];Break[],
			$MinPrecision+=4,
			$MinPrecision+=4
		];
		On[Precision::mnprec,Accuracy::mnprec];
		pflag=True;
		If[++pcount==10,Print["Excessive increase in Precision, Abort"];Abort[]];
		sol=RadialCFRem[n,s,m,a,Alm,\[Omega],Nrcf];
	];
	If[pflag,Print["Set $MinPrecision = ",$MinPrecision]];
	RedRd\[Omega]r=Re[\[CapitalDelta]sol]/\[CapitalDelta]\[Omega]r;
	ImdRd\[Omega]r=Im[\[CapitalDelta]sol]/\[CapitalDelta]\[Omega]r;
	ImdRd\[Omega]i=RedRd\[Omega]r;
	RedRd\[Omega]i=-ImdRd\[Omega]r;
	If[Abs[sol[[1]]]==0,
		Return[{Solve[{\[Delta]\[Omega]r==0,\[Delta]\[Omega]i==0},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol,{{RedRd\[Omega]r,RedRd\[Omega]i},{ImdRd\[Omega]r,ImdRd\[Omega]i}}}];
	]; (* Avoid error case *)
	If[radialdebug>4,
		Print[RedRd\[Omega]r,"\[Delta]\[Omega]r + ",RedRd\[Omega]i,"\[Delta]\[Omega]i==",-Re[sol[[1]]]];
		Print[ImdRd\[Omega]r,"\[Delta]\[Omega]r + ",ImdRd\[Omega]i,"\[Delta]\[Omega]i==",-Im[sol[[1]]]]
	];
	If[Accuracy[RedRd\[Omega]r]<=0,RedRd\[Omega]r=0];
	If[Accuracy[RedRd\[Omega]i]<=0,RedRd\[Omega]i=0];
	If[Accuracy[ImdRd\[Omega]r]<=0,ImdRd\[Omega]r=0];
	If[Accuracy[ImdRd\[Omega]i]<=0,ImdRd\[Omega]i=0];
	If[RedRd\[Omega]r==0 && RedRd\[Omega]i==0,
		Return[{Solve[{\[Delta]\[Omega]r==2 10^(\[Epsilon]+2),\[Delta]\[Omega]i==-3 10^(\[Epsilon]+2)},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol,Null}];
	]; (* Avoid error case *)
	If[ImdRd\[Omega]r==0 && ImdRd\[Omega]i==0,
		Return[{Solve[{\[Delta]\[Omega]r==-2 10^(\[Epsilon]+2),\[Delta]\[Omega]i==3 10^(\[Epsilon]+2)},{\[Delta]\[Omega]r,\[Delta]\[Omega]i}][[1]],sol,Null}];
	]; (* Avoid error case *)
	WorkPrec = 2Max[$MinPrecision,Quiet[Precision[sol[[1]]]]];
	\[Delta]\[Omega]=Solve[{RedRd\[Omega]r \[Delta]\[Omega]r + RedRd\[Omega]i \[Delta]\[Omega]i==-Re[sol[[1]]], ImdRd\[Omega]r \[Delta]\[Omega]r + ImdRd\[Omega]i \[Delta]\[Omega]i==-Im[sol[[1]]]},{\[Delta]\[Omega]r,\[Delta]\[Omega]i},WorkingPrecision->WorkPrec][[1]];
	If[radialdebug>3,Print["RadialLentzStep : ",\[Delta]\[Omega]," : ",sol]];
	If[radialdebug>5,Print["Precision - \[Delta]\[Omega] : ",Precision[\[Delta]\[Omega][[1]]]," " ,Precision[\[Delta]\[Omega][[2]]]]];
	{\[Delta]\[Omega],sol,{{RedRd\[Omega]r,RedRd\[Omega]i},{ImdRd\[Omega]r,ImdRd\[Omega]i}}}
]


Options[RadialLentzRoot]=Union[Options[RadialLentzStep],{RadialRelax->1,JacobianStep->-10}];


RadialLentzRoot[n_Integer,s_Integer,m_Integer,
				a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ,
				Nrcf_Integer,\[Epsilon]_Integer,
				Radius_Real|Radius_Rational|Radius_Integer,
				opts:OptionsPattern[]]:= 
Module[{sol1,\[Omega]root=\[Omega],\[Delta]\[Omega]1,\[Delta]\[Omega]2,Ninv,iteration=0,slow=False,convcount=0,
		\[Omega]step=10^(OptionValue[JacobianStep]),
		radialdebug=OptionValue[RadialDebug],
		radialrelax=Rationalize[OptionValue[RadialRelax]]},
	sol1=RadialLentzStep[n,s,m,a,Alm,\[Omega],\[Omega]step,Nrcf,\[Epsilon],FilterRules[{opts},Options[RadialLentzStep]]];
	\[Delta]\[Omega]1 = sol1[[1,1,2]]+I sol1[[1,2,2]];
	If[Not[NumberQ[\[Delta]\[Omega]1]],Print["RadialLentzStep failed, returning ",sol1];Abort[]];
	Ninv=n;
	If[radialdebug>0,Print["\[Delta]\[Omega]= ",\[Delta]\[Omega]1," root= ",Abs[sol1[[2,1]]]]];
	While[Abs[\[Delta]\[Omega]1]>10^\[Epsilon],
		If[++iteration > 40,
			Return[{{False,slow},{\[Omega]root,Ninv,sol1[[2,2]],\[Epsilon],Abs[\[Delta]\[Omega]1]}}]
		];
		\[Delta]\[Omega]2=\[Delta]\[Omega]1;
		\[Delta]\[Omega]1=If[Abs[\[Delta]\[Omega]1]>Radius,Radius \[Delta]\[Omega]1/Abs[\[Delta]\[Omega]1],\[Delta]\[Omega]1];
		\[Omega]root+= radialrelax \[Delta]\[Omega]1; (* NOTE: UNDER RELAXATION *)
		If[radialdebug>1,Print["\[Omega]= ",\[Omega]root]];
		sol1=RadialLentzStep[Ninv,s,m,a,Alm,\[Omega]root,\[Omega]step,Nrcf,\[Epsilon],FilterRules[{opts},Options[RadialLentzStep]]];
		\[Delta]\[Omega]1 = sol1[[1,1,2]]+I sol1[[1,2,2]];
		If[Not[NumberQ[\[Delta]\[Omega]1]],Print["RadialLentzStep failed, returning ",sol1];Abort[]];
		If[radialdebug>0,Print["\[Delta]\[Omega]= ",\[Delta]\[Omega]1," root= ",Abs[sol1[[2,1]]]]];
		If[Abs[\[Delta]\[Omega]1]/Abs[\[Delta]\[Omega]2]>1/2,
		If[++convcount>5,slow=True],convcount=0,slow=False];
		If[radialdebug>2,Print["Conv.Rate: ",Abs[\[Delta]\[Omega]2]/Abs[\[Delta]\[Omega]1]," : ",slow]]; 
	];
	\[Delta]\[Omega]1=If[Abs[\[Delta]\[Omega]1]>Radius,Radius \[Delta]\[Omega]1/Abs[\[Delta]\[Omega]1],\[Delta]\[Omega]1];
	\[Omega]root+= radialrelax \[Delta]\[Omega]1; (* NOTE: UNDER RELAXATION *)
	{{True,slow},{\[Omega]root,Ninv,Nrcf,\[Epsilon],Abs[\[Delta]\[Omega]1]}}
]


Options[RadialLentzRoot2]=Union[Options[RadialLentzStep2D],
							{RadialRelax->1,JacobianStep->-10,Root\[Epsilon]->Null[]}];


RadialLentzRoot2[n_Integer,s_Integer,m_Integer,
				a_Rational|a_Integer,Alm_?NumberQ,\[Omega]_?NumberQ,
				Nrcf_Integer,Nm_Integer,\[Epsilon]_Integer,
				Radius_Real|Radius_Rational|Radius_Integer,
				opts:OptionsPattern[]]:= 
Module[{sol1,\[Omega]root=\[Omega],\[Delta]\[Omega]1,\[Delta]\[Omega]2,Almc,Ninv,iteration=0,slow=False,convcount=0,
		refliptotal=0,reflipcount=0,lastrealsign=0,
		largeroot=0,largerootratio=0,largerootcount=0,largeroottotal=0,
		jacobianmatrix=Null,
		\[Omega]step=10^(OptionValue[JacobianStep]),\[Epsilon]root=OptionValue[Root\[Epsilon]],
		radialdebug=OptionValue[RadialDebug],
		radialrelax=Rationalize[OptionValue[RadialRelax]]},
	If[Not[NumberQ[\[Epsilon]root]],\[Epsilon]root=\[Epsilon]];
	lastrealsign=Sign[Re[\[Omega]]];
	Almc=AngularSpectralRoot[s,m,a*\[Omega]root,Alm,Nm][[1]];
	sol1=RadialLentzStep2D[n,s,m,a,Almc,\[Omega]root,\[Omega]step,Nrcf,Nm,\[Epsilon],FilterRules[{opts},Options[RadialLentzStep2D]]];
	\[Delta]\[Omega]1 = sol1[[1,1,2]]+I sol1[[1,2,2]];
	If[Not[NumberQ[\[Delta]\[Omega]1]],Print["RadialLentzStep failed, returning ",sol1];Abort[]];
	Ninv=n;
	If[radialdebug>0,Print["\[Delta]\[Omega]= ",\[Delta]\[Omega]1," root= ",Abs[sol1[[2,1]]]]];
	jacobianmatrix=sol1[[3]];
	While[Abs[\[Delta]\[Omega]1]>10^\[Epsilon] || Abs[sol1[[2,1]]]>10^\[Epsilon]root,
		If[++iteration > 40,
			Return[{{False,slow,jacobianmatrix},{\[Omega]root,Ninv,sol1[[2,2]],\[Epsilon],Abs[\[Delta]\[Omega]1]}}]
		];
		\[Delta]\[Omega]2=\[Delta]\[Omega]1;
		\[Delta]\[Omega]1=If[Abs[\[Delta]\[Omega]1]>Radius,Radius \[Delta]\[Omega]1/Abs[\[Delta]\[Omega]1],\[Delta]\[Omega]1];
		\[Omega]root+= radialrelax \[Delta]\[Omega]1; (* NOTE: UNDER RELAXATION *)
		If[Abs[sol1[[2,1]]]>10,
			largerootratio=largeroot;
			largeroot=Abs[sol1[[2,1]]];
			largerootratio/=largeroot;
			++largeroottotal;
			If[1/4<largerootratio<4,
				If[++largerootcount>5, (* root remains large and doesn't change much *)
					Return[{{False,True,jacobianmatrix}}];
				],
				If[largeroottotal>20, (* root is large too often *)
					Return[{{False,True,jacobianmatrix}}];
				];
				largerootcount=0;
			],
			largeroot=0;
			largerootcount=0;
		];
		If[Sign[Re[\[Omega]root]]!=lastrealsign,
			lastrealsign*=-1;
			++refliptotal;
			If[++reflipcount>5, (* real part of \[Omega] is repeatedly flipping sign *)
				Return[{{False,True,jacobianmatrix}}];
			],
			If[refliptotal>20,  (* real part of \[Omega] is flipping sign intermitently *)
				Return[{{False,True,jacobianmatrix}}];
			];
			reflipcount=0;
		];
		If[radialdebug>1,Print["\[Omega]= ",\[Omega]root]];
		Almc=AngularSpectralRoot[s,m,a*\[Omega]root,Alm,Nm][[1]];
		sol1=RadialLentzStep2D[Ninv,s,m,a,Almc,\[Omega]root,\[Omega]step,Nrcf,Nm,\[Epsilon],FilterRules[{opts},Options[RadialLentzStep2D]]];
		\[Delta]\[Omega]1 = sol1[[1,1,2]]+I sol1[[1,2,2]];
		If[Not[NumberQ[\[Delta]\[Omega]1]],Print["RadialLentzStep failed, returning ",sol1];Abort[]];
		If[radialdebug>0,Print["\[Delta]\[Omega]= ",\[Delta]\[Omega]1," root= ",Abs[sol1[[2,1]]]]];
		If[Abs[\[Delta]\[Omega]1]/Abs[\[Delta]\[Omega]2]>1/2,
		If[++convcount>5,slow=True],convcount=0,slow=False];
		If[radialdebug>2,Print["Conv.Rate: ",Abs[\[Delta]\[Omega]2]/Abs[\[Delta]\[Omega]1]," : ",slow]]; 
		If[Head[sol1[[3]]]==List,jacobianmatrix=sol1[[3]]];
	];
	\[Delta]\[Omega]1=If[Abs[\[Delta]\[Omega]1]>Radius,Radius \[Delta]\[Omega]1/Abs[\[Delta]\[Omega]1],\[Delta]\[Omega]1];
	\[Omega]root+= radialrelax \[Delta]\[Omega]1; (* NOTE: UNDER RELAXATION *)
	{{True,slow,jacobianmatrix},{\[Omega]root,Ninv,Nrcf,\[Epsilon],Abs[\[Delta]\[Omega]1]}}
]


If[!QNMDebug,Protect[RadialLentzStep,RadialLentzStep2D,RadialLentzRoot,RadialLentzRoot2]];


(* ::Subsection::Closed:: *)
(*Evaluate nth inversion of the Radial Equation' s continued fraction equation*)


(* ::Subsubsection::Closed:: *)
(*Approximation of remainder at the Nmax element of the continued fraction*)


RadialCFRemainder[s_Integer,m_Integer,a_Rational|a_Integer,
				  Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:=
Module[{C12tmp,C1,C2,C3,C4,C5,Rem,Err},
	C12tmp=-4I Sqrt[1-a^2]\[Omega];
	C1 = Sqrt[C12tmp];
	If[Re[C1]<0,C1=-C1]; (* Note: we are finding -Rem, so we want C1>0. Correspond to u1<0 *)
	C2=(3-4s+8I(\[Omega]+Sqrt[1-a^2]\[Omega]))/4;
	C3=(3+16Alm+16s+16s^2-32a m \[Omega]-128\[Omega]^2+80a^2\[Omega]^2-128Sqrt[1-a^2]\[Omega]^2-32I(-5Sqrt[1-a^2]\[Omega]+4Sqrt[1-a^2]s \[Omega]))/(32C1);
	C4=-(32\[Omega](9+4Alm+2s(-5+4s)+4a \[Omega](-2m+a \[Omega]))+(I(3+16Alm+96a m \[Omega]+16(s^2+(32-51a^2+32Sqrt[1-a^2])\[Omega]^2+s(1+16a \[Omega](-m+2a \[Omega])))))/Sqrt[1-a^2])/(256\[Omega]);
	C5=-(1/(8192(-1+a^2)C1 \[Omega]))I(-256Sqrt[1-a^2]Alm^2+Sqrt[1-a^2](63+288s+32s^2-512s^3-256s^4)+64I(21+100s+48s^2-64s^3+I a Sqrt[1-a^2]m(9+48s+112s^2)+a^2(-21-100s-48s^2+64s^3))\[Omega]+32(320I a m(-3+4s)-320I a^3m(-3+4s)-8(-3+225Sqrt[1-a^2]-16(1+19Sqrt[1-a^2])s+16(-1+3Sqrt[1-a^2])s^2)+a^2(3(-8+591Sqrt[1-a^2])+224Sqrt[1-a^2]m^2-16(8+133Sqrt[1-a^2])s+16(-8+59Sqrt[1-a^2])s^2))\[Omega]^2-1024I(-8I a(1+Sqrt[1-a^2])m-I a^3(-8+15Sqrt[1-a^2])m-5a^4(-7+4s)+8(1+Sqrt[1-a^2])(3+4s)-a^2(59+24Sqrt[1-a^2]+4(3+8Sqrt[1-a^2])s))\[Omega]^3+256(-128(1+Sqrt[1-a^2])+a^4(-80+7Sqrt[1-a^2])+16a^2(13+9Sqrt[1-a^2]))\[Omega]^4-32Alm(Sqrt[1-a^2](-9+16s+16s^2)-32I(7+3I a Sqrt[1-a^2]m-4s+a^2(-7+4s))\[Omega]-16(8(1+Sqrt[1-a^2])+a^2(-8+11Sqrt[1-a^2]))\[Omega]^2));
	Rem=-1+C1/Sqrt[Nmax]+C2/Nmax+C3/(Nmax^(3/2))+C4/Nmax^2+C5/(Nmax^(5/2));
	Err=Abs[(C5/(Nmax^(5/2)))/Rem];
	{Rem,Err} (* returns the negative of the remainder term *)
]


(* ::Subsubsection::Closed:: *)
(*Simple "bottom-up" evaluation at the Nmax element with no remainder approximation*)


RadialCF[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
		 Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:= 
Module[{i,Ri},
	If[n>0,
		For[{i=0;Ri=\[Beta]r[0,s,m,a,Alm,\[Omega]]},i<n,++i,
			Ri=\[Beta]r[i+1,s,m,a,Alm,\[Omega]]-\[Alpha]r[i,s,m,a,Alm,\[Omega]]\[Gamma]r[i+1,s,m,a,Alm,\[Omega]]/Ri
		]; Ri,
		\[Beta]r[0,s,m,a,Alm,\[Omega]],
		\[Omega]
	]
	-If[Nmax>0,
		For[{i=Nmax;Ri=\[Alpha]r[Nmax,s,m,a,Alm,\[Omega]];},i>n,{Ri=\[Alpha]r[i-1,s,m,a,Alm,\[Omega]]\[Gamma]r[i,s,m,a,Alm,\[Omega]]/(\[Beta]r[i,s,m,a,Alm,\[Omega]]-Ri);--i;}
		];
		If[Nmax>=n,Ri,0],
		0
	]
]


(* ::Subsubsection::Closed:: *)
(*"Lentz's Method" evaluation of continued fraction to accuracy of 10^\[Epsilon]*)


RadialCFLentz[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			  Alm_?NumberQ,\[Omega]_?NumberQ,\[Epsilon]_Integer]:=
Module[{f,i,Ri,C,D,\[CapitalDelta],\[Alpha],\[Beta],j},
	f=If[n>0,
		For[{i=0;Ri=\[Beta]r[0,s,m,a,Alm,\[Omega]]},i<n,++i,
			\[Beta]=\[Beta]r[i+1,s,m,a,Alm,\[Omega]];
			Ri=\[Beta]-\[Alpha]r[i,s,m,a,Alm,\[Omega]]\[Gamma]r[i+1,s,m,a,Alm,\[Omega]]/Ri
		];Ri,
		\[Beta]r[0,s,m,a,Alm,\[Omega]]
	];
	C=f;D=0;\[CapitalDelta]=10\[Epsilon];i=n;
	j=0;
	(*Print["Precision - C:",Precision[C]," D:",Precision[D]];
	Print["j=",j," C=",C," D=",D," \[CapitalDelta]=",\[CapitalDelta]];*)
	While[(Abs[Abs[\[CapitalDelta]]-1]>=\[Epsilon])||(++j<300),
		(*j++;*)
		\[Beta]=\[Beta]r[++i,s,m,a,Alm,\[Omega]];\[Alpha]=-\[Alpha]r[i-1,s,m,a,Alm,\[Omega]]\[Gamma]r[i,s,m,a,Alm,\[Omega]];
		D=\[Beta]+\[Alpha] D; If[Abs[D]==0,D=10^(-100);Print["Lentz D=0 : ",i]];
		C=\[Beta]+\[Alpha]/C; If[Abs[C]==0,C=10^(-100);Print["Lentz C=0 : ",i]];
		(*Print["Precision - C:",Precision[C]," D:",Precision[D]];*)
		D=1/D;
		\[CapitalDelta]=C D;
		f = f \[CapitalDelta];
		(*If[Precision[f]<20,
		If[Mod[j, 10]==0,Print["j : ",j," \[CapitalDelta] : ",\[CapitalDelta]," : C= ",C," : D= ",D];
		Print["Precision - C:",Precision[C]," D:",Precision[D]];
		Print["Precision - f:",Precision[f]," \[CapitalDelta]:",Precision[\[CapitalDelta]]]];
		];*)
		(*f[Mod[j,1000]\[Equal]0,Abort[]];*)
	];
	(* Print["f = ",f," (",Precision[f],") i= ",i]; *)
	{f,i}
]


(* ::Subsubsection::Closed:: *)
(*"Bottom-up" evaluation at the Nmax element with remainder approximation*)


RadialCFRem[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:= 
Module[{Rem,i,func,t},
	If[n>Nmax,Print["inversion greater than CF depth"];Abort[]];
	Rem=If[Nmax>0,RadialCFRemainder[s,m,a,Alm,\[Omega],Nmax][[1]],-1];
(*Print["Start RadialCFRem, Remainder = ",Rem];*)
	func=Simplify[{\[Alpha]r[i,s,m,a,Alm,\[Omega]],\[Beta]r[i,s,m,a,Alm,\[Omega]],\[Gamma]r[i,s,m,a,Alm,\[Omega]]}];
	t=Table[func,{i,0,Nmax}];
	{Evaluate[t[[n+1,2]]+t[[n+1,3]]Fold[-#2[[1]]/(#2[[2]]+#2[[3]]#1)&,0,Take[t,n]]]+
		t[[n+1,1]]Fold[-#2[[3]]/(#2[[2]]+#2[[1]]#1)&,-Rem,Reverse[Drop[t,n+1]]],Nmax}
]


RadialCFRem2[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:= 
Module[{Rem,i,func,t},
	If[n>Nmax,Print["inversion greater than CF depth"];Abort[]];
	Rem=If[Nmax>0,RadialCFRemainder[s,m,a,Alm,\[Omega],Nmax][[1]],-1];
(*Print["Start RadialCFRem, Remainder = ",Rem];*)
	func = Simplify[{\[Alpha]r[i,s,m,a,Alm,\[Omega]],\[Beta]r[i,s,m,a,Alm,\[Omega]]}/\[Gamma]r[i,s,m,a,Alm,\[Omega]]];
	t=Table[func,{i,0,Nmax}];
	{\[Gamma]r[n,s,m,a,Alm,\[Omega]](t[[n+1,2]]+Fold[-#2[[1]]/(#2[[2]]+#1)&,0,Take[t,n]]+
		t[[n+1,1]]Fold[-1/(#2[[2]]+#2[[1]]#1)&,-Rem,Reverse[Drop[t,n+1]]]),Nmax}
]


RadialCFRem1[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
			Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:= 
Module[{Rem,i,\[CapitalGamma]1,\[CapitalGamma]2,\[CapitalGamma],cf,const,func,Ri},
	If[n>Nmax,Print["inversion greater than CF depth"];Abort[]];
	Rem=If[Nmax>0,RadialCFRemainder[s,m,a,Alm,\[Omega],Nmax][[1]],-1];
(*Print["Start RadialCFRem, Remainder = ",Rem];*)
	If[Nmax<2000000,
		func:=Evaluate[Simplify[-\[Alpha]r[n+2-i,s,m,a,Alm,\[Omega]]/\[Gamma]r[n+2-i,s,m,a,Alm,\[Omega]]]];
		\[CapitalGamma]1=Join[{1},Table[func,{i,2,n+1}]];
(*Print["\[CapitalGamma]1 : ",\[CapitalGamma]1];*)
		(*\[CapitalGamma]1=Join[{1},Table[-Global`\[Alpha][n+2-i]/Global`\[Gamma][n+2-i],{i,2,n+1}]];*)
		\[CapitalGamma]2[i_Integer]:=\[CapitalGamma]2[i]=\[CapitalGamma]1[[i]]/\[CapitalGamma]2[i-1];\[CapitalGamma]2[1]=1;
		\[CapitalGamma]=Join[{1},Table[\[CapitalGamma]1[[i]]/\[CapitalGamma]2[i-1],{i,2,n+1}]];
(*Print["\[CapitalGamma] : ",\[CapitalGamma]];*)
		Clear[\[CapitalGamma]1,\[CapitalGamma]2];
		func:=Evaluate[Simplify[\[Beta]r[n+1-i,s,m,a,Alm,\[Omega]]/\[Alpha]r[n+1-i,s,m,a,Alm,\[Omega]]]];
		cf=Join[{Piecewise[{{(func/.i->1),n<Nmax}},(func/.i->n+1-Nmax)-Rem]},Table[func \[CapitalGamma][[i]],{i,2,n+1}]];
		(*cf=Join[{Piecewise[{{Global`\[Beta][n]/Global`\[Alpha][n],n<Nmax}},Global`\[Beta][Nmax]/Global`\[Alpha][Nmax]-Global`r[Nmax]]},Table[Global`\[Beta][n+1-i]/Global`\[Alpha][n+1-i]\[CapitalGamma][[i]],{i,2,n+1}]];*)
		Clear[\[CapitalGamma]];
(*Print["cf1 : ",cf];*)
		const=FromContinuedFraction[cf];
		Clear[cf];
		func:=Evaluate[Simplify[-\[Alpha]r[n+i-1,s,m,a,Alm,\[Omega]]/\[Gamma]r[n+i-1,s,m,a,Alm,\[Omega]]]];
		\[CapitalGamma]1=Join[{1},Table[func,{i,2,Nmax-n+1}]];
		(*\[CapitalGamma]1=Join[{1},Table[-Global`\[Alpha][n+i-1]/Global`\[Gamma][n+i-1],{i,2,Nmax-n+1}]];*)
		\[CapitalGamma]2[i_Integer]:=\[CapitalGamma]2[i]=\[CapitalGamma]1[[i]]/\[CapitalGamma]2[i-1];\[CapitalGamma]2[1]=1;
		\[CapitalGamma]=Join[{1},Table[\[CapitalGamma]1[[i]]/\[CapitalGamma]2[i-1],{i,2,Nmax-n+1}]];
		Clear[\[CapitalGamma]1,\[CapitalGamma]2];
		func:=Evaluate[Simplify[\[Beta]r[n+i-1,s,m,a,Alm,\[Omega]]/\[Alpha]r[n+i-1,s,m,a,Alm,\[Omega]]]];
		cf=Join[{const},
				Table[func \[CapitalGamma][[i]],{i,2,Nmax-n}],
				Piecewise[{{{((func/.i->Nmax-n+1)-Rem)\[CapitalGamma][[Nmax-n+1]]},n<Nmax}},{}]];
		(*cf=Join[{const},Table[Global`\[Beta][n+i-1]/Global`\[Alpha][n+i-1]\[CapitalGamma][[i]],{i,2,Nmax-n}],Piecewise[{{{(Global`\[Beta][Nmax]/Global`\[Alpha][Nmax]-Global`r[Nmax])\[CapitalGamma][[Nmax-n+1]]},n<Nmax}},{}]];*)
		Clear[\[CapitalGamma]];
(*Print["cf2 : ",cf];*)
		{FromContinuedFraction[cf]\[Alpha]r[n,s,m,a,Alm,\[Omega]],Nmax},
		(* Use "original" method for large Nmax *)
		{If[n>0,
			For[{i=0;Ri=\[Beta]r[0,s,m,a,Alm,\[Omega]]},i<n,++i,
				Ri=\[Beta]r[i+1,s,m,a,Alm,\[Omega]]-\[Alpha]r[i,s,m,a,Alm,\[Omega]]\[Gamma]r[i+1,s,m,a,Alm,\[Omega]]/Ri
			];Ri,
			\[Beta]r[0,s,m,a,Alm,\[Omega]],
			\[Omega]
		]
		-If[Nmax>0,
			For[{i=Nmax;Ri=\[Alpha]r[Nmax,s,m,a,Alm,\[Omega]]Rem;},i>n,{Ri=\[Alpha]r[i-1,s,m,a,Alm,\[Omega]]\[Gamma]r[i,s,m,a,Alm,\[Omega]]/(\[Beta]r[i,s,m,a,Alm,\[Omega]]-Ri);--i;}
			];
			If[Nmax>=n,Ri,0],
			0
		],Nmax}
	]
]


(* ::Subsubsection::Closed:: *)
(*Test convergence of continued fraction approximation*)


TestRadialCFConvergence[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
						Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:= 
Module[{N2,Rem,CFval},
	N2=IntegerPart[3*Nmax/2];
	Rem=RadialCFRemainder[s,m,a,Alm,\[Omega],Nmax];
	CFval=RadialCFRem[n,s,m,a,Alm,\[Omega],Nmax][[1]];
	{Abs[CFval-RadialCFRem[n,s,m,a,Alm,\[Omega],N2][[1]]],CFval,Rem}
]


TestRadialCFConvergence2[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
						Alm_?NumberQ,\[Omega]_?NumberQ,Nmax_Integer]:= 
Module[{N1,N2,Rem,CFval,CFval1,CFval2,cfpow},
	N1=IntegerPart[2*Nmax/3];
	N2=IntegerPart[3*Nmax/2];
	Rem=RadialCFRemainder[s,m,a,Alm,\[Omega],Nmax];
	CFval=RadialCFRem[n,s,m,a,Alm,\[Omega],Nmax][[1]];
	CFval1=RadialCFRem[n,s,m,a,Alm,\[Omega],N1][[1]];
	CFval2=RadialCFRem[n,s,m,a,Alm,\[Omega],N2][[1]];
	cfpow=(Log10[Abs[CFval-CFval2]]-Log10[Abs[CFval1-CFval2]])/Log10[3/2];
	{Abs[CFval-CFval2],CFval,Rem,cfpow}
]


Options[TestRadialCFConvergence3]=Options[RadialLentzRoot2];


TestRadialCFConvergence3[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
						Alm_?NumberQ,\[Omega]_?NumberQ,Nrcf_Integer,jacobian_,\[Epsilon]_Integer,Nrcfmin_Integer,
						Nm_Integer,\[Alpha]_Real|\[Alpha]_Rational|\[Alpha]_Integer,opts:OptionsPattern[]]:= 
Module[{N1,N2,Rem,CFval,CFval1,CFval2,cfpow,newNrcf,saveNrcf,diff,diffh,diffl,sol,cfpowcut=-2},
	If[Head[jacobian]==List,Null[],Null[],
		Print["WARNING: Jacobian not set when Nradialnew needs computation!"];
		Return[{Max[Nrcfmin,Ceiling[3/2 Nrcf]],0,0,0,Null[]}]
	];
	N1=Ceiling[2*Nrcf/3];
	N2=Ceiling[3*Nrcf/2];
	Rem=RadialCFRemainder[s,m,a,Alm,\[Omega],Nrcf];
	CFval=RadialCFRem[n,s,m,a,Alm,\[Omega],Nrcf][[1]];
	CFval1=RadialCFRem[n,s,m,a,Alm,\[Omega],N1][[1]];
	CFval2=RadialCFRem[n,s,m,a,Alm,\[Omega],N2][[1]];
	diff=Abs[CFval-CFval2];
	diffh=Abs[CFval1-CFval2];
	If[diff==0,
(*Print["Debug 1: diff==0"];*)
		If[-\[Epsilon]<=IntegerPart[Accuracy[diff]-(Log10[Det[jacobian]]/2)],
			Return[{Max[Nrcfmin,Ceiling[1/2 Nrcf]],diff,CFval,Rem,Null[]}],
			Print["WARNING: \[CapitalDelta]CF=0 testing RCF Depth with Acc : 10^(-",
					Ceiling[Accuracy[diff]-(Log10[Det[jacobian]]/2)],")"];
			Return[{Nrcf,diff,CFval,Rem,Null[]}];
		];
	];
	cfpow=(Log10[diff]-Log10[diffh])/Log10[3/2];
(*Print["Debug 2: cfpow = ",cfpow];*)
	If[cfpow>-1/2,(* Untrusted Sloap *)
		Print["Debug 3: cfpow>-1/2 (diff = ",diff,", diffh = ",diffh];
		Return[{Max[Nrcfmin,Ceiling[3/2 Nrcf]],diff,CFval,Rem,Null[]}]
	];
	newNrcf=Max[Nrcfmin,Ceiling[Nrcf/10],Ceiling[Nrcf (Sqrt[Det[jacobian]]10^\[Epsilon]/diff)^(1/cfpow)]];
(*Print["Debug 4: newNrcf = ",newNrcf," jacobian : ",Sqrt[Det[jacobian]]," diff : ",diff];*)
	If[newNrcf>=2Nrcf/3 || (cfpow<cfpowcut && newNrcf>Max[Nrcfmin,Ceiling[Nrcf/10]]),
		(*Print["Debug 5"];*)
		Return[{newNrcf,diff,CFval,Rem,cfpow}]
	];
	newNrcf=Max[Nrcfmin,N1];
	sol=RadialLentzRoot2[n,s,m,a,Alm,\[Omega],newNrcf,Nm,\[Epsilon],10^(-3),
							RadialRelax->\[Alpha],FilterRules[{opts},Options[RadialLentzRoot2]]];
	If[sol[[1,1]] && !sol[[1,2]],
		If[Log10[Abs[\[Omega]-sol[[2,1]]]]<\[Epsilon],newNrcf=Max[Nrcfmin,N1],newNrcf=Nrcf,newNrcf=Nrcf],
		newNrcf=Nrcf,newNrcf=Nrcf;
	];
(*Print["Debug 9: sol = ",sol];*)
Print["Warning: resetting newNrcf = ",newNrcf];
	{newNrcf,diff,CFval,Rem,cfpow}
]


If[!QNMDebug,Protect[RadialCFRemainder,RadialCF,RadialCFLentz,RadialCFRem,
					TestRadialCFConvergence,TestRadialCFConvergence2]];


(* ::Subsection::Closed:: *)
(*Mathematica Routine for finding roots of Radial Equation*)


RadialMode[n_Integer,s_Integer,m_Integer,a_Rational|a_Integer,
		   Alm_?NumberQ,\[Omega]guess_?NumberQ,Nmax_Integer]:=
Module[{\[Omega]},
	Reap[\[Omega]/.FindRoot[RadialCF[n,s,m,a,Alm,\[Omega],Nmax]==0,{\[Omega],\[Omega]guess},
		EvaluationMonitor:>Sow[{\[Omega],RadialCF[n,s,m,a,Alm,\[Omega],Nmax]}]]]
]


If[!QNMDebug,Protect[RadialMode]];


(* ::Section:: *)
(*Kerr QNM methods*)


(* ::Subsection::Closed:: *)
(*Iterative simultaneous solution of radial & angular Teukolsky equations*)


Options[Set\[CapitalDelta]a]={Min\[CapitalDelta]alevel->1,Max\[CapitalDelta]alevel->4,\[Omega]poly->Null[],Max\[CapitalDelta]\[Phi]->0.0005,Max\[CapitalDelta]\[Omega]->0.01};


Options[QNMSolution]=Union[{SolutionDebug->0,NoNeg\[Omega]->False,RadialCFMinDepth->300,QNMPrecision->24,
							SolutionSlow->10,SolutionOscillate->10,SolutionIter->50,
							RadialCFDigits->8,RCFPower->Null[]},Options[RadialLentzRoot2]];


QNMSolution[n_Integer,s_Integer,l_Integer,m_Integer,
			a_Rational|a_Integer,\[Omega]g_?NumberQ,Almg_?NumberQ,\[Epsilon]_Integer,
			relax_Real|relax_Rational|relax_Integer,
			Nrcf_Integer,Nm_Integer,
			\[Omega]0_?NumberQ,Alm0_?NumberQ,
			rl_Real|rl_Rational|rl_Integer,
			rt_Real|rt_Rational|rt_Integer,
			opts:OptionsPattern[]]:=
Module[{c,old\[Omega],oldAlm,radialsol,angularsol,lmin,lmax,Nradial,Nmatrix,
		\[CapitalDelta]\[Omega],\[CapitalDelta]\[Omega]2,iteration=0,nrpow,rcfpow,NradFlag=False,jacobianmatrix=Null,invJacobian,
		converged=False,count,expconv,rcferr,\[Alpha],inversion,i,invcount,slowcount,
		radialfail=0,slowcount2=0,oscillate=0,\[Epsilon]2=\[Epsilon],Nradialnew,rcfpower=0,testNrcf,
		nonegfreq=OptionValue[NoNeg\[Omega]],
		slowval=OptionValue[SolutionSlow],oscval=OptionValue[SolutionOscillate],
		iterval=OptionValue[SolutionIter],jacobianstep=OptionValue[JacobianStep],
		solutiondebug=OptionValue[SolutionDebug],RCFmin=OptionValue[RadialCFMinDepth],
		precision=OptionValue[QNMPrecision],RCFdigits=OptionValue[RadialCFDigits],
		RCFpowoverride=OptionValue[RCFPower]},

	lmin = Max[Abs[m],Abs[s]];
	lmax = Max[l+Ceiling[Nm/2],lmin+Nm-1];
	Nmatrix=lmax-lmin+1;
	Nradialnew=Nradial=Nrcf;
	If[solutiondebug>1,Print["Initial Nradial : ",Nradial]];
	old\[Omega]=\[Omega]g;oldAlm=Almg;
	inversion=n;
	count=0;
	\[Alpha]=relax; (* under-relaxation parameter *)
	If[solutiondebug>0,Print["\[Omega]g = ",\[Omega]g," : Almg = ",Almg]];
	If[solutiondebug>3,Print["\[Omega]0 = ",\[Omega]0," : Alm0 = ",Alm0]];
	c=a \[Omega]g;
	angularsol=AngularSpectralRoot[s,m,c,oldAlm,Nmatrix];
	expconv=Max[Take[Abs[angularsol[[3]]],-2]];
	While[expconv >= 10^\[Epsilon]2 && ++count<5,  (* Make sure Spectral resolution is good before doing lots of work *)
		angularsol=AngularSpectralRoot[s,m,c,oldAlm,++Nmatrix];
		If[solutiondebug>4,Print["Expconv = ",expconv," : N = ",Nmatrix]];
		expconv=Max[Take[Abs[angularsol[[3]]],-2]];
	];
	count=0;
	While[Not[converged],
		If[(++iteration>iterval && \[Epsilon]2==\[Epsilon] && \[Alpha]==1) || (\[Alpha]<1 && iteration > (3/5)iterval/\[Alpha]),
			If[\[Alpha]>0.05,
				\[Alpha]*=3/5;iteration=0;
				Print["Reduce Under-relaxation parameter to \[Alpha] = ",N[\[Alpha]]],
				Print["Too Many Iterations"];
				(*Print["a=",a," \[Omega]=",old\[Omega]," Alm=",oldAlm]*);
				Return[{False}];
			];
		];
		radialsol = RadialLentzRoot2[inversion,s,m,a,oldAlm,old\[Omega],Nradial,Nmatrix,\[Epsilon]2,10^(-3),
										RadialRelax->\[Alpha],FilterRules[{opts},Options[RadialLentzRoot2]]];
		If[Not[radialsol[[1,1]]] && Length[radialsol]==1, (* Re[\[Omega]] flipping sign, solution fails *)
			Return[{False,\[Alpha],Nradialnew,{a,Null[],Null[]}}]
		];
		invcount=0;
		slowcount=0;
		If[Head[radialsol[[1,3]]]==List,jacobianmatrix=radialsol[[1,3]]];
		While[Not[radialsol[[1,1]]]&&Not[radialsol[[1,2]]], 
			(* Radial solution not converged, but convergence not slow *)
			If[++slowcount>=20,Break[]];
			(* Print["Looping with non-slow radial convergence"]; *)
			radialsol = RadialLentzRoot2[inversion,s,m,a,oldAlm,radialsol[[2,1]],Nradial,Nmatrix,\[Epsilon]2,10^(-3),
											RadialRelax->\[Alpha],FilterRules[{opts}, Options[RadialLentzRoot2]]];
			If[Not[radialsol[[1,1]]] && Length[radialsol]==1, (* Re[\[Omega]] flipping sign, solution fails *)
				Return[{False,\[Alpha],Nradialnew,{a,Null[],Null[]}}]
			];
			If[Head[radialsol[[1,3]]]==List,jacobianmatrix=radialsol[[1,3]]];
		];
		If[Not[radialsol[[1,1]]],++radialfail;(*Print["WARNING: RadialLentzRoot failed to converged"]*)];
		If[Not[NumberQ[radialsol[[2,1]]]],Print["Failure: Solution of RadialLentzRoot not a number"];Abort[]];
		c=a radialsol[[2,1]];
		angularsol=AngularSpectralRoot[s,m,c,oldAlm,Nmatrix];
		\[CapitalDelta]\[Omega]2=\[CapitalDelta]\[Omega];
		\[CapitalDelta]\[Omega]=radialsol[[2,5]];
		If[TrueQ[nonegfreq] && Re[radialsol[[2,1]]]<0,
				radialsol[[2,1]]=-Conjugate[radialsol[[2,1]]];
				angularsol[[1]]=Conjugate[angularsol[[1]]];
		];
		old\[Omega]=radialsol[[2,1]];oldAlm=angularsol[[1]];
		If[solutiondebug>2,Print["a = ",a," \[Omega] = ",old\[Omega]," alm = ",oldAlm]];
		If[solutiondebug>1||NumberQ[RCFpowoverride],Print["\[CapitalDelta]\[Omega] = ",\[CapitalDelta]\[Omega]]];
		If[\[CapitalDelta]\[Omega]>0,
			If[Abs[\[CapitalDelta]\[Omega]2]/Abs[\[CapitalDelta]\[Omega]]<1,
				If[++slowcount2>slowval,
					(* Try under-relaxation if solution is slow *)
					slowcount2=0;radialfail=0;
					\[Alpha]*=3/5;iteration=0;
					Print["Persistent slow convergence: \[Alpha]= ",N[\[Alpha]]];
					If[\[Alpha]<0.05,Return[{False}]]
				], 
				If[slowcount2 >0,slowcount2=0;
					(* Try lower precision if persisitent oscillations *)
					If[++oscillate>oscval,
						Print["Persistent oscillations: Abort"];Abort[]
						(*
						Print["Persistent oscillations: \[Epsilon]= ",\[Epsilon]2+1];
						oscillate=0;radialfail=0;
						If[\[Epsilon]2-\[Epsilon]<2,++\[Epsilon]2;iteration/=2,Return[{False}]]
						*)
					]
				]
			]
		];
		If[radialfail>10,
			(*
			radialfail=0;slowcount2=0;oscillate=0;
			Print["Persistent failure of RadialLentzRoot: \[Epsilon]= ",\[Epsilon]2+1];
			If[\[Epsilon]2-\[Epsilon]<2,++\[Epsilon]2;iteration/=2,Return[{False}]]
			*)
			Print["Persistent failure of RadialLentzRoot: Abort"];Abort[]
		];
		If[\[CapitalDelta]\[Omega]<10^\[Epsilon]2,
			If[count==0,
				rcfpower=0;
				rcferr=TestRadialCFConvergence3[inversion,s,m,a,oldAlm(1+10^(\[Epsilon]+4)),old\[Omega](1-10^(\[Epsilon]+4)),
						Nradial,jacobianmatrix,\[Epsilon]2,RCFmin,Nmatrix,\[Alpha],FilterRules[{opts},Options[TestRadialCFConvergence3]]];
				Nradialnew=rcferr[[1]];rcfpower=rcferr[[5]];
				If[solutiondebug>5,
					Print["Jacobian matrix : ",jacobianmatrix];
					Print["rcferr : ",rcferr];
				];
(*
				If[rcferr[[4]]>=0 || rcferr[[4]]<0,
					If[rcferr[[1]]>=0 && rcferr[[4]]<=-1/2 && Head[jacobianmatrix]==List,
						rcfpower=rcferr[[4]];
						Nradialnew=Max[IntegerPart[Nradial (Sqrt[Det[jacobianmatrix]]10^\[Epsilon]2/rcferr[[1]])^(1/rcfpower)],RCFmin];
						,
Print["*** rcferr[[1]]<=0 when Nradialnew needs computation!"];
(*Print["#### rcferr = ",rcferr];
Print["$$$$ jacobianmatric = ",jacobianmatrix];*)
						Nradialnew=IntegerPart[(3/2)Nradial];
						,
Print["*** Jacobian not set when Nradialnew needs computation!"];
						Nradialnew=IntegerPart[(3/2)Nradial];
					],
					Null[],
					If[rcferr[[1]]\[Equal]0,
						If[-\[Epsilon]2\[LessEqual]IntegerPart[Accuracy[rcferr[[1]]]-(Log10[Det[jacobianmatrix]]/2)],
							Nradialnew=Max[IntegerPart[Nradial/2],RCFmin];
							rcfpower=Null[],
							Print["WARNING: \[CapitalDelta]CF=0 testing RCF Depth with Acc : 10^(-",
									IntegerPart[Accuracy[rcferr[[1]]]-(Log10[Det[jacobianmatrix]]/2)],")"];
						],
						Print["Unknown error testing RCF Depth"];Abort[];
					];
				];
*)
				If[solutiondebug>1,Print["Nradial : ",Nradial," ; Nradialnew : ",Nradialnew]];
				If[Nradialnew>Nradial,
					If[Nradialnew>(11/10)Nradial,
						NradFlag=True;count=0;iteration/=2;slowcount2=0;oscillate=0;
						If[solutiondebug>4,Print["Increase Nradial to ",Nradialnew]];
					];
					Nradial=Nradialnew;
					,
					If[solutiondebug>4 && Nradialnew<Nradial,Print["Decrease Nradial to ",Nradialnew]];
				];
			];

			expconv=Max[Take[Abs[angularsol[[3]]],-2]];
			(*expconv=0;*)
			If[expconv < 10^\[Epsilon]2,
				If[!NradFlag,++count];
				If[count>1&&radialsol[[1,1]],converged=True],
				If[solutiondebug>4,Print["Expconv = ",expconv," : count = ",count," : N = ",Nmatrix+1]];
				Nmatrix+=1;count=0;iteration/=2;slowcount2=0;oscillate=0
			];
			NradFlag=False,
			count=0
		];
	];
	If[rl!=0 &&rt!=0,
		If[Not[SolutionWindow[\[Omega]0,\[Omega]g,old\[Omega],rl,rt,True]],
			Print["\[Omega] = ",old\[Omega]," outside solution window"];
			Return[{False,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix]}],angularsol}}]
		];
		If[Not[SolutionWindow[Alm0,Almg,oldAlm,rl,rt,True]],
			Print["Alm outside solution window"];
			Return[{False,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix]}],angularsol}}]
		];
	];
(*If[Nradialnew==Nrcf,Print["WARNING: Nradialnew not reset from Nrcf"]];*)
	{True,\[Alpha],Nradialnew,{a,Join[radialsol[[2]],{Nradialnew,rcfpower,Det[jacobianmatrix],$MinPrecision}],angularsol}}
]


(* ::Subsection::Closed:: *)
(*Basic sequencer for increasing spin*)


Options[KerrQNMSequence]=Union[{QNMSpinWeight->Null[],QNMaStart->0,QNMGuess->0,
								QNMPrecision->24,
								SolutionRelax->1,RadialCFDepth->1,
								SolutionWindowl->1/2,SolutionWindowt->1/3},
								Options[Set\[CapitalDelta]a],Options[QNMSolution]];


KerrQNMSequence[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,KerrSEQ,SeqStatus,context,
		QNMguess,inversion,\[Omega],Alm,\[Omega]try,Almtry,\[Omega]save,Almsave,QNMsol,a,
		Nrcf,Nm=4,order,rl=0,rt=0,save1,save2,save3,NKQNM=0,
		\[CapitalDelta]a=10^(-3),new\[CapitalDelta]a,lasta,acount,\[CapitalDelta]alevel=1,
		\[CapitalDelta]alist={10^(-3),10^(-4),10^(-5),10^(-6),10^(-7),10^(-8),10^(-9)},\[CapitalDelta]aincflag=False,
		\[CapitalDelta]aredcount=0,c1,c2,c3,ExtrapInfo={},levelcount,al,\[Epsilon]a,
		qnmastart=OptionValue[QNMaStart],guess=OptionValue[QNMGuess],
		precision=OptionValue[QNMPrecision],
		solwinl=OptionValue[SolutionWindowl],solwint=OptionValue[SolutionWindowt],
		relax,srelax=Rationalize[OptionValue[SolutionRelax]],rcfdepth=OptionValue[RadialCFDepth]},
		$MinPrecision=precision; (* Sets the minimum precision for entire calculation *)
	SpinWeightTable:=Switch[s,
						-2,Global`KerrQNM,
						-1,Global`KerrQNMe,
						 0,Global`KerrQNMs,
						 _,Print["Invalid QNMSpinWeight"];Abort[]
					];
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n]
					];
	SeqStatus=If[Head[KerrSEQ]==List,If[Length[KerrSEQ]>0,True,False,False],False,False];
	relax=srelax;
	If[SeqStatus,
		(* Sequence exists, extend *)
		NKQNM=Length[KerrSEQ];
		Print["KerrQNM[",l,",",m,",",n,"] sequence exists with ",NKQNM," entries"];
		lasta=Round[10^12KerrSEQ[[NKQNM,1]]]10^(-12);a=lasta;
		acount=1;
		inversion=If[Head[n]==Integer,n,Null[],n[[1]]];
		\[Omega]=SetPrecision[KerrSEQ[[NKQNM,2,1]],precision];Alm=SetPrecision[KerrSEQ[[NKQNM,3,1]],precision];
		Nm=KerrSEQ[[NKQNM,3,2]];
		Nrcf=If[Length[KerrSEQ[[NKQNM,2]]]>=6,KerrSEQ[[NKQNM,2,6]],KerrSEQ[[NKQNM,2,3]]];
		If[rcfdepth>300,Nrcf=IntegerPart[rcfdepth]];
		If[rcfdepth<1 && rcfdepth>0,Nrcf=Max[300,IntegerPart[Nrcf*Rationalize[rcfdepth]]]];
		levelcount=1;
		\[Omega]save={\[Omega],0,0}; Almsave={Alm,0,0};order=0;
		If[NKQNM>1,
			levelcount=1;
			\[CapitalDelta]a=KerrSEQ[[NKQNM,1]]-KerrSEQ[[NKQNM-1,1]];
			\[CapitalDelta]alevel = Position[Chop[\[CapitalDelta]alist-\[CapitalDelta]a,10^(-10)],0][[1,1]];
			\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]];
			\[Omega]save[[2]]=SetPrecision[KerrSEQ[[NKQNM-1,2,1]],precision];
			Almsave[[2]]=SetPrecision[KerrSEQ[[NKQNM-1,3,1]],precision];
			\[Omega]=2\[Omega]save[[1]]-\[Omega]save[[2]];
			Alm=2Almsave[[1]]-Almsave[[2]];
			++order;
		];
		If[NKQNM>2,
			\[CapitalDelta]a=KerrSEQ[[NKQNM,1]]-KerrSEQ[[NKQNM-1,1]];
			\[CapitalDelta]alevel = Position[Chop[\[CapitalDelta]alist-\[CapitalDelta]a,10^(-10)],0][[1,1]];
			\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]];
			ExtrapInfo=RestartExtrapolation[s,l,m,n,\[CapitalDelta]alevel,\[CapitalDelta]alist,precision];
			(* Print["Restart: ",ExtrapInfo]; *)
			If[ExtrapInfo[[\[CapitalDelta]alevel,1]]==0,
				If[\[CapitalDelta]alevel>1&&ExtrapInfo[[\[CapitalDelta]alevel-1,2]]==2,
					If[Increase\[CapitalDelta]a[a,\[CapitalDelta]alevel,ExtrapInfo,\[CapitalDelta]alist,FilterRules[{opts},Options[Increase\[CapitalDelta]a]]],
						\[CapitalDelta]a=\[CapitalDelta]alist[[--\[CapitalDelta]alevel]];
						\[CapitalDelta]aincflag=True;
						Print["*** Increasing \[CapitalDelta]a *** "];
					];
				];
			];
			(* Print["Using Level ",\[CapitalDelta]alevel]; *)
			levelcount=ExtrapInfo[[\[CapitalDelta]alevel,1]];
			order = ExtrapInfo[[\[CapitalDelta]alevel,2]];
			\[Omega]save=ExtrapInfo[[\[CapitalDelta]alevel,3]];Almsave=ExtrapInfo[[\[CapitalDelta]alevel,4]];
			If[Not[\[CapitalDelta]aincflag],ExtrapInfo=Take[ExtrapInfo,\[CapitalDelta]alevel-1]];
			If[order==0,\[Omega]=\[Omega]save[[1]];Alm=Almsave[[1]];rl=0;rt=0,
				If[order==1,
					\[Omega]=2\[Omega]save[[1]]-\[Omega]save[[2]];
					Alm=2Almsave[[1]]-Almsave[[2]];
					rl=0;rt=0,
					\[Omega]=3\[Omega]save[[1]]-3\[Omega]save[[2]]+\[Omega]save[[3]];
					Alm=3Almsave[[1]]-3Almsave[[2]]+Almsave[[3]];
					rl=solwinl;rt=solwint;
				];
			];
			If[Head[guess]==List,
				\[Omega]=guess[[1]];
				Alm=guess[[2]];
				If[Length[guess]>=3,Nrcf=guess[[3]]];
				If[Length[guess]==4,Nm=guess[[4]]];
				Print["Guesses set: ",\[Omega]," : ",Alm," : ",Nrcf," : ",Nm];
			];
			If[order>1,
				new\[CapitalDelta]a=Set\[CapitalDelta]a[a,\[CapitalDelta]alevel,\[Omega]save,\[CapitalDelta]alist,FilterRules[{opts},Options[Set\[CapitalDelta]a]]];
				If[new\[CapitalDelta]a[[1]],
					Print["*** Reducing \[CapitalDelta]a ***"];
					\[CapitalDelta]a=new\[CapitalDelta]a[[2]];\[CapitalDelta]alevel=new\[CapitalDelta]a[[3]];
					If[Not[\[CapitalDelta]aincflag],
						AppendTo[ExtrapInfo,{Mod[levelcount,10],Min[order,2],\[Omega]save,Almsave}];
						(* Print["New Level: ",ExtrapInfo]; *)
						\[CapitalDelta]aredcount+=1;
						levelcount=0;order=1;
						\[Omega]=(231\[Omega]save[[1]]-42\[Omega]save[[2]]+11\[Omega]save[[3]])/200;
						Alm=(231Almsave[[1]]-42Almsave[[2]]+11Almsave[[3]])/200,
						(* reuse info *)
						\[CapitalDelta]aincflag=False;
						levelcount=ExtrapInfo[[\[CapitalDelta]alevel,1]];
						order = ExtrapInfo[[\[CapitalDelta]alevel,2]];
						\[Omega]save=ExtrapInfo[[\[CapitalDelta]alevel,3]];Almsave=ExtrapInfo[[\[CapitalDelta]alevel,4]];
						If[order==0,\[Omega]=\[Omega]save[[1]];Alm=Almsave[[1]];rl=0;rt=0,
							If[order==1,
								\[Omega]=2\[Omega]save[[1]]-\[Omega]save[[2]];
								Alm=2Almsave[[1]]-Almsave[[2]];
								rl=0;rt=0,
								\[Omega]=3\[Omega]save[[1]]-3\[Omega]save[[2]]+\[Omega]save[[3]];
								Alm=3Almsave[[1]]-3Almsave[[2]]+Almsave[[3]];
								rl=solwinl;rt=solwint;
							];
						];
						If[Head[OptionValue[QNMGuess]]==List,
							\[Omega]=OptionValue[QNMGuess][[1]];
							Alm=OptionValue[QNMGuess][[2]];
							If[Length[guess]>=3,Nrcf=guess[[3]]];
							If[Length[guess]==4,Nm=guess[[4]]];
							Print["Guesses set: ",\[Omega]," : ",Alm," : ",Nrcf," : ",Nm];
						];
					];
				]; 
			];
		],
		(* Sequence does not exist, start it *)
		If[Head[KerrSEQ]==SpinWeightTable || Length[KerrSEQ]==0,
			Print["Starting KerrQNM[",l,",",m,",",n,"] sequence"];
			\[CapitalDelta]a=\[CapitalDelta]alist[[1]];lasta=0;
			order=0;levelcount=-1;inversion=If[Head[n]==Integer,n,Null[],n[[1]]];
			a=0;
			\[CapitalDelta]alevel=Min[OptionValue[Min\[CapitalDelta]alevel],OptionValue[Max\[CapitalDelta]alevel]];
			\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]];
			If[Head[qnmastart]==List,
			(* For cases that cannot start at a=0 *)
				acount=IntegerPart[qnmastart[[1]]/\[CapitalDelta]a];
				\[Omega]=SetPrecision[qnmastart[[2]],precision];
				Alm=SetPrecision[qnmastart[[3]],precision];
				If[Length[qnmastart]==4,Nm=qnmastart[[4]]],
				Null[],
				acount=0;
				QNMguess=If[Head[n]==Integer,
							SetPrecision[SchQNMguess[l,n],precision],
							Null[],
							SetPrecision[SchQNMguess[l,n[[1]]],precision]
							];
				\[Omega]=SetPrecision[QNMguess[[1]],precision];
				Alm = l(l+1)-s(s+1);
			];
			\[Omega]save={\[Omega],0,0}; Almsave={Alm,0,0};
			Nrcf=300;
			If[rcfdepth>300,Nrcf=IntegerPart[rcfdepth]];
			Switch[s,
				   -2,Global`KerrQNM[l,m,n]={},
				   -1,Global`KerrQNMe[l,m,n]={},
					0,Global`KerrQNMs[l,m,n]={}
					],
			Print["Error determining status of KerrQNM[",l,",",m,",",n,"] sequence: Abort"];
			Abort[],
			Print["Error determining status of KerrQNM[",l,",",m,",",n,"] sequence: Abort"];
			Abort[];
		];
	];
	While[a+\[CapitalDelta]a<1,
		(* Print["lasta = ",lasta," acount = ",acount," \[CapitalDelta]a = ",\[CapitalDelta]a]; *)
		a=lasta+acount*\[CapitalDelta]a;
		QNMsol = QNMSolution[inversion,s,l,m,a,SetPrecision[\[Omega],precision],SetPrecision[Alm,precision],\[Epsilon],relax,Nrcf,Nm,\[Omega]save[[1]],Almsave[[1]],rl,rt, FilterRules[{opts}, Options[QNMSolution]]];
		If[Not[QNMsol[[1]]],
			save1=If[Length[QNMsol]>1,QNMsol[[4,2,1]],0];
			Print["Solution failed, Restart a=",a];
			(* \[Omega]window*=2;Almwindow*=2; *)
			\[Omega]try=(12\[Omega]-2\[Omega]save[[1]])/10;
			Almtry=(12Alm-2Almsave[[1]])/10;
			(*Print["Try \[Omega]=",\[Omega]try," Alm=",Almtry];*)
			QNMsol = QNMSolution[inversion,s,l,m,a,SetPrecision[\[Omega]try,precision],SetPrecision[Almtry,precision],\[Epsilon],relax,Nrcf,Nm,\[Omega]save[[1]],Almsave[[1]],rl,rt FilterRules[{opts},Options[QNMSolution]]];
			If[Not[QNMsol[[1]]],
				save2=If[Length[QNMsol]>1,QNMsol[[4,2,1]],-1];
				Print["Solution failed, Restart a=",a];
				\[Omega]try=(8\[Omega]+2\[Omega]save[[1]])/10;
				Almtry=(8Alm+2Almsave[[1]])/10;
				(*Print["Try \[Omega]=",\[Omega]try," Alm=",Almtry];*)
				QNMsol = QNMSolution[inversion,s,l,m,a,SetPrecision[\[Omega]try,precision],SetPrecision[Almtry,precision],\[Epsilon],relax,Nrcf,Nm,\[Omega]save[[1]],Almsave[[1]],rl,rt FilterRules[{opts}, Options[QNMSolution]]];
				If[Not[QNMsol[[1]]] && \[CapitalDelta]alevel==Length[\[CapitalDelta]alist],
					save3=If[Length[QNMsol]>1,QNMsol[[4,2,1]],-2];
					Print["save1 = ",save1];
					Print["save2 = ",save2," 1-2 : ",Abs[save1-save2]];
					Print["save3 = ",save3," 1-3 : ",Abs[save1-save3]," 2-3 : ",Abs[save2-save3]];
					If[Abs[save1-save2]<10^(\[Epsilon]+2)&&Abs[save1-save3]<10^(\[Epsilon]+2)&&Abs[save2-save3]<10^(\[Epsilon]+2),
						QNMsol={True,QNMsol[[2]],QNMsol[[3]],QNMsol[[4]]};
					];
				];
			];
		];
		If[QNMsol[[1]],
			If[\[CapitalDelta]aincflag,
				ExtrapInfo=Take[ExtrapInfo,\[CapitalDelta]alevel-1];
				\[CapitalDelta]aincflag=False;
			];
			\[CapitalDelta]aredcount=0;
			levelcount+=1;
			Print["QNMsol a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
			(* Print["levelcount = ",levelcount]; *)
			relax=Min[srelax,5/3QNMsol[[2]]];
			\[Omega]save=RotateRight[\[Omega]save]; Almsave=RotateRight[Almsave];
			\[Omega]save[[1]]=QNMsol[[4,2,1]];Almsave[[1]]=QNMsol[[4,3,1]];
			If[levelcount>0 && Mod[levelcount,10]==0 && \[CapitalDelta]alevel>1,
				levelcount=0;
				For[al=\[CapitalDelta]alevel-1,al>=1,--al,
					ExtrapInfo[[al,1]]+=1;
					ExtrapInfo[[al,2]]=Min[ExtrapInfo[[al,2]]+1,2];
					ExtrapInfo[[al,3]]=RotateRight[ExtrapInfo[[al,3]]];
					ExtrapInfo[[al,4]]=RotateRight[ExtrapInfo[[al,4]]];
					ExtrapInfo[[al,3,1]]=QNMsol[[4,2,1]];
					ExtrapInfo[[al,4,1]]=QNMsol[[4,3,1]];
					If[Mod[ExtrapInfo[[al,1]],10]==0,
						ExtrapInfo[[al,1]]=0,
						Break[];
					];
				];
				(* Print["10 step: ",ExtrapInfo]; *)
				If[levelcount==0 ,
					If[Increase\[CapitalDelta]a[a,\[CapitalDelta]alevel,ExtrapInfo,\[CapitalDelta]alist,FilterRules[{opts},Options[Increase\[CapitalDelta]a]]],
						\[CapitalDelta]aincflag=True;
						Print["*** Increasing \[CapitalDelta]a ***"];AppendTo[ExtrapInfo,{Mod[levelcount,10],Min[order,2],\[Omega]save,Almsave}];
						\[CapitalDelta]a=\[CapitalDelta]alist[[--\[CapitalDelta]alevel]];
						Print["    \[CapitalDelta]a = ",\[CapitalDelta]a," \[CapitalDelta]alevel = ",\[CapitalDelta]alevel];
						lasta=a;
						acount=0;
						levelcount=ExtrapInfo[[\[CapitalDelta]alevel,1]];
						order = ExtrapInfo[[\[CapitalDelta]alevel,2]];
						\[Omega]save=ExtrapInfo[[\[CapitalDelta]alevel,3]];Almsave=ExtrapInfo[[\[CapitalDelta]alevel,4]]
						(*ExtrapInfo=Take[ExtrapInfo,\[CapitalDelta]alevel-1];*)
						(* Print["Old Level: ",ExtrapInfo] *)
					];
				];
			];
			If[order==0,
				\[Omega]=\[Omega]save[[1]];
				Alm=Almsave[[1]];
				rl=0;rt=0;
			];
			If[order==1,
				\[Omega]=2\[Omega]save[[1]]-\[Omega]save[[2]];
				Alm=2Almsave[[1]]-Almsave[[2]];
				rl=0;rt=0;
			];
			If[order>=2,
				\[Omega]=3\[Omega]save[[1]]-3\[Omega]save[[2]]+\[Omega]save[[3]];
				Alm=3Almsave[[1]]-3Almsave[[2]]+Almsave[[3]];
				rl=solwinl;rt=solwint;
				new\[CapitalDelta]a=Set\[CapitalDelta]a[a,\[CapitalDelta]alevel,\[Omega]save,\[CapitalDelta]alist,FilterRules[{opts},Options[Set\[CapitalDelta]a]]];
				If[new\[CapitalDelta]a[[1]],
					Print["*** Reducing \[CapitalDelta]a ***"];
					AppendTo[ExtrapInfo,{Mod[levelcount,10],Min[order,2],\[Omega]save,Almsave}];
					(* Print["New Level: ",ExtrapInfo]; *)
					\[CapitalDelta]aredcount+=1;
					\[CapitalDelta]a=new\[CapitalDelta]a[[2]];
					\[CapitalDelta]alevel=new\[CapitalDelta]a[[3]];
					lasta=a;
					acount=0;
					levelcount=0;
					order=0;
					\[Omega]=(231\[Omega]save[[1]]-42\[Omega]save[[2]]+11\[Omega]save[[3]])/200;
					Alm=(231Almsave[[1]]-42Almsave[[2]]+11Almsave[[3]])/200;
				];
			];
			++order;
			Nm=QNMsol[[4,3,2]];
			Nrcf=QNMsol[[3]];
			(* Nrcf=Max[300,Nrcf 4/5]; speed up solution? *)
			Switch[s,
				   -2,AppendTo[Global`KerrQNM[l,m,n], QNMsol[[4]]],
				   -1,AppendTo[Global`KerrQNMe[l,m,n], QNMsol[[4]]],
					0,AppendTo[Global`KerrQNMs[l,m,n], QNMsol[[4]]]
				  ];
			++acount,
			(* No solution found, try decreasing a *)
			Print["No solution found, try decreasing a"];
			new\[CapitalDelta]a=Decrease\[CapitalDelta]a[\[CapitalDelta]alevel,\[CapitalDelta]alist,FilterRules[{opts},Options[Decrease\[CapitalDelta]a]]];
			If[Not[\[CapitalDelta]aincflag],
				If[new\[CapitalDelta]a[[1]],
					Print["*** Reducing \[CapitalDelta]a ***"];
					AppendTo[ExtrapInfo,{Mod[levelcount,10],Min[order,2],\[Omega]save,Almsave}];
					(* Print["New Level: ",ExtrapInfo]; *)
					\[CapitalDelta]aredcount+=1;
					lasta+=(acount-1)*\[CapitalDelta]a;
					acount=1;
					levelcount=0;
					order=1;
					\[CapitalDelta]a=new\[CapitalDelta]a[[2]];
					\[CapitalDelta]alevel=new\[CapitalDelta]a[[3]];
					c1=(1+10^(-\[CapitalDelta]aredcount)) (2+10^(-\[CapitalDelta]aredcount))/2;
					c2=-10^(-\[CapitalDelta]aredcount) (2+10^(-\[CapitalDelta]aredcount));
					c3=10^(-\[CapitalDelta]aredcount) (1+10^(-\[CapitalDelta]aredcount))/2;
					\[Omega]=c1 \[Omega]save[[1]]+c2 \[Omega]save[[2]]+c3 \[Omega]save[[3]];
					Alm=c1 Almsave[[1]]+c2 Almsave[[2]]+c3 Almsave[[3]], (* caution *)
					Print["Sequence ABORTED: l=",l," m=",m," n=",n," at a=",a];
					Return[]
				],
				(* Use prior to last increase *)
				\[CapitalDelta]a=new\[CapitalDelta]a[[2]];
				\[CapitalDelta]alevel=new\[CapitalDelta]a[[3]];
				levelcount=ExtrapInfo[[\[CapitalDelta]alevel,1]];
				order = ExtrapInfo[[\[CapitalDelta]alevel,2]];
				\[Omega]save=ExtrapInfo[[\[CapitalDelta]alevel,3]];Almsave=ExtrapInfo[[\[CapitalDelta]alevel,4]];
				ExtrapInfo=Take[ExtrapInfo,\[CapitalDelta]alevel-1];
				\[CapitalDelta]aincflag=False;
				If[order==0,\[Omega]=\[Omega]save[[1]];Alm=Almsave[[1]];rl=0;rt=0,
					If[order==1,
						\[Omega]=2\[Omega]save[[1]]-\[Omega]save[[2]];
						Alm=2Almsave[[1]]-Almsave[[2]];
						rl=0;rt=0,
						\[Omega]=3\[Omega]save[[1]]-3\[Omega]save[[2]]+\[Omega]save[[3]];
						Alm=3Almsave[[1]]-3Almsave[[2]]+Almsave[[3]];
						rl=solwinl;rt=solwint;
					];
				];
			];
		];
	];
]


(* ::Subsection::Closed:: *)
(*Reverse sequencer for decreasing spin*)


Options[KerrQNMSequenceReverse]=Union[{QNMSpinWeight->Null[],QNMGuess->0,
									SolutionRelax->1,RadialCFDepth->1,QNMPrecision->24},
									Options[QNMSolution]];


KerrQNMSequenceReverse[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,refine_Symbol:False,
						opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,KerrSEQ,
		inversion,Nrcf,Nm,i,j,NKQNM,a,\[CapitalDelta]a,\[CapitalDelta]alevel=1,
		\[CapitalDelta]alist={10^(-3),10^(-4),10^(-5),10^(-6),10^(-7),10^(-8),10^(-9)},
		\[Omega]save,Almsave,order,\[Omega],Alm,sol,
		precision=OptionValue[QNMPrecision],relax,srelax=Rationalize[OptionValue[SolutionRelax]],
		rcfdepth=OptionValue[RadialCFDepth],guess=OptionValue[QNMGuess]},
	If[!(refine \[Element] Booleans),Print[refine," is not a Boolean."];Abort[]];
	$MinPrecision=precision;
	SpinWeightTable:=Switch[s,
						-2,Global`KerrQNM,
						-1,Global`KerrQNMe,
						 0,Global`KerrQNMs,
						 _,Print["Invalid QNMSpinWeight"];Abort[]
					];
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n]
					];
	relax=srelax;
	If[Head[KerrSEQ]==List,
		NKQNM=Length[KerrSEQ];
		(* Check to see if sequence starts as assumed *)
		If[NKQNM<3,Print["Sequence must have at least 3 modes."];Abort[]];
		If[(KerrSEQ[[3,1]]-KerrSEQ[[2,1]])!=(KerrSEQ[[2,1]]-KerrSEQ[[1,1]]),
			Print["Sequence must have the same step size in a between first 3 modes."];
			Abort[]
		];
		inversion=If[Head[n]==Integer,n,Null[],n[[1]]];
		Nm=KerrSEQ[[1,3,2]];
		Nrcf=If[Length[KerrSEQ[[1,2]]]>=6,KerrSEQ[[1,2,6]],KerrSEQ[[1,2,3]]];
		If[rcfdepth>300,Nrcf=IntegerPart[rcfdepth]];
		If[rcfdepth<1 && rcfdepth>0,Nrcf=Max[300,IntegerPart[Nrcf*Rationalize[rcfdepth]]]];
		a=Round[10^12KerrSEQ[[1,1]]]10^(-12);
		\[CapitalDelta]a=KerrSEQ[[2,1]]-KerrSEQ[[1,1]];
		\[CapitalDelta]alevel = Position[Chop[\[CapitalDelta]alist-\[CapitalDelta]a,10^(-12)],0][[1,1]];
		If[refine,++\[CapitalDelta]alevel];
		\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]];
		If[refine,Print["Refine \[CapitalDelta]a to ",\[CapitalDelta]a]];
		\[Omega]save={KerrSEQ[[1,2,1]],KerrSEQ[[2,2,1]],KerrSEQ[[3,2,1]]};
		Almsave={KerrSEQ[[1,3,1]],KerrSEQ[[2,3,1]],KerrSEQ[[3,3,1]]};
		\[Omega]=3\[Omega]save[[1]]-3\[Omega]save[[2]]+\[Omega]save[[3]];
		Alm=Alm=3Almsave[[1]]-3Almsave[[2]]+Almsave[[3]];
		order=2;
		If[Head[guess]==List,
			\[Omega]=guess[[1]];
			Alm=guess[[2]];
			If[Length[guess]>=3,Nrcf=guess[[3]]];
			If[Length[guess]==4,Nm=guess[[4]]];
			Print["Guesses set: ",\[Omega]," : ",Alm," : ",Nrcf," : ",Nm];
		];
		If[refine,
			If[a-\[CapitalDelta]a<=0,Print["Two refinements needed to keep a>0; Abort"];Return[]];
			\[Omega]=(231\[Omega]save[[1]]-42\[Omega]save[[2]]+11\[Omega]save[[3]])/200;
			Alm=(231Almsave[[1]]-42Almsave[[2]]+11Almsave[[3]])/200;
			order=1;
			If[Head[guess]==List,
				\[Omega]=guess[[1]];
				Alm=guess[[2]];
				If[Length[guess]>=3,Nrcf=guess[[3]]];
				If[Length[guess]==4,Nm=guess[[4]]];
				Print["Guesses set: ",\[Omega]," : ",Alm," : ",Nrcf," : ",Nm];
			];
		];
		If[a-\[CapitalDelta]a<=0,
			If[++\[CapitalDelta]alevel>Length[\[CapitalDelta]alist],Print["No more points to add."];Return[]];
			\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]];
			If[a-\[CapitalDelta]a<=0,Print["Two refinements needed to keep a>0; Abort"];Return[]];
			\[Omega]=(231\[Omega]save[[1]]-42\[Omega]save[[2]]+11\[Omega]save[[3]])/200;
			Alm=(231Almsave[[1]]-42Almsave[[2]]+11Almsave[[3]])/200;
			order=1;
			If[Head[guess]==List,
				\[Omega]=guess[[1]];
				Alm=guess[[2]];
				If[Length[guess]>=3,Nrcf=guess[[3]]];
				If[Length[guess]==4,Nm=guess[[4]]];
				Print["Guesses set: ",\[Omega]," : ",Alm," : ",Nrcf," : ",Nm];
			];
		];
		While[a-\[CapitalDelta]a>0,
			a-=\[CapitalDelta]a;
			If[refine,Print["Guess a=",SetPrecision[a,10]," \[Omega]=",\[Omega]," Alm=",Alm]];
			sol=QNMSolution[inversion,s,l,m,a,\[Omega],Alm,\[Epsilon],relax,Nrcf,Nm,0,0,0,0,FilterRules[{opts}, Options[QNMSolution]]];
			If[sol[[1]],
				Print["QNMsol a=",Block[{$MinPrecision=0},N[sol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[sol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[sol[[4,3,1]],MachinePrecision]];
				\[Omega]save=RotateRight[\[Omega]save]; Almsave=RotateRight[Almsave];
				\[Omega]save[[1]]=sol[[4,2,1]];Almsave[[1]]=sol[[4,3,1]];
				Nm=sol[[4,3,2]];
				Nrcf=sol[[3]];
				relax=Min[srelax,5/3 sol[[2]]];
				Switch[s,
					   -2,PrependTo[Global`KerrQNM[l,m,n], sol[[4]]],
					   -1,PrependTo[Global`KerrQNMe[l,m,n], sol[[4]]],
						0,PrependTo[Global`KerrQNMs[l,m,n], sol[[4]]]
					  ];
				If[order==1,
					\[Omega]=2\[Omega]save[[1]]-\[Omega]save[[2]];
					Alm=2Almsave[[1]]-Almsave[[2]];
					++order,
					\[Omega]=3\[Omega]save[[1]]-3\[Omega]save[[2]]+\[Omega]save[[3]];
					Alm=3Almsave[[1]]-3Almsave[[2]]+Almsave[[3]];
				];
				If[a-\[CapitalDelta]a<=0,
					If[++\[CapitalDelta]alevel>Length[\[CapitalDelta]alist],Return[]];
					\[Omega]=(231\[Omega]save[[1]]-42\[Omega]save[[2]]+11\[Omega]save[[3]])/200;
					Alm=(231Almsave[[1]]-42Almsave[[2]]+11Almsave[[3]])/200;
					order=1;
					\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]]
				],
				Print["Solution failed for a=",SetPrecision[a,10]];
				Return[];
			];
		];
	];
]


(* ::Subsection::Closed:: *)
(*Insertion sequencer for increasing spin*)


Options[KerrQNMInsert]=Union[{QNMSpinWeight->Null[],InterpOrder->3,NInsert->9,
									SolutionRelax->1,RadialCFDepth->1,QNMPrecision->24},
									Options[QNMSolution]];


KerrQNMInsert[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
					a0_Real|a0_Rational|a0_Integer,
					opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,KerrSEQ,
		inversion,
		Nrcf,Nm,i,j,NKQNM,a,\[CapitalDelta]a,InsertIndex,aend,
		ReFitData\[Omega],ImFitData\[Omega],ReFitDataAlm,ImFitDataAlm,
		ReFit\[Omega],ImFit\[Omega],ReFitAlm,ImFitAlm,
		\[Omega],Alm,sol,tmpSEQ={},intrporder=OptionValue[InterpOrder],nInsert=OptionValue[NInsert],
		precision=OptionValue[QNMPrecision],relax,srelax=Rationalize[OptionValue[SolutionRelax]],
		rcfdepth=OptionValue[RadialCFDepth]},
	relax=srelax;
	$MinPrecision=precision;
	inversion=If[Head[n]==Integer,n,Null[],n[[1]]];
	SpinWeightTable:=Switch[s,
						-2,Global`KerrQNM,
						-1,Global`KerrQNMe,
						 0,Global`KerrQNMs,
						 _,Print["Invalid QNMSpinWeight"];Abort[]
					];
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n]
					];
	If[Head[KerrSEQ]==List,
		NKQNM=Length[KerrSEQ];
		For[InsertIndex=1,InsertIndex<Length[KerrSEQ],++InsertIndex,
			If[KerrSEQ[[InsertIndex,1]]==a0,Break[]];
		];
		If[InsertIndex==Length[KerrSEQ],
			Print["Cannot insert at end of List"];Return[];
		];
Print["Insert after index : ",InsertIndex];
		Nm=KerrSEQ[[InsertIndex,3,2]];
		Nrcf=If[Length[KerrSEQ[[InsertIndex,2]]]>=6,KerrSEQ[[InsertIndex,2,6]],KerrSEQ[[InsertIndex,2,3]]];
		If[rcfdepth>300,Nrcf=IntegerPart[rcfdepth]];
		If[rcfdepth<1 && rcfdepth>0,Nrcf=Max[300,IntegerPart[Nrcf*Rationalize[rcfdepth]]]];
		ReFitData\[Omega]=Table[{KerrSEQ[[i,1]],Re[KerrSEQ[[i,2,1]]]},{i,1,Length[KerrSEQ]}];
		ImFitData\[Omega]=Table[{KerrSEQ[[i,1]],Im[KerrSEQ[[i,2,1]]]},{i,1,Length[KerrSEQ]}];
		ReFitDataAlm=Table[{KerrSEQ[[i,1]],Re[KerrSEQ[[i,3,1]]]},{i,1,Length[KerrSEQ]}];
		ImFitDataAlm=Table[{KerrSEQ[[i,1]],Im[KerrSEQ[[i,3,1]]]},{i,1,Length[KerrSEQ]}];
		ReFit\[Omega]=Interpolation[ReFitData\[Omega],InterpolationOrder->intrporder,Method->"Spline"];
		ImFit\[Omega]=Interpolation[ImFitData\[Omega],InterpolationOrder->intrporder,Method->"Spline"];
		ReFitAlm=Interpolation[ReFitDataAlm,InterpolationOrder->intrporder,Method->"Spline"];
		ImFitAlm=Interpolation[ImFitDataAlm,InterpolationOrder->intrporder,Method->"Spline"];

		a=Round[10^14KerrSEQ[[InsertIndex,1]]]10^(-14);
		aend=Round[10^14KerrSEQ[[InsertIndex+1,1]]]10^(-14);
		\[CapitalDelta]a=(aend-a)/(nInsert+1);
Print["Insert between a=",a,"\[Rule]",aend," by ",\[CapitalDelta]a];
		While[a+\[CapitalDelta]a<aend,
			a+=\[CapitalDelta]a;
			\[Omega]=SetPrecision[ReFit\[Omega][a]+I ImFit\[Omega][a],precision];
			Alm=SetPrecision[ReFitAlm[a]+I ImFitAlm[a],precision];
			Print["Guess a=",Block[{$MinPrecision=0},N[a,{Infinity,20}]]," \[Omega]=",\[Omega]," Alm=",Alm]; (* *)
			sol=QNMSolution[inversion,s,l,m,a,\[Omega],Alm,\[Epsilon],relax,Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
			If[sol[[1]],
				Print["QNMsol a=",Block[{$MinPrecision=0},N[sol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[sol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[sol[[4,3,1]],MachinePrecision]];
				Nm=sol[[4,3,2]];
				Nrcf=sol[[3]];
				relax=Min[srelax,5/3 sol[[2]]];
				AppendTo[tmpSEQ, sol[[4]]],
				Print["Solution failed for a=",SetPrecision[a,10]];
				Return[];
			];
		];
		KerrSEQ=Take[KerrSEQ,InsertIndex]~Join~tmpSEQ~Join~Drop[KerrSEQ,InsertIndex];
		Switch[s,
			   -2,Global`KerrQNM[l,m,n]=KerrSEQ,
			   -1,Global`KerrQNMe[l,m,n]=KerrSEQ,
				0,Global`KerrQNMs[l,m,n]=KerrSEQ
				];
	];
]


(* ::Subsection::Closed:: *)
(*Accumulation-point sequencer for increasing spin*)


Options[KerrQNMAccumulation]=Union[{QNMSpinWeight->Null[],
									OTskip->{},OTmultiple->{},
									SolutionRelax->1,RadialCFDepth->1,QNMPrecision->24},
									Options[QNMSolution]];


KerrQNMAccumulation[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
					overtones_List,Nv_Integer,
					a0_Real|a0_Rational|a0_Integer,
					opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,KerrSEQ,
		skip=OptionValue[OTskip],multiple=OptionValue[OTmultiple],
		fitlists,inversion,nmode,nf,nprime,
		Nrcf,Nm,i,j,NKQNM,a,\[CapitalDelta]a,\[CapitalDelta]alevel=1,
		\[CapitalDelta]alist={10^(-3),10^(-4),10^(-5),10^(-6),10^(-7),10^(-8)(*,10^(-9)*)},
		qnmend,qnmRe\[Omega],qnmIm\[Omega],qnmReAlm,qnmImAlm,
		ReList\[Omega]={},ImList\[Omega]={},ReListAlm={},ImListAlm={},
		ReFitData\[Omega],ImFitData\[Omega],ReFitDataAlm,ImFitDataAlm,Refit\[Omega],Imfit\[Omega],RefitAlm,ImfitAlm,
		nfit,\[Epsilon]fit,\[Omega],Alm,sol,
		precision=OptionValue[QNMPrecision],relax,srelax=Rationalize[OptionValue[SolutionRelax]],
		rcfdepth=OptionValue[RadialCFDepth]},
	relax=srelax;
	$MinPrecision=precision;
	fitlists=OvertoneLists[n,overtones,FilterRules[{opts},Options[OvertoneLists]]];
	nmode=fitlists[[1]];
	inversion=If[Head[nmode]==Integer,nmode,Null[],nmode[[1]]];
	nfit=fitlists[[2]];
	nf=fitlists[[3]];
	nprime=fitlists[[4]];
(*Print["nmode : ",nmode];
Print["inversion : ",inversion];
Print["nfit : ",nfit];
Print["nf : ",nf];
Print["nprime : ",nprime];*)
	SpinWeightTable:=Switch[s,
						-2,Global`KerrQNM,
						-1,Global`KerrQNMe,
						 0,Global`KerrQNMs,
						 _,Print["Invalid QNMSpinWeight"];Abort[]
					];
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,nmode],
					-1,Global`KerrQNMe[l,m,nmode],
					 0,Global`KerrQNMs[l,m,nmode]
					];
	If[Head[KerrSEQ]==List,
		NKQNM=Length[KerrSEQ];
		Nm=KerrSEQ[[NKQNM,3,2]];
		Nrcf=If[Length[KerrSEQ[[NKQNM,2]]]>=6,KerrSEQ[[NKQNM,2,6]],KerrSEQ[[NKQNM,2,3]]];
		If[rcfdepth>300,Nrcf=IntegerPart[rcfdepth]];
		If[rcfdepth<1 && rcfdepth>0,Nrcf=Max[300,IntegerPart[Nrcf*Rationalize[rcfdepth]]]];
		a=Round[10^12KerrSEQ[[NKQNM,1]]]10^(-12);
		\[CapitalDelta]a=KerrSEQ[[NKQNM,1]]-KerrSEQ[[NKQNM-1,1]];
		\[CapitalDelta]alevel = Position[Chop[\[CapitalDelta]alist-\[CapitalDelta]a,10^(-12)],0][[1,1]];
		\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]];
		While[a+\[CapitalDelta]a>=1,
			If[++\[CapitalDelta]alevel>Length[\[CapitalDelta]alist],Print["No more points to add."];Return[]];
			\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]];
		];
		For[i=1,i<=Length[nf],++i,
			qnmend=Take[SpinWeightTable[l,m,nf[[i]]],-Nv];
			qnmRe\[Omega]={};
			qnmIm\[Omega]={};
			qnmReAlm={};
			qnmImAlm={};
			For[j=1,j<=Nv,++j,
				If[1-qnmend[[j,1]]<=a0,
					AppendTo[qnmRe\[Omega],{1-qnmend[[j,1]],Re[qnmend[[j,2,1]]]}];
					AppendTo[qnmIm\[Omega],{1-qnmend[[j,1]],Im[qnmend[[j,2,1]]]}];
					AppendTo[qnmReAlm,{1-qnmend[[j,1]],Re[qnmend[[j,3,1]]]}];
					AppendTo[qnmImAlm,{1-qnmend[[j,1]],Im[qnmend[[j,3,1]]]}];
				];
			];
			AppendTo[ReList\[Omega],qnmRe\[Omega]];
			AppendTo[ImList\[Omega],qnmIm\[Omega]];
			AppendTo[ReListAlm,qnmReAlm];
			AppendTo[ImListAlm,qnmImAlm];
			If[skip !=0,If[nprime[[i]]>skip,nprime[[i]]-=1]];
		];
		ReFitData\[Omega]=Flatten[Table[{nprime[[i]],ReList\[Omega][[i,j,1]],ReList\[Omega][[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ReList\[Omega][[i]]]}],1];
		ImFitData\[Omega]=Flatten[Table[{nprime[[i]],ImList\[Omega][[i,j,1]],ImList\[Omega][[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ImList\[Omega][[i]]]}],1];
		ReFitDataAlm=Flatten[Table[{nprime[[i]],ReListAlm[[i,j,1]],ReListAlm[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ReListAlm[[i]]]}],1];
		ImFitDataAlm=Flatten[Table[{nprime[[i]],ImListAlm[[i,j,1]],ImListAlm[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ImListAlm[[i]]]}],1];
		If[nf[[1]]==0,
			(* 0 fit*)
			Refit\[Omega]=NonlinearModelFit[ReFitData\[Omega],m/2-\[Alpha]1 Sqrt[\[Epsilon]1/2]+(\[Alpha]2+\[Alpha]3 n1)\[Epsilon]1,{\[Alpha]1,\[Alpha]2,\[Alpha]3},{n1,\[Epsilon]1}];
			Imfit\[Omega]=NonlinearModelFit[ImFitData\[Omega],-(n1+1/2)(Sqrt[\[Epsilon]1/2]-\[Alpha]4 \[Epsilon]1),{\[Alpha]4},{n1,\[Epsilon]1}];
			RefitAlm=NonlinearModelFit[ReFitDataAlm,l(l+1)-2+\[Beta]1+\[Beta]2 Sqrt[\[Epsilon]1/2]+(\[Beta]3+\[Beta]4 n1+\[Beta]5 n1^2)\[Epsilon]1,{\[Beta]1,\[Beta]2,\[Beta]3,\[Beta]4,\[Beta]5},{n1,\[Epsilon]1}];
			ImfitAlm=NonlinearModelFit[ImFitDataAlm,(n1+1/2)(\[Beta]6 Sqrt[\[Epsilon]1/2]+\[Beta]7 \[Epsilon]1),{\[Beta]6,\[Beta]7},{n1,\[Epsilon]1}],
			(* Normal fit*)
			Refit\[Omega]=NonlinearModelFit[ReFitData\[Omega],m/2+(\[Alpha]1+\[Alpha]2 n1)\[Epsilon]1,{\[Alpha]1,\[Alpha]2},{n1,\[Epsilon]1}];
			Imfit\[Omega]=NonlinearModelFit[ImFitData\[Omega],-1(\[Alpha]3+n1+1/2)Sqrt[\[Epsilon]1/2]+(\[Alpha]4+\[Alpha]5 n1)\[Epsilon]1,{\[Alpha]3,\[Alpha]4,\[Alpha]5},{n1,\[Epsilon]1}];
			RefitAlm=NonlinearModelFit[ReFitDataAlm,l(l+1)-2+\[Beta]1+(\[Beta]2+\[Beta]3 n1+\[Beta]4 n1^2)\[Epsilon]1,{\[Beta]1,\[Beta]2,\[Beta]3,\[Beta]4},{n1,\[Epsilon]1}];
			ImfitAlm=NonlinearModelFit[ImFitDataAlm,(\[Beta]5+n1+1/2)\[Beta]6 Sqrt[\[Epsilon]1/2]+(\[Beta]7+\[Beta]8 n1)\[Epsilon]1,{\[Beta]5,\[Beta]6,\[Beta]7,\[Beta]8},{n1,\[Epsilon]1}]
		];
		(*Print["Re(\[Omega])= "Normal[Refit\[Omega]]];
		Print["Im(\[Omega])= "Normal[Imfit\[Omega]]];
		Print["Re(Alm)= "Normal[RefitAlm]];
		Print["Im(Alm)= "Normal[ImfitAlm]];
		Off[SetPrecision::"precsm"];
		Off[N::"precsm"];
		Print["Re(\[Omega]) :",Refit\[Omega]["ParameterConfidenceIntervalTable"],"  Im(\[Omega]) :",Imfit\[Omega]["ParameterConfidenceIntervalTable"]];
		Print["Re(Alm) :",RefitAlm["ParameterConfidenceIntervalTable"],"  Im(Alm) :",ImfitAlm["ParameterConfidenceIntervalTable"]];
		On[SetPrecision::"precsm"];
		On[N::"precsm"];*)
		While[a+\[CapitalDelta]a<1,
			a+=\[CapitalDelta]a;
			\[Epsilon]fit=(1-a);
			\[Omega]=SetPrecision[Refit\[Omega][nfit,\[Epsilon]fit]+I Imfit\[Omega][nfit,\[Epsilon]fit],precision];
			Alm=SetPrecision[RefitAlm[nfit,\[Epsilon]fit]+I ImfitAlm[nfit,\[Epsilon]fit],precision];
			(* Print["Guess a=",SetPrecision[a,10]," \[Omega]=",\[Omega]," Alm=",Alm]; *)
			sol=QNMSolution[inversion,s,l,m,a,\[Omega],Alm,\[Epsilon],relax,Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
			If[sol[[1]],
				Print["QNMsol a=",Block[{$MinPrecision=0},N[sol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[sol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[sol[[4,3,1]],MachinePrecision]];
				Nm=sol[[4,3,2]];
				Nrcf=sol[[3]];
				relax=Min[srelax,5/3 sol[[2]]];
				Switch[s,
					   -2,AppendTo[Global`KerrQNM[l,m,nmode], sol[[4]]],
					   -1,AppendTo[Global`KerrQNMe[l,m,nmode], sol[[4]]],
						0,AppendTo[Global`KerrQNMs[l,m,nmode], sol[[4]]]
						];
				If[a+\[CapitalDelta]a>=1,
					If[++\[CapitalDelta]alevel>Length[\[CapitalDelta]alist],Return[]];
					\[CapitalDelta]a=\[CapitalDelta]alist[[\[CapitalDelta]alevel]]
				],
				Print["Solution failed for a=",SetPrecision[a,10]];
				Return[];
			];
		];
	];
]


(* ::Subsection:: *)
(*Adaptive Bisection sequencer*)


Options[AdaptCheck3]=Union[{Minblevel->0,Maxblevel->20,CurvatureRatio->1/2,Max\[CapitalDelta]\[Omega]->0.01,ExtrapolationOrder->2},Options[QNMSolution]];


Options[KerrQNMSequenceB]=Union[{QNMSpinWeight->Null[],QNMaStart->0,QNMGuess->0,
								SeqDirection->Forward,Maximala\[Epsilon]->10,
								SolutionRelax->1,RadialCFDepth->1,
								SolutionWindowl->1/2,SolutionWindowt->1/3},
								Options[AdaptCheck3]];
Options[KerrQNMSequenceBwork]=Options[KerrQNMSequenceB];


KerrQNMSequenceB[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{QNMSavePrecision=$MinPrecision,saneopts},
	(* saneopts ensures options set via SetOptions[KerQNMSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrQNMSequenceB],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrQNMSequenceBwork[l,m,n,\[Epsilon],FilterRules[saneopts,Options[KerrQNMSequenceB]]],
				$MinPrecision=QNMSavePrecision;Abort[]];
	$MinPrecision=QNMSavePrecision;
]


KerrQNMSequenceBwork[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]max_Integer,
					opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,KerrSEQ,KerrSEQret,AC3ret,SeqStatus,context,
		QNMguess,inversion,\[Omega],Alm,\[Omega]try,Almtry,QNMsol,a,index0=0,index0p,index0m,
		\[Epsilon]=\[Epsilon]max,Nrcf,Nm=4,rl=0,rt=0,NKQNM=0,edat0,edat,ef,iv,afit,
		maxmimuma=1,blevel=0,\[CapitalDelta]a=10^(-3),\[CapitalDelta]a2,\[CapitalDelta]a3,dir=1,\[CapitalDelta]aincflag=False,\[CapitalDelta]aincstep=0,blevelsave,
		\[Omega]0,\[Omega]p,\[Omega]m,Alm0,Almp,Almm,\[Omega]w=0,Almw=0,precisionsave,precisioncount=0,
		Minb=OptionValue[Minblevel],Maxb=OptionValue[Maxblevel],
		forward=If[OptionValue[SeqDirection]==Backward,False,True,True],
		qnmastart=OptionValue[QNMaStart],guess=OptionValue[QNMGuess],
		precision=OptionValue[QNMPrecision],extraporder=OptionValue[ExtrapolationOrder],
		solwinl=OptionValue[SolutionWindowl],solwint=OptionValue[SolutionWindowt],
		relax,srelax=Rationalize[OptionValue[SolutionRelax]],RCFmin=OptionValue[RadialCFMinDepth],
		rcfdepth=OptionValue[RadialCFDepth]},

	If[precision!=$MinPrecision,Print["Set $MinPrecision = ",precision]];
	$MinPrecision=precisionsave=precision; (* Sets the minimum precision for entire calculation *)
	maxmimuma=If[OptionValue[Maximala\[Epsilon]] \[Element] Integers,
				1-(2^-OptionValue[Maximala\[Epsilon]])/1000,
				If[OptionValue[Maximala\[Epsilon]] \[Element] Booleans && OptionValue[Maximala\[Epsilon]]==False,
					1,
					Print["Invalide value for Maximala\[Epsilon]"];Abort[]
				]];
	SpinWeightTable:=Switch[s,
						-2,Global`KerrQNM,
						-1,Global`KerrQNMe,
						 0,Global`KerrQNMs,
						 _,Print["Invalid QNMSpinWeight"];Abort[]
					];
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n]
					];
	SeqStatus=If[Head[KerrSEQ]==List,If[Length[KerrSEQ]>0,True,False,False],False,False];
	dir=If[forward,1,-1];
	relax=srelax;
	blevel=Minb;
	If[SeqStatus,
(*Print["Untested section of code! 0"];*)
		(* Sequence exists, extend *)
		NKQNM=Length[KerrSEQ];
		Print["KerrQNM[",l,",",m,",",n,"] sequence exists with ",NKQNM," entries"];
		If[forward,
			\[Epsilon]=Min[\[Epsilon]max,KerrSEQ[[NKQNM,2,4]]];
			If[Length[KerrSEQ[[NKQNM,2]]]>=9,$MinPrecision=KerrSEQ[[NKQNM,2,9]],Print["Set $MinPrecision = ",$MinPrecision]],
			\[Epsilon]=Min[\[Epsilon]max,KerrSEQ[[1,2,4]]];
			If[Length[KerrSEQ[[1,2]]]>=9,$MinPrecision=KerrSEQ[[1,2,9]],Print["Set $MinPrecision = ",$MinPrecision]]
		];
		a=If[forward,KerrSEQ[[NKQNM,1]],KerrSEQ[[1,1]]];
		inversion=If[Head[n]==Integer,n,Null[],n[[1]]];
		If[forward,
			\[Omega]=SetPrecision[KerrSEQ[[NKQNM,2,1]],Max[precision,$MinPrecision]];
			Alm=SetPrecision[KerrSEQ[[NKQNM,3,1]],Max[precision,$MinPrecision]],
			\[Omega]=SetPrecision[KerrSEQ[[1,2,1]],Max[precision,$MinPrecision]];
			Alm=SetPrecision[KerrSEQ[[1,3,1]],Max[precision,$MinPrecision]]
		];
		Nm=If[forward,KerrSEQ[[NKQNM,3,2]],KerrSEQ[[1,3,2]]];
		Nrcf=If[forward,
			If[Length[KerrSEQ[[NKQNM,2]]]>=6,KerrSEQ[[NKQNM,2,6]],KerrSEQ[[NKQNM,2,3]]],
			If[Length[KerrSEQ[[1,2]]]>=6,KerrSEQ[[1,2,6]],KerrSEQ[[1,2,3]]]
		];
		If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
		If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
		Nrcf=Max[Nrcf,RCFmin];
		If[NKQNM>1,
			If[Head[guess]==List,
Print["Untested section of code! 1"];
				\[Omega]=guess[[1]];
				Alm=guess[[2]];
				If[Length[guess]>=3,Nrcf=guess[[3]]];
				If[Length[guess]==4,Nm=guess[[4]]];
				Print["Guesses set: ",\[Omega]," : ",Alm," : ",Nrcf," : ",Nm];
			];
			If[NKQNM==2,
Print["Untested section of code! 2"];
				If[forward,
					\[Omega]=2KerrSEQ[[2,2,1]]-KerrSEQ[[1,2,1]];Alm=2KerrSEQ[[2,3,1]]-KerrSEQ[[1,3,1]],
					\[Omega]=2KerrSEQ[[1,2,1]]-KerrSEQ[[2,2,1]];Alm=2KerrSEQ[[1,3,1]]-KerrSEQ[[2,3,1]]
				]
			];
			If[NKQNM>2,
(*Print["Untested section of code! 3"];*)
				index0=If[forward,NKQNM-1,2];
				\[CapitalDelta]a=If[forward,KerrSEQ[[NKQNM,1]]-KerrSEQ[[index0,1]],KerrSEQ[[index0,1]]-KerrSEQ[[1,1]]];
				blevelsave=blevel=Round[-(3+Log10[\[CapitalDelta]a])/Log10[2]];
				\[CapitalDelta]a2=If[forward,KerrSEQ[[index0,1]]-KerrSEQ[[index0-1,1]],KerrSEQ[[index0+1,1]]-KerrSEQ[[index0,1]]];

				If[\[CapitalDelta]a > 2\[CapitalDelta]a2,Print["Sequence has unusual stepsizes, Aborting"];Abort[]];
				\[CapitalDelta]a3=2^(-Maxb)/1000;
				If[a+dir*\[CapitalDelta]a3>maxmimuma,Return[]];
				If[\[CapitalDelta]a>=\[CapitalDelta]a2,
					If[\[CapitalDelta]a==2\[CapitalDelta]a2,
						{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
						AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,forward,True,False,FilterRules[{opts},Options[AdaptCheck3]]],
						{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
						AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,forward,False,False,FilterRules[{opts},Options[AdaptCheck3]]]
					];
					Switch[s,
						   -2,Global`KerrQNM[l,m,n]=KerrSEQret,
						   -1,Global`KerrQNMe[l,m,n]=KerrSEQret,
							0,Global`KerrQNMs[l,m,n]=KerrSEQret
						  ];
				];
				(*If[\[CapitalDelta]a==2\[CapitalDelta]a2,\[CapitalDelta]aincflag=True];
				{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
					AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,forward,\[CapitalDelta]aincflag,False,FilterRules[{opts},Options[AdaptCheck3]]];
				Switch[s,
					   -2,Global`KerrQNM[l,m,n]=KerrSEQret,
					   -1,Global`KerrQNMe[l,m,n]=KerrSEQret,
						0,Global`KerrQNMs[l,m,n]=KerrSEQret
					  ];*)
				NKQNM=Length[KerrSEQ];
				index0=If[forward,NKQNM-1,2];
				\[CapitalDelta]a=2^(-blevel)/1000;
				If[blevel>blevelsave,Print["Decreasing \[CapitalDelta]a, blevel = ",blevel]];
				If[blevel<blevelsave,Print["Increasing \[CapitalDelta]a, blevel = ",blevel]];
				index0p=If[forward,index0+1,index0+\[CapitalDelta]aincstep];
				index0m=If[forward,index0-\[CapitalDelta]aincstep,index0-1];
				\[Omega]0=SetPrecision[KerrSEQ[[index0,2,1]],Max[precision,$MinPrecision]];
				\[Omega]p=SetPrecision[KerrSEQ[[index0p,2,1]],Max[precision,$MinPrecision]];
				\[Omega]m=SetPrecision[KerrSEQ[[index0m,2,1]],Max[precision,$MinPrecision]];
				Alm0=SetPrecision[KerrSEQ[[index0,3,1]],Max[precision,$MinPrecision]];
				Almp=SetPrecision[KerrSEQ[[index0p,3,1]],Max[precision,$MinPrecision]];
				Almm=SetPrecision[KerrSEQ[[index0m,3,1]],Max[precision,$MinPrecision]];
				\[Omega]w=If[forward,\[Omega]p,\[Omega]m];
				Almw=If[forward,Almp,Almm];
				rl=solwinl;rt=solwint;
				If[\[CapitalDelta]aincflag,
(*Print["Untested section of code! 4"];*)
					(* Do not reset \[CapitalDelta]aincflag, it is needed for next call to AdaptCheck3 *)
					If[forward,
						\[Omega]=6\[Omega]p-8\[Omega]0+3\[Omega]m;Alm=6Almp-8Alm0+3Almm,
						\[Omega]=6\[Omega]m-8\[Omega]0+3\[Omega]p;Alm=6Almm-8Alm0+3Almp
					],
(*Print["Untested section of code! 4a"];*)
					If[forward,
						\[Omega]=3(\[Omega]p-\[Omega]0)+\[Omega]m;Alm=3(Almp-Alm0)+Almm,
						\[Omega]=3(\[Omega]m-\[Omega]0)+\[Omega]p;Alm=3(Almm-Alm0)+Almp
					];
					If[extraporder==Accumulate,
						If[!forward,
							Print["Cannot use Accumulation extrapolation with backward sequencing"];
							Abort[]
						];
						edat0=SetPrecision[Take[KerrSEQ,-10],Max[precision,$MinPrecision]];
						edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
						\[Omega]=afit[1-(KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a)];
						edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						\[Omega]+=I afit[1-(KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a)];
						,Null[],
						If[extraporder>2,
							edat0=Take[KerrSEQ,If[forward,-(extraporder+1),extraporder+1]];
							edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
							ef[iv_]=InterpolatingPolynomial[edat,iv];
							\[Omega]=ef[If[forward,KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]];
						];
					];
				];
			];
		],
		(* Sequence does not exist, start it *)
		If[Head[KerrSEQ]==SpinWeightTable || Length[KerrSEQ]==0,
(*Print["Untested section of code! 5"];*)
			Print["Starting KerrQNM[",l,",",m,",",n,"] sequence"];
			a=0;
			inversion=If[Head[n]==Integer,n,Null[],n[[1]]];
			If[Head[qnmastart]==List,
(*Print["Untested section of code! 6"];*)
			(* For cases that cannot start at a=0 *)
				a=qnmastart[[1]];
				\[Omega]=SetPrecision[qnmastart[[2]],Max[precision,$MinPrecision]];
				Alm=SetPrecision[qnmastart[[3]],Max[precision,$MinPrecision]];
				If[Length[qnmastart]==4,Nm=qnmastart[[4]]],
				Null[],
(*Print["Untested section of code! 6a"];*)
				QNMguess=If[Head[n]==Integer,
							SetPrecision[SchQNMguess[l,n],Max[precision,$MinPrecision]],
							Null[],
							SetPrecision[SchQNMguess[l,n[[1]]],Max[precision,$MinPrecision]]
							];
				\[Omega]=SetPrecision[QNMguess[[1]],Max[precision,$MinPrecision]];
				Alm = l(l+1)-s(s+1);
			];
			a=a-dir*\[CapitalDelta]a; (* offset to "previous" a *)
			Nrcf=RCFmin;
			If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
(*Print["Starting with ",\[Omega]," : ",Alm," | a : ",a+dir*\[CapitalDelta]a," | Nrcf : ",Nrcf];*)
			Switch[s,
				   -2,Global`KerrQNM[l,m,n]={},
				   -1,Global`KerrQNMe[l,m,n]={},
					0,Global`KerrQNMs[l,m,n]={}
					],
			Print["Error determining status of KerrQNM[",l,",",m,",",n,"] sequence: Abort"];
			Abort[],
			Print["Error determining status of KerrQNM[",l,",",m,",",n,"] sequence: Abort"];
			Abort[];
		];
	];
	While[0 <= a+dir*\[CapitalDelta]a <= maxmimuma,
		(* Try lowering precision every so often *)
		If[$MinPrecision<precisionsave,precisionsave=$MinPrecision;precisioncount=0];
		If[$MinPrecision>precisionsave,
			precisioncount=0;precisionsave=$MinPrecision,
			If[$MinPrecision==precisionsave && ++precisioncount>10,
				$MinPrecision=Max[precision,precisionsave-=4];precisioncount=0]
		];
		(* Print["lasta = ",Block[{$MinPrecision=0},N[a,{Infinity,20}]]," \[CapitalDelta]a = ",Block[{$MinPrecision=0},N[dir*\[CapitalDelta]a,{Infinity,20}]]]; *)
		QNMsol = QNMSolution[inversion,s,l,m,a+dir*\[CapitalDelta]a,SetPrecision[\[Omega],Max[precision,$MinPrecision]],SetPrecision[Alm,Max[precision,$MinPrecision]],\[Epsilon],relax,Nrcf,Nm,\[Omega]w,Almw,rl,rt,FilterRules[{opts},Options[QNMSolution]]];
		If[QNMsol[[1]],
(*Print["Untested section of code! 7"];*)
			(* Solution found, save solution check sequence smoothness *)
			a=a+dir*\[CapitalDelta]a;
			Print["QNMsol a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
			(* Print["levelcount = ",levelcount]; *)
			relax=Min[srelax,5/3QNMsol[[2]]];
			Nm=QNMsol[[4,3,2]];
			Nrcf=QNMsol[[3]];
			Nrcf=Max[Nrcf,RCFmin];
			(* Nrcf=Max[RCFmin,Nrcf 4/5]; speed up solution? *)
			If[forward,
				Switch[s,
					   -2,AppendTo[Global`KerrQNM[l,m,n], QNMsol[[4]]],
					   -1,AppendTo[Global`KerrQNMe[l,m,n], QNMsol[[4]]],
						0,AppendTo[Global`KerrQNMs[l,m,n], QNMsol[[4]]]
					  ];
				,
				Switch[s,
					   -2,PrependTo[Global`KerrQNM[l,m,n], QNMsol[[4]]],
					   -1,PrependTo[Global`KerrQNMe[l,m,n], QNMsol[[4]]],
						0,PrependTo[Global`KerrQNMs[l,m,n], QNMsol[[4]]]
					  ];
			];
			NKQNM=Length[KerrSEQ];
			If[NKQNM==1,\[Omega]=KerrSEQ[[1,2,1]];Alm=KerrSEQ[[1,3,1]]];
			If[NKQNM>1,index0=If[forward,NKQNM-1,2]];
			If[NKQNM==2,
(*Print["Untested section of code! 8"];*)
				If[forward,
					\[Omega]=2KerrSEQ[[2,2,1]]-KerrSEQ[[1,2,1]];Alm=2KerrSEQ[[2,3,1]]-KerrSEQ[[1,3,1]],
					\[Omega]=2KerrSEQ[[1,2,1]]-KerrSEQ[[2,2,1]];Alm=2KerrSEQ[[1,3,1]]-KerrSEQ[[2,3,1]]
				]
			];
			If[NKQNM>2,
(*Print["Untested section of code! 9"];*)
				blevelsave=blevel;
				{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
					AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,forward,\[CapitalDelta]aincflag,False,FilterRules[{opts},Options[AdaptCheck3]]];
				Switch[s,
					   -2,Global`KerrQNM[l,m,n]=KerrSEQret,
					   -1,Global`KerrQNMe[l,m,n]=KerrSEQret,
						0,Global`KerrQNMs[l,m,n]=KerrSEQret
					  ];
				NKQNM=Length[KerrSEQ];
				index0=If[forward,NKQNM-1,2];
(*Print["Ret from AC3: blevel = ",blevel," \[CapitalDelta]aincflag = ",\[CapitalDelta]aincflag];*)
				\[CapitalDelta]a=2^(-blevel)/1000;
				index0p=If[forward,index0+1,index0+\[CapitalDelta]aincstep];
				index0m=If[forward,index0-\[CapitalDelta]aincstep,index0-1];
				\[Omega]0=KerrSEQ[[index0,2,1]];
				\[Omega]p=KerrSEQ[[index0p,2,1]];
				\[Omega]m=KerrSEQ[[index0m,2,1]];
				Alm0=KerrSEQ[[index0,3,1]];
				Almp=KerrSEQ[[index0p,3,1]];
				Almm=KerrSEQ[[index0m,3,1]];
				\[Omega]w=If[forward,\[Omega]p,\[Omega]m];
				Almw=If[forward,Almp,Almm];
				rl=solwinl;rt=solwint;
				If[\[CapitalDelta]aincflag,
(*Print["Untested section of code! 9a"];*)
					(* Do not reset \[CapitalDelta]aincflag, it is needed for next call to AdaptCheck3 *)
					Print["Increasing \[CapitalDelta]a, blevel = ",blevel];
					If[forward,
						\[Omega]=6\[Omega]p-8\[Omega]0+3\[Omega]m;Alm=6Almp-8Alm0+3Almm,
						\[Omega]=6\[Omega]m-8\[Omega]0+3\[Omega]p;Alm=6Almm-8Alm0+3Almp
					],
(*Print["Untested section of code! 9b"];*)
					If[blevel>blevelsave,Print["Decreasing \[CapitalDelta]a, blevel = ",blevel]];
					If[forward,
						\[Omega]=3(\[Omega]p-\[Omega]0)+\[Omega]m;Alm=3(Almp-Alm0)+Almm,
						\[Omega]=3(\[Omega]m-\[Omega]0)+\[Omega]p;Alm=3(Almm-Alm0)+Almp
					];
					If[extraporder==Accumulate,
						If[!forward,
							Print["Cannot use Accumulation extrapolation with backward sequencing"];
							Abort[]
						];
						edat0=SetPrecision[Take[KerrSEQ,-10],Max[precision,$MinPrecision]];
						edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
						\[Omega]=afit[1-(KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a)];
						edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
						afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
													{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
						\[Omega]+=I afit[1-(KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a)];
						,Null[],
						If[extraporder>2,
							edat0=Take[KerrSEQ,If[forward,-(extraporder+1),extraporder+1]];
							edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
							ef[iv_]=InterpolatingPolynomial[edat,iv];
							\[Omega]=ef[If[forward,KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]];
						];
					];
				];
			],
			(* No solution found, try decreasing a *)
(*Print["Untested section of code! 10"];*)
			blevelsave= ++blevel;
			If[blevel > Maxb,Return[]];
			If[NKQNM==0,Print["No solution found at a = ",Block[{$MinPrecision=0},N[a+dir*\[CapitalDelta]a,{Infinity,20}]]];Abort[]];
			Print["No solution found at a = ",Block[{$MinPrecision=0},N[a+dir*\[CapitalDelta]a,{Infinity,20}]],", try decreasing \[CapitalDelta]a, blevel = ",blevel];
			If[NKQNM>2,
				index0p=If[forward,index0+1,index0+\[CapitalDelta]aincstep];
				index0m=If[forward,index0-\[CapitalDelta]aincstep,index0-1],
Print["NKQNM <= 2 : index0 = ",index0];
				If[NKQNM==2,
					index0=If[forward,1,2];
					index0p=If[forward,2,1];
					index0m=If[forward,1,2],
					index0p=index0m=index0=1;
				];
			];
			\[Omega]0=KerrSEQ[[index0,2,1]];
			\[Omega]p=KerrSEQ[[index0p,2,1]];
			\[Omega]m=KerrSEQ[[index0m,2,1]];
			Alm0=KerrSEQ[[index0,3,1]];
			Almp=KerrSEQ[[index0p,3,1]];
			Almm=KerrSEQ[[index0m,3,1]];
			\[Omega]w=If[forward,\[Omega]p,\[Omega]m];
			Almw=If[forward,Almp,Almm];
			rl=solwinl;rt=solwint;
			If[\[CapitalDelta]aincflag,
			(* \[CapitalDelta]a increase caused failure, reset and try again *)
(*Print["Untested section of code! 10a"];*)
				\[CapitalDelta]a=2^(-blevel)/1000;
				\[CapitalDelta]aincflag=False; (* \[CapitalDelta]a caused failure, reset and try again *)
				If[forward,
					\[Omega]=3(\[Omega]p-\[Omega]0)+\[Omega]m;Alm=3(Almp-Alm0)+Almm,
					\[Omega]=3(\[Omega]m-\[Omega]0)+\[Omega]p;Alm=3(Almm-Alm0)+Almp
				],
			(* failure not due to increase in \[CapitalDelta]a, refine previous step *)
				If[forward,
					\[Omega] = (6\[Omega]0+3\[Omega]p-\[Omega]m)/8; Alm = (6Alm0+3Almp-Almm)/8,
					\[Omega] = (6\[Omega]0+3\[Omega]m-\[Omega]p)/8; Alm = (6Alm0+3Almm-Almp)/8;
				];
				(*Print["Fail call QNMSol, a= ",Block[{$MinPrecision=0},N[KerrSEQ[[index0,1]]+dir \[CapitalDelta]a/2,{Infinity,20}]]];*)
				If[forward,
					Nrcf=Max[If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]],
								If[Length[KerrSEQ[[index0p,2]]]>=6,KerrSEQ[[index0p,2,6]],KerrSEQ[[index0p,2,3]]]],
					Nrcf=Max[If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]],
								If[Length[KerrSEQ[[index0m,2]]]>=6,KerrSEQ[[index0m,2,6]],KerrSEQ[[index0m,2,3]]]]
				];
				Nrcf=Max[Nrcf,RCFmin];
				Nm=KerrSEQ[[index0,3,2]];
				QNMsol=QNMSolution[inversion,s,l,m,
									KerrSEQ[[index0,1]]+dir \[CapitalDelta]a/2,
									SetPrecision[\[Omega],Max[precision,$MinPrecision]],
									SetPrecision[Alm,Max[precision,$MinPrecision]],\[Epsilon],relax,
									Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
				If[Not[QNMsol[[1]]],Print["a+/- solution failed."];Abort[]];
(*Print["Untested section of code! 11"];*)
				Print["QNMsol+/- a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
				If[forward,index0=index0+1];
				Switch[s,
					   -2,Global`KerrQNM[l,m,n] =Insert[KerrSEQ,QNMsol[[4]],index0],
					   -1,Global`KerrQNMe[l,m,n]=Insert[KerrSEQ,QNMsol[[4]],index0],
						0,Global`KerrQNMs[l,m,n]=Insert[KerrSEQ,QNMsol[[4]],index0]
					  ];
				blevelsave=blevel;
				{KerrSEQret,blevel,\[CapitalDelta]aincflag,\[CapitalDelta]aincstep,\[Epsilon]}=
					AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,forward,False,False,Minblevel->Max[blevel,OptionValue[Minblevel]],FilterRules[{opts},Options[AdaptCheck3]]];
				Switch[s,
					   -2,Global`KerrQNM[l,m,n]=KerrSEQret,
					   -1,Global`KerrQNMe[l,m,n]=KerrSEQret,
						0,Global`KerrQNMs[l,m,n]=KerrSEQret
					  ];
				NKQNM=Length[KerrSEQ];
				index0=If[forward,NKQNM-1,2];
				\[CapitalDelta]a=2^(-blevel)/1000;
				a=KerrSEQ[[index0,1]]+dir*\[CapitalDelta]a;
				If[\[CapitalDelta]aincflag,\[CapitalDelta]aincflag=False,blevel=blevelsave]; (* Don't allow \[CapitalDelta]a to increase *)
				index0p=If[forward,index0+1,index0+\[CapitalDelta]aincstep];
				index0m=If[forward,index0-\[CapitalDelta]aincstep,index0-1];
				\[Omega]0=KerrSEQ[[index0,2,1]];
				\[Omega]p=KerrSEQ[[index0p,2,1]];
				\[Omega]m=KerrSEQ[[index0m,2,1]];
				Alm0=KerrSEQ[[index0,3,1]];
				Almp=KerrSEQ[[index0p,3,1]];
				Almm=KerrSEQ[[index0m,3,1]];
				\[Omega]w=If[forward,\[Omega]p,\[Omega]m];
				Almw=If[forward,Almp,Almm];
				rl=solwinl;rt=solwint;
				If[blevel>blevelsave,Print["Further decreasing \[CapitalDelta]a, blevel = ",blevel]];
				If[forward,
					\[Omega]=3(\[Omega]p-\[Omega]0)+\[Omega]m;Alm=3(Almp-Alm0)+Almm,
					\[Omega]=3(\[Omega]m-\[Omega]0)+\[Omega]p;Alm=3(Almm-Alm0)+Almp
				];
				If[extraporder==Accumulate,
					If[!forward,
						Print["Cannot use Accumulation extrapolation with backward sequencing"];
						Abort[]
					];
					edat0=SetPrecision[Take[KerrSEQ,-10],Max[precision,$MinPrecision]];
					edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
					afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
												{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
					(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
					\[Omega]=afit[1-(KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a)];
					edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
					afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
												{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
					\[Omega]+=I afit[1-(KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a)];
					,Null[],
					If[extraporder>2,
						edat0=Take[KerrSEQ,If[forward,-(extraporder+1),extraporder+1]];
						edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
						ef[iv_]=InterpolatingPolynomial[edat,iv];
						\[Omega]=ef[If[forward,KerrSEQ[[NKQNM,1]]+\[CapitalDelta]a,KerrSEQ[[1,1]]-\[CapitalDelta]a]];
					];
				]
			],
			(* Unknown Failure *)
				Print["Invalid call to QNMSolution"];Abort[];
		];
	];
]


AdaptCheck3[KerrTMP_List,inversion_Integer,s_Integer,l_Integer,m_Integer,\[Epsilon]max_Integer,
			relax_Real|relax_Rational|relax_Integer,
			index0_Integer,blevel_Integer,
			forward_/;forward \[Element] Booleans,
			incflag_/;incflag \[Element] Booleans,
			recursflag_/;recursflag \[Element] Booleans,
			opts:OptionsPattern[]]:=
Module[{KerrSEQ=KerrTMP,AC3ret,ind0,index0p=index0+1,index0m=index0-1,blevelp=blevel,blevelm=blevel,
		a0,\[Epsilon]=\[Epsilon]max,\[Epsilon]2,\[Epsilon]p,\[Epsilon]m,\[Epsilon]min=-13,Nrcf,Nm,\[CapitalDelta]a,\[CapitalDelta]ap,\[CapitalDelta]am,\[Omega]m,\[Omega]0,\[Omega]p,curvrat,\[CapitalDelta]\[Omega],d\[Omega],dd\[Omega],Almm,Alm0,Almp,\[Omega]g,Almg,
		QNMsol,incres=False,incstep=1,inctmp,
		edat,edat0,ef,iv,afit,extraporder=OptionValue[ExtrapolationOrder],
		Minb=OptionValue[Minblevel],Maxb=OptionValue[Maxblevel],RCFmin=OptionValue[RadialCFMinDepth],
		maxcurvrat=OptionValue[CurvatureRatio],max\[CapitalDelta]\[Omega]=Abs[OptionValue[Max\[CapitalDelta]\[Omega]]],precision=OptionValue[QNMPrecision]},
	(*If[blevel>Maxb || blevel<Minb,Print["\[CapitalDelta]a level out of bounds."];Abort[]];*)
	a0=KerrSEQ[[index0,1]];
	\[CapitalDelta]a=2^(-blevel)/1000;
	If[Not[forward] && a0==\[CapitalDelta]a,Return[{KerrSEQ,blevel,incres,incstep,\[Epsilon]}]]; (* Don't step past a=0 *)
	Nm=KerrSEQ[[index0,3,2]];
	If[incflag,
(*Print["Untested section of code! 12"];*)
		If[forward,
			While[KerrSEQ[[index0,1]]-KerrSEQ[[index0m,1]]<\[CapitalDelta]a,--index0m;++incstep],
			While[KerrSEQ[[index0p,1]]-KerrSEQ[[index0,1]]<\[CapitalDelta]a,++index0p;++incstep]
		]
	];
	If[KerrSEQ[[index0,1]]-KerrSEQ[[index0m,1]]!=\[CapitalDelta]a,Print["Incorrect \[CapitalDelta]a(-)"];Abort[]];
	If[KerrSEQ[[index0p,1]]-KerrSEQ[[index0,1]]!=\[CapitalDelta]a,Print["Incorrect \[CapitalDelta]a(+)"];Abort[]];
	\[Omega]0=KerrSEQ[[index0,2,1]];
	\[Omega]p=KerrSEQ[[index0p,2,1]];
	\[Omega]m=KerrSEQ[[index0m,2,1]];
	Alm0=KerrSEQ[[index0,3,1]];
	Almp=KerrSEQ[[index0p,3,1]];
	Almm=KerrSEQ[[index0m,3,1]];
	index0p=index0m=index0; (* reset to central index *)
	blevelp=blevelm=blevel;
	d\[Omega]=\[Omega]p-\[Omega]m; dd\[Omega]=\[Omega]p-2\[Omega]0+\[Omega]m;
	curvrat=4Sqrt[Abs[d\[Omega]]^2Abs[dd\[Omega]]^2-(Re[d\[Omega]]Re[dd\[Omega]]+Im[d\[Omega]]Im[dd\[Omega]])^2]/(Abs[d\[Omega]]^2);
	\[CapitalDelta]\[Omega]=Abs[\[Omega]0-If[forward,\[Omega]p,\[Omega]m]];
	\[Epsilon]p=\[Epsilon]m=\[Epsilon]=Max[Min[\[Epsilon]max,Floor[Log10[Abs[\[CapitalDelta]\[Omega]]]-2.5]],\[Epsilon]min];
	\[Epsilon]2=Max[Min[\[Epsilon]max,Floor[Log10[Abs[\[CapitalDelta]\[Omega]/2]]-2.5]],\[Epsilon]min];
	If[blevel < Minb || 
		(Maxb > blevel && (curvrat > maxcurvrat || \[CapitalDelta]\[Omega] > max\[CapitalDelta]\[Omega] || 
		(forward && a0+2\[CapitalDelta]a>=1) || (Not[forward] && a0-2\[CapitalDelta]a<0))),
(*Print["Untested section of code! 13"];*)
		If[incflag,incstep=1]; (* shouldn't have increased step size *)
		(* Increase resolution *)
		If[forward || (!forward && !incflag),
(*Print["Untested section of code! 14"];*)
			(* Compute solution at a+ *)
			\[Omega]g = (6\[Omega]0+3\[Omega]p-\[Omega]m)/8;
			Almg = (6Alm0+3Almp-Almm)/8;
			If[extraporder==Accumulate,
				If[!forward,
					Print["Cannot use Accumulation extrapolation with backward sequencing"];
					Abort[]
				];
				edat0=SetPrecision[Take[KerrSEQ,{index0p-9,index0p}],Max[precision,$MinPrecision]];
				edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
				afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
											{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
				(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
				\[Omega]g=afit[1-(KerrSEQ[[index0,1]]+\[CapitalDelta]a/2)];
				edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
				afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
											{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
				\[Omega]g+=I afit[1-(KerrSEQ[[index0,1]]+\[CapitalDelta]a/2)];
				Null[]
				,Null[],
				If[extraporder>2,
						edat0=Take[KerrSEQ,If[forward,{index0p-extraporder,index0p},{index0m,index0m+extraporder}]];
						edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
						ef[iv_]=InterpolatingPolynomial[edat,iv];
						\[Omega]g=ef[KerrSEQ[[index0,1]]+\[CapitalDelta]a/2];
				];
			];
			Nrcf=Max[If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]],
						If[Length[KerrSEQ[[index0p,2]]]>=6,KerrSEQ[[index0p,2,6]],KerrSEQ[[index0p,2,3]]]];
			Nrcf=Max[Nrcf,RCFmin];
			QNMsol=QNMSolution[inversion,s,l,m,a0+\[CapitalDelta]a/2,
								SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
								SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon]2,relax,
								Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
			If[QNMsol[[1]],(* valid solution *)
(*Print["Untested section of code! 15"];*)
				Print["QNMsol+ a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
				index0p=index0+1;
				blevelp=blevel+1;
				KerrSEQ=Insert[KerrSEQ,QNMsol[[4]],index0p];
				AC3ret=AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon]2,relax,index0p,blevelp,forward,False,True,FilterRules[{opts},Options[AdaptCheck3]]];
				KerrSEQ=AC3ret[[1]];
				blevelp=AC3ret[[2]];
				\[Epsilon]p=AC3ret[[5]];
				,(* invalid solution *)
				Print["a+ solution failed."];
				Abort[];
			];
		];
		If[!forward || (forward && !incflag),
(*Print["Untested section of code! 16"];*)
			(* Compute solution at a- *)
			\[Omega]g = (6\[Omega]0+3\[Omega]m-\[Omega]p)/8;
			Almg = (6Alm0+3Almm-Almp)/8;
			If[extraporder==Accumulate,
				If[!forward,
					Print["Cannot use Accumulation extrapolation with backward sequencing"];
					Abort[]
				];
				edat0=SetPrecision[Take[KerrSEQ,{index0p-9,index0p}],Max[precision,$MinPrecision]];
				edat=Table[{1-edat0[[i,1]],Re[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
				afit=NonlinearModelFit[edat,m/2+\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
											{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
				(*afit=NonlinearModelFit[edat,m/2+\[Beta] eps,{\[Beta]},eps];*)
				\[Omega]g=afit[1-(KerrSEQ[[index0,1]]-\[CapitalDelta]a/2)];
				edat=Table[{1-edat0[[i,1]],Im[edat0[[i,2,1]]]},{i,1,Length[edat0]}];
				afit=NonlinearModelFit[edat,\[Alpha] Sqrt[eps]+\[Beta] eps+\[Gamma] eps^(3/2)+\[Delta] eps^2+\[Zeta] eps^(5/2)+\[Eta] eps^3,
											{\[Alpha],\[Beta],\[Gamma],\[Delta],\[Zeta],\[Eta]},eps];
				\[Omega]g+=I afit[1-(KerrSEQ[[index0,1]]-\[CapitalDelta]a/2)];
				Null[]
				,Null[],
				If[extraporder>2,
						edat0=Take[KerrSEQ,If[forward,{index0p-extraporder,index0p},{index0m,index0m+extraporder}]];
						edat=Table[{edat0[[i,1]],edat0[[i,2,1]]},{i,1,Length[edat0]}];
						ef[iv_]=InterpolatingPolynomial[edat,iv];
						\[Omega]g=ef[KerrSEQ[[index0,1]]-\[CapitalDelta]a/2];
				];
			];
			Nrcf=Max[If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]],
						If[Length[KerrSEQ[[index0m,2]]]>=6,KerrSEQ[[index0m,2,6]],KerrSEQ[[index0m,2,3]]]];
			Nrcf=Max[Nrcf,RCFmin];
			QNMsol=QNMSolution[inversion,s,l,m,a0-\[CapitalDelta]a/2,
								SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
								SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon]2,relax,
								Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
			If[QNMsol[[1]],(* valid solution *)
(*Print["Untested section of code! 17"];*)
				Print["QNMsol- a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
				blevelm=blevel+1;
				KerrSEQ=Insert[KerrSEQ,QNMsol[[4]],index0];
				AC3ret=AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon]2,relax,index0m,blevelm,forward,False,True,FilterRules[{opts},Options[AdaptCheck3]]];
				KerrSEQ=AC3ret[[1]];
				blevelm=AC3ret[[2]];
				\[Epsilon]m=AC3ret[[5]];
				,(* invalid solution *)
				Print["a- solution failed."];
				Abort[];
			];
		];
		(* get \[CapitalDelta]a left and right of index0.  If off by more than 2x, refine once to help smooth near rapid changes *)
		ind0=index0;
		While[KerrSEQ[[ind0,1]]>a0,--ind0];
		While[KerrSEQ[[ind0,1]]<a0,++ind0];
		\[CapitalDelta]ap = KerrSEQ[[ind0+1,1]]-KerrSEQ[[ind0,1]];
		\[CapitalDelta]am = KerrSEQ[[ind0,1]]-KerrSEQ[[ind0-1,1]];
		While[\[CapitalDelta]ap>2\[CapitalDelta]am || \[CapitalDelta]am>2\[CapitalDelta]ap,
			If[\[CapitalDelta]ap>2\[CapitalDelta]am,
				\[Omega]g = (KerrSEQ[[ind0+1,2,1]]+KerrSEQ[[ind0,2,1]])/2;
				Almg = (KerrSEQ[[ind0+1,3,1]]+KerrSEQ[[ind0,3,1]])/2;
				Nrcf=Max[If[Length[KerrSEQ[[ind0,2]]]>=6,KerrSEQ[[ind0,2,6]],KerrSEQ[[ind0,2,3]]],
							If[Length[KerrSEQ[[ind0+1,2]]]>=6,KerrSEQ[[ind0+1,2,6]],KerrSEQ[[ind0+1,2,3]]]];
				Nrcf=Max[Nrcf,RCFmin];
				\[Epsilon]2=Max[Min[\[Epsilon]max,Floor[Log10[Abs[(KerrSEQ[[ind0+1,2,1]]+KerrSEQ[[ind0,2,1]])/2]]-2.5]],\[Epsilon]min];
				QNMsol=QNMSolution[inversion,s,l,m,a0+\[CapitalDelta]ap/2,
									SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
									SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon]2,relax,
									Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
				If[QNMsol[[1]],
					Print["QNMsol++ a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
					KerrSEQ=Insert[KerrSEQ,QNMsol[[4]],ind0+1];
					,(* invalid solution *)
					Print["a++ solution failed."];
					Abort[];
				];
			,If[\[CapitalDelta]am>2\[CapitalDelta]ap,
				\[Omega]g = (KerrSEQ[[ind0,2,1]]+KerrSEQ[[ind0-1,2,1]])/2;
				Almg = (KerrSEQ[[ind0,3,1]]+KerrSEQ[[ind0-1,3,1]])/2;
				Nrcf=Max[If[Length[KerrSEQ[[ind0,2]]]>=6,KerrSEQ[[ind0,2,6]],KerrSEQ[[ind0,2,3]]],
							If[Length[KerrSEQ[[ind0-1,2]]]>=6,KerrSEQ[[ind0-1,2,6]],KerrSEQ[[ind0-1,2,3]]]];
				Nrcf=Max[Nrcf,RCFmin];
				\[Epsilon]2=Max[Min[\[Epsilon]max,Floor[Log10[Abs[(KerrSEQ[[ind0,2,1]]-KerrSEQ[[ind0-1,2,1]])/2]]-2.5]],\[Epsilon]min];
				QNMsol=QNMSolution[inversion,s,l,m,a0-\[CapitalDelta]am/2,
									SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
									SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon]2,relax,
									Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
				If[QNMsol[[1]],
					Print["QNMsol-- a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
					KerrSEQ=Insert[KerrSEQ,QNMsol[[4]],ind0++]; (* must increment to keep ind0 at same a *)
					,(* invalid solution *)
					Print["a-- solution failed."];
					Abort[];
				];
			]];
			\[CapitalDelta]ap = KerrSEQ[[ind0+1,1]]-KerrSEQ[[ind0,1]];
			\[CapitalDelta]am = KerrSEQ[[ind0,1]]-KerrSEQ[[ind0-1,1]];
		];
		,
		If[Not[recursflag] && 
			(blevel > Maxb || (Minb < blevel && curvrat < maxcurvrat/2 && \[CapitalDelta]\[Omega] < max\[CapitalDelta]\[Omega]/2
								&& ((forward && a0+3\[CapitalDelta]a<1) || (Not[forward] && a0-3\[CapitalDelta]a>=0)))),
(*Print["Untested section of code! 18"];*)
			(* Reduce resoltuion *)
			If[Mod[1000(a0+\[CapitalDelta]a),2000\[CapitalDelta]a]==0, (* Make sure we can end up one sensible values of a *)
				blevelm=blevelp=blevel-1;
				incres=True;
			]
		]
	];
(*
	KerrSEQ: full sequence list [must be saved at some point]
	blevel[p,m]: blevel to use for next extrapolated step
	incres: true if next step will have larger step size
	incstep: for a successful increas in step size, offset to [+/-] data
*)
(*Print["Untested section of code! 19"];*)
	If[forward,
		{KerrSEQ,blevelp,incres,incstep,\[Epsilon]p},
		{KerrSEQ,blevelm,incres,incstep,\[Epsilon]m}
	]
]


Options[KerrQNMRefineSequenceB]=Union[{QNMSpinWeight->Null[],Index->False,
								Refinement->All,RefinementAction->None,ForceRefinement->False,
								RefinementPlot->SeqLevel,LimitRefinement->None,
								SolutionRelax->1,RadialCFDepth->1,RadialCFMaxGuess->20000000},
								Options[AdaptCheck3],Options[ListLinePlot]];
Options[KerrQNMRefineSequenceBwork]=Options[KerrQNMRefineSequenceB];


KerrQNMRefineSequenceB[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]_Integer,
				opts:OptionsPattern[]]:=
Module[{QNMSavePrecision=$MinPrecision,saneopts},
	(* saneopts ensures options set via SetOptions[KerQNMRefineSequenceB,...] are used *)
	saneopts=Flatten[Union[{opts},FilterRules[Options[KerrQNMRefineSequenceB],Except[Flatten[{opts}]]]]];
	CheckAbort[KerrQNMRefineSequenceBwork[l,m,n,\[Epsilon],saneopts],
				$MinPrecision=QNMSavePrecision;Abort[]];
	$MinPrecision=QNMSavePrecision;
]


KerrQNMRefineSequenceBwork[l_Integer,m_Integer,n_Integer|n_List,\[Epsilon]max_Integer,
				opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,KerrSEQ,
		NKQNM,index0,index0m,index0p,alow,ahigh,refinechange=False,
		i,plotdata,plotdata1,\[Omega]dat,d\[Omega],dd\[Omega],ind0,a0,\[CapitalDelta]ap,\[CapitalDelta]am,\[Omega]g,Almg,QNMsol,inversion,\[Epsilon]=\[Epsilon]max,Nrcf,Nm,
		KerrSEQret,dummy,blevel,forward,incflag,limitlist={},ll,inc,dec,last,re\[Omega],width,
		indexmin,indexmax,offset,oldNrcf,newNrcf,oldCf,newCf,rcfmin,ref\[Epsilon],
		useindex=OptionValue[Index],
		precision=OptionValue[QNMPrecision],refinement=OptionValue[Refinement],
		action=OptionValue[RefinementAction],forcerefine=OptionValue[ForceRefinement],
		plottype=OptionValue[RefinementPlot],limitrefine=OptionValue[LimitRefinement],
		Minb=OptionValue[Minblevel],Maxb=OptionValue[Maxblevel],
		relax=Rationalize[OptionValue[SolutionRelax]],
		RCFmin=OptionValue[RadialCFMinDepth],RCFmax=OptionValue[RadialCFMaxGuess],
		rcfdepth=OptionValue[RadialCFDepth]},
	KerrQNMRefineSequenceB::Refinement="The value of Refinement (`1`) is not an integer index, real value for a, a list of either specifying a range, or ALL.";
	KerrQNMRefineSequenceB::index="using index `1` instead of `2`";
	KerrQNMRefineSequenceB::value="using a = `1` instead of `2`";
	KerrQNMRefineSequenceB::range="using range `1` instead of `2`";
	KerrQNMRefineSequenceB::list="`1` must be a 2 element list of either integers or reals";
	KerrQNMRefineSequenceB::invalidlist="`1` is an invalid range of elements";
	KerrQNMRefineSequenceB::sequence="Sequence has `1` elements; must have at least 3 to refine";
	KerrQNMRefineSequenceB::badplot="Invalid RefinementPlot `1` given";
	KerrQNMRefineSequenceB::badaction="Invalid RefinementAction `1` given";
	KerrQNMRefineSequenceB::limits="Invalid LimitRefinement `1` given";

	SpinWeightTable:=Switch[s,
						-2,Global`KerrQNM,
						-1,Global`KerrQNMe,
						 0,Global`KerrQNMs,
						 _,Print["Invalid QNMSpinWeight"];Abort[]
					];
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n]
					];
	NKQNM=Length[KerrSEQ];
	If[NKQNM<3,Message[KerrQNMRefineSequenceB::sequence,NKQNM];Return[]];
	If[action==RefineAccuracy || action==RefinePrecision || action==Update || action==None,
		indexmin=1;indexmax=NKQNM;offset=0
		,Null[],
		indexmin=2;indexmax=NKQNM-1;offset=1
	];
(* Parse Refinement option to set range of values to refine *)
	Switch[refinement
		,_Symbol,
			Switch[refinement
				,All,
					index0m=1;
					index0p=NKQNM;
					alow=KerrSEQ[[1,1]];
					ahigh=KerrSEQ[[NKQNM,1]];
					index0=indexmax;
				,_,Message[KerrQNMRefineSequenceB::Refinement,refinement];Return[]
			];
		,_Integer,
			index0=refinement;
			If[index0<0,index0+=NKQNM+1];
			If[index0<indexmin,refinechange=True;index0=indexmin];
			If[index0>indexmax,refinechange=True;index0=indexmax];
			If[refinechange,If[refinement>=0,
				Message[KerrQNMRefineSequenceB::index,index0,refinement],
				Message[KerrQNMRefineSequenceB::index,index0-NKQNM-1,refinement]
			]];
			index0m=index0-offset;
			index0p=index0+offset;
			alow=KerrSEQ[[index0m,1]];ahigh=KerrSEQ[[index0p,1]];
		,_Rational|_Real,
			index0=1;
			While[KerrSEQ[[index0,1]]<refinement && index0<NKQNM,++index0];
			If[index0<indexmin,refinechange=True;index0=indexmin];
			If[index0>indexmax,refinechange=True;index0=indexmax];
			If[KerrSEQ[[index0,1]]!=refinement,
				refinechange=True;
				If[index0<NKQNM-1 &&
					(Abs[KerrSEQ[[index0,1]]-refinement]>Abs[KerrSEQ[[index0+1,1]]-refinement]),
						++index0];
			];
			If[refinechange,Message[KerrQNMRefineSequenceB::value,Block[{$MinPrecision=0},N[KerrSEQ[[index0,1]],{Infinity,20}]],refinement]];
			index0m=index0-offset;
			index0p=index0+offset;
			alow=KerrSEQ[[index0m,1]];ahigh=KerrSEQ[[index0p,1]];
		,_List,
			If[Length[refinement]!=2,Message[KerrQNMRefineSequenceB::list,refinement];Return[]];
			Switch[refinement[[1]]
				,_Integer,
					index0m=refinement[[1]];index0p=refinement[[2]];
					If[Head[index0p]==Integer,Null,Null,Message[KerrQNMRefineSequenceB::list,refinement];Return[]];
					If[index0m<0,index0m+=NKQNM+1];
					If[index0p<0,index0p+=NKQNM+1];
					If[index0m>index0p,Message[KerrQNMRefineSequenceB::invalidlist,refinement];Return[]];
					If[index0m<1,refinechange=True;index0m=1];
					If[index0p>NKQNM,refinechange=True;index0p=NKQNM];
					If[index0m>index0p-offset,refinechange=True;index0m=index0p-offset];
					If[index0p<index0m+offset,refinechange=True;index0p=index0m+offset];
					If[index0m==NKQNM,index0m=indexmax-offset];
					If[refinechange,Message[KerrQNMRefineSequenceB::range,
						{If[refinement[[1]]<0,index0m-NKQNM-1,index0m],
						 If[refinement[[2]]<0,index0p-NKQNM-1,index0p]},refinement]];
					alow=KerrSEQ[[index0m,1]];
					ahigh=KerrSEQ[[index0p,1]];
					index0=index0p-1;
				,_Rational|_Real,
					If[Head[refinement[[2]]]==Real||Head[refinement[[2]]]==Rational,Null,Null,Message[KerrQNMRefineSequenceB::list,refinement];Return[]];
					If[refinement[[1]]>refinement[[2]],Message[KerrQNMRefineSequenceB::invalidlist,refinement];Return[]];
					index0m=1;
					While[KerrSEQ[[index0m,1]]<refinement[[1]] && index0m<NKQNM,++index0m];
					If[KerrSEQ[[index0m,1]]!=refinement[[1]],
						refinechange=True;
						If[index0m<NKQNM &&
							(Abs[KerrSEQ[[index0m,1]]-refinement[[1]]]>Abs[KerrSEQ[[index0m+1,1]]-refinement[[1]]]),
								++index0m];
					];
					index0p=1;
					While[KerrSEQ[[index0p,1]]<refinement[[2]] && index0p<NKQNM,++index0p];
					If[KerrSEQ[[index0p,1]]!=refinement[[2]],
						refinechange=True;
						If[index0p<NKQNM &&
							(Abs[KerrSEQ[[index0p,1]]-refinement[[2]]]>Abs[KerrSEQ[[index0p+1,1]]-refinement[[2]]]),
								++index0p];
					];
					If[index0p>NKQNM,refinechange=True;index0p=NKQNM];
					If[index0m>index0p-1-offset,refinechange=True;index0m=index0p-1-offset;];
					If[index0p<index0m+1+offset,refinechange=True;index0p=index0m+1+offset;];
					alow=KerrSEQ[[index0m,1]];
					ahigh=KerrSEQ[[index0p,1]];
					If[refinechange,Message[KerrQNMRefineSequenceB::range,Block[{$MinPrecision=0},N[{alow,ahigh},{Infinity,20}]],refinement]];
					index0=index0p-1;
				,_,
					Message[KerrQNMRefineSequenceB::list,refinement];Return[]
			];
		,_,Message[KerrQNMRefineSequenceB::Refinement,refinement];Return[]
	];
	Switch[limitrefine
		,None,
			limitlist={{index0m,index0p}};
		,Minima,
			inc=False;dec=False;
			last=Re[KerrSEQ[[index0m,2,1]]];
			For[i=index0m+1,i<=index0p,++i,
				re\[Omega]=Re[KerrSEQ[[i,2,1]]];
				If[re\[Omega]>=last,
					If[dec && Not[inc],PrependTo[limitlist,{i-1,i-1}]];
					If[Not[inc],inc=True,dec=False],
					(*re\[Omega]<last*)
					If[Not[dec],dec=True,inc=False]
				];
				last=re\[Omega];
			];
		,_List,
			Switch[limitrefine[[1]]
				,Minima,
					Switch[limitrefine[[2]]
						,_Integer,
							width=limitrefine[[2]];
						,_,
							Message[KerrQNMRefineSequenceB::limits,limitrefine];Return[]
					];
					inc=False;dec=False;
					last=Re[KerrSEQ[[index0m,2,1]]];
					For[i=index0m+1,i<=index0p,++i,
						re\[Omega]=Re[KerrSEQ[[i,2,1]]];
						If[re\[Omega]>=last,
							If[dec && Not[inc],PrependTo[limitlist,{i-1-width,i-1+width}]];
							If[Not[inc],inc=True,dec=False],
							(*re\[Omega]<last*)
							If[Not[dec],dec=True,inc=False]
						];
						last=re\[Omega];
					];
					If[Length[limitrefine]==3,
						rcfmin=limitrefine[[3]];
						ll=limitlist;limitlist={};
						For[i=1,i<=Length[ll],++i,
							PrependTo[limitlist,Table[j,{j,ll[[i,1]],ll[[i,2]]}]];
						];
						ll=Flatten[limitlist];limitlist={};
						For[i=1,i<=Length[ll],++i,
							If[Length[KerrSEQ[[ll[[i]],2]]]>=6 && KerrSEQ[[ll[[i]],2,6]]<=rcfmin,
								PrependTo[limitlist,{ll[[i]],ll[[i]]}];
							];
						];
					];
				,RadialCFMinDepth,
					Switch[limitrefine[[2]]
						,_Integer,
							rcfmin=limitrefine[[2]];
						,_,
							Message[KerrQNMRefineSequenceB::limits,limitrefine];Return[]
					];
					For[i=index0m,i<=index0p,++i,
						If[Length[KerrSEQ[[i,2]]]>=6 && KerrSEQ[[i,2,6]]<=rcfmin,
							PrependTo[limitlist,{i,i}];
						];
					];
				,_,
					Message[KerrQNMRefineSequenceB::limits,limitrefine];Return[]
			];
		,_,
			Message[KerrQNMRefineSequenceB::limits,limitrefine];Return[]
	];
	Switch[plottype
		,None,
			plotdata=None;
		,SeqLevel,
			plotdata=Table[{If[useindex,i,KerrSEQ[[i,1]]],Round[-(3+Log10[KerrSEQ[[i+1,1]]-KerrSEQ[[i,1]]])/Log10[2]]},{i,Min[NKQNM-1,index0m],Min[NKQNM-1,index0p]}];
		,RadialCFLevel,
			plotdata=Table[{If[useindex,i,KerrSEQ[[i,1]]],If[Length[KerrSEQ[[i,2]]]>=6,KerrSEQ[[i,2,6]],KerrSEQ[[i,2,3]],0]},{i,Min[NKQNM-1,index0m],Min[NKQNM-1,index0p]}];
		,AccuracyLevel,
			plotdata=Table[{If[useindex,i,KerrSEQ[[i,1]]],KerrSEQ[[i,2,4]]},{i,Max[1,index0m],Min[NKQNM,index0p]}];
		,PrecisionLevel,
			plotdata=Table[{If[useindex,i,KerrSEQ[[i,1]]],MyPrecision[KerrSEQ[[i,2,1]]]},{i,Max[1,index0m],Min[NKQNM,index0p]}];
		,StepRatio,
			plotdata1=Table[{KerrSEQ[[i,1]],i},{i,index0m,index0p}];
			plotdata=RotateLeft[plotdata1,1]-plotdata1;
			plotdata=RotateRight[plotdata,1]/plotdata;
			plotdata=Table[{If[useindex,plotdata1[[i,2]],plotdata1[[i,1]]],plotdata[[i,1]]},{i,2,Length[plotdata1]-1}];
		,CurveRatio,
			plotdata1=Table[{KerrSEQ[[i,1]],i},{i,index0m,index0p}];
			plotdata=\[Omega]dat=Table[KerrSEQ[[i,2,1]],{i,index0m,index0p}];
			For[i=2,i<=Length[plotdata1]-1,++i,
				If[plotdata1[[i+1,1]]-plotdata1[[i,1]]==plotdata1[[i,1]]-plotdata1[[i-1,1]] || i==2 || i==Length[plotdata1]-1,
					d\[Omega]=\[Omega]dat[[i+1]]-\[Omega]dat[[i-1]];
					dd\[Omega]=\[Omega]dat[[i+1]]-2\[Omega]dat[[i]]+\[Omega]dat[[i-1]];
				,If[plotdata1[[i+1,1]]-plotdata1[[i,1]]>plotdata1[[i,1]]-plotdata1[[i-1,1]],
					d\[Omega]=\[Omega]dat[[i+1]]-\[Omega]dat[[i-2]];
					dd\[Omega]=\[Omega]dat[[i+1]]-2\[Omega]dat[[i]]+\[Omega]dat[[i-2]];
				,
					d\[Omega]=\[Omega]dat[[i+2]]-\[Omega]dat[[i-1]];
					dd\[Omega]=\[Omega]dat[[i+2]]-2\[Omega]dat[[i]]+\[Omega]dat[[i-1]];
				]];
				plotdata[[i]]=4Sqrt[Abs[d\[Omega]]^2Abs[dd\[Omega]]^2-(Re[d\[Omega]]Re[dd\[Omega]]+Im[d\[Omega]]Im[dd\[Omega]])^2]/(Abs[d\[Omega]]^2);
			];
			plotdata=Table[{If[useindex,plotdata1[[i,2]],plotdata1[[i,1]]],plotdata[[i]]},{i,2,Length[plotdata1]-1}];
		,_,Message[KerrQNMRefineSequenceB::badplot,plottype];Return[]
	];
	If[plottype==SeqLevel || plottype==RadialCFLevel || plottype==AccuracyLevel || plottype==PrecisionLevel || plottype==StepRatio || plottype==CurveRatio,
		Print[ListLinePlot[plotdata,FilterRules[{opts},Options[ListLinePlot]]]];
	];
	Switch[action
		,None,Null
		,RefineAccuracy,
			If[forcerefine && precision!=$MinPrecision,Print["Set $MinPrecision = ",precision]];
			$MinPrecision=precision;
			For[index0=index0m,index0<=index0p,++index0,
				While[index0>limitlist[[-1,2]],
					limitlist=Drop[limitlist,-1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0p+1;Continue[]];
				If[(\[Epsilon]<KerrSEQ[[index0,2,4]] || forcerefine) && (limitlist[[-1,1]]<=index0<=limitlist[[-1,2]]),
Print["RefineAcc at : ",index0];
					inversion=KerrSEQ[[index0,2,2]];
					oldNrcf=If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]];
(*Print["Old Nrcf : ",oldNrcf];*)
					Nrcf=If[Length[KerrSEQ[[index0,2]]]>=8,
						If[NumberQ[KerrSEQ[[index0,2,7]]],
							IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/KerrSEQ[[index0,2,7]])],
							oldNrcf],
						IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/3)]];
(*Print["New Nrcf : ",Nrcf];*)
					If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
					If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
					Nrcf=Max[Nrcf,RCFmin];
					If[Nrcf>RCFmax,
						Print["Warning: Computed radial CF depth (",Nrcf,") too large."];
						Print["         Setting CF depth to ",RCFmax];
						Nrcf=RCFmax;
					];
					Nm=KerrSEQ[[index0,3,2]];
					a0=KerrSEQ[[index0,1]];
					\[Omega]g=KerrSEQ[[index0,2,1]];
					Almg=KerrSEQ[[index0,3,1]];
					If[!forcerefine,
						$MinPrecision=Max[precision,
							If[Length[KerrSEQ[[index0,2]]]>=9,KerrSEQ[[index0,2,9]],IntegerPart[MyPrecision[\[Omega]g]]]]
					];
					QNMsol=QNMSolution[inversion,s,l,m,a0,
										SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
										SetPrecision[Almg,Max[precision,$MinPrecision]],\[Epsilon],relax,
										Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
					If[QNMsol[[1]],
						Print["QNMsol a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision],
							" Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision],
							"  |\[CapitalDelta]\[Omega]| = ",SetPrecision[Abs[\[Omega]g-QNMsol[[4,2,1]]],MachinePrecision]];
						(*
						oldCf=RadialCFRem[inversion,s,m,a0,Almg,\[Omega]g,oldNrcf];
						newNrcf=QNMsol[[4,2,3]];(* Must be the value used, not the "best" *)
						newCf=RadialCFRem[inversion,s,m,a0,Almg,\[Omega]g,newNrcf];
						Print["Prior Accuracy : ",KerrSEQ[[index0,2,4]]," Pred|\[CapitalDelta]\[Omega]| = ",1/Sqrt[Det[QNMsol[[5]]]]Abs[newCf[[1]]-oldCf[[1]]] ];
						*)
						Switch[s,
							   -2,Global`KerrQNM[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]],
							   -1,Global`KerrQNMe[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]],
								0,Global`KerrQNMs[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]]
							  ];
						,(* invalid solution *)
						Print["Solution failed at index ",index0];
					];
				];
			];
		,RefinePrecision,
			If[precision!=$MinPrecision,Print["Set $MinPrecision = ",precision]];
			$MinPrecision=precision;
			For[index0=index0m,index0<=index0p,++index0,
				While[index0>limitlist[[-1,2]],
					limitlist=Drop[limitlist,-1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0p+1;Continue[]];
				If[(precision>MyPrecision[KerrSEQ[[index0,2,1]]] || forcerefine) && (limitlist[[-1,1]]<=index0<=limitlist[[-1,2]]),
					inversion=KerrSEQ[[index0,2,2]];
					oldNrcf=If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]];
					Nrcf=If[Length[KerrSEQ[[index0,2]]]>=8,
						If[NumberQ[KerrSEQ[[index0,2,7]]],
							IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/KerrSEQ[[index0,2,7]])],
							oldNrcf],
						IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/3)]];
					If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
					If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
					Nrcf=Max[Nrcf,RCFmin];
					If[Nrcf>RCFmax,
						Print["Warning: Computed radial CF depth (",Nrcf,") too large."];
						Print["         Setting CF depth to ",RCFmax];
						Nrcf=RCFmax;
					];
					Nm=KerrSEQ[[index0,3,2]];
					a0=KerrSEQ[[index0,1]];
					\[Omega]g=KerrSEQ[[index0,2,1]];
					Almg=KerrSEQ[[index0,3,1]];
					ref\[Epsilon]=If[forcerefine,\[Epsilon],Min[\[Epsilon],KerrSEQ[[index0,2,4]]]];
					If[!forcerefine,
						$MinPrecision=Max[precision,
							If[Length[KerrSEQ[[index0,2]]]>=9,KerrSEQ[[index0,2,9]],IntegerPart[MyPrecision[\[Omega]g]]]]
					];
					QNMsol=QNMSolution[inversion,s,l,m,a0,
										SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
										SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
										Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
					If[QNMsol[[1]],
						Print["QNMsol a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision],
							" Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision],
							"  |\[CapitalDelta]\[Omega]| = ",SetPrecision[Abs[\[Omega]g-QNMsol[[4,2,1]]],MachinePrecision]];
						Switch[s,
							   -2,Global`KerrQNM[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]],
							   -1,Global`KerrQNMe[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]],
								0,Global`KerrQNMs[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]]
							  ];
						,(* invalid solution *)
						Print["Solution failed at index ",index0];
					];
				];
			];
		,RefineAdapt,
			If[precision!=$MinPrecision,Print["Set $MinPrecision = ",precision]];
			$MinPrecision=precision;
			incflag=False;
			While[index0>index0m,
				While[index0<limitlist[[1,1]],
					limitlist=Drop[limitlist,1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0m;Continue[]];
				If[index0>limitlist[[1,2]],--index0;Continue[]];
				inversion=KerrSEQ[[index0,2,2]];
				Nrcf=If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]];
				If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
				If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
				Nrcf=Max[Nrcf,RCFmin];
				If[Nrcf>RCFmax,
					Print["Warning: Computed radial CF depth (",Nrcf,") too large."];
					Print["         Setting CF depth to ",RCFmax];
					Nrcf=RCFmax;
				];
				Nm=KerrSEQ[[index0,3,2]];
				blevel=Round[-(3+Log10[KerrSEQ[[index0+1,1]]-KerrSEQ[[index0,1]]])/Log10[2]];
				forward=True;
				incflag=False;
				If[KerrSEQ[[index0+1,1]]-KerrSEQ[[index0,1]]>KerrSEQ[[index0,1]]-KerrSEQ[[index0-1,1]],
					incflag=True,
					forward=False;
					blevel=Round[-(3+Log10[KerrSEQ[[index0,1]]-KerrSEQ[[index0-1,1]]])/Log10[2]];
					If[KerrSEQ[[index0+1,1]]-KerrSEQ[[index0,1]]<KerrSEQ[[index0,1]]-KerrSEQ[[index0-1,1]],
						incflag=True;
					];
				];
				\[Epsilon]=Min[\[Epsilon]max,KerrSEQ[[index0,2,4]]];
				$MinPrecision=If[Length[KerrSEQ[[index0,2]]]>=9,KerrSEQ[[index0,2,9]],IntegerPart[MyPrecision[KerrSEQ[[index0,2,1]]]]];
(*Print["At a = ",N[KerrSEQ[[index0,1]]]," index0 = ",index0," forward : ",forward," incflag : ",incflag];*)
				{KerrSEQret,blevel,dummy,dummy,\[Epsilon]}=
					AdaptCheck3[KerrSEQ,inversion,s,l,m,\[Epsilon],relax,index0,blevel,forward,incflag,False,FilterRules[{opts},Options[AdaptCheck3]]];
				Switch[s,
					   -2,Global`KerrQNM[l,m,n]=KerrSEQret,
					   -1,Global`KerrQNMe[l,m,n]=KerrSEQret,
						0,Global`KerrQNMs[l,m,n]=KerrSEQret
					  ];
				incflag=True;
				--index0;
			];
		,FixAdapt,
			If[precision!=$MinPrecision,Print["Set $MinPrecision = ",precision]];
			$MinPrecision=precision;
			plotdata1=Table[KerrSEQ[[i,1]],{i,index0m,index0p}];
			plotdata=RotateLeft[plotdata1,1]-plotdata1;
			plotdata=RotateRight[plotdata,1]/plotdata;
			For[i=Length[plotdata1]-1,i>1,--i,
				ind0=index0m+i-1;
				While[ind0<limitlist[[1,1]],
					limitlist=Drop[limitlist,1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,i=1;Continue[]];
				If[ind0>limitlist[[1,2]],Continue[]];
				If[Not[1/2<=plotdata[[i]]<=2],
					a0=plotdata1[[i]];
					inversion=KerrSEQ[[ind0,2,2]];
					Nrcf=If[Length[KerrSEQ[[ind0,2]]]>=6,KerrSEQ[[ind0,2,6]],KerrSEQ[[ind0,2,3]]];
					If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
					If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
					Nrcf=Max[Nrcf,RCFmin];
					If[Nrcf>RCFmax,
						Print["Warning: Computed radial CF depth (",Nrcf,") too large."];
						Print["         Setting CF depth to ",RCFmax];
						Nrcf=RCFmax;
					];
					Nm=KerrSEQ[[ind0,3,2]];
					Print["Adaptation error at a = ",N[plotdata1[[i]]]," : ratio = ",plotdata[[i]]];
					\[CapitalDelta]ap=plotdata1[[i+1]]-plotdata1[[i]];
					\[CapitalDelta]am=plotdata1[[i]]-plotdata1[[i-1]];
					While[\[CapitalDelta]ap>2\[CapitalDelta]am || \[CapitalDelta]am>2\[CapitalDelta]ap,
						If[\[CapitalDelta]ap>2\[CapitalDelta]am,
							\[Omega]g=(KerrSEQ[[ind0+1,2,1]]+KerrSEQ[[ind0,2,1]])/2;
							Almg = (KerrSEQ[[ind0+1,3,1]]+KerrSEQ[[ind0,3,1]])/2;
							ref\[Epsilon]=If[forcerefine,\[Epsilon],Min[\[Epsilon],Max[KerrSEQ[[index0,2,4]],KerrSEQ[[index0+1,2,4]]]]];
							$MinPrecision=If[Length[KerrSEQ[[ind0,2]]]>=9,KerrSEQ[[ind0,2,9]],IntegerPart[MyPrecision[KerrSEQ[[ind0,2,1]]]]];
							QNMsol=QNMSolution[inversion,s,l,m,a0+\[CapitalDelta]ap/2,
												SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
												SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
											Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
							If[QNMsol[[1]],
								Print["QNMsol++ a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
								Switch[s,
								   -2,Global`KerrQNM[l,m,n] =Insert[KerrSEQ,QNMsol[[4]],ind0+1],
								   -1,Global`KerrQNMe[l,m,n]=Insert[KerrSEQ,QNMsol[[4]],ind0+1],
									0,Global`KerrQNMs[l,m,n]=Insert[KerrSEQ,QNMsol[[4]],ind0+1]
								  ];
								,(* invalid solution *)
								Print["a++ solution failed."];
								Abort[];
							];
						,If[\[CapitalDelta]am>2\[CapitalDelta]ap,
							\[Omega]g = (KerrSEQ[[ind0,2,1]]+KerrSEQ[[ind0-1,2,1]])/2;
							Almg = (KerrSEQ[[ind0,3,1]]+KerrSEQ[[ind0-1,3,1]])/2;
							ref\[Epsilon]=If[forcerefine,\[Epsilon],Min[\[Epsilon],Max[KerrSEQ[[index0,2,4]],KerrSEQ[[index0-1,2,4]]]]];
							$MinPrecision=If[Length[KerrSEQ[[ind0,2]]]>=9,KerrSEQ[[ind0,2,9]],IntegerPart[MyPrecision[KerrSEQ[[ind0,2,1]]]]];
							QNMsol=QNMSolution[inversion,s,l,m,a0-\[CapitalDelta]am/2,
												SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
												SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
												Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
							If[QNMsol[[1]],
								Print["QNMsol-- a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision]," Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision]];
								Switch[s,
								   -2,Global`KerrQNM[l,m,n] =Insert[KerrSEQ,QNMsol[[4]],ind0],
								   -1,Global`KerrQNMe[l,m,n]=Insert[KerrSEQ,QNMsol[[4]],ind0],
									0,Global`KerrQNMs[l,m,n]=Insert[KerrSEQ,QNMsol[[4]],ind0]
								  ];
								++ind0; (* must increment to keep ind0 at same a *)
								,(* invalid solution *)
								Print["a-- solution failed."];
								Abort[];
							];
						]];
						\[CapitalDelta]ap = KerrSEQ[[ind0+1,1]]-KerrSEQ[[ind0,1]];
						\[CapitalDelta]am = KerrSEQ[[ind0,1]]-KerrSEQ[[ind0-1,1]];
					];
				];
			];
		,RemoveLevels,
			While[index0>index0m,
				While[index0<limitlist[[1,1]],
					limitlist=Drop[limitlist,1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0m;Continue[]];
				If[index0>limitlist[[1,2]],--index0;Continue[]];
				If[Mod[1000KerrSEQ[[index0,1]],2^(-Maxb)]!=0,
					Switch[s,
					   -2,Global`KerrQNM[l,m,n] =Drop[KerrSEQ,{index0}],
					   -1,Global`KerrQNMe[l,m,n]=Drop[KerrSEQ,{index0}],
						0,Global`KerrQNMs[l,m,n]=Drop[KerrSEQ,{index0}]
					  ];
				];
				--index0;
			];
		,Update,
			If[precision!=$MinPrecision,Print["Set $MinPrecision = ",precision]];
			$MinPrecision=precision;
			For[index0=index0m,index0<=index0p,++index0,
				While[index0>limitlist[[-1,2]],
					limitlist=Drop[limitlist,-1];
					If[Length[limitlist]==0,Break[]];
				];
				If[Length[limitlist]==0,index0=index0p+1;Continue[]];
				If[(Length[KerrSEQ[[index0,2]]]<9 || forcerefine) && (limitlist[[-1,1]]<=index0<=limitlist[[-1,2]]),
					inversion=KerrSEQ[[index0,2,2]];
					oldNrcf=If[Length[KerrSEQ[[index0,2]]]>=6,KerrSEQ[[index0,2,6]],KerrSEQ[[index0,2,3]]];
					Nrcf=If[Length[KerrSEQ[[index0,2]]]>=8,
						If[NumberQ[KerrSEQ[[index0,2,7]]],
							IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/KerrSEQ[[index0,2,7]])],
							oldNrcf],
						IntegerPart[oldNrcf 10^((\[Epsilon]-KerrSEQ[[index0,2,4]])/3)]];
					If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
					If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
					Nrcf=Max[Nrcf,RCFmin];
					If[Nrcf>RCFmax,
						Print["Warning: Computed radial CF depth (",Nrcf,") too large."];
						Print["         Setting CF depth to ",RCFmax];
						Nrcf=RCFmax;
					];
					Nm=KerrSEQ[[index0,3,2]];
					a0=KerrSEQ[[index0,1]];
					\[Omega]g=KerrSEQ[[index0,2,1]];
					Almg=KerrSEQ[[index0,3,1]];
					ref\[Epsilon]=If[forcerefine,\[Epsilon],Min[\[Epsilon],KerrSEQ[[index0,2,4]]]];
					If[!forcerefine,
						$MinPrecision=Max[precision,
							If[Length[KerrSEQ[[index0,2]]]>=9,KerrSEQ[[index0,2,9]],IntegerPart[MyPrecision[\[Omega]g]]]]
					];
					QNMsol=QNMSolution[inversion,s,l,m,a0,
										SetPrecision[\[Omega]g,Max[precision,$MinPrecision]],
										SetPrecision[Almg,Max[precision,$MinPrecision]],ref\[Epsilon],relax,
										Nrcf,Nm,0,0,0,0,FilterRules[{opts},Options[QNMSolution]]];
					If[QNMsol[[1]],
						Print["QNMsol a=",Block[{$MinPrecision=0},N[QNMsol[[4,1]],{Infinity,20}]]," \[Omega]=",SetPrecision[QNMsol[[4,2,1]],MachinePrecision],
							" Alm=",SetPrecision[QNMsol[[4,3,1]],MachinePrecision],
							"  |\[CapitalDelta]\[Omega]| = ",SetPrecision[Abs[\[Omega]g-QNMsol[[4,2,1]]],MachinePrecision]];
						Switch[s,
							   -2,Global`KerrQNM[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]],
							   -1,Global`KerrQNMe[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]],
								0,Global`KerrQNMs[l,m,n]=ReplacePart[KerrSEQ,index0->QNMsol[[4]]]
							  ];
						,(* invalid solution *)
						Print["Solution failed at index ",index0];
					];
				];
			];
		,_,Message[KerrQNMRefineSequenceB::badaction,action];Return[]
	];
]




(* ::Subsection::Closed:: *)
(*Stepsize and solution validation routines*)


SolutionWindow[F0_?NumberQ,Fg_NumberQ,Fs_NumberQ,
				rl_Real|rl_Rational|rl_Integer,
				rt_Real|rt_Rational|rt_Integer,
				plot_/;plot \[Element] Booleans]:=
Module[{Fn,Ff,Fp,Fm,Fsv,tp,tm,rad},
(*
F0: Prior solution;
Fg: Guessed solution;
Fs: Numerical solution, should be close to guess;
rl: Fraction of distance between F0 and Fg for length of error wedge (lateral ratio);
rt: Fraction of distance between F0 and Fg for width of error wedge (transverse ratio);
*)
	If[rt==0||rl==0,Return[True]];
	Fn=Fg-rl(Fg-F0);Ff=Fg+rl(Fg-F0);
	Fp=Ff-F0+I rt(Fg-F0);Fm=Ff-F0-I rt(Fg-F0);
	Fsv=Fs-F0;Fn=Abs[Fn-F0];Ff=Abs[Ff-F0];
	tp=Re[Fsv]Im[Fp]-Im[Fsv]Re[Fp];
	tm=Re[Fsv]Im[Fm]-Im[Fsv]Re[Fm];
	Fsv=Abs[Fsv];
	If[Fsv<Fn,(* Near fail *)Null,
	If[Fsv>Ff,(* Far fail *)Null,
	If[tp<0,(* Lateral fail + *)Null,
	If[tm>0,(* Lateral fail - *)Null,
	Return[True]]]]];
	If[plot,
		rad=Abs[Fg-F0]/50;
		Print[Show[Graphics[{Red,Circle[{Re[F0],Im[F0]},rad]}],Graphics[{Blue,Circle[{Re[Fg],Im[Fg]},rad]}],Graphics[{Black,Circle[{Re[Fp+F0],Im[Fp+F0]},rad]}],Graphics[{Green,Circle[{Re[Fm+F0],Im[Fm+F0]},rad]}],Graphics[{Black,Disk[{Re[Fs],Im[Fs]},rad]}]]]
	];
	Return[False];
]


Set\[CapitalDelta]a[a_Rational|a_Integer,\[CapitalDelta]alevel_Integer,\[Omega]save_List,\[CapitalDelta]alist_List,
		OptionsPattern[]]:=
Module[{\[CapitalDelta]\[Omega],\[CapitalDelta]\[Phi],\[Omega],\[Omega]p,\[Omega]pp,\[Omega]pdpp,\[Omega]pow,Ni,\[Omega]s1,\[Omega]s2,\[Omega]s3,\[Omega]si,
		a1,a2,a3,dat,fit,io,curv,\[CapitalDelta]abar,
		min\[CapitalDelta]alevel=OptionValue[Min\[CapitalDelta]alevel],
		max\[CapitalDelta]alevel=OptionValue[Max\[CapitalDelta]alevel],
		max\[CapitalDelta]\[Phi]=Abs[OptionValue[Max\[CapitalDelta]\[Phi]]],max\[CapitalDelta]\[Omega]=Abs[OptionValue[Max\[CapitalDelta]\[Omega]]],
		opoly=OptionValue[\[Omega]poly]},
	If[\[CapitalDelta]alevel>=Min[Length[\[CapitalDelta]alist],max\[CapitalDelta]alevel],Return[{False}]];
	If[a+\[CapitalDelta]alist[[\[CapitalDelta]alevel]]>=1,Return[{True,\[CapitalDelta]alist[[\[CapitalDelta]alevel+1]],\[CapitalDelta]alevel+1}]];
	If[NumberQ[opoly] && \[Omega]save[[3]]!=0,
		a1=a;a2=a-\[CapitalDelta]alist[[\[CapitalDelta]alevel]];a3=a-2\[CapitalDelta]alist[[\[CapitalDelta]alevel]];
		dat={{-Im[\[Omega]save[[1]]],Re[\[Omega]save[[1]]]},
			{-Im[\[Omega]save[[2]]],Re[\[Omega]save[[2]]]},
			{-Im[\[Omega]save[[3]]],Re[\[Omega]save[[3]]]}};
		fit=Fit[dat,{1,io,io^2},io];
(*Print["fit \[Omega] = ",fit];
Print[Plot[fit,{io,dat[[3,1]],2dat[[1,1]]-dat[[3,1]]}]];*)
		curv=2 Coefficient[fit,io^2];
		If[Abs[Re[\[Omega]save[[1]]]]<opoly && curv>0,
			\[CapitalDelta]abar=-Coefficient[fit,io]/curv;
			dat={{-Im[\[Omega]save[[1]]],a1},
				{-Im[\[Omega]save[[2]]],a2}};
			fit=Fit[dat,{1,io},io];
(*Print["fit a = ",fit];
Print[Show[Plot[fit,{io,-Im[\[Omega]save[[3]]],2dat[[1,1]]+Im[\[Omega]save[[3]]]}],ListPlot[dat]]];*)
			\[CapitalDelta]abar=fit/.io->\[CapitalDelta]abar;
			\[CapitalDelta]abar=\[CapitalDelta]abar-a;
			Print["     \[Omega] poly ratio : ",\[CapitalDelta]abar/\[CapitalDelta]alist[[\[CapitalDelta]alevel]]];
			If[(\[CapitalDelta]abar>0 && Abs[\[CapitalDelta]abar]/\[CapitalDelta]alist[[\[CapitalDelta]alevel]]<2)||
					(\[CapitalDelta]abar<0 && Abs[\[CapitalDelta]abar]/\[CapitalDelta]alist[[\[CapitalDelta]alevel]]<0.75),
				Return[{True,\[CapitalDelta]alist[[\[CapitalDelta]alevel+1]],\[CapitalDelta]alevel+1}]
			];
		]
	];
	\[Omega]s1=\[Omega]save[[1]];
	\[Omega]s2=\[Omega]save[[2]];
	\[Omega]s3=\[Omega]save[[3]];
(*
	\[CapitalDelta]\[Omega]={\[Omega]s2-\[Omega]s1,\[Omega]s3-\[Omega]s2};
	\[CapitalDelta]\[Phi]= Abs[Arg[\[CapitalDelta]\[Omega][[2]]]-Arg[\[CapitalDelta]\[Omega][[1]]]];
	\[CapitalDelta]\[Phi]=Abs[\[CapitalDelta]\[Omega][[2]]]Min[\[CapitalDelta]\[Phi],2\[Pi]-\[CapitalDelta]\[Phi]];
	(*\[Omega]p=Abs[\[Omega]s3-\[Omega]s1]/2;
	\[Omega]pp=Abs[\[Omega]s3-2\[Omega]s2+\[Omega]s1];
	\[Omega]pdpp=Re[\[Omega]s3-\[Omega]s1]/2 Re[\[Omega]s3-2\[Omega]s2+\[Omega]s1]+Im[\[Omega]s3-\[Omega]s1]/2 Im[\[Omega]s3-2\[Omega]s2+\[Omega]s1];
	\[CapitalDelta]\[Phi]=Sqrt[(\[Omega]p \[Omega]pp)^2-\[Omega]pdpp^2]/\[Omega]p^2;*)
Print["Old \[CapitalDelta]\[Phi]=",\[CapitalDelta]\[Phi]," Abs[\[CapitalDelta]\[Omega][[1]]]=",Abs[\[CapitalDelta]\[Omega][[1]]]];
*)
	\[Omega]si=3\[Omega]s1-3\[Omega]s2+\[Omega]s3;
(*Print["Extrapolated \[Omega] ",\[Omega]si];*)
	\[Omega]s3=\[Omega]s2;\[Omega]s2=\[Omega]s1;\[Omega]s1=\[Omega]si;
	\[CapitalDelta]\[Omega]={\[Omega]s2-\[Omega]s1,\[Omega]s3-\[Omega]s2};
	\[CapitalDelta]\[Phi]= Abs[Arg[\[CapitalDelta]\[Omega][[2]]]-Arg[\[CapitalDelta]\[Omega][[1]]]];
	\[CapitalDelta]\[Phi]=Abs[\[CapitalDelta]\[Omega][[2]]]Min[\[CapitalDelta]\[Phi],2\[Pi]-\[CapitalDelta]\[Phi]];
(*	\[Omega]p=Abs[\[Omega]s3-\[Omega]s1]/2;
	\[Omega]pp=Abs[\[Omega]s3-2\[Omega]s2+\[Omega]s1];
	\[Omega]pdpp=Re[\[Omega]s3-\[Omega]s1]/2 Re[\[Omega]s3-2\[Omega]s2+\[Omega]s1]+Im[\[Omega]s3-\[Omega]s1]/2 Im[\[Omega]s3-2\[Omega]s2+\[Omega]s1];
	\[CapitalDelta]\[Phi]=Sqrt[(\[Omega]p \[Omega]pp)^2-\[Omega]pdpp^2]/\[Omega]p^2;*)
(*Print["testing level ",\[CapitalDelta]alevel];
Print["New \[CapitalDelta]\[Phi]=",\[CapitalDelta]\[Phi]," Abs[\[CapitalDelta]\[Omega][[1]]]=",Abs[\[CapitalDelta]\[Omega][[1]]]];*)
	(* Force higher res *) 
	If[\[CapitalDelta]alevel<min\[CapitalDelta]alevel,Return[{True,\[CapitalDelta]alist[[\[CapitalDelta]alevel+1]],\[CapitalDelta]alevel+1}]];
	If[(\[CapitalDelta]\[Phi]>max\[CapitalDelta]\[Phi] || Abs[\[CapitalDelta]\[Omega][[1]]]>max\[CapitalDelta]\[Omega]),
		Return[{True,\[CapitalDelta]alist[[\[CapitalDelta]alevel+1]],\[CapitalDelta]alevel+1}];
	];
	Return[{False}];
]


Options[Decrease\[CapitalDelta]a]=Options[Set\[CapitalDelta]a];


Decrease\[CapitalDelta]a[\[CapitalDelta]alevel_Integer,\[CapitalDelta]alist_List,OptionsPattern[]]:=
Module[{max\[CapitalDelta]alevel=OptionValue[Max\[CapitalDelta]alevel]},
	If[\[CapitalDelta]alevel<Min[Length[\[CapitalDelta]alist],max\[CapitalDelta]alevel],
		Return[{True,\[CapitalDelta]alist[[\[CapitalDelta]alevel+1]],\[CapitalDelta]alevel+1}];
	];
	Return[{False}];
]


Options[Increase\[CapitalDelta]a]=Options[Set\[CapitalDelta]a];


Increase\[CapitalDelta]a[a_Rational|a_Integer,\[CapitalDelta]alevel_Integer,ExtrapInfo_List,\[CapitalDelta]alist_List,
			opts:OptionsPattern[]]:=
Module[{min\[CapitalDelta]alevel=OptionValue[Min\[CapitalDelta]alevel]},
	(* Force higher res *)
	If[\[CapitalDelta]alevel<=min\[CapitalDelta]alevel,Return[{False}]];
	If[\[CapitalDelta]alevel>1,
		If[a+\[CapitalDelta]alist[[\[CapitalDelta]alevel-1]]>=1,Return[False]];
		If[Set\[CapitalDelta]a[a,\[CapitalDelta]alevel-1,ExtrapInfo[[\[CapitalDelta]alevel-1,3]],\[CapitalDelta]alist,FilterRules[{opts},Options[Set\[CapitalDelta]a]]][[1]],Return[False]];
		Return[True];
	];
	Return[False];
]


RestartExtrapolation[s_Integer,l_Integer,m_Integer,n_Integer|n_List,
					\[CapitalDelta]alevel_Integer,\[CapitalDelta]alist_List,precision_]:=
Module[{KerrSEQ,NKQNM,lasta,a,i,pos,alist,ExtrapInfo={}},
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n],
					_,Print["Invalid QNMSpinWeight"];Abort[]
					];
	NKQNM=Length[KerrSEQ];
	alist = Table[KerrSEQ[[i,1]],{i,1,NKQNM}];
	lasta=alist[[NKQNM]];
	a=Round[lasta/\[CapitalDelta]alist[[\[CapitalDelta]alevel]]];
	For[i=\[CapitalDelta]alevel,i>=1,--i,
		PrependTo[ExtrapInfo,{Mod[a,10],a \[CapitalDelta]alist[[i]],0,0}];a=IntegerPart[a/10];
	];
	For[i=1,i<=\[CapitalDelta]alevel,++i,
		a=ExtrapInfo[[i,2]];
		pos=Position[Chop[alist-a,10^(-12)],0];
		If[Length[pos]==0,pos={{NKQNM}}]; (* fix anomalies *)
		ExtrapInfo[[i,2]]=0;
		ExtrapInfo[[i,3]]={SetPrecision[KerrSEQ[[pos[[1,1]],2,1]],precision],0,0};
		ExtrapInfo[[i,4]]={SetPrecision[KerrSEQ[[pos[[1,1]],3,1]],precision],0,0};
		If[lasta>=\[CapitalDelta]alist[[i]],
			pos=Position[Chop[alist-a+\[CapitalDelta]alist[[i]],10^(-12)],0];
			If[Length[pos]==1,
				ExtrapInfo[[i,2]]=1;
				ExtrapInfo[[i,3,2]]=SetPrecision[KerrSEQ[[pos[[1,1]],2,1]],precision];
				ExtrapInfo[[i,4,2]]=SetPrecision[KerrSEQ[[pos[[1,1]],3,1]],precision];
			];
			If[lasta>=2\[CapitalDelta]alist[[i]],
				pos=Position[Chop[alist-a+2\[CapitalDelta]alist[[i]],10^(-12)],0];
				If[Length[pos]==1,
					ExtrapInfo[[i,2]]=2;
					ExtrapInfo[[i,3,3]]=SetPrecision[KerrSEQ[[pos[[1,1]],2,1]],precision];
					ExtrapInfo[[i,4,3]]=SetPrecision[KerrSEQ[[pos[[1,1]],3,1]],precision];
				];
			];
		];
	];
	Return[ExtrapInfo];
]


(* ::Subsection::Closed:: *)
(*Utility routines*)


MyPrecision[x_?NumberQ]:=Module[{saveprecision,returnprecision},
	saveprecision=$MinPrecision;
	$MinPrecision=0;
	returnprecision=Precision[x];
	$MinPrecision=saveprecision;
	Return[returnprecision];
]


Options[VerifyExpansion]={QNMSpinWeight->Null[]};


VerifyExpansion[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],KerrSEQ,lmin,lindex,coefs,Na,i,maxcoef,lmax,norm},
	lmin = Max[Abs[m],2];
	lindex = l-lmin+1;
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n],
					_,Print["Invalid QNMSpinWeight"];Abort[]
					];
	Na=Length[KerrSEQ];
	If[Na == 0,Print["No data"];Return[]];
	For[i=1,i<=Na,++i,
		coefs = Chop[KerrSEQ[[i,3,3]],10^(-13)];
		maxcoef=Max[Abs[coefs]];
		lmax = Position[Abs[coefs],maxcoef];
		norm = Chop[coefs.Conjugate[coefs],10^(-13)];
		If[lmax!=lindex || norm!=1,
			Print["a= ",KerrSEQ[[i,1]]," : ",lmax-lindex," : ",norm];
			Return[];
		];
	];
	Print["Expansion Coefficients correct"];
]


Options[RadialCFErrorEst]={QNMSpinWeight->Null[],QNMPrecision->24};
Options[RadialCFErrorEstwork]=Options[RadialCFErrorEst];



RadialCFErrorEst[inversion_Integer,m_Integer,a_Rational|a_Integer,
				\[Omega]_?NumberQ,Alm_?NumberQ,\[Epsilon]_Integer,Nrcf_List,opts:OptionsPattern[]]:=
Module[{QNMSavePrecision=$MinPrecision},
	CheckAbort[RadialCFErrorEstwork[inversion,m,a,\[Omega],Alm,\[Epsilon],Nrcf,FilterRules[{opts},Options[RadialCFErrorEstwork]]],
				$MinPrecision=QNMSavePrecision;Abort[]];
	$MinPrecision=QNMSavePrecision;
]


RadialCFErrorEstwork[inversion_Integer,m_Integer,a_Rational|a_Integer,
				\[Omega]_?NumberQ,Alm_?NumberQ,\[Epsilon]_Integer,Nrcf_List,opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],precision=OptionValue[QNMPrecision],i,rcf={},rcferr={}},
	$MinPrecision=precision;
	For[i=1,i<=Length[Nrcf],++i,
		Print["Nrcf = ",Nrcf[[i]]," [Remainder : ",
			Abs[RadialCFRemainder[s,m,a,Alm(1+10^(\[Epsilon]+4)),\[Omega](1-10^(\[Epsilon]+4)),Nrcf[[i]]]],"]"];
		AppendTo[rcf,RadialCFRem[inversion,s,m,a,Alm(1+10^(\[Epsilon]+4)),\[Omega](1-10^(\[Epsilon]+4)),Nrcf[[i]]]];
		Print["Sol = ",RadialCFRem[inversion,s,m,a,Alm,\[Omega],Nrcf[[i]]]];(**)
	];
	For[i=1,i<Length[Nrcf],++i,
		AppendTo[rcferr,{rcf[[i,2]],Abs[rcf[[i,1]]-rcf[[-1,1]]]}];
	];
	Print[ListLogLogPlot[rcferr,PlotRange->All]];
	Print[CoefficientList[LinearModelFit[Log10[rcferr],n,n]["BestFit"],n]];
	Print[CoefficientList[LinearModelFit[Take[Log10[rcferr],2],n,n]["BestFit"],n]];
	Print[CoefficientList[LinearModelFit[Take[Log10[rcferr],-2],n,n]["BestFit"],n]];
	rcferr
]


RadialCFErrorEst[l_Integer,m_Integer,n_Integer|n_List,index_Integer,
				Nrcf_List,opts:OptionsPattern[]]:=
Module[{QNMSavePrecision=$MinPrecision},
	CheckAbort[RadialCFErrorEstwork[l,m,n,index,Nrcf,FilterRules[{opts},Options[RadialCFErrorEstwork]]],
				$MinPrecision=QNMSavePrecision;Abort[]];
	$MinPrecision=QNMSavePrecision;
]


RadialCFErrorEstwork[l_Integer,m_Integer,n_Integer|n_List,index_Integer,
				Nrcf_List,opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],precision=OptionValue[QNMPrecision],KerrSEQ,inversion,a,\[Omega],Alm,\[Epsilon]},
	$MinPrecision=precision;
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n],
					_,Print["Invalid QNMSpinWeight"];Abort[]
					];
	inversion=KerrSEQ[[index,2,2]];
	a=KerrSEQ[[index,1]];
	\[Omega]=KerrSEQ[[index,2,1]];
	Alm=KerrSEQ[[index,3,1]];
	\[Epsilon]=KerrSEQ[[index,2,4]];
	RadialCFErrorEst[inversion,m,a,\[Omega],Alm,\[Epsilon],Nrcf,FilterRules[{opts},Options[RadialCFErrorEst]]]
]


Options[KerrQNMMakeMultiplet]={QNMSpinWeight->Null[],OTmultiple->0};


KerrQNMMakeMultiplet[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],nm=OptionValue[OTmultiple],KerrSEQ},
	KerrSEQ=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n],
					_,Print["Invalid QNMSpinWeight"];Abort[]
					];
	If[Head[n]==Integer,
		Switch[s,
				-2,Global`KerrQNM[l,m,{n,nm}]=KerrSEQ;Global`KerrQNM[l,m,n]=.,
				-1,Global`KerrQNMe[l,m,{n,nm}]=KerrSEQ;Global`KerrQNMe[l,m,n]=.,
				 0,Global`KerrQNMs[l,m,{n,nm}]=KerrSEQ;Global`KerrQNMs[l,m,n]=.
				],
		Null[],
		If[n[[2]]!=nm,
			Switch[s,
					-2,Global`KerrQNM[l,m,{n[[1]],nm}]=KerrSEQ;Global`KerrQNM[l,m,n]=.,
					-1,Global`KerrQNMe[l,m,{n[[1]],nm}]=KerrSEQ;Global`KerrQNMe[l,m,n]=.,
					 0,Global`KerrQNMs[l,m,{n[[1]],nm}]=KerrSEQ;Global`KerrQNMs[l,m,n]=.
					]
		]
	];
]


Options[ShortenQNMSequence]={QNMSpinWeight->Null[]};


ShortenQNMSequence[l_Integer,m_Integer,n_Integer|n_List,N_Integer,OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],KerrSEQ,SeqStatus,na},
	KerrSEQ:=Switch[s,
					-2,Global`KerrQNM[l,m,n],
					-1,Global`KerrQNMe[l,m,n],
					 0,Global`KerrQNMs[l,m,n],
					_,Print["Invalid QNMSpinWeight"];Abort[]
					];
	SeqStatus=If[Head[KerrSEQ]==List,If[Length[KerrSEQ]>0,True,False,False],False,False];
	If[!SeqStatus,Print["KerrQNM[",l,",",m,",",n,"] does not exist."];Return[]];
	na=Length[KerrSEQ];
	Print["Original Length of KerrQNM[",l,",",m,",",n,"]=",na];
	Print["Removing last ",N," elements"];
	Switch[s,
		   -2,Global`KerrQNM[l,m,n]=Take[KerrSEQ,na-N],
		   -1,Global`KerrQNMe[l,m,n]=Take[KerrSEQ,na-N],
		    0,Global`KerrQNMs[l,m,n]=Take[KerrSEQ,na-N]
		  ];
]


MergeQNMSequence[seq1_List,seq2_List]:=
Module[{mergedseq={},ind1=1,ind2=1,
		N1=Length[seq1],N2=Length[seq2]},
	While[ind1<=N1 || ind2<=N2,
		If[ind1>N1,
			AppendTo[mergedseq,seq2[[ind2++]]]
			Continue[];
		];
		If[ind2>N2,
			AppendTo[mergedseq,seq1[[ind1++]]]
			Continue[];
		];
		If[seq1[[ind1,1]]<seq2[[ind2,1]],
			AppendTo[mergedseq,seq1[[ind1++]]],
			If[seq1[[ind1,1]]>seq2[[ind2,1]],
				AppendTo[mergedseq,seq2[[ind2++]]],
				If[seq1[[ind1,1]]==seq2[[ind2,1]],
					If[seq1[[ind1,2,4]]<seq2[[ind2,2,4]],
						AppendTo[mergedseq,seq1[[ind1]]],
						If[seq1[[ind1,2,4]]>seq2[[ind2,2,4]],
							AppendTo[mergedseq,seq2[[ind2]]],
							If[MyPrecision[seq1[[ind1,2,1]]]>=MyPrecision[seq2[[ind2,2,1]]],
								AppendTo[mergedseq,seq1[[ind1]]],
								AppendTo[mergedseq,seq2[[ind2]]]
							];
						];	
					];
					++ind1;++ind2;
				];
			];
		];
	];
	Return[mergedseq];
]


SetSpinWeight[s_Integer]:=
Module[{},
	Switch[s,
		   -2,Null[],
		   -1,Null[],
		    0,Null[],
			_,Print["Invalid Spin Weight : ",s];Abort[]
		  ];
	SetOptions[KerrQNMSequenceB,QNMSpinWeight->s];
	SetOptions[KerrQNMRefineSequenceB,QNMSpinWeight->s];
	SetOptions[KerrQNMSequenceB,QNMSpinWeight->s];
	SetOptions[KerrQNMRefineSequenceB,QNMSpinWeight->s];
	SetOptions[KerrQNMSequence,QNMSpinWeight->s];
	SetOptions[KerrQNMSequenceReverse,QNMSpinWeight->s];
	SetOptions[KerrQNMInsert,QNMSpinWeight->s];
	SetOptions[KerrQNMAccumulation,QNMSpinWeight->s];
	SetOptions[VerifyExpansion,QNMSpinWeight->s];
	SetOptions[RadialCFErrorEstwork,QNMSpinWeight->s];
	SetOptions[KerrQNMMakeMultiplet,QNMSpinWeight->s];
	SetOptions[ShortenQNMSequence,QNMSpinWeight->s];
	SetOptions[SchQNMguess,QNMSpinWeight->s];
	SetOptions[SchwarzschildQNM,QNMSpinWeight->s];
	SetOptions[QNMPlotSch,QNMSpinWeight->s];
	SetOptions[SchwarzschildOmega,QNMSpinWeight->s];
	SetOptions[KerrOmegaList,QNMSpinWeight->s];
	SetOptions[KerrOmegaListS,QNMSpinWeight->s];
	SetOptions[KerrAList,QNMSpinWeight->s];
	SetOptions[KerrAListS,QNMSpinWeight->s];
	SetOptions[QNMPlotOmega,QNMSpinWeight->s];
	SetOptions[QNMPlotA,QNMSpinWeight->s];
	SetOptions[QNMPlotOmegaTones,QNMSpinWeight->s];
	SetOptions[QNMPlotATones,QNMSpinWeight->s];
	SetOptions[QNMPlotAccumulation\[Omega],QNMSpinWeight->s];
	SetOptions[QNMPlotAccumulationAlm,QNMSpinWeight->s];
	Print["All KerrQNM routines set for Spin-Weight s = ",s];
]


If[!QNMDebug,Protect[QNMSolution,KerrQNMSequenceB,KerrQNMRefineSequenceB,
					KerrQNMSequenceBwork,KerrQNMRefineSequenceBwork,AdaptCheck3,
					KerrQNMSequence,KerrQNMSequenceReverse,
					KerrQNMAccumulation,SolutionWindow,Set\[CapitalDelta]a,Decrease\[CapitalDelta]a,
					Increase\[CapitalDelta]a,RestartExtrapolation,VerifyExpansion,
					ShortenQNMSequence,SetSpinWeight]];




(* ::Section::Closed:: *)
(*Initial Guesses*)


Options[SchQNMguess]={QNMSpinWeight->Null[]};


SchQNMguess[l_Integer,n_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,guess},
	SpinWeightTable:=Switch[s,
				   -2,Global`SchQNMTable,
				   -1,Global`SchQNMeTable,
					0,Global`SchQNMsTable,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	If[n<SpinWeightTable[l],
		(* Initial guess is in Table *)
		SpinWeightTable[l,n],
		(* Interpolate initial Guess *)
		If[n<=0,Print[SpinWeightTable," not properly initialized"];Abort[]];
		guess=2 SchQNMguess[l,n-1,FilterRules[{opts},Options[SchQNMguess]]]-SchQNMguess[l,n-2,FilterRules[{opts},Options[SchQNMguess]]];
		guess[[1]]=Abs[Re[guess[[1]]]]+I Im[guess[[1]]];
		guess[[2]]=n; (* Use preferred inversion *)
		guess,
		(* No guesses yest for l *)
		If[l<=s,Print["Problem with ",SpinWeightTable];Abort[]];
		guess=2 SchQNMguess[l-1,0,FilterRules[{opts},Options[SchQNMguess]]]-SchQNMguess[l-2,0,FilterRules[{opts},Options[SchQNMguess]]];
		guess[[1]]=Abs[Re[guess[[1]]]]+I Im[guess[[1]]];
		guess[[2]]=n; (* Use preferred inversion *)
		guess[[3]]=0;guess[[4]]=0;guess[[5]]=0;
		Switch[s,
			   -2,Global`SchQNMTable[l,0]=guess,
			   -1,Global`SchQNMeTable[l,0]=guess,
				0,Global`SchQNMsTable[l,0]=guess
			   ];
		guess=2 SchQNMguess[l-1,1,FilterRules[{opts},Options[SchQNMguess]]]-SchQNMguess[l-2,1,FilterRules[{opts},Options[SchQNMguess]]];
		guess[[1]]=Abs[Re[guess[[1]]]]+I Im[guess[[1]]];
		guess[[2]]=n; (* Use preferred inversion *)
		guess[[3]]=0;guess[[4]]=0;guess[[5]]=0;
		Switch[s,
			   -2,Global`SchQNMTable[l,1]=guess;Global`SchQNMTable[l]=2,
			   -1,Global`SchQNMeTable[l,1]=guess;Global`SchQNMeTable[l]=2,
				0,Global`SchQNMsTable[l,1]=guess;Global`SchQNMsTable[l]=2
			   ];
		SchQNMguess[l,n,FilterRules[{opts},Options[SchQNMguess]]]
	]
]


Options[SchwarzschildQNM]=Union[{QNMSpinWeight->Null[],SchDebug->0,SchAnSol->False,
								RadialCFMinDepth->300,RadialCFDepth->1,QNMPrecision->24},
								Options[TestRadialCFConvergence3]];


SchwarzschildQNM[l_Integer,n_Integer,Nmax_Integer,opts:OptionsPattern[]] := 
Module[{s=OptionValue[QNMSpinWeight],rcfdepth=OptionValue[RadialCFDepth],guess, sol,Nrcf=Nmax,
		precision=OptionValue[QNMPrecision]},
	$MinPrecision = precision;
	guess = SchQNMguess[l, n];
	Print["guess = ", guess];
	If[rcfdepth>Nrcf,Nrcf=IntegerPart[rcfdepth]];
	sol=RadialMode[guess[[2]],s,0,0,l(l+1)-s(s+1),guess[[1]],Nrcf];
	Print[sol];
	{sol[[1]], n, Nmax, 0, 0}
]


SchwarzschildQNM[l_Integer,n_Integer,opts:OptionsPattern[]] :=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,debug=OptionValue[SchDebug],
		analytic=OptionValue[SchAnSol],rcfdepth=OptionValue[RadialCFDepth],
		Nrcf=300,testinv0,testinv1,Ninv,rcferr,rcfpow,nrpow,jacobianmatrix,Nradialnew,
		notconverged,guess,sol0,sol1,\[Epsilon]=-14,
		RCFmin=OptionValue[RadialCFMinDepth],
		precision=OptionValue[QNMPrecision]},
	$MinPrecision = precision;
	SpinWeightTable:=Switch[s,
				   -2,Global`SchQNMTable,
				   -1,Global`SchQNMeTable,
					0,Global`SchQNMsTable,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	If[SpinWeightTable[l]<n, 
		SchwarzschildQNM[l,n-1,FilterRules[{opts},Options[SchwarzschildQNM]]],
		Null[],
		(* No initial guesses yet for l *)
		SchQNMguess[l,n,FilterRules[{opts},Options[RadialLentzRoot]]];
	];
	If[rcfdepth>RCFmin,Nrcf=IntegerPart[rcfdepth]];
	If[rcfdepth<1 && rcfdepth>0,Nrcf=IntegerPart[Nrcf*Rationalize[rcfdepth]]];
	Nrcf=Max[Nrcf,RCFmin];
	If[SpinWeightTable[l]>=n,
		Print["Computing (l=", l, ",n=", n, ")"];
		If[Head[SpinWeightTable[l, n]]==SpinWeightTable,
			Print["No prior guess"];
			If[analytic,
				sol1=SchwarzschildQNM[l, n, Nrcf,FilterRules[{opts},Options[SchwarzschildQNM]]];
				Switch[s,
					   -2,Global`SchQNMTable[l,n]=sol1;Global`SchQNMTable[l]=n+1,
					   -1,Global`SchQNMeTable[l,n]=sol1;Global`SchQNMeTable[l]=n+1,
						0,Global`SchQNMsTable[l,n]=sol1;Global`SchQNMsTable[l]=n+1
					   ],
				sol1=SchQNMguess[l, n];
			];
			Ninv = n, (* Use preferred inversion *)
			(*False Case*)
			sol1 = SchQNMguess[l, n];
			Ninv = n, (* Use preferred inversion *)
			(*Unevaluated Case*)
			sol1 = SchQNMguess[l, n];
			Ninv = n; (* Use preferred inversion *)
		];
		notconverged = True;
		While[notconverged,
			sol1=RadialLentzRoot2[Ninv,s,0,0,l(l+1)-s(s+1),SetPrecision[sol1[[1]],precision],Nrcf,l+2,\[Epsilon],10^(-3),FilterRules[{opts},Options[RadialLentzRoot2]]];
			If[sol1[[1, 1]],
				jacobianmatrix=sol1[[1,3]];
				sol1=sol1[[2]];
				rcferr=TestRadialCFConvergence3[Ninv,s,0,0,l(l+1)-s(s+1),SetPrecision[sol1[[1]],precision],Nrcf,jacobianmatrix,\[Epsilon],RCFmin,l+2,1,FilterRules[{opts},Options[TestRadialCFConvergence3]]];
				Nradialnew=rcferr[[1]];
				If[debug>0,Print["RadialConverg : ",rcferr]];
				If[(Nradialnew>Nrcf),
					Nrcf=Nradialnew;
					If[debug>0,Print["Increase Nradial to ",Nrcf]];
					,
					notconverged=False;
				],
				Return[];
			];
		];
		Switch[s,
			   -2,Global`SchQNMTable[l,n]=sol1;Global`SchQNMTable[l]=Max[SpinWeightTable[l],n+1],
			   -1,Global`SchQNMeTable[l,n]=sol1;Global`SchQNMeTable[l]=Max[SpinWeightTable[l],n+1],
				0,Global`SchQNMsTable[l,n]=sol1;Global`SchQNMsTable[l]=Max[SpinWeightTable[l],n+1]
			   ];
		Print[sol1];
	];
]


If[!QNMDebug,Protect[SchQNMguess,SchwarzschildQNM]];


(* ::Subsection::Closed:: *)
(*Gravitational (s=-2)*)


(* Table for use in initializing initial-guess table

SchQNMTable[2]=2;
SchQNMTable[2,0]={0.3736716844180418`-0.08896231568893573` \[ImaginaryI],0,433,-14,2.7652268320563298`*^-17};
SchQNMTable[2,1]={0.34671099687916546`-0.2739148752912331` \[ImaginaryI],0,639,-14,2.680837030854711`*^-15};

SchQNMTable[3]=2;
SchQNMTable[3,0]={0.5994432884374897`-0.09270304794494767` \[ImaginaryI],0,252,-14,3.128874914131178`*^-16};
SchQNMTable[3,1]={0.5826438030332974`-0.2812981134350426` \[ImaginaryI],0,328,-14,2.495504809490339`*^-15};

Save["SchQNMTable.dat",SchQNMTable];

*)


(* ::Subsection:: *)
(*Electromagnetic (s=-1)*)


(* ::Subsection:: *)
(*Scalar (s = 0)*)


(* ::Section::Closed:: *)
(*Graphics*)


Options[QNMPlotSch]=Union[{QNMSpinWeight->Null[]},Options[ListPlot]];


QNMPlotSch[l_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,plist,mlist},
	SpinWeightTable:=Switch[s,
				   -2,Global`SchQNMTable,
				   -1,Global`SchQNMeTable,
					0,Global`SchQNMsTable,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	plist=Table[{Re[SpinWeightTable[l,n][[1]]],-Im[SpinWeightTable[l,n][[1]]]},{n,0,SpinWeightTable[l]-1}];
	mlist=Table[{-Re[SpinWeightTable[l,n][[1]]],-Im[SpinWeightTable[l,n][[1]]]},{n,0,SpinWeightTable[l]-1}];
	ListLinePlot[{plist,mlist},PlotMarkers->Automatic,FilterRules[{opts},Options[ListLinePlot]]]
]


Options[SchwarzschildOmega]={QNMSpinWeight->Null[]};


SchwarzschildOmega[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],KerrSEQ,not},
	If[Head[n]==Integer,not=n,Null[],not=n[[1]]];
	If[l==2 && not==8,Return[{0,2}]]; (* Special Cases *)
	If[l==3 && not==40,Return[{0,10}]]; (* Special Cases *)
	KerrSEQ=Switch[s,
				   -2,Global`KerrQNM[l,m,n],
				   -1,Global`KerrQNMe[l,m,n],
					0,Global`KerrQNMs[l,m,n],
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				  ];
	{Re[KerrSEQ[[1,2,1]]],-Im[KerrSEQ[[1,2,1]]]}
]


Options[KerrOmegaList]={QNMSpinWeight->Null[]};


KerrOmegaList[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:= 
Module[{s=OptionValue[QNMSpinWeight],KerrSEQ,Na},
	KerrSEQ=Switch[s,
				   -2,Global`KerrQNM[l,m,n],
				   -1,Global`KerrQNMe[l,m,n],
					0,Global`KerrQNMs[l,m,n],
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				  ];
	Na = Length[KerrSEQ];
	Table[{Re[KerrSEQ[[i,2,1]]],-Im[KerrSEQ[[i,2,1]]]},{i,1,Na}]
]


Options[KerrOmegaListS]={QNMSpinWeight->Null[]};


KerrOmegaListS[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:= 
Module[{s=OptionValue[QNMSpinWeight],KerrSEQ,Na,Nend,i,Slist={}},
	KerrSEQ=Switch[s,
				   -2,Global`KerrQNM[l,m,n],
				   -1,Global`KerrQNMe[l,m,n],
					0,Global`KerrQNMs[l,m,n],
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				  ];
	Na = Length[KerrSEQ];
	Nend = If[KerrSEQ[[Na,1]]<999999/1000000,Na,Na-1];
	For[i=1,i<=Nend,++i,
		If[Mod[KerrSEQ[[i,1]],1/20]==0,
			AppendTo[Slist,{Re[KerrSEQ[[i,2,1]]],-Im[KerrSEQ[[i,2,1]]]}]
		];
	];
	If[KerrSEQ[[Na,1]]>=999999/1000000,
		Append[Slist,{Re[KerrSEQ[[Na,2,1]]],-Im[KerrSEQ[[Na,2,1]]]}],
		Slist
	]
]


Options[KerrAList]={QNMSpinWeight->Null[]};


KerrAList[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:= 
Module[{s=OptionValue[QNMSpinWeight],KerrSEQ,Na},
	KerrSEQ=Switch[s,
				   -2,Global`KerrQNM[l,m,n],
				   -1,Global`KerrQNMe[l,m,n],
					0,Global`KerrQNMs[l,m,n],
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	Na = Length[KerrSEQ];
	Table[{Re[KerrSEQ[[i,3,1]]],Im[KerrSEQ[[i,3,1]]]},{i,1,Na}]
]


Options[KerrAListS]={QNMSpinWeight->Null[]};


KerrAListS[l_Integer,m_Integer,n_Integer|n_List,OptionsPattern[]]:= 
Module[{s=OptionValue[QNMSpinWeight],KerrSEQ,Na,Nend,i,Slist={}},
	KerrSEQ=Switch[s,
				   -2,Global`KerrQNM[l,m,n],
				   -1,Global`KerrQNMe[l,m,n],
					0,Global`KerrQNMs[l,m,n],
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	Na = Length[KerrSEQ];
	Nend = If[KerrSEQ[[Na,1]]<999999/1000000,Na,Na-1];
	For[i=1,i<=Nend,++i,
		If[Mod[KerrSEQ[[i,1]],1/20]==0,
			AppendTo[Slist,{Re[KerrSEQ[[i,3,1]]],Im[KerrSEQ[[i,3,1]]]}]
		];
	];
	If[KerrSEQ[[Na,1]]>=999999/1000000,
		Append[Slist,{Re[KerrSEQ[[Na,3,1]]],Im[KerrSEQ[[Na,3,1]]]}],
		Slist
	]
]


Options[QNMPlotOmega]=Union[{QNMSpinWeight->Null[],OTmultiple->{}},
							Options[ListLinePlot],Options[ListPlot]];


QNMPlotOmega[l_Integer,n_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],multiple=OptionValue[OTmultiple],
		legend=OptionValue[PlotLegends],autolegend,
		SpinWeightTable,KerrSEQ,
		mmodes={},multints,i,pos,m,linelist,pointlist,mainplot},
	SpinWeightTable:=Switch[s,
				   -2,Global`KerrQNM,
				   -1,Global`KerrQNMe,
					0,Global`KerrQNMs,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
	For[m=-l,m<=l,++m,
		If[MemberQ[multints,m],
			(* m,n is a multplet *)
			pos=Flatten[Position[multiple,{m,_}]][[1]];
			For[i=0,i<multiple[[pos,2]],++i,
				KerrSEQ:=SpinWeightTable[l,m,{n,i}];
				If[Head[KerrSEQ]==List,AppendTo[mmodes,{m,{n,i}}]]
			],
			(* Not a multplet *)
			KerrSEQ:=SpinWeightTable[l,m,n];
			If[Head[KerrSEQ]==List,AppendTo[mmodes,{m,n}]]
		]
	];
	autolegend=Table[If[Head[mmodes[[i,2]]]==List,Subscript[mmodes[[i,1]], mmodes[[i,2,2]]],Null,mmodes[[i,1]]],{i,1,Length[mmodes]}];
	If[legend==Automatic,legend=autolegend];
	If[Head[legend]==Placed,If[legend[[1]]==Automatic,legend=Placed[autolegend,legend[[2]]]]];
	linelist=KerrOmegaList[l,#[[1]],#[[2]],FilterRules[{opts},Options[KerrOmegaList]]]&/@  mmodes;
	pointlist=KerrOmegaListS[l,#[[1]],#[[2]],FilterRules[{opts},Options[KerrOmegaListS]]]&/@  mmodes;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,PlotLegends->legend,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic]],
		Show[mainplot]
	]
]


QNMPlotOmega[l_Integer,m_Integer,n_Integer|n_List,opts:OptionsPattern[]]:=
Module[{mmodes={},linelist,pointlist,mainplot},
	mmodes={m};
	linelist=KerrOmegaList[l,#,n,FilterRules[{opts},Options[KerrOmegaList]]]&/@  mmodes;
	pointlist=KerrOmegaListS[l,#,n,FilterRules[{opts},Options[KerrOmegaListS]]]&/@  mmodes;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic]],
		Show[mainplot]
	]
]


Options[QNMPlotA]=Union[{QNMSpinWeight->Null[],OTmultiple->{}},
						Options[ListLinePlot],Options[ListPlot]];


QNMPlotA[l_Integer,n_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],multiple=OptionValue[OTmultiple],
		legend=OptionValue[PlotLegends],autolegend,
		SpinWeightTable,KerrSEQ,
		mmodes={},multints,i,pos,m,linelist,pointlist,mainplot},
	SpinWeightTable:=Switch[s,
				   -2,Global`KerrQNM,
				   -1,Global`KerrQNMe,
					0,Global`KerrQNMs,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
	For[m=-l,m<=l,++m,
		If[MemberQ[multints,m],
			(* m,n is a multplet *)
			pos=Flatten[Position[multiple,{m,_}]][[1]];
			For[i=0,i<multiple[[pos,2]],++i,
				KerrSEQ:=SpinWeightTable[l,m,{n,i}];
				If[Head[KerrSEQ]==List,AppendTo[mmodes,{m,{n,i}}]]
			],
			(* Not a multplet *)
			KerrSEQ:=SpinWeightTable[l,m,n];
			If[Head[KerrSEQ]==List,AppendTo[mmodes,{m,n}]]
		]
	];
	autolegend=Table[If[Head[mmodes[[i,2]]]==List,Subscript[mmodes[[i,1]], mmodes[[i,2,2]]],Null,mmodes[[i,1]]],{i,1,Length[mmodes]}];
	If[legend==Automatic,legend=autolegend];
	If[Head[legend]==Placed,If[legend[[1]]==Automatic,legend=Placed[autolegend,legend[[2]]]]];
	linelist=KerrAList[l,#[[1]],#[[2]],FilterRules[{opts},Options[KerrAList]]]&/@  mmodes;
	pointlist=KerrAListS[l,#[[1]],#[[2]],FilterRules[{opts},Options[KerrAListS]]]&/@  mmodes;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,PlotLegends->legend,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic]],
		Show[mainplot]
	]
]


QNMPlotA[l_Integer,m_Integer,n_Integer|n_List,opts:OptionsPattern[]]:=
Module[{mmodes={},linelist,pointlist,mainplot},
	mmodes={m};
	linelist=KerrAList[l,#,n,FilterRules[{opts},Options[KerrAList]]]&/@  mmodes;
	pointlist=KerrAListS[l,#,n,FilterRules[{opts},Options[KerrAListS]]]&/@  mmodes;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic]],
		Show[mainplot]
	]
]


Options[QNMPlotOmegaTones]=Union[{QNMSpinWeight->Null[],OTmultiple->{}},
									Options[ListLinePlot],Options[ListPlot]];


QNMPlotOmegaTones[l_Integer,m_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],multiple=OptionValue[OTmultiple],
		legend=OptionValue[PlotLegends],autolegend,
		SpinWeightTable,KerrSEQ,
		ntones={},schtones={},n,multints,pos,i,j,linelist,pointlist,mainplot,Schlist,amin,amini},
	SpinWeightTable:=Switch[s,
				   -2,Global`KerrQNM,
				   -1,Global`KerrQNMe,
					0,Global`KerrQNMs,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
	For[n=0,n<=100,++n,
		If[MemberQ[multints,n],
			pos=Flatten[Position[multiple,{n,_}]][[1]];
			For[i=0,i<multiple[[pos,2]],++i,
				KerrSEQ:=SpinWeightTable[l,m,{n,i}];
				If[Head[KerrSEQ]==List,AppendTo[ntones,{n,i}]]
			]
			,
			KerrSEQ:=SpinWeightTable[l,m,n];
			If[Head[KerrSEQ]==List,AppendTo[ntones,n]]
		]
	];
	schtones=ntones;
	For[i=1,i<=Length[multints],++i,
		pos=Position[schtones,{multints[[i]],_}];
		amini=1;
		amin=SpinWeightTable[l,m,schtones[[pos[[1,1]]]]][[1,1]];
		For[j=2,j<=Length[pos],++j,
			If[SpinWeightTable[l,m,schtones[[pos[[j,1]]]]][[1,1]]<amin,
				amini=j;
				amin=SpinWeightTable[l,m,schtones[[pos[[j,1]]]]][[1,1]];
			];
		];
		pos=Drop[pos,amini];
		schtones=Delete[schtones,pos];
	];
	autolegend=Table[If[Head[ntones[[i]]]==List,Subscript[ntones[[i,1]], ntones[[i,2]]],Null,ntones[[i]]],{i,1,Length[ntones]}];
	If[legend==Automatic,legend=autolegend];
	If[Head[legend]==Placed,If[legend[[1]]==Automatic,legend=Placed[autolegend,legend[[2]]]]];
	linelist=KerrOmegaList[l,m,#,FilterRules[{opts},Options[KerrOmegaList]]]&/@ ntones;
	pointlist=KerrOmegaListS[l,m,#,FilterRules[{opts},Options[KerrOmegaListS]]]&/@ ntones;
	Schlist=SchwarzschildOmega[l,m,#,FilterRules[{opts},Options[SchwarzschildOmega]]]&/@ schtones;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,PlotLegends->legend,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic],
			ListLinePlot[Schlist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers,PlotStyle}]],PlotStyle->{Gray,Dashed},PlotRange->All],ImageSize->800],
		Show[mainplot,
			ListLinePlot[Schlist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers,PlotStyle}]],PlotStyle->{Gray,Dashed},PlotRange->All],ImageSize->800]
	]
]


Options[QNMPlotATones]=Union[{QNMSpinWeight->Null[],OTmultiple->{}},
								Options[ListLinePlot],Options[ListPlot]];


QNMPlotATones[l_Integer,m_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],multiple=OptionValue[OTmultiple],
		legend=OptionValue[PlotLegends],autolegend,
		SpinWeightTable,KerrSEQ,
		ntones={},n,multints,pos,i,linelist,pointlist,mainplot},
	SpinWeightTable:=Switch[s,
				   -2,Global`KerrQNM,
				   -1,Global`KerrQNMe,
					0,Global`KerrQNMs,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
	For[n=0,n<=100,++n,
		If[MemberQ[multints,n],
			pos=Flatten[Position[multiple,{n,_}]][[1]];
			For[i=0,i<multiple[[pos,2]],++i,
				KerrSEQ:=SpinWeightTable[l,m,{n,i}];
				If[Head[KerrSEQ]==List,AppendTo[ntones,{n,i}]]
			]
			,
			KerrSEQ:=SpinWeightTable[l,m,n];
			If[Head[KerrSEQ]==List,AppendTo[ntones,n]]
		]
	];
	autolegend=Table[If[Head[ntones[[i]]]==List,Subscript[ntones[[i,1]], ntones[[i,2]]],Null,ntones[[i]]],{i,1,Length[ntones]}];
	If[legend==Automatic,legend=autolegend];
	If[Head[legend]==Placed,If[legend[[1]]==Automatic,legend=Placed[autolegend,legend[[2]]]]];
	linelist=KerrAList[l,m,#,FilterRules[{opts},Options[KerrAList]]]&/@ ntones;
	pointlist=KerrAListS[l,m,#,FilterRules[{opts},Options[KerrAListS]]]&/@ ntones;
	mainplot=ListLinePlot[linelist,FilterRules[FilterRules[{opts},Options[ListLinePlot]],Except[{PlotLegends,PlotMarkers}]],PlotRange->All];
	If[Length[pointlist[[1]]]>0,
		Show[mainplot,
			ListPlot[pointlist,PlotLegends->legend,FilterRules[{opts},Options[ListPlot]],PlotRange->All,PlotMarkers->Automatic],ImageSize->800],
		Show[mainplot,ImageSize->800]
	]
]


(* eg: OTskip\[Rule]{5,8,{8,1},9}  and  OTmultiple\[Rule]{{8,2}}, assuming n=8 has two multiplets *)
(* can also replace 8 with {8,0}                                                        *)
Options[OvertoneLists]={OTskip->{},OTmultiple->{}};


OvertoneLists[n_Integer|n_List,overtones_List,OptionsPattern[]]:=
Module[{skip=DeleteDuplicates[OptionValue[OTskip]],
		multiple=DeleteDuplicates[OptionValue[OTmultiple]],otsort,
		fitot=n,fitotsave=n,fitind,nf=DeleteDuplicates[overtones],nfus,intnf={},
		i,multints,posnf,possk,pos,totlist,nprime,nrm,intlist},
	otsort[a_,b_]:=If[Head[a]==Head[b],OrderedQ[{a,b}],Null[],
					If[Head[a]==List,a[[1]]<=b,Null[],
					a<=b[[1]]]];
	nfus=nf;
(*
Print["Orig. nf = ",nf];
Print["Orig. nfus = ",nfus];
Print["Orig. multiple = ",multiple];
*)
	multints=Sort[DeleteDuplicates[Table[multiple[[i,1]],{i,Length[multiple]}]]];
	skip=Sort[skip,otsort];
	multiple=Sort[multiple,otsort];
	nf=Sort[nf,otsort];
(*
Print["1 : multints = ",multints];
Print["1 : skip = ",skip];
Print["1 : multiple = ",multiple];
Print["1 : nf = ",nf];
Print["1 : nfus = ",nfus];
*)
	(* Abort if any multiple overtone has length one in OTmultiple *)
	For[i=1,i<=Length[multiple],++i,
		If[multiple[[i,2]]<2,Print["Otmultiple ",multiple[[i]]," is not a multiple."];Abort[]]
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
(*
Print["2 : multints = ",multints];
Print["2 : skip = ",skip];
Print["2 : multiple = ",multiple];
Print["2 : nf = ",nf];
*)
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
(*
Print["3 : multints = ",multints];
Print["3 : skip = ",skip];
Print["3 : multiple = ",multiple];
Print["3 : nf = ",nf];
*)
	(* Get list of all overtones used in fit, removing specified overtones *)
	nrm=Flatten[Position[nf,#] & /@ skip,1];
	nf=Delete[nf,nrm];
	nrm=Flatten[Position[nfus,#] & /@ skip,1];
	nfus=Delete[nfus,nrm];
(*
Print["4 : multints = ",multints];
Print["4 : skip = ",skip];
Print["4 : multiple = ",multiple];
Print["4 : nf = ",nf];
Print["4 : nfus = ",nfus];
*)
	(* Get list of all remaining non-multiple overtones in list *)
	nrm=Drop[Position[Head /@ nf,List],1];
	intlist=Delete[nf,nrm];
	(* Get list of all possible integer overtones *)
	totlist=Table[i,{i,Min[intlist,multints],Max[intlist,multints]}];
	nrm=Flatten[Position[totlist,#] & /@ multints,1];
	totlist=Delete[totlist,nrm];
(*
Print["5 : intlist = ",intlist];
Print["5 : totlist = ",totlist];
*)
	(* Add list of all possible multiple overtones *)
	For[i=1,i<=Length[multiple],++i,
		totlist=Join[totlist,Table[{multiple[[i,1]],j},{j,0,multiple[[i,2]]-1}]];
	];
	(* remove specified overtones *)
	nrm=Flatten[Position[totlist,#] & /@ skip,1];
	totlist=Sort[Delete[totlist,nrm],otsort];
	nprime=Table[Position[totlist,nf[[i]]][[1,1]]-1,{i,1,Length[nf]}];
(*
Print["6 : totlist = ",totlist];
Print["6 :  nprime= ",nprime];
*)
	(* return integer overtones in nf *)
	For[i=1,i<=Length[intnf],++i,
		pos=Flatten[Position[nf,{intnf[[i]],0}]];
		If[Length[pos]>0,nf[[pos[[1]]]]=intnf[[i]]];
	];
	pos=Flatten[Position[totlist,fitot]];
	If[Length[pos]>0,fitind=pos[[1]]-1,fitind=-1];
(*
Print["7 : nf = ",nf];
Print["7 : nfus = ",nfus];
Print["7 :  nprime= ",nprime];
*)
	{fitotsave,fitind,nfus,nprime}
]


Options[QNMPlotAccumulation\[Omega]]={QNMSpinWeight->Null[],OTskip->{},OTmultiple->{}};


QNMPlotAccumulation\[Omega][l_Integer,m_Integer,n_List,Nv_Integer,
					a0_Real|a0_Rational|a0_Integer,
					opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,KerrSEQ,
		nf,nprime,fitlists,i,j,qnmend,qnmRe,qnmIm,ReList={},ImList={},
		Refit,ReLOF,ReParams,ReFitData,Imfit,ImLOF,ImParams,ImFitData,plot,label},
	SpinWeightTable:=Switch[s,
				   -2,Global`KerrQNM,
				   -1,Global`KerrQNMe,
					0,Global`KerrQNMs,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	fitlists=OvertoneLists[-1,n,FilterRules[{opts},Options[OvertoneLists]]];
	nf=fitlists[[3]];
	nprime=fitlists[[4]];
	For[i=1,i<=Length[nf],++i,
		KerrSEQ:=SpinWeightTable[l,m,nf[[i]]];
		qnmend=Take[KerrSEQ,-Nv];
		qnmRe={};
		qnmIm={};
		For[j=1,j<=Nv,++j,
			If[1-qnmend[[j,1]]<=a0,
				AppendTo[qnmRe,{1-qnmend[[j,1]],Re[qnmend[[j,2,1]]]}];
				AppendTo[qnmIm,{1-qnmend[[j,1]],Im[qnmend[[j,2,1]]]}];
			];
		];
		AppendTo[ReList,qnmRe];
		AppendTo[ImList,qnmIm];
	];
	ReFitData=Flatten[Table[{nprime[[i]],ReList[[i,j,1]],ReList[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ReList[[i]]]}],1];
	ImFitData=Flatten[Table[{nprime[[i]],ImList[[i,j,1]],ImList[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ImList[[i]]]}],1];
	If[nf[[1]]==0,
		Refit=NonlinearModelFit[ReFitData,m/2-\[Alpha]1 Sqrt[\[Epsilon]1/2]+(\[Alpha]2+\[Alpha]3 n1)\[Epsilon]1,{\[Alpha]1,\[Alpha]2,\[Alpha]3},{n1,\[Epsilon]1}];
		Imfit=NonlinearModelFit[ImFitData,-(n1+1/2)(Sqrt[\[Epsilon]1/2]-\[Alpha]4 \[Epsilon]1),{\[Alpha]4},{n1,\[Epsilon]1}],
		Refit=NonlinearModelFit[ReFitData,m/2+(\[Alpha]1+\[Alpha]2 n1)\[Epsilon]1,{\[Alpha]1,\[Alpha]2},{n1,\[Epsilon]1}];
		Imfit=NonlinearModelFit[ImFitData,-(\[Alpha]3+n1+1/2)Sqrt[\[Epsilon]1/2]+(\[Alpha]4+\[Alpha]5 n1)\[Epsilon]1,{\[Alpha]3,\[Alpha]4,\[Alpha]5},{n1,\[Epsilon]1}]
	];
	(* Modify nf for multiplet labels *)
	For[i=1,i<=Length[nf],++i,
		If[Length[nf[[i]]]==2,nf[[i]]=Subscript[nf[[i,1]], nf[[i,2]]]];
	];
	label=DisplayForm[RowBox[{"l=",l," m=",m," n=",nf}]];
	ReParams=Refit["BestFitParameters"];
	ImParams=Imfit["BestFitParameters"];
	If[nf[[1]]==0,
		ReLOF[nf_,\[Epsilon]f_]:=m/2-\[Alpha]1 Sqrt[\[Epsilon]f/2]/.ReParams;
		ImLOF[nf_,\[Epsilon]f_]:=-(nf+1/2)Sqrt[\[Epsilon]f/2]/.ImParams,
		ReLOF[nf_,\[Epsilon]f_]:=m/2+(\[Alpha]1+\[Alpha]2 nf)\[Epsilon]f;
		ImLOF[nf_,\[Epsilon]f_]:=-1(\[Alpha]3+nf+1/2)Sqrt[\[Epsilon]f/2]/.ImParams
	];
	ReList=Map[Function[x,x-{0,m/2}],ReList,{2}];
	plot=Show[ListPlot[ReList,PlotRange->{{0,a0},All},PlotMarkers->Automatic,AxesOrigin->{Automatic,0},TicksStyle->Directive[14],AxesLabel->{Style["1-\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],Style["Re(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))"<>ToString[NumberForm[N[-m/2],3]],16]},PlotLabel->Style[label,16]],Plot[ReLOF[#,\[Epsilon]]&/@nprime-m/2,{\[Epsilon],0,a0},PlotStyle->{Red,Dashed}],Plot[Refit[#,\[Epsilon]]&/@nprime-m/2,{\[Epsilon],0,a0}],ImageSize->800];
	Print[plot];
	plot=Show[ListPlot[ImList,PlotRange->{{0,a0},All},PlotMarkers->Automatic,AxesOrigin->{Automatic,0},TicksStyle->Directive[14],AxesLabel->{Style["1-\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],Style["Im(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))",16]},PlotLabel->Style[label,16]],Plot[ImLOF[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0},PlotStyle->{Red,Dashed}],Plot[Imfit[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0}],ImageSize->800];
	Print[plot];
	Print["Re(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))= "Normal[Refit]];
	Print["Im(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\))= "Normal[Imfit]];
	Off[SetPrecision::"precsm"];
	Off[N::"precsm"];
	Print["\n  Re(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\)) :\n",Refit["ParameterConfidenceIntervalTable"],"\n  Im(\!\(\*OverscriptBox[\(\[Omega]\), \(_\)]\)) :\n",Imfit["ParameterConfidenceIntervalTable"]];
	On[SetPrecision::"precsm"];
	On[N::"precsm"];
]


Options[QNMPlotAccumulationAlm]={QNMSpinWeight->Null[],OTskip->{},OTmultiple->{}};


QNMPlotAccumulationAlm[l_Integer,m_Integer,n_List,Nv_Integer,
					a0_Real|a0_Rational|a0_Integer,
					opts:OptionsPattern[]]:=
Module[{s=OptionValue[QNMSpinWeight],SpinWeightTable,KerrSEQ,
		nf=n,nprime,fitlists,i,j,qnmend,qnmRe,qnmIm,ReList={},ImList={},
		Refit,ReLOF,ReParams,ReFitData,Imfit,ImLOF,ImParams,ImFitData,
		A0,\[Delta],intercept,plot,label},
	SpinWeightTable:=Switch[s,
				   -2,Global`KerrQNM,
				   -1,Global`KerrQNMe,
					0,Global`KerrQNMs,
				   _,Print["Invalid QNMSpinWeight"];Abort[]
				   ];
	fitlists=OvertoneLists[-1,n,FilterRules[{opts},Options[OvertoneLists]]];
	nf=fitlists[[3]];
	nprime=fitlists[[4]];
	For[i=1,i<=Length[nf],++i,
		KerrSEQ:=SpinWeightTable[l,m,nf[[i]]];
		qnmend=Take[KerrSEQ,-Nv];
		qnmRe={};
		qnmIm={};
		For[j=1,j<=Nv,++j,
			If[1-qnmend[[j,1]]<=a0,
				AppendTo[qnmRe,{1-qnmend[[j,1]],Re[qnmend[[j,3,1]]]}];
				AppendTo[qnmIm,{1-qnmend[[j,1]],Im[qnmend[[j,3,1]]]}];
			];
		];
		AppendTo[ReList,qnmRe];
		AppendTo[ImList,qnmIm];
	];
	ReFitData=Flatten[Table[{nprime[[i]],ReList[[i,j,1]],ReList[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ReList[[i]]]}],1];
	ImFitData=Flatten[Table[{nprime[[i]],ImList[[i,j,1]],ImList[[i,j,2]]},{i,1,Length[nprime]},{j,1,Length[ImList[[i]]]}],1];
	If[nf[[1]]==0,
		Refit=NonlinearModelFit[ReFitData,l(l+1)-s(s+1)+\[Beta]1+\[Beta]2 Sqrt[\[Epsilon]1/2]+(\[Beta]3+\[Beta]4 n1+\[Beta]5 n1^2)\[Epsilon]1,{\[Beta]1,\[Beta]2,\[Beta]3,\[Beta]4,\[Beta]5},{n1,\[Epsilon]1}];
		Imfit=NonlinearModelFit[ImFitData,(n1+1/2)(\[Beta]6 Sqrt[\[Epsilon]1/2]+\[Beta]7 \[Epsilon]1),{\[Beta]6,\[Beta]7},{n1,\[Epsilon]1}],
		Refit=NonlinearModelFit[ReFitData,l(l+1)-s(s+1)+\[Beta]1+(\[Beta]2+\[Beta]3 n1+\[Beta]4 n1^2)\[Epsilon]1,{\[Beta]1,\[Beta]2,\[Beta]3,\[Beta]4},{n1,\[Epsilon]1}];
		Imfit=NonlinearModelFit[ImFitData,(\[Beta]5+n1+1/2)\[Beta]6 Sqrt[\[Epsilon]1/2]+(\[Beta]7+\[Beta]8 n1)\[Epsilon]1,{\[Beta]5,\[Beta]6,\[Beta]7,\[Beta]8},{n1,\[Epsilon]1}]
	];
	A0=Refit["BestFitParameters"];
	intercept=l(l+1)-2+\[Beta]1/.A0;
	\[Delta]=Sqrt[7/4 m^2 - (-2+1/2)^2-(l(l+1)-2+\[Beta]1)/.A0];
	Print["\[Delta] = ",\[Delta]];
	(* Modify nf for multiplet labels *)
	For[i=1,i<=Length[nf],++i,
		If[Length[nf[[i]]]==2,nf[[i]]=Subscript[nf[[i,1]], nf[[i,2]]]];
	];
	label=DisplayForm[RowBox[{"l=",l," m=",m," n=",nf}]];
	ReParams=Refit["BestFitParameters"];
	ImParams=Imfit["BestFitParameters"];
	If[nf[[1]]==0,
		ReLOF[nf_,\[Epsilon]f_]:=l(l+1)-2+\[Beta]1+\[Beta]2 Sqrt[\[Epsilon]f/2]/.ReParams;
		ImLOF[nf_,\[Epsilon]f_]:=(nf+1/2)\[Beta]6 Sqrt[\[Epsilon]f/2]/.ImParams,
		ReLOF[nf_,\[Epsilon]f_]:=l(l+1)-2+\[Beta]1+(\[Beta]2+\[Beta]3 nf+\[Beta]4 nf^2)\[Epsilon]f;
		ImLOF[nf_,\[Epsilon]f_]:=(\[Beta]5+nf+1/2)\[Beta]6 Sqrt[\[Epsilon]f/2]/.ImParams
	];
	ReList=Map[Function[x,x-{0,intercept}],ReList,{2}];
	plot=Show[ListPlot[ReList,PlotRange->{{0,a0},All},PlotMarkers->Automatic,AxesOrigin->{Automatic,0},TicksStyle->Directive[14],AxesLabel->{Style["1-\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],Style["Re(\!\(\*SubscriptBox[\(A\), \(lm\)]\))"<>ToString[NumberForm[-intercept,3]],16]},PlotLabel->Style[label,16]],Plot[ReLOF[#,\[Epsilon]]&/@nprime-intercept,{\[Epsilon],0,a0},PlotStyle->{Red,Dashed}],Plot[Refit[#,\[Epsilon]]&/@nprime-intercept,{\[Epsilon],0,a0}],ImageSize->800];
	Print[plot];
	plot=Show[ListPlot[ImList,PlotRange->{{0,a0},All},PlotMarkers->Automatic,AxesOrigin->{Automatic,0},TicksStyle->Directive[14],AxesLabel->{Style["1-\!\(\*OverscriptBox[\(a\), \(_\)]\)",16],Style["Im(\!\(\*SubscriptBox[\(A\), \(lm\)]\))",16]},PlotLabel->Style[label,16]],Plot[ImLOF[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0},PlotStyle->{Red,Dashed}],Plot[Imfit[#,\[Epsilon]]&/@nprime,{\[Epsilon],0,a0}],ImageSize->800];
	Print[plot];
	Print["Re(\!\(\*SubscriptBox[\(A\), \(lm\)]\))= "Normal[Refit]];
	Print["Im(\!\(\*SubscriptBox[\(A\), \(lm\)]\))= "Normal[Imfit]];
	Off[SetPrecision::"precsm"];
	Off[N::"precsm"];
	Print["\n  Re(\!\(\*SubscriptBox[\(A\), \(lm\)]\)) :\n",Refit["ParameterConfidenceIntervalTable"],"\n  Im(\!\(\*SubscriptBox[\(A\), \(lm\)]\)) :\n",Imfit["ParameterConfidenceIntervalTable"]];
	On[SetPrecision::"precsm"];
	On[N::"precsm"];
]


QNMColor[i_Integer]:=ColorData[1,i];
QNMMark[i_Integer]:={Global`\[FilledCircle],Global`\[FilledSmallSquare],Global`\[FilledDiamond],Global`\[FilledUpTriangle],Global`\[FilledDownTriangle],
					Global`\[EmptyCircle],Global`\[EmptySquare],Global`\[EmptyDiamond],Global`\[EmptyUpTriangle],Global`\[EmptyDownTriangle]}[[Mod[i-1,10]+1]];


If[!QNMDebug,Protect[QNMPlotSch,SchwarzschildOmega,KerrOmegaList,KerrOmegaListS,
					KerrAList,KerrAListS,QNMPlotOmega,QNMPlotA,QNMPlotOmegaTones,
					QNMPlotATones,OvertoneLists,QNMPlotAccumulation\[Omega],
					QNMPlotAccumulationAlm,QNMColor,QNMMark]];


(* ::Section::Closed:: *)
(*End of KerrQNM Package*)


End[] (* `Private` *)


EndPackage[]
