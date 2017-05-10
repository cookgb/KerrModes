(* ::Package:: *)

(* ::Title:: *)
(*QuasiNormal Modes of Kerr*)


(* ::Section::Closed:: *)
<<<<<<< HEAD
(*Begin KerrQNM Package*)
||||||| merged common ancestors
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


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*Sequencers*)


KerrQNMSequenceB::usage=""


KerrQNMRefineSequenceB::usage=""


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

=======
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


KerrQNMSequenceB::usage=
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
	"in the range 0,1,...,(Nmult-1), with 'Nmult' the number of sequences with the same overtone "


KerrQNMRefineSequenceB::usage=
"KerrQNMSequence[l,m,n,\[Epsilon]] computes a sequence of Quasi-Normal Mode solutions "<>
	"for overtone n of mode (l,m).  The solutions are computed to an absolute "<>
	"accuracty of \!\(\*SuperscriptBox[\(10\), \(\[Epsilon]\)]\).  The sequence is "<>
	"parameterized by increasing values of the dimensionless angular momentum "<>
	"'a' starting at a=0 (or the largest value of 'a' already computed) up to "<>
	"(but not including) a=1.\n \n "


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

>>>>>>> 8ae13c6e42a08f8b15e7c1e7c9aa8d23e59edb44


BeginPackage["KerrQNM`",{"KerrModes`"}]


Unprotect[KerrQNMDebug];
KerrQNMDebug=False; (* Set this to True to allow reloading of Package with changes *)
If[KerrQNMDebug,Unprotect["KerrQNM`*"];Unprotect["KerrQNM`Private`*"]];
Protect[KerrQNMDebug];


(* ::Section::Closed:: *)
(*Documentation of External Functions in KerrModes Namespace*)


SetSpinWeight::usage=
	"SetSpinWeight[s] sets the value of the spin-weight used in all subsequent "<>
	"QNM computations:\n"<>
	"\t s=-2 : Gravitational perturbations\n"<>
	"\t s=-1 : Electro-Magnetic perturbations\n"<>
	"\t s= 0 : Scalar perturbations."


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


Clear[rp,rm,\[Sigma]p,C0,C1,C2,C3,C4,D0,D1,D2,D3,D4];
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
(*Set SpinWeight and Data-Variable Names*)


SetSpinWeight[s_Integer]:=
Module[{},
	SetSpinWeight::spinweight="Invalid QNM Spin Weight : `1`";
	Switch[s,
		   -2,modeName:=Global`KerrQNM; SchTable:=Global`SchQNMTable,
		   -1,modeName:=Global`KerrQNMe; SchTable:=Global`SchQNMeTable,
		    0,modeName:=Global`KerrQNMs; SchTable:=Global`SchQNMsTable,
			_,Message[SetSpinWeight::spinweight,s];Abort[]
		  ];
	SetOptions[KerrQNM`SchwarzschildQNM,SpinWeight->s];
	Print["All KerrMode routines (QNM) set for Spin-Weight s = ",s];
]


If[!KerrModeDebug,Protect[SetSpinWeight]];


End[] (* KerrModes`Private` *)


(* ::Section::Closed:: *)
(*Documentation of External Functions in KerrQNM Namespace*)


SchwarzschildQNM::usage=
	"SchwarzschildQNM[l,n] computes the Quasi-Normal Mode solutions for overtone n of mode l.  "<>
	"The mode is computed to an accuracy of \!\(\*SuperscriptBox[\(10\), \(-14\)]\).  For given 'l', if solutions with overtones "<>
	"(n-1) and (n-2) have not been computed, then the routine is recersively called for overtone "<>
	"'(n-1).  If no solutions exist for mode l, then the first two overtones of (l-1) and (l-2) "<>
	"are used to extrapolate initial guesses for these modes.\n\n"<>
	"Options:\n"<>
	"\t SpinWeight\[Rule]Null : -2,-1,0\n"<>
	"\t\t The spin weight must be set via a call to SetSpinWeight before any KerrQNM\n"<>
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


PlotSchQNM::usage=
	"PlotSchQNM[l] plots both the \"positive\" and \"negative\" frequency QNMs.  "<>
	"By default, the gravitational modes are plotted, but the PlotSpinWeight option "<>
	"can be set to change this.\n\n"<>
	"Options:\n"<>
	"\t PlotSpinWeight->-2 : -2,-1,0\n\n"<>
	"PlotSchQNM also take all of the options available to ListPlot.\n"


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[PlotSpinWeight];


Begin["`Private`"]


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


Options[PlotSchQNM]=Union[{PlotSpinWeight->-2},Options[KerrModes`Private`PlotSchModes]];


PlotSchQNM[l_Integer,opts:OptionsPattern[]]:=
Module[{s=OptionValue[PlotSpinWeight],ptable},
	PlotSchQNM::spinweight="Invalid QNM Spin Weight : `1`";
	Switch[s,
		   -2,ptable:=Global`SchQNMTable,
		   -1,ptable:=Global`SchQNMeTable,
		    0,ptable:=Global`SchQNMsTable,
			_,Message[PlotSchQNM::spinweight,s];Abort[]
		  ];
	KerrModes`Private`PlotSchModes[l,PlotTable->ptable,
									FilterRules[{opts},Options[KerrModes`Private`PlotSchModes]]]
]


(* ::Section::Closed:: *)
(*End of KerrQNM Package*)


End[] (* `Private` *)


EndPackage[]
