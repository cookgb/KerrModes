# KerrModes
Set of *Mathematica Paclets* for computing various modes of the Kerr Geometry.
The following paclets are based directly on the original Mathematica packages used in [Cook and Zalutsky (2014)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.90.124021), [Cook and Zalutskiy (2016)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.104074), [Cook, Annichiarico, and Vickers (2019)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.024008), and subsequent works.  Please see these references for details and cite these references when publishing results based on these packages.

---
## Paclets

#### SWSpheroidal\`
Mathematica Paclet to solve the spin-weighted spheroical function equation (also known as the Teukolsky angular equation).

#### KerrModes\`
Mathematica Paclet to solve the coupled radial and angular Teukolsky equations.  (This paclet is not intended to be loaded directly, but provides all of the common functionality needed by KerrQNM\`, KerrTTML\`, and KerrTTMR\`.)

#### KerrQNM\`
Mathematica Paclet to find Quasi-Normal Modes of the Kerr geometry.

#### KerrTTML\`
Mathematica Paclet to find Left Total-Transmission Modes of the Kerr geometry.

#### KerrTTMR\`
Mathematica Paclet to find Right Total-Transmission Modes of the Kerr geometry.

---
## Support notebooks

#### CreatePaclets.nb
Mathematica notebook to compile each package into a paclet.

#### InstallPaclets.nb
Mathematica notebook to install each paclet into a local Mathematica environment.

---
## Mathematica-Style Package Documentation

All public functions in each paclet (and some private functions as well) have full Mathematica-type documentation.  Guides and tutorials are also provided within the Mathematica documentation.  **However**, the current implementation of paclets does not fully integrate in with the Mathematica's feature of "hovering over" a function and asking for help.  
1) Help can be obtained by searching the Wolfram Documentation in a Mathematica environment where the KerrModes paclets have been installed.
2) Help can also be obtained by using the ?FunctionName method for getting information about a function.  In this case, clicking on the blue circle with an "i" in the middle will bring up the full Mathematica documentation for the function.  This documentation will include links to guides and tutorials.
