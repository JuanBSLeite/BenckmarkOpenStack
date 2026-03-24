# Laura++ release notes

===
## Laura++ v3r7

24th February 2026 Thomas Latham
* Update compiler versions for dev3 and dev4 builds in CI

21st August 2025 Thomas Latham
* Fix missing includes of TMath

17th July 2025 Thomas Latham
* CMake fixes:
  - Fix missing PUBLIC keyword in `target_link_libraries`
  - Increase minimum CMake version to 3.13

17th July 2025 Thomas Latham
* Modernise LauRooFitTask to work with newer ROOT versions
* Implement toy generation in LauRooFitTask

11th July 2025 Thomas Latham
* Update GenFitKpipi histograms URL

3rd October 2024 Thomas Latham
* Doxygen fix

===
## Laura++ v3r6

5th June 2024 Thomas Latham
* Apply formatting to source files
* Enable formatting check in CI
* See https://phab.hepforge.org/T111 for more details

5th June 2024 Thomas Latham
* Revamp how parameters are passed to the fit model
  - Better handling of clones
  - Reduces code in the derived classes
* Fix minor memory leak in LauSimFitTask

11th January 2024 Thomas Latham
* Fix handling of multi-dimensional constraints in sim fit coordinator
* Store per-experiment values of constraint means in fit results ntuple

10th January 2024 Thomas Latham
* Remove all uses of deprecated ROOT ClassImp macro
* Take advantage of improvements in ROOT dictionary generation
  - Consequently bump the minimum ROOT version to 6.18
* See https://phab.hepforge.org/T213 for more details

10th January 2024 Thomas Latham
* Fix python example

14th December 2023 Thomas Latham
* Fix missing includes / fwd declarations
* Fix issues with Doxygen documentation
* Minor formatting fixes for examples
* Minor improvements to CMake config

7th December 2023 Thomas Latham
* Create clang-format (version 16) style file
* Create scripts to:
  - Apply clang-format to all source files
  - Check and fix indentation of Doxygen comment blocks
  - Fix common issues in formatting of CMake and Markdown files
* See https://phab.hepforge.org/T111 for more details

20th November 2023 Thomas Latham
* Allow increase of C++ standard via environment variable
* See https://phab.hepforge.org/D105 for more details

16th November 2023 Thomas Latham
* Add functions to fit models that return:
  - signal DP amplitudes
  - signal and background DP likelihood(s)
* Add example to plot histograms of likelihood values
* Caveat: SCF signal is not currently handled - an appropriate warning is printed
* See https://phab.hepforge.org/T227 for more details

9th August 2023 Thomas Latham
* Update CI config to use el9 for LCG dev builds

8th August 2023 Andy Morris
* Suppress 'up-to-date' messages in install step

22nd June 2023 Thomas Latham
* Added functionality to properly manage Gaussian constraints in toy generation
  - see https://phab.hepforge.org/T45

9th June 2023 Mark Whitehead
* Added functionality to include n-dimensional Gaussian constraints
  - see https://phab.hepforge.org/T214

10th February 2023 Thomas Latham
* Store the sqDP jacobian separately in weightEvents routines
  - Previously the weight was multiplied by the jacobian
  - Now the jacobian is stored in a separate branch
  - New info messages give instructions on what to do depending on how the
    input sample was generated

21st November 2022 Andy Morris
* Use Minuit to automatically determine ASq max for signal in LauIsobarDynamics
  - see https://phab.hepforge.org/D85

28th June 2022 Thomas Latham
* Improvements to Doxygen config

31st May 2022 Thomas Latham
* Add examples for easily generating toy samples uniformly in any DP or SqDP
* Very minor updates to related QuasiFlatSqDalitz examples

30th November 2021 Thomas Latham
* Update CI config to:
  - Simplify changing compiler version
  - Disallow dev build failure in scheduled pipelines

7th September 2021 Thomas Latham
* Update CI to test future ROOT versions in dev3/dev4 LCG builds

1st July 2021 Thomas Latham
* Add version of QuasiFlatSqDalitz where particle masses can be supplied on command line

28th June 2021 Thomas Latham
* Fix for compilation with ROOT 6.24

18th May 2021 Mark Whitehead
* Add new MIPW example
  - see https://phab.hepforge.org/T147

26th March 2021 Thomas Latham
* Fix bug in setting for linear interpolation in Lau1DHistPdf
  - see https://phab.hepforge.org/T142

9th March 2021 Mark Whitehead
* Update particle properties to PDG 2020
  - see https://phab.hepforge.org/T140

5th March 2021 Dan Johnson
* Alterations to allow Blatt-Weisskopf factor for parent to be created in a model with only a K-matrix
  - see https://phab.hepforge.org/D52

4th February 2021 Dan Johnson
* Extend the K-matrix implementation to handle non-zero spin
  - see https://phab.hepforge.org/T135

1st February 2021 Dan Johnson
* Allow floating of parameters in the K-matrix
  - see https://phab.hepforge.org/T59

2nd December 2020 Thomas Latham
* Fix LauFormulaPar to follow change in behaviour of TFormula
  - see https://phab.hepforge.org/T129

27th November 2020 Dan Johnson
* Allow slope of NR exponential model to vary negative - improves fit stability in low-statistics fits
  - see https://phab.hepforge.org/T128

17th September 2020 Mark Whitehead
* Begin updates to use inclusive language. Simultaneous fits now handled by Coordinator and Tasks
  - see https://phab.hepforge.org/T112

19th August 2020 Thomas Latham
* Remove explicit normalisation of numerator Blatt-Weisskopf factors
  - See https://phab.hepforge.org/T93

22nd May 2020 Thomas Latham
* Fix uninitialised variable (related to rho-omega mixing) in LauIsobarDynamics

12th December 2019 Thomas Latham & Daniel Johnson
* Fix issue with generation of events for categories with small expected yield
  - See https://phab.hepforge.org/T76

21st November 2019 Thomas Latham
* Add QuasiFlatSqDalitz example, which generates toy according to EvtGen's FLATSQDALITZ model

6th - 20th November 2019 Thomas Latham
* Adopt CMake as the build system
  - See https://phab.hepforge.org/T33

===
## Laura++ v3r5

6th August 2019 Thomas Latham
* Add some extra charmonium states to list of known resonances
* Minor modifications to code for importing parameter values from file

17th May 2019 Thomas Latham
* Fix licences and Doxygen in pole, nonresonant form factors, and rescattering line shapes
* Add journal reference for pole line shape
* Fix class names and add Doxygen to LauCalcChiSq, LauResultsExtractor, LauMergeDataFiles
* Update Doxyfile from more recent Doxygen version
* Make consistent the use of ClassDef and ClassImp everywhere
* Add functions to return the "sanitised" names of parent and daughters

8th January 2019 Thomas Latham
* Workaround for issue with splitting toy MC files when generating large number of experiments

5th December 2018 Thomas Latham
* Move sources and headers for the utilities library into the main library
  - makes examples dir executables only (in preparation for CMake build)

22nd, 23rd July and 4th August 2018 Juan Otalora
* Add pole, nonresonant form factors, and rescattering line shapes

17th April 2018 Daniel O'Hanlon
* Fix bug in rho-omega mixing fit-fraction calculation: rho fit-fraction was not normalised to the total DP rate

23rd March 2018 Thomas Latham
* Some fixes in LauRhoOmegaMix that propagate the following settings to the subcomponents:
  - spinType_
  - flipHelicity_
  - ignoreMomenta_
  - ignoreSpin_
  - ignoreBarrierScaling_

23rd February 2018 Thomas Latham
* Add section on structure of package to README file (requested by CPC technial editor)
* Add copyright notice to files in test directory

21st February 2018 Thomas Latham
* Improve comments in functions calculating helicity angles in LauKinematics

19th February 2018 Daniel O'Hanlon
* Fix bug in LauCPFitModel introduced in code for importing parameters (25/01/18)
  - When parameters are asked to be randomised via useRandomInitFitPars, those
    parameters that were already fixed in the model, in addition to those that
    are fixed on import, were being freed.

19th February 2018 Daniel O'Hanlon
* Bug fixes in rho/omega fit fractions calculation

25th January 2018 Daniel O'Hanlon
* Add feature for importing parameters into LauCPFitModels from previous fit output files.

23rd January 2018 Daniel O'Hanlon
* Calculate separate rho and omega fit-fractions for LauRhoOmegaMix, when turned on with calculateRhoOmegaFitFractions in LauIsobarDynamics.

===
## Laura++ v3r4

16th January 2018 Thomas Latham
* Update licence for all files to the Apache Software License Version 2.0

14th December 2017 Thomas Latham
* Correct LauFlatteRes to use the exact formulae from the papers that report the default parameter values for each supported resonance
  - Deals correctly now with the cases where the m0 is absorbed into the couplings
  - Now uses the correct value of g1 for the a_0(980) (needed to be squared)
  - Also sets the mass values to those reported in those papers
  - Added printout to specify what is being done in each case
* Improve the consistency between the weightEvents functions in LauSimpleFitModel and LauCPFitModel
  - Also remove the unecessary division by ASqMax in LauIsobarDynamics::getEventWeight

4th December 2017 Thomas Latham
* Allow background yields to be specified as any LauAbsRValue (i.e. LauParameter or now also LauFormulaPar)

30th November 2017 Thomas Latham
* In LauFitDataTree::readExperimentData add a check that the tree contains the branch "iExpt" and if not:
  - If the requested experiment is 0 then read all data in the tree (and print a warning)
  - Otherwise print an error message and return

29th November 2017 Thomas Latham
* Improve error messages in LauCPFitModel::weightEvents

===
## Laura++ v3r3

23rd November 2017 Thomas Latham
* Add an example written as a python script: GenFit3pi.py
* Various improvements to the other example code

23rd November 2017 Thomas Latham
* Fix bug in the LauAbsResonance constructor used by the K-matrix
  - "resonance" charge is now set to the sum of the daughter charges

1st November 2017 Thomas Latham
* Fix bug in LauFlatteRes
  - m_0 factor multiplying the width was, if useAdlerTerm_ is true, replaced by
    f_Adler rather than being multiplied by it
  - only affected case where using LauFlatteRes to describe K_0*(1430)

===
## Laura++ v3r2

18th October 2017 Thomas Latham
* Modify LauDaughters::testDPSymmetry to:
  - produce an ERROR message and exit if the symmetric DP does not have its daughters arranged correctly
  - detect also flavour-conjugate DPs and set a corresponding flag that can be retrieved via LauDaughters::gotFlavourConjugateDP
    - in this case if the daughters are sub-optimally arranged print a WARNING
* Modify LauIsobarDynamics::calcDPNormalisationScheme so that any narrow resonance regions are symmetrised appropriately if the DP is fully-symmetric, symmetric or is flavour-conjugate (and in this last case, symmetrisation of the integration has been forced via a new boolean function)

11th October 2017 Thomas Latham
* Allow the user to specify the randomiser to be used to randomise the initial values of the isobar parameters
  - simply use the new static function LauAbsCoeffSet::setRandomiser
  - if not specified it defaults to the current use of LauRandom::zeroSeedRandom

10th October 2017 Thomas Latham
* Make symmetrisation of vetoes automatic in LauVetoes for symmetric DPs
  - but add new allowed indices (4 and 5) to allow vetoes to be applied to mMin or mMax if desired
* Also apply to fully symmetric DPs
* In addition implement the more efficient calculation of the amplitude in the fully symmetric case

23rd September 2017 Thomas Latham
* Fix expressions for Covariant option for Blatt-Weisskopf momentum for spin > 1
* Make use of setSpinFormalism to set the formalism to Legendre for NR types (instead of kludging with ignoreMomentum)

15th September 2017 Thomas Latham
* Various improvements to the examples

6th September 2017 Thomas Latham
* Improve doxygen comments for event-embedding functions in fit models

5th September 2017 Thomas Latham
* Improve efficiency of Covariant spin factor calculations

31st August 2017 Thomas Latham
* Add further option for parent Blatt-Weisskopf momentum - Covariant (corresponding to Covariant spin factor)

8th August 2017 Thomas Latham
* Implement expressions for the spin factor based on covariant tensor formalism
* Make it much simpler to switch between formalism options via new functions in LauResonanceMaker:
 - setSpinFormalism, setBWType, setBWBachelorRestFrame
 - default settings are unchanged (and are equivalent to setting LauAbsResonance::Zemach_P, LauBlattWeisskopfFactor::BWPrimeBarrier, LauBlattWeisskopfFactor::ResonanceFrame, respectively)

21st July 2017 Thomas Latham
* Add note to README file about compilation issue on Ubuntu 16.04 LTS

21st July 2017 Thomas Latham
* Create public functions to update the kinematics based on one invariant mass and the corresponding helicity angle

20th June 2017 Daniel O'Hanlon
* Terminate when asked to read a non-existent K-matrix parameter file

8th June 2017 Thomas Latham
* Fix compilation error on gcc 4.9

30th May 2017 Thomas Latham
* Permit different efficiency histograms to be defined in classic/square DP and force enable calculation of square DP coords (here and for background PDFs) if required

29th May 2017 Thomas Latham
* Propagate information on EDM from the minimiser to the fit models and into the fit results ntuple

29th May 2017 Thomas Latham
* Ensure that the kinematics will calculate the square DP co-ordinates if the integration requires them

29th May 2017 Thomas Latham
* Reintegrate the RooFit-slave branch into the trunk (and delete the branch)

28th March 2017 Thomas Latham
(in branch for developing RooFit-based slave)
* Rename cacheFitData to verifyFitData in all fit models and RF-slave

28th March 2017 Daniel O'Hanlon
* Fix bug in LauCPFitModel::weightEvents

24th March 2017 Thomas Latham
(in branch for developing RooFit-based slave)
* Refactor code between fit models, master, slave and fit-object classes

22nd March 2017 Thomas Latham
(in branch for developing RooFit-based slave)
* Make the compilation of the RF-based slave class, and its associated example binary, optional

22nd March 2017 Thomas Latham
(in branch for developing RooFit-based slave)
* Complete working example of RooFit-based slave
 - Have identified some scope for refactoring of code
 - Also want to make the compilation of this class, and its associated example binary, optional

15th March 2017 Mark Whitehead
(in branch for developing time-dependent model)
* Handle event-by-event mistag probabilities

7th March 2017 Thomas Latham
* Rename the command for weighting events from "reweight" to "weight"

3rd March 2017 Daniel O'Hanlon
* Implement LauCPFitModel::weightEvents based on the same function in LauSimpleFitModel

1st March 2017 Thomas Latham
(in branch for developing RooFit-based slave)
* Fix root-cling related errors

28th February 2017 Thomas Latham
(in branch for developing RooFit-based slave)
* Start work on creation of RooFit-based slave class
 - Will allow RooFit-based fitters to plug-in straightforwardly to the simultaneous fitting (JFit) framework

20th February 2017 Thomas Latham
* Add warning messages to LauBkgndDPModel::setBkgndHisto when supplied backgroud histogram pointers are null

31st January 2017 Thomas Latham
* Audit of the code to automate the creation of the integration binning scheme.
 - Fix one bug, some tidy-ups, and reintroduce a couple of lost features:
  - the use of a single square-DP grid if a narrow resonance is found in m12
  - extension of region down to threshold if the narrow resonance is close to threshold

22nd January 2017 Daniel O'Hanlon
* Fix bug in automated integration scheme where resonances were assumed to have been added in order of ascending mass

14th December 2016 Daniel O'Hanlon
* Add several light resonances to the list in LauResonanceMaker

13th December 2016 Daniel O'Hanlon
* Automate the determination of the integration scheme for an arbitrary number of narrow resonances in m13 and m23

12th December 2016 Daniel O'Hanlon
* Efficiency saving from modifying the behaviour of LauIsobarDynamics::calculateAmplitudes, LauIsobarDynamics::resAmp and LauIsobarDynamics::incohResAmp in the case of symmetric DPs
 - Previously, there were 2N calls to LauKinematics::flipAndUpdateKinematics for N resonances, now reduced to 2

21st November 2016 John Back
* Added two K-matrix examples for creating pedagogical plots (with associated data files):
 - B3piKMatrixPlots for K-matrix elements, propagator, and pole and SVP production terms
 - B3piKMatrixMassProj for pole and SVP mass projections generated over the DP

17th November 2016 John Back
* Modifications of the K matrix implementation:
 - Changed the format of the parameter file to use keywords, with updated examples
 - The scattering coefficients are now symmetrised, unless "ScattSymmetry 0" is used
 - Allow the option to use the Adler zero suppression term for the production poles and SVPs (default = off).
   These are specified with the useProdAdler boolean in the LauKMatrixProdPole/SVP constructors
 - Added a few helper functions for obtaining coefficients and (internal) matrices
 - Added a method to return the Lorentz-invariant transition amplitude for plots/tests

1st November 2016 Wenbin Qian
* Add calculation of so-called covariant factors for spin amplitude


===
## Laura++ v3r1

9th September 2016 Thomas Latham
* Modification of LauFitNtuple to check the size of the covariance matrix wrt known number of parameters and act accordingly:
  - If it is empty, just fill a diagonal correlation matrix and issue a warning.
  - If it results from a failed first stage of a two-stage fit then the correlation matrix is padded with 1s and 0s for the parameters that were fixed in the first stage.
* Remove the feature that allows parameters to float in first stage of a two-stage fit but to be fixed in second.
* Minor fix to LauAbsCoeffSet: the names of parameters would be mangled if they were the same as the basename, e.g. A1_A would become A1_ while A1_B would be correct.

8th September 2016 Thomas Latham
* Modifications to LauResonanceInfo to allow customisation of which extra parameters are shared between charge conjugate or shared-parameter records.
  - Where the parameters are not shared they are independently created instead of being cloned.
* Modification of LauRhoOmegaMix to take advantage of the above to have the magB and phiB parameters independent.
* Addition to LauResonanceMaker of rho(770)_COPY record that is a shared record with rho(770) - this allows the above new feature to be used.
* Minor unrelated improvement to information messages in LauIsobarDynamics.

25th August 2016 John Back
* Modified LauRhoOmegaMix to allow either a RelBW or GS lineshape for the rho, as well as
  allowing the option to set the second-order denominator term to be equal to unity.
  Also fixed the bug where the spinTerm was not included in the rho-omega amplitude.
  Changed the LauAbsResonance enum from RhoOmegaMix to RhoOmegaMix_GS, RhoOmegaMix_RBW,
  RhoOmegaMix_GS_1 and RhoOmegaMix_RBW_1, where _1 sets the denominator term to unity and
  _GS or _RBW sets the appropriate rho lineshape. The omega is always a RelBW.

* Added to LauAbsResonance: ignoreSpin and ignoreBarrierScaling boolean flags, which ignore
  either the spin terms or the barrier scaling factors for the amplitude. The latter does not
  turn off the resonance barrier factor for mass-dependent widths

* Added to LauResonanceMaker: getResInfo(resonanceName), which retrieves the resonance information.
  This is used to obtain the PDG omega values for initialising LauRhoOmegaMix

* Made LauRelBreitWignerRes ignore momentum-dependent terms for the resonance width if ignoreMomenta is set.
  This is used in the LauRhoOmegaMix for the omega lineshape where its width does not depend on momentum

23rd May 2016 John Back
* Added new lineshape model for rho-omega mass mixing, LauRhoOmegaMix.

12th April 2016 Thomas Latham
* Switch to integrating in square DP when narrow resonances are found in m12.
  - The integration grid size can be specified by the user

19th January 2016 John Back
* Correct the f(m^2) factor in the denominator of the LauGounarisSakuraiRes
  lineshape to use Gamma_0 instead of Gamma(m)

14th January 2016 Thomas Latham
* Documentation improvements

7th December 2015 Thomas Latham
* Resolve bug that meant the order of resonances in LauIsobarDynamics was assumed to match with the order in which the complex coefficients are supplied to the fit model
  - The ordering of resonances is defined by LauIsobarDynamics:
    - Firstly all coherent resonances in order of addition
    - Followed by all incoherent resonances in order of addition
  - The complex coefficients are now rearranged to match that order
  - Printout of the model summary at the end of initialisation has been enhanced to indicate the ordering
  - Doxygen updated to reflect these changes

12th November 2015 Daniel Craik
* Added support for Akima splines and linear interpolation to Lau1DCubicSpline
* LauAbsModIndPartWave, LauModIndPartWaveRealImag and LauModIndPartWaveMagPhase updated to allow choice of spline interpolation method
* LauEFKLLMRes updated to use Akima splines

10th November 2015 Thomas Latham & Daniel Craik
* Add the EFKLLM form-factor model for the Kpi S-wave and an example using this lineshape
* Modify LauResonanceMaker::getResonance to use a switch for greater clarity and easier checking on missing cases

4th November 2015 Daniel Craik
* Add checks to LauIsobarDynamics::addResonance and LauIsobarDynamics::addIncohResonance to stop the wrong type of LauResonanceModel being used
  - LauAbsResonance::isIncoherentModel method added to identify incoherent models

8th September 2015 Mark Whitehead
* Add the ability to modify the error of parameters via the CoeffSet
  - setParameterError added to LauAbsCoeffSet
* Tweak the handling of initial error values passed to MINUIT (to determine initial step size) in LauMinuit


7th September 2015 Mark Whitehead
* Add the ability to Gaussian constrain parameters via the CoeffSet
  - addGaussianConstraint added to LauAbsCoeffSet

12th June 2015 Thomas Latham
* Modifications to Belle-style nonresonant models
  - LauBelleNR modified to use pure Legendre polynomials of cos(theta) in the spin term (i.e. to remove the q*p factors)
  - New form added to LauBelleSymNR (LauAbsResonance::BelleSymNRNoInter) that removes the interference term between the two DP halves
  - The new form also works with non-zero spin (warning added if non-zero spin specified for BelleSymNR and TaylorNR)

8th June 2015 Thomas Latham
* Further work on the blinding mechanism:
  - New method added LauParameter::blindParameter that activates the blinding.
  - The rest of the framework updated to use another new method LauParameter::unblindedValue in all likelihood calculations etc.
  - Example GenFitNoDP updated to include lines to optionally blind the yield parameters.

29th May 2015 Daniel Craik
* Added LauBlind class for blinding and unblinding a value with an offset based on a blinding string

26th May 2015 Daniel Craik
* Stopped LauCPFitModel passing fixed signal/background yields or
  asymmetries to Minuit to avoid hitting limit of 101 fixed parameters

22nd April 2015 Daniel Craik
* Updated MIPW classes to use Lau1DCubicSpline

19th April 2015 Daniel Craik
* Added Lau1DCubicSpline class for 1D spline interpolation

26th March 2015 Thomas Latham
* Reworked MIPW code into abstract base class and derived classes
  to allow different representations of the amplitude at each knot

31st December 2015 Daniel Craik
* Added unbinned goodness of fit tests to examples

12th January 2015 Daniel Craik
* Calculate effective masses for virtual resonances above the upper kinematic limit

10th December 2014 Daniel Craik
* Added coefficient sets to extract gamma from a simultaneous fit to CP and nonCP final states,
  such as the B0->D_CP K pi and B0->D0bar K pi Dalitz plots, as proposed in Phys. Rev. D79, 051301 (2009)
  - LauPolarGammaCPCoeffSet uses the CP parameters r, delta and gamma directly
  - LauRealImagGammaCPCoeffSet parameterises CPV as X_CP+/- and Y_CP+/-
  - LauCartesianGammaCPCoeffSet parameterises CPV as X_CP, Y_CP DeltaX_CP DeltaY_CP
  - Fixed CP parameters are not passed to the fitter so the same coefficient sets can be used for both the
    CP and nonCP Dalitz plots
  - LauPolarGammaCPCoeffSet allows for a single gamma parameter to be shared between multiple resonances
  - LauAbsCoeffSet::adjustName made virtual to allow global variables such as gamma to not receive a prefix

===
## Laura++ v3r0p1

19th June 2015 Thomas Latham
* Factor out the JFit slave code from LauAbsFitModel into a new base class LauSimFitSlave

19th June 2015 Thomas Latham
* Fix check in LauIsobarDynamics::calcDPNormalisationScheme to avoid using hardcoded number

===
## Laura++ v3r0

24th October 2014 Thomas Latham
* Fixed bug in floating of Blatt-Weisskopf barrier radii
  - The values at the pole mass were not being updated when the radii changed

21st October 2014 Daniel Craik
* Fixed bug in LauIsobarDynamics where multiple incoherent amplitudes led to nonsensical fit fractions

17th October 2014 John Back
* Added the ability to calculate the transition amplitude matrix T in LauKMatrixPropagator,
  as well as a few other minor speed-up changes and code checks. Example/PlotKMatrixTAmp.cc
  can be used to check the T amplitude variation, phase shift and inelasticity, for a given
  K matrix channel, as a function of the invariant mass squared variable s

15th October 2014 Thomas Latham
* Add methods to LauIsobarDynamics to make the integration binning more tunable by the user:
  - setNarrowResonanceThreshold - modify the value below which a resonance is considered to be narrow (defaults to 0.02 GeV/c2)
  - setIntegralBinningFactor - modify the factor by which the narrow resonance width is divided to obtain the bin size (defaults to 100)
* Print warning messages if the memory usage is likely to be very large

13th October 2014 Thomas Latham
* Modify Makefile to allow compilation with ROOT 6 (in addition to maintaining support for ROOT 5)
* Fix a few compilation errors on MacOSX 10.9

13th October 2014 Daniel Craik
* Update LauModIndPartWave to allow knots at kinematic limits to be modified
  - Add new method setKnotAmp to modify existing knots (and the knot at the upper kinematic limit which is automatically added at initialisation)
* Update doxygen for LauIsobarDynamics::addIncoherentResonance to mention that incoherent resonances must be added last

10th October 2014 Thomas Latham
* Add new method to LauResonanceMaker to set whether the radius of a given Blatt-Weisskopf category should be fixed of floated
* Modify the methods of LauResonanceMaker to set the radius value and whether it should be floated so that they work before and after the resonances have been created

9th October 2014 John Back
* Corrected the eta-eta' and 4pi phase space factors in LauKMatrixPropagator,
  which is used for the K-matrix amplitude:
  - calcEtaEtaPRho() does not include the mass difference term m_eta - m_eta'
    following the recommendation in hep-ph/0204328 and from advice from M Pennington
  - calcFourPiRho() incorporates a better parameterisation of the double integral of Eq 4 in
    hep-ph/0204328 which avoids the exponential increase for small values of s (~< 0.1)
  - More detailed comments are provided in the above two functions to explain what is
    going on and the reason for the choices made

6th October 2014 Thomas Latham
* Implement the mechanism for floating Blatt-Weisskopf barrier factor radius parameters

30th September 2014 Thomas Latham
* Fix issue in the checks on toy MC generation validity
  - in the case of exceeding max iterations it was possible to enter an infinite loop
  - the checks now detect all three possible states:
    - aSqMaxSet_ is too low (generation is biased) => increase aSqMaxSet_ value
    - aSqMaxSet_ is too high => reduce aSqMaxSet_ value to improve efficiency
    - aSqMaxSet_ is high (causing low efficiency) but cannot be lowered without biasing => increase iterationsMax_ limit
* Update resonance parameter in LauResonanceMaker to match PDG 2014
* Modify behaviour when TTree objects are saved into files to avoid having multiple cycle numbers present

29th September 2014 Daniel Craik
* Add support for incoherent resonances in the signal model
  - LauIsobarDynamics updated to include incoherent terms
  - ABC for incoherent resonances, LauAbsIncohRes, added deriving from LauAbsResonance
  - LauGaussIncohRes added to implement a Gaussian incoherent resonance, deriving from LauAbsIncohRes
  - Small changes to various other classes to enable incoherent terms
* Fixed small bug in LauMagPhaseCoeffSet which could hang if phase is exactly pi or -pi
* Added charged version of the BelleNR resonances to LauResonanceMaker
* Updated parameters in LauConstants to match PDG 2014

14th July 2014 Thomas Latham
* Add intial support for fully-symmetric final states such as B0 -> KS KS KS
  - Performs the symmetrisation of the signal model
  - Background (and efficiency) histogram classes need some work if the user wants to provide folded histograms

8th July 2014 Daniel Craik
* Add class for model-independent partial wave
  - Uses splines to produce a smooth amplitude from a set of magnitude and phase values at given invariant masses
  - The individual magnitudes and phases can be floated in the fit

16th June 2014 Thomas Latham
* Allow floating of resonance parameters in simultaneous fits

13th June 2014 Thomas Latham
* Fix bug in LauResonanceInfo cloning method, where the width parameter was given a clone of the mass

10th June 2014 Thomas Latham
* Add new function to allow sharing of resonance parameters between components that are not charged conjugates, e.g. LASS_BW and LASS_NR

9th June 2014 Thomas Latham and Daniel Craik
* Fix bug in the new integration scheme
  - Was not accounting for cases where several resonances share a floating parameter
  - Meant that the integrals and caches for that resonance were not being updated
* Introduce a change in the implementation of the helicity flip for neutral parent decays
  - Prior to this change the helicity was flipped for any neutral resonance in the decay of a neutral particle.
  - Now the flip no longer occurs in flavour-specific decays (such as Bs0 -> D0bar K- pi+ or B0 -> K+ pi- pi0) since it is only required in flavour-conjugate modes (such as B0 -> KS pi+ pi-).
  - This does not affect any physical results but it does contribute a pi phase flip to some odd-spin resonances (for example K*(892)0 in Bs0->D0barKpi).
  - Therefore results are not necessarily comparable between fits run before and after this changeset.
  - This change will first be released in v3r0.

===
## Laura++ v2r2

5th June 2014 Thomas Latham
* Fix issue in asymmetric efficiency histogram errors
  - Fluctuation of bins was incorrectly sampling - assumed area each side of peak was the same

5th June 2014 Thomas Latham
(in branch for release in v3r0)
* Introduce intelligent calculation of amplitudes during recalculation of integrals and recaching of data points
  - Floating resonance parameters is now much more efficient
* Make resonance parameters second-stage, also improves fit timing when floating them

3rd June 2014 Rafael Coutinho
* Implement generation of toy MC from fit results in fitSlave

27th May 2014 Thomas Latham
(in branch for release in v3r0)
* Complete audit of special functions
* Remove unncessary LauAbsDPDynamics base class and move all functionality into LauIsobarDynamics

20th May 2014 Daniel Craik
(in branch for release in v3r0)
* Add support for K*_0(1430) and a_0(980) to LauFlatteRes

16th-19th May 2014 Thomas Latham and Daniel Craik
(in branch for release in v3r0)
* Update all other lineshapes so that their parameters can float
  - The only resonance parameters that now cannot float are the Blatt-Weisskopf barrier factor radii

15th May 2014 Thomas Latham
(in branch for release in v3r0)
* Change the mechanism for getting floating resonance parameters into the fit
  - Moved from LauResonanceMaker to the resonances themselves
  - Lays some groundwork for improving the efficiency of recalculating the integrals
  - LauBelleNR and LauBelleSymNR lineshape parameters can now be floated

13th May 2014 Daniel Craik
(in branch for release in v3r0)
* Fix bug where illegal characters were being propagated from resonance names into TBranch names

6th May 2014 Thomas Latham
(in branch for release in v3r0)
* Provide accessors for mass and width parameters

5th May 2014 Louis Henry
* Fix compilation problem by making LauDatabasePDG destructor virtual

4th May 2014 Thomas Latham
* Provide a new argument to LauSimpleFitModel::splitSignalComponent and
  LauCPFitModel::splitSignalComponent to allow fluctuation of the bins on the
  SCF fraction histogram

29th April 2014 Thomas Latham
* Fix bug in the determination of the integration scheme
  - Nearby narrow resonances caused problems if their "zones" overlap
  - These zones are now merged together

29th April 2014 Thomas Latham
(in branch for release in v3r0)
* Some improvments to integration scheme storage

26th April 2014 Juan Otalora
(in branch for release in v3r0)
* Make integation scheme fixed after first determination
  - is stored in a new class LauDPPartialIntegralInfo
  - used on subsequent re-evaluations of the integrals

23rd April 2014 Thomas Latham
* Attempt to improve clarity of LauIsobarDynamics::addResonance function
  - the 3rd argument is now an enumeration of the various resonance models
  - removed the optional arguments regarding the change of mass, width & spin
  - the same functionality can be obtained by using the returned pointer and calling changeResonance
  - the resonance name now only affects the lineshape in the LauPolNR case
  - the BelleSymNR / TaylorNR choice is now made in the 3rd argument
  - similarly for the BelleNR / (newly introduced) PowerLawNR choice
* Add new PowerLawNR nonresonant model (part of LauBelleNR)
* All examples updated to use new interface

23rd April 2014 Thomas Latham
* Address issue of setting values of resonance parameters for all models
  - decided to do away with need to have LauIsobarDynamics know everything
  - LauIsobarDynamics::addResonance now returns a pointer to LauAbsResonance
  - parameters can be changed through LauAbsResonance::setResonanceParameter
  - LauIsobarDynamics still knows about Blatt-Weisskopf factors - better to only need to set this once
  - Update GenFit3pi example to demonstrate mechanism

22nd April 2014 Thomas Latham
* Allow Gaussian constraints to be added in simultaneous fitting
  - constraints will be ignored by the slaves
  - those added directly to fit parameters will be propogated to the master and handled there
  - those on combinations of parameters should be added in the master process

22nd April 2014 Mark Whitehead
* Update Laura to cope with resonances of spin 4 and spin 5
  - Zemach spin terms added to src/LauAbsResonance.cc
  - BW barrier factors added to src/LauRelBreitWignerRes.cc

19th April 2014 Daniel Craik
* Add LauWeightedSumEffModel which gives an efficiency model from the weighted sum of several LauEffModel objects.
* Added pABC, LauAbsEffModel, for LauEffModel and LauWeightedSumEffModel.
* Various classes updated to use pointers to LauAbsEffModel instead of LauEffModel.

15th April 2014 Daniel Craik
* Enable LauEfficiencyModel to contain several Lau2DAbsDP objects with the total efficiency calculated as the product.

10th April 2014 Mark Whitehead
* Fix an issue with the likelihood penalty term for Gaussian constraints
  - Factor two missing in the denominator
  - New penalty term is: ( x-mean )^2 / 2*(width^2)

4th April 2014 Thomas Latham
* Add storage of fit fractions that have not been efficiency corrected

===
## Laura++ v2r1

1st April 2014 Thomas Latham
* Fix issue in LauFitter that prevents compilation with g++ 4.8
  - Missing virtual destructor
  - Take opportunity to audit other special functions

31st March 2014 Mark Whitehead
(in branch for release in v2r1)
* Added an efficiency branch to the ntuple produced for toy data samples
  - Both LauSimpleFitModel and LauCPFitModel updated

28th March 2014 Daniel Craik
(in branch for release in v2r2)
* Added support for asymmetric errors to Lau2DHistDP, Lau2DSplineDP and LauEffModel.

27th March 2014 Daniel Craik
* Changed histogram classes to use seeded random number generator for
  fluctuation and raising or lowering of bins and updated doxygen.

20th March 2014 Mark Whitehead
(in branch for release in v2r1)
* Added the ability to add Gaussian contraints to LauFormulaPars of fit parameters
  - User supplies the information but the LauFormulaPar is constructed behind the scenes

18th March 2014 Thomas Latham
* Improve behaviour of toy generation from fit results

13th March 2014 Juan Otalora
(in branch for release in v3r0)
* Extended ability to float mass and width to other resonance lineshapes (Flatte, LASS and G-S)

11th March 2014 Mark Whitehead
(in branch for release in v2r1)
* Added the functionality to make LauFormulaPars usable in fits
  - Added a new class LauAbsRValue which LauParameter and LauFormularPar inherit from
  - Many files updated to accept LauAbsRValues instead of LauParameters
* Fairly major clean up of obsolete functions in LauAbsPdf
  - Only LauLinearPdf used any of them, this has now been updated

10th March 2014 Thomas Latham
(in branch for release in v3r0)
* First attempt at floating resonance parameters (work mostly from Juan)
  - Only works for RelBW lineshape
  - Can only float mass and width
  - Works nicely!
  - Still needs much work to generalise and make more efficient

===
## Laura++ v2r0

8th March 2014 Thomas Latham
* Some additional functionality for the CoeffSet classes:
  - allow the parameter values to be set (optionally setting the initial and generated values as well)
  - allow the parameters to be set to float or to be fixed in the fit
  These are needed when cloning but wanting some of the parameters to have different values and/or floating behaviour from the cloned set.
* Improve the printout of the setting of the coefficient values in the fit models and the creation of resonances
* Add LauFlatNR for the unform NR model - ends its rather strange special treatment
* Fix bug setting resAmpInt to 0 for LauPolNR
* Many other improvements to the info/warning/error printouts
* Modify GenFitBelleCPKpipi example to demonstrate cloning in action
* Add -Werror to compiler flags (treats warnings as errors)

5th March 2014 Thomas Latham
* Some improvements to LauPolNR to speed up amplitude calculation

2nd March 2014 Thomas Latham
* A number of updates to the CoeffSet classes:
  - allow specification of the basename just after construction (before being given to the fit model)
  - allow configuration of the parameter fit ranges (through static methods of base class)
  - more adaptable cloning of the parameters (e.g. can only clone phase but allow magnitude to float independently)
  - allow CP-violating parameters to be second-stage in Belle and CLEO formats
* Some improvements to the Doxygen and runtime information printouts

20th February 2014 Louis Henry
* Add LauPolNR - class for modelling the nonresonant contribution based on BaBar 3K model (arXiv:1201.5897)

6th February 2014 Thomas Latham
* Correct helicity convention information in Doxygen for LauKinematics

===
## Laura++ v1r2

5th February 2014 Thomas Latham
* Add rule to the Makefile that creates a rootmap file for the library

4th February 2014 Daniel Craik
* Fixed bug in Lau2DSplineDPPdf - normalisation was not being calculated
* Added out-of-range warning in LauBkgndDPModel and supressed excessive warnings

3rd February 2014 Mark Whitehead
(in branch for release in v2r0)
* Added a new class to allow parameters to be a function of other parameters
  - inc/LauFormulaPar.hh
  - src/LauFormulaPar.cc

28th January 2014 Daniel Craik
* Improved out-of-range efficiency warnings in LauEffModel and supressed excessive errors
* Modified LauIsobarDynamics to allow LASS parameters to be configured for LauLASSBWRes and
  LauLASSNRRes

27th January 2014 Daniel Craik
* Added spline interpolation to DP backgrounds
  - Added Lau2DSplineDPPdf which uses a spline to model a normalised PDF across a DP
  - Added pABC, Lau2DAbsDPPdf, for Lau2DHistDPPdf and Lau2DSplineDPPdf and moved common
    code in Lau2DHistDPPdf and Lau2DSplineDPPdf into ABC Lau2DAbsHistDPPdf
  - LauBkgndDPModel, modified to use Lau2DAbsDPPdf instead of Lau2DHistDPPdf
  - setBkgndSpline method added to LauBkgndDPModel to allow use of splines

22nd January 2014 Thomas Latham
* Improve some error checks and corresponding warning messages in
  LauCPFitModel::setSignalDPParameters

16th January 2014 Thomas Latham
* Add LauRealImagCPCoeffSet, which provides an (x,y), (xbar,ybar) way of
  parametrising the complex coefficients.
* Try to improve timing in the *CoeffSet classes by making the complex coeffs
  into member variables.

20th December 2013 Daniel Craik
* Added Lau2DCubicSpline which provides cubic spline interpolation of a histogram
  - Added Lau2DSplineDP which uses a spline to model variation across a DP (eg efficiency)
  - Added pABC, Lau2DAbsDP, for Lau2DHistDP and Lau2DSplineDP and moved common code
    in Lau2DHistDP and Lau2DSplineDP into ABC Lau2DAbsHistDP
  - LauEffModel, LauDPDepSumPdf and LauDPDepMapPdf modified to use Lau2DAbsDP instead of
    Lau2DHistDP
  - setEffSpline method added to LauEffModel and useSpline flag added to constructor for
    LauDPDepSumPdf to allow use of splines

18th December 2013 Mark Whitehead
(in branch for release in v2r0)
* Added functionality to include Gaussian constraints on floated
  parameters in the fit.
  The files updated are:
  - inc/LauAbsFitModel.hh
  - inc/LauParameter.hh
  - src/LauAbsFitModel.cc
  - src/LauParameter.cc

5th December 2013 Thomas Latham
* Fix small bug in GenFitKpipi example where background asymmetry parameter had
  its limits the wrong way around

4th December 2013 Daniel Craik
* Updated 2D chi-squared code to use adaptive binning.

3rd December 2013 Thomas Latham
* Generalise the Makefile in the examples directory
  - results in minor changes to the names of 3 of the binaries

3rd December 2013 Thomas Latham
(in branch for release in v2r0)
* Have the master save an ntuple with all fitter parameters and the full correlation matrix information.

29th November 2013 Thomas Latham
* Fixed bug in ResultsExtractor where the output file was not written

29th November 2013 Thomas Latham
(in branch for release in v2r0)
* Allow the slave ntuples to store the partial covariance matrices in the simultaneous fitting

26th November 2013 Thomas Latham
(in branch for release in v2r0)
* Added first version of the simultaneous fitting framework

===
## Laura++ v1r1p1

22nd November 2013 Thomas Latham
* Fixed bug in LauCPFitModel where values of q = -1 extra PDFs
  were used regardless of the event charge.

===
## Laura++ v1r1

20th November 2013 Mark Whitehead
* Changed convention for the barrier factors, swapping from p* to p.
  This seems to give more physically reasonable results.
  The files updated are:
  - src/LauGounarisSakuraiRes.cc
  - src/LauRelBreitWignerRes.cc

18th October 2013 Thomas Latham
* Fix dependency problem in Makefile

8th October 2013 Thomas Latham
* Some fixes to yield implementation
* Minor bug fix in DP background histogram class

7th October 2013 Mark Whitehead
* Update to accept the yields and yield asymmetries as LauParameters.
  All examples have been updated to reflect this change.
  This updated the following files:
  - inc/LauCPFitModel.hh
  - inc/LauSimpleFitModel.hh
  - inc/LauAbsFitModel.hh
  - src/LauCPFitModel.cc
  - src/LauSimpleFitModel.cc

* Addition of the following particles to src/LauResonanceMaker.cc
  Ds*+-, Ds0*(2317)+-, Ds2*(2573)+-, Ds1*(2700)+- and Bs*0

===
## Laura++ v1r0

13th September 2013 Thomas Latham
* Initial import of the package into HEPforge

