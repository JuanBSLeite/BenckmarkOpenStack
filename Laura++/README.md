# Laura++ README

## Introduction
This package is a C++ development of LAURA, the FORTRAN "Likelihood Analysis
Unofficial RooFitDalitz Alternative" code written by Paul Harrison for
performing Dalitz plot analyses of 3-body decays of B mesons at BaBar.

This C++ version was originally developed by John Back and is now mainly
developed by Tom Latham with continuting contributions from John Back and
others.

## Structure of the Package
The package consists of the following directories and files:

```
laura.git
├── CMakeLists.txt - the top-level CMake configuration file
├── LICENSE-2.0 - the Apache License, Version 2.0
├── README.md - this README file
├── applyFormatting.sh - script to apply formatting style to source files
├── cmake - directory containing custom CMake modules and templates
│   ├── Modules
│   │   ├── LauraCompilerFlags.cmake
│   │   └── LauraExternalDependencies.cmake
│   └── Templates
│       └── LauraConfig.cmake.in
├── doc - directory related to building Doxygen documentation
│   ├── CMakeLists.txt - CMake configuration file
│   ├── Doxyfile.in - Doxygen configuration file
│   ├── mainpage.dox.in - front page of the Doxygen documents
│   └── ReleaseNotes.md - history of commits to the package
├── examples - directory containing code examples and their associated data files and job scripts
│   ├── B3piKMNoAdler.dat
│   ├── B3piKMPoles.dat
│   ├── B3piKMSVPs.dat
│   ├── B3piKMatrixCoeff.dat
│   ├── B3piKMatrixMassProj.cc
│   ├── B3piKMatrixPlots.cc
│   ├── CMakeLists.txt
│   ├── CalcChiSq.cc
│   ├── DToKspipiKMatrixCoeff.dat
│   ├── FOCUSD3pi.dat
│   ├── FlatPhaseSpace-CustomMasses.cc
│   ├── FlatPhaseSpace.cc
│   ├── FlatSqDalitz-CustomMasses.cc
│   ├── FlatSqDalitz.cc
│   ├── GenFit3K.cc
│   ├── GenFit3KS.cc
│   ├── GenFit3pi.cc
│   ├── GenFit3pi.py.in
│   ├── GenFitBelleCPKpipi.cc
│   ├── GenFitDpipi.cc
│   ├── GenFitDs2KKpi.cc
│   ├── GenFitEFKLLM.cc
│   ├── GenFitKpipi.cc
│   ├── GenFitNoDP.cc
│   ├── GenFitNoDPMultiDim.cc
│   ├── KMatrixDto3pi.cc
│   ├── KMatrixExample.cc
│   ├── LauKMatrixCoeff.dat
│   ├── LauKMatrixCoeff2.dat
│   ├── MergeDataFiles.cc
│   ├── PlotKMatrixTAmp.cc
│   ├── PlotLikelihood.cc
│   ├── PlotResults.cc
│   ├── QuasiFlatSqDalitz-CustomMasses.cc
│   ├── QuasiFlatSqDalitz.cc
│   ├── ResultsExtractor.cc
│   ├── SimFitCoordinator.cc
│   ├── SimFitTask.cc
│   ├── SimFitTaskRooFit.cc
│   ├── chiSqInput.txt
│   ├── mixedSampleTest.cc
│   ├── point2PointTestSample.cc
│   ├── runCoordinatorRooFitTask.sh
│   ├── runCoordinatorTask.sh
│   ├── runPoint2PointTest.sh
│   └── usfactor.dat
├── fix-comments.py
├── fix-comments.sh
├── inc - directory containing the header files for the Laura++ library
│   ├── CMakeLists.txt
│   ├── Lau1DCubicSpline.hh
│   ├── Lau1DHistPdf.hh
│   ├── Lau2DAbsDP.hh
│   ├── Lau2DAbsDPPdf.hh
│   ├── Lau2DAbsHistDP.hh
│   ├── Lau2DAbsHistDPPdf.hh
│   ├── Lau2DCubicSpline.hh
│   ├── Lau2DHistDP.hh
│   ├── Lau2DHistDPPdf.hh
│   ├── Lau2DHistPdf.hh
│   ├── Lau2DSplineDP.hh
│   ├── Lau2DSplineDPPdf.hh
│   ├── LauASqMaxFinder.hh
│   ├── LauAbsBkgndDPModel.hh
│   ├── LauAbsCoeffSet.hh
│   ├── LauAbsEffModel.hh
│   ├── LauAbsFitModel.hh
│   ├── LauAbsFitter.hh
│   ├── LauAbsIncohRes.hh
│   ├── LauAbsModIndPartWave.hh
│   ├── LauAbsPdf.hh
│   ├── LauAbsRValue.hh
│   ├── LauAbsResonance.hh
│   ├── LauArgusPdf.hh
│   ├── LauAsymmCalc.hh
│   ├── LauBelleCPCoeffSet.hh
│   ├── LauBelleNR.hh
│   ├── LauBelleSymNR.hh
│   ├── LauBifurcatedGaussPdf.hh
│   ├── LauBkgndDPModel.hh
│   ├── LauBlattWeisskopfFactor.hh
│   ├── LauBlind.hh
│   ├── LauBreitWignerRes.hh
│   ├── LauCPFitModel.hh
│   ├── LauCacheData.hh
│   ├── LauCalcChiSq.hh
│   ├── LauCartesianCPCoeffSet.hh
│   ├── LauCartesianGammaCPCoeffSet.hh
│   ├── LauChebychevPdf.hh
│   ├── LauCleoCPCoeffSet.hh
│   ├── LauComplex.hh
│   ├── LauConstants.hh
│   ├── LauCruijffPdf.hh
│   ├── LauCrystalBallPdf.hh
│   ├── LauDPDepBifurGaussPdf.hh
│   ├── LauDPDepCruijffPdf.hh
│   ├── LauDPDepGaussPdf.hh
│   ├── LauDPDepMapPdf.hh
│   ├── LauDPDepSumPdf.hh
│   ├── LauDPPartialIntegralInfo.hh
│   ├── LauDabbaRes.hh
│   ├── LauDatabasePDG.hh
│   ├── LauDaughters.hh
│   ├── LauEFKLLMRes.hh
│   ├── LauEffModel.hh
│   ├── LauEmbeddedData.hh
│   ├── LauExponentialPdf.hh
│   ├── LauFitDataTree.hh
│   ├── LauFitNtuple.hh
│   ├── LauFitObject.hh
│   ├── LauFitter.hh
│   ├── LauFlatNR.hh
│   ├── LauFlatteRes.hh
│   ├── LauFormulaPar.hh
│   ├── LauGaussIncohRes.hh
│   ├── LauGaussPdf.hh
│   ├── LauGenNtuple.hh
│   ├── LauGounarisSakuraiRes.hh
│   ├── LauIntegrals.hh
│   ├── LauIsobarDynamics.hh
│   ├── LauKMatrixProdPole.hh
│   ├── LauKMatrixProdSVP.hh
│   ├── LauKMatrixPropFactory.hh
│   ├── LauKMatrixPropagator.hh
│   ├── LauKappaRes.hh
│   ├── LauKinematics.hh
│   ├── LauLASSBWRes.hh
│   ├── LauLASSNRRes.hh
│   ├── LauLASSRes.hh
│   ├── LauLinearPdf.hh
│   ├── LauMagPhaseCPCoeffSet.hh
│   ├── LauMagPhaseCoeffSet.hh
│   ├── LauMergeDataFiles.hh
│   ├── LauMinuit.hh
│   ├── LauModIndPartWaveMagPhase.hh
│   ├── LauModIndPartWaveRealImag.hh
│   ├── LauNRAmplitude.hh
│   ├── LauNovosibirskPdf.hh
│   ├── LauParamFixed.hh
│   ├── LauParameter.hh
│   ├── LauParametricStepFuncPdf.hh
│   ├── LauParticlePDG.hh
│   ├── LauPolNR.hh
│   ├── LauPolarFormFactorNR.hh
│   ├── LauPolarFormFactorSymNR.hh
│   ├── LauPolarGammaCPCoeffSet.hh
│   ├── LauPoleRes.hh
│   ├── LauPrint.hh
│   ├── LauRandom.hh
│   ├── LauRealImagCPCoeffSet.hh
│   ├── LauRealImagCoeffSet.hh
│   ├── LauRealImagGammaCPCoeffSet.hh
│   ├── LauRelBreitWignerRes.hh
│   ├── LauRescattering2Res.hh
│   ├── LauRescatteringRes.hh
│   ├── LauResonanceInfo.hh
│   ├── LauResonanceMaker.hh
│   ├── LauResultsExtractor.hh
│   ├── LauRhoOmegaMix.hh
│   ├── LauRooFitTask.hh
│   ├── LauSPlot.hh
│   ├── LauScfMap.hh
│   ├── LauSigmaRes.hh
│   ├── LauSigmoidPdf.hh
│   ├── LauSimFitCoordinator.hh
│   ├── LauSimFitTask.hh
│   ├── LauSimpleFitModel.hh
│   ├── LauString.hh
│   ├── LauSumPdf.hh
│   ├── LauTextFileParser.hh
│   ├── LauVetoes.hh
│   ├── LauWeightedSumEffModel.hh
│   └── Laura++_LinkDef.h
├── macros - directory containing some ROOT macro files for plotting results, etc.
│   ├── makeFitFile.C
│   ├── plotComplexVars.C
│   ├── plotCorrs.C
│   ├── plotCorrsFromToy.C
│   └── plotDataIsobars.C
├── src - directory containing the source files for the Laura++ library
│   ├── CMakeLists.txt
│   ├── Lau1DCubicSpline.cc
│   ├── Lau1DHistPdf.cc
│   ├── Lau2DAbsHistDP.cc
│   ├── Lau2DAbsHistDPPdf.cc
│   ├── Lau2DCubicSpline.cc
│   ├── Lau2DHistDP.cc
│   ├── Lau2DHistDPPdf.cc
│   ├── Lau2DHistPdf.cc
│   ├── Lau2DSplineDP.cc
│   ├── Lau2DSplineDPPdf.cc
│   ├── LauASqMaxFinder.cc
│   ├── LauAbsBkgndDPModel.cc
│   ├── LauAbsCoeffSet.cc
│   ├── LauAbsFitModel.cc
│   ├── LauAbsIncohRes.cc
│   ├── LauAbsModIndPartWave.cc
│   ├── LauAbsPdf.cc
│   ├── LauAbsResonance.cc
│   ├── LauArgusPdf.cc
│   ├── LauAsymmCalc.cc
│   ├── LauBelleCPCoeffSet.cc
│   ├── LauBelleNR.cc
│   ├── LauBelleSymNR.cc
│   ├── LauBifurcatedGaussPdf.cc
│   ├── LauBkgndDPModel.cc
│   ├── LauBlattWeisskopfFactor.cc
│   ├── LauBlind.cc
│   ├── LauBreitWignerRes.cc
│   ├── LauCPFitModel.cc
│   ├── LauCacheData.cc
│   ├── LauCalcChiSq.cc
│   ├── LauCartesianCPCoeffSet.cc
│   ├── LauCartesianGammaCPCoeffSet.cc
│   ├── LauChebychevPdf.cc
│   ├── LauCleoCPCoeffSet.cc
│   ├── LauComplex.cc
│   ├── LauCruijffPdf.cc
│   ├── LauCrystalBallPdf.cc
│   ├── LauDPDepBifurGaussPdf.cc
│   ├── LauDPDepCruijffPdf.cc
│   ├── LauDPDepGaussPdf.cc
│   ├── LauDPDepMapPdf.cc
│   ├── LauDPDepSumPdf.cc
│   ├── LauDPPartialIntegralInfo.cc
│   ├── LauDabbaRes.cc
│   ├── LauDatabasePDG.cc
│   ├── LauDaughters.cc
│   ├── LauEFKLLMRes.cc
│   ├── LauEffModel.cc
│   ├── LauEmbeddedData.cc
│   ├── LauExponentialPdf.cc
│   ├── LauFitDataTree.cc
│   ├── LauFitNtuple.cc
│   ├── LauFitObject.cc
│   ├── LauFitter.cc
│   ├── LauFlatNR.cc
│   ├── LauFlatteRes.cc
│   ├── LauFormulaPar.cc
│   ├── LauGaussIncohRes.cc
│   ├── LauGaussPdf.cc
│   ├── LauGenNtuple.cc
│   ├── LauGounarisSakuraiRes.cc
│   ├── LauIntegrals.cc
│   ├── LauIsobarDynamics.cc
│   ├── LauKMatrixProdPole.cc
│   ├── LauKMatrixProdSVP.cc
│   ├── LauKMatrixPropFactory.cc
│   ├── LauKMatrixPropagator.cc
│   ├── LauKappaRes.cc
│   ├── LauKinematics.cc
│   ├── LauLASSBWRes.cc
│   ├── LauLASSNRRes.cc
│   ├── LauLASSRes.cc
│   ├── LauLinearPdf.cc
│   ├── LauMagPhaseCPCoeffSet.cc
│   ├── LauMagPhaseCoeffSet.cc
│   ├── LauMergeDataFiles.cc
│   ├── LauMinuit.cc
│   ├── LauModIndPartWaveMagPhase.cc
│   ├── LauModIndPartWaveRealImag.cc
│   ├── LauNRAmplitude.cc
│   ├── LauNovosibirskPdf.cc
│   ├── LauParameter.cc
│   ├── LauParametricStepFuncPdf.cc
│   ├── LauParticlePDG.cc
│   ├── LauPolNR.cc
│   ├── LauPolarFormFactorNR.cc
│   ├── LauPolarFormFactorSymNR.cc
│   ├── LauPolarGammaCPCoeffSet.cc
│   ├── LauPoleRes.cc
│   ├── LauPrint.cc
│   ├── LauRandom.cc
│   ├── LauRealImagCPCoeffSet.cc
│   ├── LauRealImagCoeffSet.cc
│   ├── LauRealImagGammaCPCoeffSet.cc
│   ├── LauRelBreitWignerRes.cc
│   ├── LauRescattering2Res.cc
│   ├── LauRescatteringRes.cc
│   ├── LauResonanceInfo.cc
│   ├── LauResonanceMaker.cc
│   ├── LauResultsExtractor.cc
│   ├── LauRhoOmegaMix.cc
│   ├── LauRooFitTask.cc
│   ├── LauSPlot.cc
│   ├── LauScfMap.cc
│   ├── LauSigmaRes.cc
│   ├── LauSigmoidPdf.cc
│   ├── LauSimFitCoordinator.cc
│   ├── LauSimFitTask.cc
│   ├── LauSimpleFitModel.cc
│   ├── LauString.cc
│   ├── LauSumPdf.cc
│   ├── LauTextFileParser.cc
│   ├── LauVetoes.cc
│   └── LauWeightedSumEffModel.cc
└── test - directory containing code for some test executables
    ├── CMakeLists.txt
    ├── TestCovariant.cc
    ├── TestCovariant2.cc
    └── TestNewKinematicsMethods.cc
```

## Building the library

### Prerequisites
Compilation of Laura++ on a UNIX operating system requires:
* the [CMake](https://www.cmake.org) build tool
* a C++17 compatible compiler ([GCC](https://gcc.gnu.org) 8 or better, [Clang](https://clang.llvm.org) 8 or better are recommended)
* [GNU Make](https://www.gnu.org/software/make/) or [Ninja](https://ninja-build.org) (not currently tested)

The package depends only on [ROOT](https://root.cern.ch).
Before building the code, it is necessary that the ROOT package be findable by CMake, which can be achieved by doing one of the following:
- the installation prefix of the ROOT package is in the `CMAKE_PREFIX_PATH` environment variable
- the `ROOT_DIR` environment variable is set to the installation prefix of the ROOT package
- the `ROOTSYS` environment variable is set to the installation prefix of the ROOT package
- the directory containing the `root-config` program is in the `PATH` environment variable

Since ROOT depends on the [nlohmann\_json](https://json.nlohmann.me) package, this must be similarly locatable by CMake.

In order to setup an environment that satifies all of the above requirements we highly recommend use of the [LCG views](https://lcginfo.cern.ch).
For example, on a machine with the [CVMFS client](https://cvmfs.readthedocs.io/en/stable/) installed, one can do:
```
source /cvmfs/sft.cern.ch/lcg/views/setupViews.(c)sh <LCG release> <arch-os-compiler>
```
for example:
```
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_104 x86_64-el9-clang16-opt
```

### Build procedure
To build from a clone of this repository, open a terminal and change directory
into that holding this README file.
Create a build directory in which to run `cmake` and the build, and change into
it by doing the following:

```
$ mkdir ../laura.build
$ cd ../laura.build
```

It is best to install the package after building.
You can either install into an existing location or you can create an install
directory specifically for Laura++ (we assume the latter in the following
instructions):

```
$ mkdir ../laura.install
```

Run `cmake` in the build directory, pointing it to the directory (referred to
as `<source dir>` in the following) holding this README, and consequently the
top level CMake script for the project:

```
$ cmake ../<source dir> -DCMAKE_INSTALL_PREFIX=../laura.install
... system specific output ...
-- Configuring done
-- Generating done
-- Build files have been written to: ... your build dir path ...
$
```

The exact output will depend on your system, compiler and build directory
location, but you should not see any errors.
CMake will generate Makefiles by default on UNIX platforms, so to build, simply
run `make` in the build directory:

```
$ make
... verbose output ...
[100%] Built target Laura++
...
$ make install
Install the project...
... verbose output ...
$
```

Again, the exact output will be system specific, but you should see the
`Laura++` target built without error and that the installation was successful.


## Examples and documentation

Example code is included in the examples directory.

To enable the building and installation of the example executables you need to
supply the following option when running `cmake`:
```
-DLAURA_BUILD_EXAMPLES=ON
```
After building and installing, the example executables will be in the `bin`
directory in the install location.

To build the Doxygen documentation you need to supply the following option when
running `cmake`:
```
-DLAURA_BUILD_DOCS=ON
```
After building and installing, you can load the
`share/doc/Laura++/html/index.html` file (found in the install location) into
your web browser.


## Building your own code

If you wish to build your own code that uses the Laura++ library you should add the Laura++ installation area to the `CMAKE_PREFIX_PATH` environment variable.
You can then add lines like the following to your own `CMakeLists.txt` file:
```cmake
# Find Laura++ - all dependencies will automatically be searched for as well
find_package(Laura REQUIRED)

# Now build the executable
add_executable(MyExe MyExe.cc)
target_link_libraries(MyExe PRIVATE Laura::Laura++)
```


## Authors and contributors

As mentioned in the Introduction above, the authors of the package are:\
Thomas Latham\
John Back\
Paul Harrison

The authors would also like to thank the following people for their invaluable
contributions:\
Sian Morgan\
Tim Gershon\
Pablo del Amo Sanchez\
Jelena Ilic\
Eugenia Puccio\
Mark Whitehead\
Daniel Craik\
Rafael Coutinho\
Charlotte Wallace\
Juan Otalora\
Wenbin Qian\
Daniel O'Hanlon\
Daniel Johnson\
Andy Morris

Many thanks also go to numerous members of the BaBar, Belle, and LHCb
collaborations for their helpful input.

Contact: Tom Latham (laura@projects.hepforge.org)

