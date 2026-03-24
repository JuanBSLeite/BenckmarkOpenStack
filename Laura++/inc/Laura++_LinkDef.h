
/*
Copyright 2013 University of Warwick

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*
Laura++ package authors:
John Back
Paul Harrison
Thomas Latham
*/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ class Lau1DCubicSpline+;
#pragma link C++ class Lau1DHistPdf+;
#pragma link C++ class Lau2DAbsDP+;
#pragma link C++ class Lau2DAbsDPPdf+;
#pragma link C++ class Lau2DAbsHistDP+;
#pragma link C++ class Lau2DAbsHistDPPdf+;
#pragma link C++ class Lau2DCubicSpline+;
#pragma link C++ class Lau2DHistDP+;
#pragma link C++ class Lau2DHistDPPdf+;
#pragma link C++ class Lau2DHistPdf+;
#pragma link C++ class Lau2DSplineDP+;
#pragma link C++ class Lau2DSplineDPPdf+;
#pragma link C++ class LauAbsBkgndDPModel+;
#pragma link C++ class LauAbsCoeffSet+;
#pragma link C++ class LauAbsEffModel+;
#pragma link C++ class LauAbsFitter+;
#pragma link C++ class LauAbsFitModel+;
#pragma link C++ class LauAbsIncohRes+;
#pragma link C++ class LauAbsModIndPartWave+;
#pragma link C++ class LauAbsPdf+;
#pragma link C++ class LauAbsResonance+;
#pragma link C++ class LauAbsRValue+;
#pragma link C++ class LauArgusPdf+;
#pragma link C++ class LauASqMaxFinder+;
#pragma link C++ class LauAsymmCalc+;
#pragma link C++ class LauBelleCPCoeffSet+;
#pragma link C++ class LauBelleNR+;
#pragma link C++ class LauBelleSymNR+;
#pragma link C++ class LauBifurcatedGaussPdf+;
#pragma link C++ class LauBkgndDPModel+;
#pragma link C++ class LauBlattWeisskopfFactor+;
#pragma link C++ class LauBlind+;
#pragma link C++ class LauBreitWignerRes+;
#pragma link C++ class LauCacheData+;
#pragma link C++ class LauCalcChiSq+;
#pragma link C++ class LauCartesianCPCoeffSet+;
#pragma link C++ class LauCartesianGammaCPCoeffSet+;
#pragma link C++ class LauChebychevPdf+;
#pragma link C++ class LauCleoCPCoeffSet+;
#pragma link C++ class LauComplex+;
#pragma link C++ class LauCPFitModel+;
#pragma link C++ class LauCruijffPdf+;
#pragma link C++ class LauCrystalBallPdf+;
#pragma link C++ class LauDabbaRes+;
#pragma link C++ class LauDatabasePDG+;
#pragma link C++ class LauDaughters+;
#pragma link C++ class LauDPDepBifurGaussPdf+;
#pragma link C++ class LauDPDepCruijffPdf+;
#pragma link C++ class LauDPDepGaussPdf+;
#pragma link C++ class LauDPDepMapPdf+;
#pragma link C++ class LauDPDepSumPdf+;
#pragma link C++ class LauEffModel+;
#pragma link C++ class LauEFKLLMRes+;
#pragma link C++ class LauEmbeddedData+;
#pragma link C++ class LauExponentialPdf+;
#pragma link C++ class LauFitDataTree+;
#pragma link C++ class LauFitNtuple+;
#pragma link C++ class LauFitter+;
#pragma link C++ class LauFitObject+;
#pragma link C++ class LauFlatteRes+;
#pragma link C++ class LauFlatNR+;
#pragma link C++ class LauFormulaPar+;
#pragma link C++ class LauGaussIncohRes+;
#pragma link C++ class LauGaussPdf+;
#pragma link C++ class LauGenNtuple+;
#pragma link C++ class LauGounarisSakuraiRes+;
#pragma link C++ class LauIntegrals+;
#pragma link C++ class LauDPPartialIntegralInfo+;
#pragma link C++ class LauIsobarDynamics+;
#pragma link C++ class LauKappaRes+;
#pragma link C++ class LauKinematics+;
#pragma link C++ class LauKMatrixProdPole+;
#pragma link C++ class LauKMatrixProdSVP+;
#pragma link C++ class LauKMatrixPropagator+;
#pragma link C++ class LauKMatrixPropFactory+;
#pragma link C++ class LauLASSBWRes+;
#pragma link C++ class LauLASSNRRes+;
#pragma link C++ class LauLASSRes+;
#pragma link C++ class LauLinearPdf+;
#pragma link C++ class LauMagPhaseCoeffSet+;
#pragma link C++ class LauMagPhaseCPCoeffSet+;
#pragma link C++ class LauMergeDataFiles+;
#pragma link C++ class LauMinuit+;
#pragma link C++ class LauModIndPartWaveMagPhase+;
#pragma link C++ class LauModIndPartWaveRealImag+;
#pragma link C++ class LauNovosibirskPdf+;
#pragma link C++ class LauNRAmplitude+;
#pragma link C++ class LauParameter+;
#pragma link C++ class LauParametricStepFuncPdf+;
#pragma link C++ class LauParamFixed+;
#pragma link C++ class LauParticlePDG+;
#pragma link C++ class LauPolNR+;
#pragma link C++ class LauPoleRes+;
#pragma link C++ class LauPolarFormFactorNR+;
#pragma link C++ class LauPolarFormFactorSymNR+;
#pragma link C++ class LauPolarGammaCPCoeffSet+;
#pragma link C++ class LauPrint+;
#pragma link C++ class LauRealImagCoeffSet+;
#pragma link C++ class LauRealImagCPCoeffSet+;
#pragma link C++ class LauRealImagGammaCPCoeffSet+;
#pragma link C++ class LauRelBreitWignerRes+;
#pragma link C++ class LauResonanceInfo+;
#pragma link C++ class LauRescatteringRes+;
#pragma link C++ class LauRescattering2Res+;
#pragma link C++ class LauResonanceMaker+;
#pragma link C++ class LauResultsExtractor+;
#pragma link C++ class LauRhoOmegaMix+;
#ifdef DOLAUROOFITTASK
#pragma link C++ class LauRooFitTask+;
#endif
#pragma link C++ class LauScfMap+;
#pragma link C++ class LauSigmaRes+;
#pragma link C++ class LauSigmoidPdf+;
#pragma link C++ class LauSimpleFitModel+;
#pragma link C++ class LauSimFitCoordinator+;
#pragma link C++ class LauSimFitTask+;
#pragma link C++ class LauSPlot+;
#pragma link C++ class LauString+;
#pragma link C++ class LauSumPdf+;
#pragma link C++ class LauTextFileParser+;
#pragma link C++ class LauVetoes+;
#pragma link C++ class LauWeightedSumEffModel+;
#pragma link C++ namespace LauConstants+;
#pragma link C++ namespace LauRandom+;

#endif
