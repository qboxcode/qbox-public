#-------------------------------------------------------------------------------
#
# Copyright (c) 2008-2014 The Regents of the University of California
#
# This file is part of Qbox
#
# Qbox is distributed under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 2 of
# the License, or (at your option) any later version.
# See the file COPYING in the root directory of this distribution
# or <http://www.gnu.org/licenses/>.
#
#-------------------------------------------------------------------------------
include $(TARGET).mk
#------------------------------------------------------------------------------
#
EXEC=qb
TESTEXECS=testAndersonMixer testBase64Transcoder testBasis testBasisMapping \
testCGOptimizer testChargeDensity testContext testEnergyFunctional \
testFourierTransform testFunction3d testLDAFunctional testLineMinimizer \
testMatrix testSample testSampleReader testSlaterDet \
testSpecies testUnitCell testVWN testWavefunction \
testXCFunctional testXMLGFPreprocessor test_fftw test_vext \
testjacobi testjade 

OBJECTS=MPIdata.o qb.o AtomSet.o Atom.o Species.o \
#       Wavefunction.o SlaterDet.o \
#       EnergyFunctional.o SampleStepper.o \
#       Basis.o BasisMapping.o FourierTransform.o Matrix.o Context.o \
#       sinft.o spline.o UnitCell.o StructureFactor.o \
#       ChargeDensity.o UserInterface.o RunCmd.o \
#       LoadCmd.o SaveCmd.o \
#       SpeciesCmd.o SpeciesReader.o SpeciesHandler.o \
#       XCOperator.o ExchangeOperator.o Bisection.o \
	XCPotential.o LDAFunctional.o VWNFunctional.o \
#       PBEFunctional.o BLYPFunctional.o SCANFunctional.o \
#       B3LYPFunctional.o  BHandHLYPFunctional.o \
	ExponentialIntegral.o HSEFunctional.o \
#       RSHFunctional.o \
	NonLocalPotential.o SampleReader.o StructuredDocumentHandler.o \
	SampleHandler.o AtomSetHandler.o WavefunctionHandler.o \
#       XMLGFPreprocessor.o Base64Transcoder.o \
#       CPSampleStepper.o BOSampleStepper.o \
#       SDWavefunctionStepper.o \
#       MDWavefunctionStepper.o SDIonicStepper.o MDIonicStepper.o \
#       BMDIonicStepper.o \
#       PSDWavefunctionStepper.o PSDAWavefunctionStepper.o \
#       JDWavefunctionStepper.o CellStepper.o SDCellStepper.o \
	CGCellStepper.o ConfinementPotential.o Preconditioner.o \
#       release.o qbox_xmlns.o isodate.o \
#       AndersonMixer.o SDAIonicStepper.o CGIonicStepper.o \
#       ConstraintSet.o Constraint.o PositionConstraint.o DistanceConstraint.o \
#       AngleConstraint.o TorsionConstraint.o jacobi.o \
#       SampleWriter.o ComputeMLWFCmd.o MLWFTransform.o \
#       jade.o PlotCmd.o \
#       ExtForceSet.o ExtForce.o AtomicExtForce.o PairExtForce.o \
#       GlobalExtForce.o \
#       uuid_str.o sampling.o CGOptimizer.o LineMinimizer.o \
#       ElectricEnthalpy.o PartialChargeCmd.o \
#       ExternalPotential.o ResponseCmd.o \
#       Function3d.o Function3dHandler.o \
#       SpectrumCmd.o \
#       $(PLTOBJECTS)
OBJECTS=MPIdata.o qb.o AtomSet.o Atom.o Species.o \
        Basis.o sinft.o spline.o UnitCell.o StructureFactor.o \
        UserInterface.o \
        release.o qbox_xmlns.o isodate.o \
        uuid_str.o sampling.o \
        $(PLTOBJECTS)
CXXFLAGS += -DTARGET='"$(TARGET)"'
# to include VERSION info in release string, use:
# env VERSION=`git describe` make
CXXFLAGS += -DVERSION='"$(VERSION)"'
 $(EXEC):    $(OBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)

 tests: $(TESTEXECS)

 testMPIdata: testMPIdata.o MPIdata.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testLineMinimizer: testLineMinimizer.o LineMinimizer.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testCGOptimizer: testCGOptimizer.o CGOptimizer.o LineMinimizer.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testAndersonMixer: testAndersonMixer.o AndersonMixer.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSample: testSample.o AtomSet.o Atom.o Species.o \
        Wavefunction.o SlaterDet.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o Context.o \
        sinft.o spline.o UnitCell.o \
        Base64Transcoder.o Constraint.o ConstraintSet.o DistanceConstraint.o \
        AngleConstraint.o TorsionConstraint.o PositionConstraint.o \
        ExtForceSet.o ExtForce.o PairExtForce.o AtomicExtForce.o \
        GlobalExtForce.o sampling.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSampleReader: testSampleReader.o AtomSet.o Atom.o Species.o \
        Wavefunction.o SlaterDet.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o Context.o \
        sinft.o spline.o UnitCell.o \
        Base64Transcoder.o Constraint.o ConstraintSet.o DistanceConstraint.o \
        AngleConstraint.o TorsionConstraint.o PositionConstraint.o \
        ExtForceSet.o ExtForce.o PairExtForce.o AtomicExtForce.o \
        GlobalExtForce.o sampling.o \
        SampleReader.o StructuredDocumentHandler.o \
        SampleHandler.o AtomSetHandler.o WavefunctionHandler.o \
        SpeciesReader.o SpeciesHandler.o \
        XMLGFPreprocessor.o Base64Transcoder.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testChargeDensity: testChargeDensity.o ChargeDensity.o \
        Wavefunction.o SlaterDet.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o UnitCell.o \
        Context.o Base64Transcoder.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testWavefunction: testWavefunction.o Wavefunction.o SlaterDet.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o UnitCell.o \
        Context.o Base64Transcoder.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testEnergyFunctional: testEnergyFunctional.o EnergyFunctional.o Basis.o \
	BasisMapping.o SlaterDet.o Matrix.o UnitCell.o Context.o \
        FourierTransform.o \
        Wavefunction.o Species.o Atom.o AtomSet.o StructureFactor.o \
        ExtForceSet.o ExtForce.o PairExtForce.o AtomicExtForce.o \
        GlobalExtForce.o \
        ConstraintSet.o Constraint.o DistanceConstraint.o \
        AngleConstraint.o TorsionConstraint.o PositionConstraint.o \
        NonLocalPotential.o sampling.o Base64Transcoder.o \
        ChargeDensity.o sinft.o spline.o \
        XCOperator.o ExchangeOperator.o Bisection.o \
	XCPotential.o LDAFunctional.o VWNFunctional.o \
        PBEFunctional.o BLYPFunctional.o B3LYPFunctional.o  \
        BHandHLYPFunctional.o SCANFunctional.o \
	ExponentialIntegral.o HSEFunctional.o RSHFunctional.o \
        ConfinementPotential.o ElectricEnthalpy.o MLWFTransform.o \
        ExternalPotential.o Function3d.o Function3dHandler.o qbox_xmlns.o \
        jade.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSlaterDet: testSlaterDet.o SlaterDet.o FourierTransform.o \
        Basis.o BasisMapping.o UnitCell.o Matrix.o Context.o \
        Base64Transcoder.o MPIdata.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testBasisMapping: testBasisMapping.o BasisMapping.o Basis.o \
        Context.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSpecies: testSpecies.o Species.o sinft.o spline.o \
	SpeciesReader.o StructuredDocumentHandler.o SpeciesHandler.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testXCFunctional: testXCFunctional.o LDAFunctional.o PBEFunctional.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testLDAFunctional: testLDAFunctional.o LDAFunctional.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testVWN: testVWN.o VWNFunctional.o LDAFunctional.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testMatrix: testMatrix.o Matrix.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testjacobi: testjacobi.o jacobi.o Matrix.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testjade: testjade.o jade.o Matrix.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testFourierTransform: testFourierTransform.o FourierTransform.o \
        Basis.o BasisMapping.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testUnitCell: testUnitCell.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testBasis: testBasis.o Basis.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testContext: testContext.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testXMLGFPreprocessor: testXMLGFPreprocessor.o XMLGFPreprocessor.o Context.o \
        Base64Transcoder.o Matrix.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testBase64Transcoder: testBase64Transcoder.o Base64Transcoder.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 test_fftw: test_fftw.o $(PLTOBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testFunction3d: testFunction3d.o Function3d.o Function3dHandler.o \
        Base64Transcoder.o qbox_xmlns.o $(PLTOBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 test_vext: test_vext.o Function3d.o Function3dHandler.o \
        Base64Transcoder.o qbox_xmlns.o $(PLTOBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
#------------------------------------------------------------------------------
# generate dependencies in makefile: use -Y to avoid library header files
# that are likely to be different on other platforms.
 depend :
	makedepend -Y -D$(PLT) $(PLTFLAGS) *.cpp
#------------------------------------------------------------------------------
#  Cleanup object files
 clean :
	rm -f *.o
	rm -f $(EXEC) $(TESTEXECS)
#------------------------------------------------------------------------------
 ctags :
	ctags -o tags *.[Ch]
#------------------------------------------------------------------------------
# DO NOT DELETE

AndersonMixer.o: AndersonMixer.h blas.h
AngleConstraint.o: AngleConstraint.h Constraint.h D3vector.h AtomSet.h Atom.h
AngleConstraint.o: UnitCell.h D3tensor.h blas.h Species.h
Atom.o: Atom.h D3vector.h
AtomSet.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h Species.h
AtomSet.o: NameOf.h MPIdata.h sampling.h
AtomSetHandler.o: AtomSetHandler.h StructureHandler.h D3vector.h AtomSet.h
AtomSetHandler.o: Atom.h UnitCell.h D3tensor.h blas.h Species.h
AtomSetHandler.o: SpeciesHandler.h SpeciesReader.h StrX.h SampleReader.h
AtomicExtForce.o: AtomicExtForce.h ExtForce.h D3vector.h AtomSet.h Atom.h
AtomicExtForce.o: UnitCell.h D3tensor.h blas.h Species.h
B3LYPFunctional.o: B3LYPFunctional.h XCFunctional.h BLYPFunctional.h
B3LYPFunctional.o: VWNFunctional.h
BHandHLYPFunctional.o: BHandHLYPFunctional.h XCFunctional.h BLYPFunctional.h
BLYPFunctional.o: BLYPFunctional.h XCFunctional.h
BMDIonicStepper.o: BMDIonicStepper.h IonicStepper.h Sample.h MPIdata.h
BMDIonicStepper.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
BMDIonicStepper.o: Control.h Species.h
BOSampleStepper.o: BOSampleStepper.h SampleStepper.h Timer.h
BOSampleStepper.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
BOSampleStepper.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
BOSampleStepper.o: UnitCell.h SlaterDet.h Basis.h Sample.h MPIdata.h
BOSampleStepper.o: AtomSet.h Atom.h D3tensor.h blas.h Control.h
BOSampleStepper.o: ChargeDensity.h WavefunctionStepper.h
BOSampleStepper.o: SDWavefunctionStepper.h PSDWavefunctionStepper.h
BOSampleStepper.o: PSDAWavefunctionStepper.h JDWavefunctionStepper.h
BOSampleStepper.o: UserInterface.h Preconditioner.h SDIonicStepper.h
BOSampleStepper.o: IonicStepper.h Species.h SDAIonicStepper.h LineMinimizer.h
BOSampleStepper.o: CGIonicStepper.h CGOptimizer.h MDIonicStepper.h
BOSampleStepper.o: BMDIonicStepper.h SDCellStepper.h CellStepper.h
BOSampleStepper.o: CGCellStepper.h AndersonMixer.h MLWFTransform.h
BOSampleStepper.o: BasisMapping.h
Base64Transcoder.o: Base64Transcoder.h
Basis.o: Basis.h D3vector.h UnitCell.h
BasisMapping.o: Basis.h D3vector.h UnitCell.h Context.h blacs.h
BasisMapping.o: BasisMapping.h blas.h
Bisection.o: Bisection.h Context.h blacs.h SlaterDet.h Basis.h D3vector.h
Bisection.o: UnitCell.h Matrix.h Timer.h jade.h FourierTransform.h
Bisection.o: BasisMapping.h
CGCellStepper.o: CGCellStepper.h CellStepper.h Sample.h MPIdata.h AtomSet.h
CGCellStepper.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h Control.h
CGCellStepper.o: CGOptimizer.h LineMinimizer.h
CGIonicStepper.o: CGIonicStepper.h IonicStepper.h Sample.h MPIdata.h
CGIonicStepper.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
CGIonicStepper.o: Control.h Species.h CGOptimizer.h LineMinimizer.h
CGOptimizer.o: CGOptimizer.h LineMinimizer.h blas.h
CPSampleStepper.o: CPSampleStepper.h SampleStepper.h Timer.h
CPSampleStepper.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
CPSampleStepper.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
CPSampleStepper.o: UnitCell.h SlaterDet.h Basis.h Sample.h MPIdata.h
CPSampleStepper.o: AtomSet.h Atom.h D3tensor.h blas.h Control.h
CPSampleStepper.o: ChargeDensity.h MDWavefunctionStepper.h
CPSampleStepper.o: WavefunctionStepper.h MDIonicStepper.h IonicStepper.h
CPSampleStepper.o: Species.h SDCellStepper.h CellStepper.h CGCellStepper.h
CPSampleStepper.o: CGOptimizer.h LineMinimizer.h
CellStepper.o: CellStepper.h Sample.h MPIdata.h AtomSet.h Atom.h D3vector.h
CellStepper.o: UnitCell.h D3tensor.h blas.h Control.h
ChargeDensity.o: ChargeDensity.h Timer.h Context.h blacs.h D3vector.h Basis.h
ChargeDensity.o: UnitCell.h Wavefunction.h FourierTransform.h BasisMapping.h
ChargeDensity.o: SlaterDet.h Matrix.h blas.h
ComputeMLWFCmd.o: ComputeMLWFCmd.h UserInterface.h Sample.h MPIdata.h
ComputeMLWFCmd.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
ComputeMLWFCmd.o: Control.h MLWFTransform.h BasisMapping.h Context.h blacs.h
ComputeMLWFCmd.o: SlaterDet.h Basis.h Matrix.h Timer.h
ConfinementPotential.o: ConfinementPotential.h Basis.h D3vector.h UnitCell.h
Constraint.o: Constraint.h
ConstraintSet.o: ConstraintSet.h PositionConstraint.h Constraint.h
ConstraintSet.o: DistanceConstraint.h AngleConstraint.h D3vector.h
ConstraintSet.o: TorsionConstraint.h Atom.h AtomSet.h UnitCell.h D3tensor.h
ConstraintSet.o: blas.h Context.h blacs.h
Context.o: Context.h blacs.h
DistanceConstraint.o: DistanceConstraint.h Constraint.h AtomSet.h Atom.h
DistanceConstraint.o: D3vector.h UnitCell.h D3tensor.h blas.h Species.h
ElectricEnthalpy.o: Timer.h Context.h blacs.h Matrix.h Sample.h MPIdata.h
ElectricEnthalpy.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
ElectricEnthalpy.o: Control.h ElectricEnthalpy.h Wavefunction.h SlaterDet.h
ElectricEnthalpy.o: Basis.h MLWFTransform.h BasisMapping.h ChargeDensity.h
ElectricEnthalpy.o: FourierTransform.h
EnergyFunctional.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
EnergyFunctional.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
EnergyFunctional.o: UnitCell.h SlaterDet.h Basis.h Timer.h Sample.h MPIdata.h
EnergyFunctional.o: AtomSet.h Atom.h D3tensor.h blas.h Control.h Species.h
EnergyFunctional.o: ChargeDensity.h FourierTransform.h BasisMapping.h
EnergyFunctional.o: XCOperator.h NonLocalPotential.h ConfinementPotential.h
EnergyFunctional.o: ExternalPotential.h
ExchangeOperator.o: VectorLess.h ExchangeOperator.h Sample.h MPIdata.h
ExchangeOperator.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
ExchangeOperator.o: Control.h SlaterDet.h Context.h blacs.h Basis.h Matrix.h
ExchangeOperator.o: Timer.h FourierTransform.h BasisMapping.h Bisection.h
ExponentialIntegral.o: ExponentialIntegral.h
ExtForce.o: ExtForce.h D3vector.h
ExtForceSet.o: ExtForceSet.h AtomicExtForce.h ExtForce.h D3vector.h
ExtForceSet.o: PairExtForce.h GlobalExtForce.h Atom.h AtomSet.h UnitCell.h
ExtForceSet.o: D3tensor.h blas.h Context.h blacs.h
ExternalPotential.o: Basis.h D3vector.h UnitCell.h ExternalPotential.h
ExternalPotential.o: Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h blas.h
ExternalPotential.o: Control.h ChargeDensity.h Timer.h Context.h blacs.h
ExternalPotential.o: FourierTransform.h BasisMapping.h Function3d.h
ExternalPotential.o: Base64Transcoder.h
FourierTransform.o: FourierTransform.h Timer.h BasisMapping.h Basis.h
FourierTransform.o: D3vector.h UnitCell.h blas.h
Function3d.o: Function3d.h D3vector.h Function3dHandler.h StrX.h
Function3d.o: Base64Transcoder.h qbox_xmlns.h Timer.h
Function3dHandler.o: Function3d.h D3vector.h Function3dHandler.h StrX.h
Function3dHandler.o: Base64Transcoder.h Timer.h
GlobalExtForce.o: GlobalExtForce.h ExtForce.h D3vector.h AtomSet.h Atom.h
GlobalExtForce.o: UnitCell.h D3tensor.h blas.h Species.h
HSEFunctional.o: HSEFunctional.h XCFunctional.h ExponentialIntegral.h
JDWavefunctionStepper.o: JDWavefunctionStepper.h WavefunctionStepper.h
JDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
JDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
JDWavefunctionStepper.o: EnergyFunctional.h StructureFactor.h
JDWavefunctionStepper.o: ElectricEnthalpy.h Sample.h MPIdata.h AtomSet.h
JDWavefunctionStepper.o: Atom.h D3tensor.h blas.h Control.h Preconditioner.h
LDAFunctional.o: LDAFunctional.h XCFunctional.h
LineMinimizer.o: LineMinimizer.h
LoadCmd.o: LoadCmd.h UserInterface.h Sample.h MPIdata.h AtomSet.h Atom.h
LoadCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h Control.h SampleReader.h
MDIonicStepper.o: MDIonicStepper.h IonicStepper.h Sample.h MPIdata.h
MDIonicStepper.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
MDIonicStepper.o: Control.h Species.h sampling.h
MDWavefunctionStepper.o: MDWavefunctionStepper.h WavefunctionStepper.h
MDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
MDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
MDWavefunctionStepper.o: Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h
MDWavefunctionStepper.o: blas.h Control.h
MLWFTransform.o: MLWFTransform.h D3vector.h BasisMapping.h Basis.h UnitCell.h
MLWFTransform.o: SlaterDet.h Context.h blacs.h Matrix.h Timer.h jade.h blas.h
MPIdata.o: MPIdata.h
Matrix.o: Context.h blacs.h Matrix.h
NonLocalPotential.o: NonLocalPotential.h AtomSet.h Atom.h D3vector.h
NonLocalPotential.o: UnitCell.h D3tensor.h blas.h Basis.h SlaterDet.h
NonLocalPotential.o: Context.h blacs.h Matrix.h Timer.h Species.h
PBEFunctional.o: PBEFunctional.h XCFunctional.h
PSDAWavefunctionStepper.o: PSDAWavefunctionStepper.h WavefunctionStepper.h
PSDAWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
PSDAWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
PSDAWavefunctionStepper.o: Preconditioner.h
PSDWavefunctionStepper.o: PSDWavefunctionStepper.h WavefunctionStepper.h
PSDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
PSDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
PSDWavefunctionStepper.o: Preconditioner.h
PairExtForce.o: PairExtForce.h ExtForce.h D3vector.h AtomSet.h Atom.h
PairExtForce.o: UnitCell.h D3tensor.h blas.h Species.h
PartialChargeCmd.o: PartialChargeCmd.h UserInterface.h Sample.h MPIdata.h
PartialChargeCmd.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
PartialChargeCmd.o: Control.h Context.h blacs.h Basis.h ChargeDensity.h
PartialChargeCmd.o: Timer.h
PlotCmd.o: PlotCmd.h UserInterface.h Sample.h MPIdata.h AtomSet.h Atom.h
PlotCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h Control.h isodate.h
PlotCmd.o: release.h Context.h blacs.h SampleReader.h Basis.h
PlotCmd.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h Matrix.h
PlotCmd.o: Wavefunction.h SlaterDet.h Timer.h FourierTransform.h
PlotCmd.o: BasisMapping.h Species.h ChargeDensity.h
PositionConstraint.o: PositionConstraint.h Constraint.h AtomSet.h Atom.h
PositionConstraint.o: D3vector.h UnitCell.h D3tensor.h blas.h Species.h
Preconditioner.o: Preconditioner.h Basis.h D3vector.h UnitCell.h
Preconditioner.o: Wavefunction.h EnergyFunctional.h StructureFactor.h
Preconditioner.o: ElectricEnthalpy.h Matrix.h Context.h blacs.h SlaterDet.h
Preconditioner.o: Timer.h Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h
Preconditioner.o: blas.h Control.h ConfinementPotential.h
RSHFunctional.o: RSHFunctional.h XCFunctional.h ExponentialIntegral.h
ResponseCmd.o: BOSampleStepper.h SampleStepper.h Timer.h EnergyFunctional.h
ResponseCmd.o: StructureFactor.h ElectricEnthalpy.h Matrix.h Context.h
ResponseCmd.o: blacs.h D3vector.h Wavefunction.h UnitCell.h SlaterDet.h
ResponseCmd.o: Basis.h Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h blas.h
ResponseCmd.o: Control.h ChargeDensity.h ExternalPotential.h
ResponseCmd.o: FourierTransform.h BasisMapping.h ResponseCmd.h
ResponseCmd.o: UserInterface.h release.h isodate.h Species.h Function3d.h
RunCmd.o: RunCmd.h UserInterface.h BOSampleStepper.h SampleStepper.h Timer.h
RunCmd.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h Matrix.h
RunCmd.o: Context.h blacs.h D3vector.h Wavefunction.h UnitCell.h SlaterDet.h
RunCmd.o: Basis.h Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h blas.h
RunCmd.o: Control.h ChargeDensity.h CPSampleStepper.h
SCANFunctional.o: SCANFunctional.h XCFunctional.h
SDAIonicStepper.o: SDAIonicStepper.h IonicStepper.h Sample.h MPIdata.h
SDAIonicStepper.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
SDAIonicStepper.o: Control.h Species.h LineMinimizer.h
SDCellStepper.o: SDCellStepper.h CellStepper.h Sample.h MPIdata.h AtomSet.h
SDCellStepper.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h Control.h
SDIonicStepper.o: SDIonicStepper.h IonicStepper.h Sample.h MPIdata.h
SDIonicStepper.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
SDIonicStepper.o: Control.h Species.h
SDWavefunctionStepper.o: SDWavefunctionStepper.h WavefunctionStepper.h
SDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
SDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
SDWavefunctionStepper.o: Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h
SDWavefunctionStepper.o: blas.h Control.h
SampleHandler.o: SampleHandler.h StructureHandler.h Sample.h MPIdata.h
SampleHandler.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
SampleHandler.o: Control.h AtomSetHandler.h WavefunctionHandler.h
SampleHandler.o: Wavefunction.h SlaterDet.h Context.h blacs.h Basis.h
SampleHandler.o: Matrix.h Timer.h StrX.h
SampleReader.o: Sample.h MPIdata.h AtomSet.h Atom.h D3vector.h UnitCell.h
SampleReader.o: D3tensor.h blas.h Control.h SampleReader.h
SampleReader.o: XMLGFPreprocessor.h Matrix.h Context.h blacs.h Timer.h
SampleReader.o: SampleHandler.h StructureHandler.h
SampleReader.o: StructuredDocumentHandler.h StrX.h
SampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h StructureFactor.h
SampleStepper.o: ElectricEnthalpy.h Matrix.h Context.h blacs.h D3vector.h
SampleStepper.o: Wavefunction.h UnitCell.h SlaterDet.h Basis.h Sample.h
SampleStepper.o: MPIdata.h AtomSet.h Atom.h D3tensor.h blas.h Control.h
SampleStepper.o: ChargeDensity.h Species.h
SampleWriter.o: SampleWriter.h Context.h blacs.h Sample.h MPIdata.h AtomSet.h
SampleWriter.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h Control.h
SampleWriter.o: qbox_xmlns.h Timer.h SharedFilePtr.h
SaveCmd.o: SaveCmd.h UserInterface.h Sample.h MPIdata.h AtomSet.h Atom.h
SaveCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h Control.h SampleWriter.h
SaveCmd.o: Context.h blacs.h isodate.h release.h
SlaterDet.o: SlaterDet.h Context.h blacs.h Basis.h D3vector.h UnitCell.h
SlaterDet.o: Matrix.h Timer.h FourierTransform.h BasisMapping.h blas.h
SlaterDet.o: Base64Transcoder.h SharedFilePtr.h MPIdata.h
Species.o: Species.h spline.h sinft.h
SpeciesCmd.o: SpeciesCmd.h UserInterface.h Sample.h MPIdata.h AtomSet.h
SpeciesCmd.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h Control.h
SpeciesCmd.o: SpeciesReader.h Species.h
SpeciesHandler.o: SpeciesHandler.h StructureHandler.h Species.h StrX.h
SpeciesReader.o: Species.h SpeciesReader.h StructuredDocumentHandler.h StrX.h
SpeciesReader.o: StructureHandler.h SpeciesHandler.h
SpectrumCmd.o: SpectrumCmd.h UserInterface.h Context.h blacs.h
SpectrumCmd.o: ChargeDensity.h Timer.h D3vector.h EnergyFunctional.h
SpectrumCmd.o: StructureFactor.h ElectricEnthalpy.h Matrix.h Wavefunction.h
SpectrumCmd.o: UnitCell.h SlaterDet.h Basis.h Sample.h MPIdata.h AtomSet.h
SpectrumCmd.o: Atom.h D3tensor.h blas.h Control.h MLWFTransform.h
SpectrumCmd.o: BasisMapping.h
StructureFactor.o: StructureFactor.h Basis.h D3vector.h UnitCell.h
StructuredDocumentHandler.o: StructuredDocumentHandler.h StrX.h
StructuredDocumentHandler.o: StructureHandler.h
TorsionConstraint.o: TorsionConstraint.h Constraint.h D3vector.h AtomSet.h
TorsionConstraint.o: Atom.h UnitCell.h D3tensor.h blas.h Species.h
UnitCell.o: UnitCell.h D3vector.h
UserInterface.o: UserInterface.h qbox_xmlns.h
VWNFunctional.o: VWNFunctional.h XCFunctional.h
Wavefunction.o: Wavefunction.h D3vector.h UnitCell.h SlaterDet.h Context.h
Wavefunction.o: blacs.h Basis.h Matrix.h Timer.h FourierTransform.h
Wavefunction.o: BasisMapping.h jacobi.h SharedFilePtr.h
WavefunctionHandler.o: WavefunctionHandler.h StructureHandler.h UnitCell.h
WavefunctionHandler.o: D3vector.h Wavefunction.h SlaterDet.h Context.h
WavefunctionHandler.o: blacs.h Basis.h Matrix.h Timer.h FourierTransform.h
WavefunctionHandler.o: BasisMapping.h SampleReader.h StrX.h
XCOperator.o: XCOperator.h Sample.h MPIdata.h AtomSet.h Atom.h D3vector.h
XCOperator.o: UnitCell.h D3tensor.h blas.h Control.h ChargeDensity.h Timer.h
XCOperator.o: Context.h blacs.h XCPotential.h ExchangeOperator.h SlaterDet.h
XCOperator.o: Basis.h Matrix.h FourierTransform.h BasisMapping.h
XCOperator.o: HSEFunctional.h XCFunctional.h RSHFunctional.h
XCPotential.o: XCPotential.h ChargeDensity.h Timer.h Context.h blacs.h
XCPotential.o: D3vector.h LDAFunctional.h XCFunctional.h VWNFunctional.h
XCPotential.o: PBEFunctional.h BLYPFunctional.h HSEFunctional.h
XCPotential.o: RSHFunctional.h B3LYPFunctional.h BHandHLYPFunctional.h
XCPotential.o: SCANFunctional.h Sample.h MPIdata.h AtomSet.h Atom.h
XCPotential.o: UnitCell.h D3tensor.h blas.h Control.h Wavefunction.h
XCPotential.o: SlaterDet.h Basis.h Matrix.h FourierTransform.h BasisMapping.h
XMLGFPreprocessor.o: Timer.h Context.h blacs.h Base64Transcoder.h Matrix.h
XMLGFPreprocessor.o: XMLGFPreprocessor.h
isodate.o: isodate.h
jacobi.o: blacs.h Context.h Matrix.h blas.h
jade.o: blacs.h Context.h Matrix.h blas.h Timer.h
qb.o: isodate.h release.h qbox_xmlns.h uuid_str.h MPIdata.h Timer.h
qb.o: UserInterface.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
qb.o: D3tensor.h blas.h Control.h QuitCmd.h
qbox_xmlns.o: qbox_xmlns.h
release.o: release.h
sinft.o: sinft.h
spline.o: spline.h
testAndersonMixer.o: AndersonMixer.h
testBase64Transcoder.o: Base64Transcoder.h
testBasis.o: Basis.h D3vector.h UnitCell.h
testBasisMapping.o: Context.h blacs.h Basis.h D3vector.h UnitCell.h
testBasisMapping.o: BasisMapping.h Timer.h
testCGOptimizer.o: CGOptimizer.h LineMinimizer.h
testChargeDensity.o: Context.h blacs.h Wavefunction.h D3vector.h UnitCell.h
testChargeDensity.o: ChargeDensity.h Timer.h SlaterDet.h Basis.h Matrix.h
testChargeDensity.o: FourierTransform.h BasisMapping.h
testContext.o: Context.h blacs.h
testEnergyFunctional.o: Context.h blacs.h Sample.h MPIdata.h AtomSet.h Atom.h
testEnergyFunctional.o: D3vector.h UnitCell.h D3tensor.h blas.h Control.h
testEnergyFunctional.o: Wavefunction.h ChargeDensity.h Timer.h
testEnergyFunctional.o: EnergyFunctional.h StructureFactor.h
testEnergyFunctional.o: ElectricEnthalpy.h Matrix.h SlaterDet.h Basis.h
testFourierTransform.o: Basis.h D3vector.h UnitCell.h FourierTransform.h
testFourierTransform.o: Timer.h BasisMapping.h
testFunction3d.o: Function3d.h D3vector.h
testLDAFunctional.o: LDAFunctional.h XCFunctional.h
testLineMinimizer.o: LineMinimizer.h
testMPIdata.o: MPIdata.h
testMatrix.o: Timer.h Context.h blacs.h Matrix.h
testSample.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h UnitCell.h
testSample.o: Matrix.h Timer.h Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h
testSample.o: blas.h Control.h
testSampleReader.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h
testSampleReader.o: UnitCell.h Matrix.h Timer.h Sample.h MPIdata.h AtomSet.h
testSampleReader.o: Atom.h D3tensor.h blas.h Control.h SampleReader.h
testSlaterDet.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h UnitCell.h
testSlaterDet.o: Matrix.h Timer.h FourierTransform.h BasisMapping.h MPIdata.h
testSpecies.o: Species.h SpeciesReader.h
testUnitCell.o: UnitCell.h D3vector.h
testVWN.o: VWNFunctional.h XCFunctional.h LDAFunctional.h
testWavefunction.o: Context.h blacs.h Wavefunction.h D3vector.h UnitCell.h
testWavefunction.o: SlaterDet.h Basis.h Matrix.h Timer.h
testXCFunctional.o: LDAFunctional.h XCFunctional.h PBEFunctional.h Timer.h
testXMLGFPreprocessor.o: Context.h blacs.h Matrix.h XMLGFPreprocessor.h
test_vext.o: Function3d.h D3vector.h
testjacobi.o: Timer.h Context.h blacs.h Matrix.h jacobi.h
testjade.o: Timer.h Context.h blacs.h Matrix.h jade.h
uuid_str.o: uuid_str.h
