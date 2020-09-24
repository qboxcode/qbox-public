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
        Wavefunction.o SlaterDet.o \
        EnergyFunctional.o SampleStepper.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o Context.o \
        sinft.o spline.o UnitCell.o StructureFactor.o \
        ChargeDensity.o UserInterface.o RunCmd.o \
        SpeciesCmd.o SpeciesReader.o SpeciesHandler.o \
        XCOperator.o ExchangeOperator.o Bisection.o \
	XCPotential.o LDAFunctional.o VWNFunctional.o \
        PBEFunctional.o BLYPFunctional.o SCANFunctional.o \
        B3LYPFunctional.o  BHandHLYPFunctional.o \
        NonLocalPotential.o \
        SDWavefunctionStepper.o PSDWavefunctionStepper.o \
        PSDAWavefunctionStepper.o JDWavefunctionStepper.o \
        release.o qbox_xmlns.o isodate.o \
        ConstraintSet.o Constraint.o PositionConstraint.o DistanceConstraint.o \
        AngleConstraint.o TorsionConstraint.o \
        ExtForceSet.o ExtForce.o AtomicExtForce.o PairExtForce.o \
        GlobalExtForce.o \
        uuid_str.o sampling.o \
        LoadCmd.o SaveCmd.o \
        ExponentialIntegral.o HSEFunctional.o \
        RSHFunctional.o \
        SampleReader.o StructuredDocumentHandler.o \
        SampleHandler.o AtomSetHandler.o WavefunctionHandler.o \
        XMLGFPreprocessor.o Base64Transcoder.o \
        CPSampleStepper.o BOSampleStepper.o \
        MDWavefunctionStepper.o SDIonicStepper.o MDIonicStepper.o \
        BMDIonicStepper.o \
        CellStepper.o SDCellStepper.o \
        CGCellStepper.o ConfinementPotential.o Preconditioner.o \
        AndersonMixer.o SDAIonicStepper.o CGIonicStepper.o \
        jacobi.o jade.o \
        SampleWriter.o ComputeMLWFCmd.o MLWFTransform.o \
        PlotCmd.o \
        sampling.o CGOptimizer.o LineMinimizer.o \
        ElectricEnthalpy.o PartialChargeCmd.o \
        ExternalPotential.o ResponseCmd.o \
        Function3d.o Function3dHandler.o \
        SpectrumCmd.o \
        cout0.o \
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
        Context.o Base64Transcoder.o MPIdata.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testWavefunction: testWavefunction.o Wavefunction.o SlaterDet.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o UnitCell.o \
        Context.o Base64Transcoder.o MPIdata.o
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
        jade.o MPIdata.o
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
BMDIonicStepper.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
BMDIonicStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
BMDIonicStepper.o: Wavefunction.h Control.h Species.h
BOSampleStepper.o: BOSampleStepper.h SampleStepper.h Timer.h
BOSampleStepper.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
BOSampleStepper.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
BOSampleStepper.o: UnitCell.h SlaterDet.h Basis.h Sample.h MPIdata.h
BOSampleStepper.o: AtomSet.h Atom.h D3tensor.h blas.h ConstraintSet.h
BOSampleStepper.o: ExtForceSet.h Control.h ChargeDensity.h
BOSampleStepper.o: WavefunctionStepper.h SDWavefunctionStepper.h
BOSampleStepper.o: PSDWavefunctionStepper.h PSDAWavefunctionStepper.h
BOSampleStepper.o: JDWavefunctionStepper.h UserInterface.h Preconditioner.h
BOSampleStepper.o: SDIonicStepper.h IonicStepper.h Species.h
BOSampleStepper.o: SDAIonicStepper.h LineMinimizer.h CGIonicStepper.h
BOSampleStepper.o: CGOptimizer.h MDIonicStepper.h BMDIonicStepper.h
BOSampleStepper.o: SDCellStepper.h CellStepper.h CGCellStepper.h
BOSampleStepper.o: AndersonMixer.h MLWFTransform.h BasisMapping.h cout0.h
Base64Transcoder.o: Base64Transcoder.h
Basis.o: Basis.h D3vector.h UnitCell.h
BasisMapping.o: Basis.h D3vector.h UnitCell.h Context.h blacs.h
BasisMapping.o: BasisMapping.h blas.h
Bisection.o: Bisection.h Context.h blacs.h SlaterDet.h Basis.h D3vector.h
Bisection.o: UnitCell.h Matrix.h Timer.h jade.h FourierTransform.h
Bisection.o: BasisMapping.h
CGCellStepper.o: MPIdata.h CGCellStepper.h CellStepper.h Sample.h Context.h
CGCellStepper.o: blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h
CGCellStepper.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
CGCellStepper.o: Control.h CGOptimizer.h LineMinimizer.h
CGIonicStepper.o: MPIdata.h CGIonicStepper.h IonicStepper.h Sample.h
CGIonicStepper.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
CGIonicStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
CGIonicStepper.o: Wavefunction.h Control.h Species.h CGOptimizer.h
CGIonicStepper.o: LineMinimizer.h
CGOptimizer.o: CGOptimizer.h LineMinimizer.h blas.h
CPSampleStepper.o: MPIdata.h CPSampleStepper.h SampleStepper.h Timer.h
CPSampleStepper.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
CPSampleStepper.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
CPSampleStepper.o: UnitCell.h SlaterDet.h Basis.h Sample.h AtomSet.h Atom.h
CPSampleStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Control.h
CPSampleStepper.o: ChargeDensity.h MDWavefunctionStepper.h
CPSampleStepper.o: WavefunctionStepper.h MDIonicStepper.h IonicStepper.h
CPSampleStepper.o: Species.h SDCellStepper.h CellStepper.h CGCellStepper.h
CPSampleStepper.o: CGOptimizer.h LineMinimizer.h
CellStepper.o: CellStepper.h Sample.h MPIdata.h Context.h blacs.h AtomSet.h
CellStepper.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
CellStepper.o: ExtForceSet.h Wavefunction.h Control.h
ChargeDensity.o: MPIdata.h ChargeDensity.h Timer.h Context.h blacs.h
ChargeDensity.o: D3vector.h Basis.h UnitCell.h Wavefunction.h
ChargeDensity.o: FourierTransform.h BasisMapping.h SlaterDet.h Matrix.h
ChargeDensity.o: blas.h
ComputeMLWFCmd.o: ComputeMLWFCmd.h UserInterface.h Sample.h MPIdata.h
ComputeMLWFCmd.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
ComputeMLWFCmd.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
ComputeMLWFCmd.o: Wavefunction.h Control.h MLWFTransform.h BasisMapping.h
ComputeMLWFCmd.o: SlaterDet.h Basis.h Matrix.h Timer.h cout0.h
ConfinementPotential.o: ConfinementPotential.h Basis.h D3vector.h UnitCell.h
Constraint.o: Constraint.h
ConstraintSet.o: MPIdata.h ConstraintSet.h PositionConstraint.h Constraint.h
ConstraintSet.o: DistanceConstraint.h AngleConstraint.h D3vector.h
ConstraintSet.o: TorsionConstraint.h Atom.h AtomSet.h UnitCell.h D3tensor.h
ConstraintSet.o: blas.h
Context.o: Context.h blacs.h
DistanceConstraint.o: DistanceConstraint.h Constraint.h AtomSet.h Atom.h
DistanceConstraint.o: D3vector.h UnitCell.h D3tensor.h blas.h Species.h
ElectricEnthalpy.o: Timer.h Context.h blacs.h Matrix.h Sample.h MPIdata.h
ElectricEnthalpy.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
ElectricEnthalpy.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ElectricEnthalpy.o: ElectricEnthalpy.h SlaterDet.h Basis.h MLWFTransform.h
ElectricEnthalpy.o: BasisMapping.h ChargeDensity.h FourierTransform.h
EnergyFunctional.o: MPIdata.h EnergyFunctional.h StructureFactor.h
EnergyFunctional.o: ElectricEnthalpy.h Matrix.h Context.h blacs.h D3vector.h
EnergyFunctional.o: Wavefunction.h UnitCell.h SlaterDet.h Basis.h Timer.h
EnergyFunctional.o: Sample.h AtomSet.h Atom.h D3tensor.h blas.h
EnergyFunctional.o: ConstraintSet.h ExtForceSet.h Control.h Species.h
EnergyFunctional.o: ChargeDensity.h FourierTransform.h BasisMapping.h
EnergyFunctional.o: XCOperator.h NonLocalPotential.h ConfinementPotential.h
EnergyFunctional.o: ExternalPotential.h
ExchangeOperator.o: VectorLess.h ExchangeOperator.h Sample.h MPIdata.h
ExchangeOperator.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
ExchangeOperator.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
ExchangeOperator.o: Wavefunction.h Control.h SlaterDet.h Basis.h Matrix.h
ExchangeOperator.o: Timer.h FourierTransform.h BasisMapping.h Bisection.h
ExponentialIntegral.o: ExponentialIntegral.h
ExtForce.o: ExtForce.h D3vector.h
ExtForceSet.o: MPIdata.h ExtForceSet.h AtomicExtForce.h ExtForce.h D3vector.h
ExtForceSet.o: PairExtForce.h GlobalExtForce.h Atom.h AtomSet.h UnitCell.h
ExtForceSet.o: D3tensor.h blas.h
ExternalPotential.o: Basis.h D3vector.h UnitCell.h ExternalPotential.h
ExternalPotential.o: Sample.h MPIdata.h Context.h blacs.h AtomSet.h Atom.h
ExternalPotential.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
ExternalPotential.o: Wavefunction.h Control.h ChargeDensity.h Timer.h
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
JDWavefunctionStepper.o: Atom.h D3tensor.h blas.h ConstraintSet.h
JDWavefunctionStepper.o: ExtForceSet.h Control.h Preconditioner.h
LDAFunctional.o: LDAFunctional.h XCFunctional.h
LineMinimizer.o: LineMinimizer.h
LoadCmd.o: LoadCmd.h UserInterface.h Sample.h MPIdata.h Context.h blacs.h
LoadCmd.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
LoadCmd.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
LoadCmd.o: SampleReader.h
MDIonicStepper.o: MPIdata.h MDIonicStepper.h IonicStepper.h Sample.h
MDIonicStepper.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
MDIonicStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
MDIonicStepper.o: Wavefunction.h Control.h Species.h sampling.h
MDWavefunctionStepper.o: MDWavefunctionStepper.h WavefunctionStepper.h
MDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
MDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
MDWavefunctionStepper.o: Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h
MDWavefunctionStepper.o: blas.h ConstraintSet.h ExtForceSet.h Control.h
MLWFTransform.o: MLWFTransform.h D3vector.h BasisMapping.h Basis.h UnitCell.h
MLWFTransform.o: SlaterDet.h Context.h blacs.h Matrix.h Timer.h jade.h blas.h
MPIdata.o: MPIdata.h
Matrix.o: Context.h blacs.h Matrix.h
NonLocalPotential.o: NonLocalPotential.h AtomSet.h Atom.h D3vector.h
NonLocalPotential.o: UnitCell.h D3tensor.h blas.h Basis.h SlaterDet.h
NonLocalPotential.o: Context.h blacs.h Matrix.h Timer.h Species.h MPIdata.h
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
PartialChargeCmd.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
PartialChargeCmd.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
PartialChargeCmd.o: Wavefunction.h Control.h Basis.h ChargeDensity.h Timer.h
PlotCmd.o: PlotCmd.h UserInterface.h Sample.h MPIdata.h Context.h blacs.h
PlotCmd.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
PlotCmd.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h isodate.h
PlotCmd.o: release.h SampleReader.h Basis.h EnergyFunctional.h
PlotCmd.o: StructureFactor.h ElectricEnthalpy.h Matrix.h SlaterDet.h Timer.h
PlotCmd.o: FourierTransform.h BasisMapping.h Species.h ChargeDensity.h
PositionConstraint.o: PositionConstraint.h Constraint.h AtomSet.h Atom.h
PositionConstraint.o: D3vector.h UnitCell.h D3tensor.h blas.h Species.h
Preconditioner.o: Preconditioner.h Basis.h D3vector.h UnitCell.h
Preconditioner.o: Wavefunction.h EnergyFunctional.h StructureFactor.h
Preconditioner.o: ElectricEnthalpy.h Matrix.h Context.h blacs.h SlaterDet.h
Preconditioner.o: Timer.h Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h
Preconditioner.o: blas.h ConstraintSet.h ExtForceSet.h Control.h
Preconditioner.o: ConfinementPotential.h
RSHFunctional.o: RSHFunctional.h XCFunctional.h ExponentialIntegral.h
ResponseCmd.o: MPIdata.h BOSampleStepper.h SampleStepper.h Timer.h
ResponseCmd.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
ResponseCmd.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
ResponseCmd.o: UnitCell.h SlaterDet.h Basis.h Sample.h AtomSet.h Atom.h
ResponseCmd.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Control.h
ResponseCmd.o: ChargeDensity.h ExternalPotential.h FourierTransform.h
ResponseCmd.o: BasisMapping.h ResponseCmd.h UserInterface.h release.h
ResponseCmd.o: isodate.h Species.h Function3d.h
RunCmd.o: RunCmd.h UserInterface.h BOSampleStepper.h SampleStepper.h Timer.h
RunCmd.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h Matrix.h
RunCmd.o: Context.h blacs.h D3vector.h Wavefunction.h UnitCell.h SlaterDet.h
RunCmd.o: Basis.h Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h blas.h
RunCmd.o: ConstraintSet.h ExtForceSet.h Control.h ChargeDensity.h
RunCmd.o: CPSampleStepper.h
SCANFunctional.o: SCANFunctional.h XCFunctional.h
SDAIonicStepper.o: MPIdata.h SDAIonicStepper.h IonicStepper.h Sample.h
SDAIonicStepper.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
SDAIonicStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
SDAIonicStepper.o: Wavefunction.h Control.h Species.h LineMinimizer.h
SDCellStepper.o: MPIdata.h SDCellStepper.h CellStepper.h Sample.h Context.h
SDCellStepper.o: blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h
SDCellStepper.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
SDCellStepper.o: Control.h
SDIonicStepper.o: SDIonicStepper.h IonicStepper.h Sample.h MPIdata.h
SDIonicStepper.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
SDIonicStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
SDIonicStepper.o: Wavefunction.h Control.h Species.h
SDWavefunctionStepper.o: SDWavefunctionStepper.h WavefunctionStepper.h
SDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
SDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
SDWavefunctionStepper.o: Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h
SDWavefunctionStepper.o: blas.h ConstraintSet.h ExtForceSet.h Control.h
SampleHandler.o: SampleHandler.h StructureHandler.h Sample.h MPIdata.h
SampleHandler.o: Context.h blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h
SampleHandler.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
SampleHandler.o: Wavefunction.h Control.h AtomSetHandler.h
SampleHandler.o: WavefunctionHandler.h SlaterDet.h Basis.h Matrix.h Timer.h
SampleHandler.o: StrX.h
SampleReader.o: Sample.h MPIdata.h Context.h blacs.h AtomSet.h Atom.h
SampleReader.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SampleReader.o: ExtForceSet.h Wavefunction.h Control.h SampleReader.h
SampleReader.o: XMLGFPreprocessor.h Matrix.h Timer.h SampleHandler.h
SampleReader.o: StructureHandler.h StructuredDocumentHandler.h StrX.h
SampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h StructureFactor.h
SampleStepper.o: ElectricEnthalpy.h Matrix.h Context.h blacs.h D3vector.h
SampleStepper.o: Wavefunction.h UnitCell.h SlaterDet.h Basis.h Sample.h
SampleStepper.o: MPIdata.h AtomSet.h Atom.h D3tensor.h blas.h ConstraintSet.h
SampleStepper.o: ExtForceSet.h Control.h ChargeDensity.h Species.h
SampleWriter.o: MPIdata.h SampleWriter.h Context.h blacs.h Sample.h AtomSet.h
SampleWriter.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
SampleWriter.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
SampleWriter.o: qbox_xmlns.h Timer.h SharedFilePtr.h
SaveCmd.o: SaveCmd.h UserInterface.h Sample.h MPIdata.h Context.h blacs.h
SaveCmd.o: AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
SaveCmd.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
SaveCmd.o: SampleWriter.h isodate.h release.h
SlaterDet.o: SlaterDet.h Context.h blacs.h Basis.h D3vector.h UnitCell.h
SlaterDet.o: Matrix.h Timer.h FourierTransform.h BasisMapping.h blas.h
SlaterDet.o: Base64Transcoder.h SharedFilePtr.h MPIdata.h
Species.o: Species.h spline.h sinft.h
SpeciesCmd.o: SpeciesCmd.h UserInterface.h Sample.h MPIdata.h Context.h
SpeciesCmd.o: blacs.h AtomSet.h Atom.h D3vector.h UnitCell.h D3tensor.h
SpeciesCmd.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
SpeciesCmd.o: SpeciesReader.h Species.h
SpeciesHandler.o: SpeciesHandler.h StructureHandler.h Species.h StrX.h
SpeciesReader.o: Species.h SpeciesReader.h StructuredDocumentHandler.h StrX.h
SpeciesReader.o: StructureHandler.h SpeciesHandler.h
SpectrumCmd.o: SpectrumCmd.h UserInterface.h Context.h blacs.h
SpectrumCmd.o: ChargeDensity.h Timer.h D3vector.h EnergyFunctional.h
SpectrumCmd.o: StructureFactor.h ElectricEnthalpy.h Matrix.h Wavefunction.h
SpectrumCmd.o: UnitCell.h SlaterDet.h Basis.h Sample.h MPIdata.h AtomSet.h
SpectrumCmd.o: Atom.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
SpectrumCmd.o: Control.h MLWFTransform.h BasisMapping.h
StructureFactor.o: StructureFactor.h Basis.h D3vector.h UnitCell.h
StructuredDocumentHandler.o: StructuredDocumentHandler.h StrX.h
StructuredDocumentHandler.o: StructureHandler.h
TorsionConstraint.o: TorsionConstraint.h Constraint.h D3vector.h AtomSet.h
TorsionConstraint.o: Atom.h UnitCell.h D3tensor.h blas.h Species.h
UnitCell.o: UnitCell.h D3vector.h
UserInterface.o: UserInterface.h qbox_xmlns.h
VWNFunctional.o: VWNFunctional.h XCFunctional.h
Wavefunction.o: MPIdata.h Wavefunction.h D3vector.h UnitCell.h SlaterDet.h
Wavefunction.o: Context.h blacs.h Basis.h Matrix.h Timer.h FourierTransform.h
Wavefunction.o: BasisMapping.h jacobi.h SharedFilePtr.h cout0.h
WavefunctionHandler.o: WavefunctionHandler.h StructureHandler.h UnitCell.h
WavefunctionHandler.o: D3vector.h Wavefunction.h SlaterDet.h Context.h
WavefunctionHandler.o: blacs.h Basis.h Matrix.h Timer.h FourierTransform.h
WavefunctionHandler.o: BasisMapping.h SampleReader.h MPIdata.h StrX.h
XCOperator.o: XCOperator.h Sample.h MPIdata.h Context.h blacs.h AtomSet.h
XCOperator.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
XCOperator.o: ExtForceSet.h Wavefunction.h Control.h ChargeDensity.h Timer.h
XCOperator.o: XCPotential.h ExchangeOperator.h SlaterDet.h Basis.h Matrix.h
XCOperator.o: FourierTransform.h BasisMapping.h HSEFunctional.h
XCOperator.o: XCFunctional.h RSHFunctional.h
XCPotential.o: XCPotential.h ChargeDensity.h Timer.h Context.h blacs.h
XCPotential.o: D3vector.h LDAFunctional.h XCFunctional.h VWNFunctional.h
XCPotential.o: PBEFunctional.h BLYPFunctional.h HSEFunctional.h
XCPotential.o: RSHFunctional.h B3LYPFunctional.h BHandHLYPFunctional.h
XCPotential.o: SCANFunctional.h Sample.h MPIdata.h AtomSet.h Atom.h
XCPotential.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
XCPotential.o: Wavefunction.h Control.h SlaterDet.h Basis.h Matrix.h
XCPotential.o: FourierTransform.h BasisMapping.h
XMLGFPreprocessor.o: Timer.h Context.h blacs.h Base64Transcoder.h Matrix.h
XMLGFPreprocessor.o: XMLGFPreprocessor.h
cout0.o: MPIdata.h cout0.h
isodate.o: isodate.h
jacobi.o: blacs.h Context.h Matrix.h blas.h
jade.o: blacs.h Context.h Matrix.h blas.h Timer.h
qb.o: isodate.h release.h qbox_xmlns.h uuid_str.h MPIdata.h Timer.h
qb.o: UserInterface.h Sample.h Context.h blacs.h AtomSet.h Atom.h D3vector.h
qb.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
qb.o: Wavefunction.h Control.h AngleCmd.h AtomCmd.h ComputeMLWFCmd.h
qb.o: MLWFTransform.h BasisMapping.h ConstraintCmd.h DistanceCmd.h
qb.o: ExtForceCmd.h FoldInWsCmd.h HelpCmd.h KpointCmd.h ListAtomsCmd.h
qb.o: ListSpeciesCmd.h LoadCmd.h MoveCmd.h PartialChargeCmd.h PlotCmd.h
qb.o: PrintCmd.h QuitCmd.h RandomizeRCmd.h RandomizeVCmd.h RandomizeWfCmd.h
qb.o: ResetRotationCmd.h ResetVcmCmd.h RescaleVCmd.h ResponseCmd.h RseedCmd.h
qb.o: RunCmd.h SaveCmd.h SetCmd.h SetVelocityCmd.h SpeciesCmd.h SpectrumCmd.h
qb.o: StatusCmd.h ChargeDensity.h FourierTransform.h StrainCmd.h TorsionCmd.h
qb.o: BisectionCmd.h Bisection.h SlaterDet.h Basis.h Matrix.h AlphaPBE0.h
qb.o: AlphaRSH.h AtomsDyn.h BetaRSH.h BlHF.h BtHF.h Cell.h CellDyn.h
qb.o: CellLock.h CellMass.h ChargeMixCoeff.h ChargeMixNdim.h ChargeMixRcut.h
qb.o: Debug.h Dspin.h Ecut.h Ecutprec.h Ecuts.h Efield.h ForceTol.h
qb.o: Polarization.h Emass.h ExtStress.h FermiTemp.h IterCmd.h
qb.o: IterCmdPeriod.h Dt.h MuRSH.h Nempty.h NetCharge.h Nspin.h Occ.h
qb.o: RefCell.h ScfTol.h Stress.h StressTol.h Thermostat.h ThTemp.h ThTime.h
qb.o: ThWidth.h Vext.h ExternalPotential.h WfDiag.h WfDyn.h Xc.h
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
testChargeDensity.o: MPIdata.h Context.h blacs.h Wavefunction.h D3vector.h
testChargeDensity.o: UnitCell.h ChargeDensity.h Timer.h SlaterDet.h Basis.h
testChargeDensity.o: Matrix.h FourierTransform.h BasisMapping.h
testContext.o: Context.h blacs.h
testEnergyFunctional.o: MPIdata.h Sample.h Context.h blacs.h AtomSet.h Atom.h
testEnergyFunctional.o: D3vector.h UnitCell.h D3tensor.h blas.h
testEnergyFunctional.o: ConstraintSet.h ExtForceSet.h Wavefunction.h
testEnergyFunctional.o: Control.h ChargeDensity.h Timer.h EnergyFunctional.h
testEnergyFunctional.o: StructureFactor.h ElectricEnthalpy.h Matrix.h
testEnergyFunctional.o: SlaterDet.h Basis.h
testFourierTransform.o: Basis.h D3vector.h UnitCell.h FourierTransform.h
testFourierTransform.o: Timer.h BasisMapping.h
testFunction3d.o: Function3d.h D3vector.h
testLDAFunctional.o: LDAFunctional.h XCFunctional.h
testLineMinimizer.o: LineMinimizer.h
testMPIdata.o: MPIdata.h
testMatrix.o: Timer.h Context.h blacs.h Matrix.h
testSample.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h UnitCell.h
testSample.o: Matrix.h Timer.h Sample.h MPIdata.h AtomSet.h Atom.h D3tensor.h
testSample.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
testSampleReader.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h
testSampleReader.o: UnitCell.h Matrix.h Timer.h Sample.h MPIdata.h AtomSet.h
testSampleReader.o: Atom.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
testSampleReader.o: Wavefunction.h Control.h SampleReader.h
testSlaterDet.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h UnitCell.h
testSlaterDet.o: Matrix.h Timer.h FourierTransform.h BasisMapping.h MPIdata.h
testSpecies.o: Species.h SpeciesReader.h
testUnitCell.o: UnitCell.h D3vector.h
testVWN.o: VWNFunctional.h XCFunctional.h LDAFunctional.h
testWavefunction.o: MPIdata.h Context.h blacs.h Wavefunction.h D3vector.h
testWavefunction.o: UnitCell.h SlaterDet.h Basis.h Matrix.h Timer.h
testXCFunctional.o: LDAFunctional.h XCFunctional.h PBEFunctional.h Timer.h
testXMLGFPreprocessor.o: Context.h blacs.h Matrix.h XMLGFPreprocessor.h
test_vext.o: Function3d.h D3vector.h
testjacobi.o: Timer.h Context.h blacs.h Matrix.h jacobi.h
testjade.o: Timer.h Context.h blacs.h Matrix.h jade.h
uuid_str.o: uuid_str.h
