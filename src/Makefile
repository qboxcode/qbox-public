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
TESTEXECS=testMatrix testBasis \
testSlaterDet testEnergyFunctional testSample \
testChargeDensity  testFourierTransform testSpecies testContext \
testFunction3d

OBJECTS=qb.o AtomSet.o Atom.o Species.o \
        Wavefunction.o SlaterDet.o \
        EnergyFunctional.o SampleStepper.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o Context.o \
        sinft.o spline.o UnitCell.o StructureFactor.o \
        ChargeDensity.o UserInterface.o RunCmd.o \
        LoadCmd.o SaveCmd.o \
        SpeciesCmd.o SpeciesReader.o SpeciesHandler.o \
        XCOperator.o ExchangeOperator.o Bisection.o \
	XCPotential.o LDAFunctional.o VWNFunctional.o \
        PBEFunctional.o BLYPFunctional.o SCANFunctional.o \
        B3LYPFunctional.o  BHandHLYPFunctional.o \
	ExponentialIntegral.o HSEFunctional.o \
        RSHFunctional.o \
	NonLocalPotential.o SampleReader.o StructuredDocumentHandler.o \
	SampleHandler.o AtomSetHandler.o WavefunctionHandler.o \
        XMLGFPreprocessor.o Base64Transcoder.o \
        CPSampleStepper.o BOSampleStepper.o \
        SDWavefunctionStepper.o \
        MDWavefunctionStepper.o SDIonicStepper.o MDIonicStepper.o \
        BMDIonicStepper.o \
        PSDWavefunctionStepper.o PSDAWavefunctionStepper.o \
        JDWavefunctionStepper.o CellStepper.o SDCellStepper.o \
	CGCellStepper.o ConfinementPotential.o Preconditioner.o \
        release.o qbox_xmlns.o isodate.o \
        AndersonMixer.o SDAIonicStepper.o CGIonicStepper.o \
        ConstraintSet.o Constraint.o PositionConstraint.o DistanceConstraint.o \
        AngleConstraint.o TorsionConstraint.o jacobi.o \
        SampleWriter.o ComputeMLWFCmd.o MLWFTransform.o \
        jade.o PlotCmd.o \
        ExtForceSet.o ExtForce.o AtomicExtForce.o PairExtForce.o \
        GlobalExtForce.o \
        uuid_str.o sampling.o CGOptimizer.o LineMinimizer.o \
        ElectricEnthalpy.o PartialChargeCmd.o \
        ExternalPotential.o ResponseCmd.o \
        Function3d.o Function3dHandler.o \
        SpectrumCmd.o \
        $(PLTOBJECTS)
CXXFLAGS += -DTARGET='"$(TARGET)"'
# to include VERSION info in release string, use:
# env VERSION=`git describe` make
CXXFLAGS += -DVERSION='"$(VERSION)"'
 $(EXEC):    $(OBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 lib:   $(OBJECTS)
	ar cr libqb.a $^
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
        BHandHLYPFunctional.o \
	ExponentialIntegral.o HSEFunctional.o RSHFunctional.o \
        ConfinementPotential.o ElectricEnthalpy.o MLWFTransform.o \
        jade.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSlaterDet: testSlaterDet.o SlaterDet.o FourierTransform.o \
        Basis.o BasisMapping.o UnitCell.o Matrix.o Context.o \
        Base64Transcoder.o
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
	makedepend -Y -D$(PLT) $(PLTFLAGS) *.[cCh]
#makedepend -D$(PLT) -I/usr/include/CC *.[cC]
#------------------------------------------------------------------------------
#  Cleanup object files
 clean :
	rm -f *.o
	rm -f $(EXEC) $(TESTEXECS)
	rm -rf ti_files
	rm -f libqb.a
#------------------------------------------------------------------------------
 ctags :
	ctags -o tags *.[Ch]
#------------------------------------------------------------------------------
 html :
	enscript -Ecpp --color -Whtml --toc -pqbsrc.html *.h *.C
#------------------------------------------------------------------------------
 ps :
	enscript -Ecpp --toc -2rGj -pqb.ps *.h *.C Makefile *.mk
#------------------------------------------------------------------------------
 tar :
	tar czf qbsrc.tgz *.[Ch] Makefile *.mk notes
#------------------------------------------------------------------------------
# DO NOT DELETE

AlphaPBE0.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
AlphaPBE0.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
AlphaPBE0.o: Wavefunction.h Control.h
AlphaRSH.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
AlphaRSH.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
AlphaRSH.o: Control.h
AndersonMixer.o: AndersonMixer.h blas.h
AngleCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
AngleCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
AngleCmd.o: ExtForceSet.h Wavefunction.h Control.h
AngleConstraint.o: AngleConstraint.h Constraint.h D3vector.h AtomSet.h
AngleConstraint.o: Context.h blacs.h Atom.h UnitCell.h D3tensor.h blas.h
AngleConstraint.o: Species.h
AngleConstraint.o: Constraint.h D3vector.h
Atom.o: Atom.h D3vector.h
Atom.o: D3vector.h
AtomCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
AtomCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
AtomCmd.o: ExtForceSet.h Wavefunction.h Control.h
AtomSet.o: AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
AtomSet.o: D3tensor.h blas.h Species.h NameOf.h sampling.h
AtomSet.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
AtomSetHandler.o: AtomSetHandler.h StructureHandler.h D3vector.h AtomSet.h
AtomSetHandler.o: Context.h blacs.h Atom.h UnitCell.h D3tensor.h blas.h
AtomSetHandler.o: Species.h SpeciesHandler.h SpeciesReader.h StrX.h
AtomSetHandler.o: SampleReader.h
AtomSetHandler.o: StructureHandler.h D3vector.h
AtomicExtForce.o: AtomicExtForce.h ExtForce.h D3vector.h AtomSet.h Context.h
AtomicExtForce.o: blacs.h Atom.h UnitCell.h D3tensor.h blas.h Species.h
AtomicExtForce.o: ExtForce.h D3vector.h
AtomsDyn.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
AtomsDyn.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
AtomsDyn.o: Control.h
B3LYPFunctional.o: B3LYPFunctional.h XCFunctional.h BLYPFunctional.h
B3LYPFunctional.o: VWNFunctional.h
B3LYPFunctional.o: XCFunctional.h
BHandHLYPFunctional.o: BHandHLYPFunctional.h XCFunctional.h BLYPFunctional.h
BHandHLYPFunctional.o: XCFunctional.h
BLYPFunctional.o: BLYPFunctional.h XCFunctional.h
BLYPFunctional.o: XCFunctional.h
BMDIonicStepper.o: BMDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
BMDIonicStepper.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
BMDIonicStepper.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
BMDIonicStepper.o: Control.h Species.h
BMDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h blacs.h Atom.h
BMDIonicStepper.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
BMDIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
BOSampleStepper.o: BOSampleStepper.h SampleStepper.h Timer.h
BOSampleStepper.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
BOSampleStepper.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
BOSampleStepper.o: UnitCell.h SlaterDet.h Basis.h Sample.h AtomSet.h Atom.h
BOSampleStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Control.h
BOSampleStepper.o: ChargeDensity.h WavefunctionStepper.h
BOSampleStepper.o: SDWavefunctionStepper.h PSDWavefunctionStepper.h
BOSampleStepper.o: PSDAWavefunctionStepper.h JDWavefunctionStepper.h
BOSampleStepper.o: UserInterface.h Preconditioner.h SDIonicStepper.h
BOSampleStepper.o: IonicStepper.h Species.h SDAIonicStepper.h LineMinimizer.h
BOSampleStepper.o: CGIonicStepper.h CGOptimizer.h MDIonicStepper.h
BOSampleStepper.o: BMDIonicStepper.h SDCellStepper.h CellStepper.h
BOSampleStepper.o: CGCellStepper.h AndersonMixer.h MLWFTransform.h
BOSampleStepper.o: BasisMapping.h
BOSampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h
BOSampleStepper.o: StructureFactor.h ElectricEnthalpy.h Matrix.h Context.h
BOSampleStepper.o: blacs.h D3vector.h Wavefunction.h UnitCell.h SlaterDet.h
BOSampleStepper.o: Basis.h Sample.h AtomSet.h Atom.h D3tensor.h blas.h
BOSampleStepper.o: ConstraintSet.h ExtForceSet.h Control.h ChargeDensity.h
Base64Transcoder.o: Base64Transcoder.h
Basis.o: Basis.h D3vector.h UnitCell.h
Basis.o: D3vector.h UnitCell.h
BasisMapping.o: Basis.h D3vector.h UnitCell.h Context.h blacs.h
BasisMapping.o: BasisMapping.h
BetaRSH.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
BetaRSH.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
BetaRSH.o: Control.h
Bisection.o: Bisection.h Context.h blacs.h SlaterDet.h Basis.h D3vector.h
Bisection.o: UnitCell.h Matrix.h Timer.h jade.h FourierTransform.h
Bisection.o: BasisMapping.h
Bisection.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h UnitCell.h
Bisection.o: Matrix.h Timer.h
BisectionCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
BisectionCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
BisectionCmd.o: ExtForceSet.h Wavefunction.h Control.h Bisection.h
BisectionCmd.o: SlaterDet.h Basis.h Matrix.h Timer.h
BlHF.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
BlHF.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
BlHF.o: Control.h
BtHF.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
BtHF.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
BtHF.o: Control.h
CGCellStepper.o: CGCellStepper.h CellStepper.h Sample.h AtomSet.h Context.h
CGCellStepper.o: blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
CGCellStepper.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
CGCellStepper.o: CGOptimizer.h LineMinimizer.h
CGCellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h blacs.h Atom.h
CGCellStepper.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
CGCellStepper.o: ExtForceSet.h Wavefunction.h Control.h CGOptimizer.h
CGCellStepper.o: LineMinimizer.h
CGIonicStepper.o: CGIonicStepper.h IonicStepper.h Sample.h AtomSet.h
CGIonicStepper.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
CGIonicStepper.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
CGIonicStepper.o: Control.h Species.h CGOptimizer.h LineMinimizer.h
CGIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h blacs.h Atom.h
CGIonicStepper.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
CGIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
CGIonicStepper.o: CGOptimizer.h LineMinimizer.h
CGOptimizer.o: CGOptimizer.h LineMinimizer.h blas.h
CGOptimizer.o: LineMinimizer.h
CPSampleStepper.o: CPSampleStepper.h SampleStepper.h Timer.h
CPSampleStepper.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
CPSampleStepper.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
CPSampleStepper.o: UnitCell.h SlaterDet.h Basis.h Sample.h AtomSet.h Atom.h
CPSampleStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Control.h
CPSampleStepper.o: ChargeDensity.h MDWavefunctionStepper.h
CPSampleStepper.o: WavefunctionStepper.h MDIonicStepper.h IonicStepper.h
CPSampleStepper.o: Species.h SDCellStepper.h CellStepper.h CGCellStepper.h
CPSampleStepper.o: CGOptimizer.h LineMinimizer.h
CPSampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h
CPSampleStepper.o: StructureFactor.h ElectricEnthalpy.h Matrix.h Context.h
CPSampleStepper.o: blacs.h D3vector.h Wavefunction.h UnitCell.h SlaterDet.h
CPSampleStepper.o: Basis.h Sample.h AtomSet.h Atom.h D3tensor.h blas.h
CPSampleStepper.o: ConstraintSet.h ExtForceSet.h Control.h ChargeDensity.h
Cell.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Cell.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Cell.o: Control.h
CellDyn.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
CellDyn.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
CellDyn.o: Control.h SlaterDet.h Basis.h Matrix.h Timer.h
CellLock.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
CellLock.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
CellLock.o: Control.h
CellMass.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
CellMass.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
CellMass.o: Control.h
CellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h blacs.h Atom.h
CellStepper.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
CellStepper.o: ExtForceSet.h Wavefunction.h Control.h
CellStepper.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
CellStepper.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
CellStepper.o: Wavefunction.h Control.h
ChargeDensity.o: ChargeDensity.h Timer.h Context.h blacs.h D3vector.h Basis.h
ChargeDensity.o: UnitCell.h Wavefunction.h FourierTransform.h BasisMapping.h
ChargeDensity.o: SlaterDet.h Matrix.h blas.h
ChargeDensity.o: Timer.h Context.h blacs.h D3vector.h
ChargeMixCoeff.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
ChargeMixCoeff.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
ChargeMixCoeff.o: Wavefunction.h Control.h
ChargeMixNdim.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
ChargeMixNdim.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
ChargeMixNdim.o: Wavefunction.h Control.h
ChargeMixRcut.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
ChargeMixRcut.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
ChargeMixRcut.o: Wavefunction.h Control.h
ComputeMLWFCmd.o: ComputeMLWFCmd.h UserInterface.h Sample.h AtomSet.h
ComputeMLWFCmd.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
ComputeMLWFCmd.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ComputeMLWFCmd.o: Control.h MLWFTransform.h BasisMapping.h SlaterDet.h
ComputeMLWFCmd.o: Basis.h Matrix.h Timer.h
ComputeMLWFCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
ComputeMLWFCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
ComputeMLWFCmd.o: ExtForceSet.h Wavefunction.h Control.h MLWFTransform.h
ComputeMLWFCmd.o: BasisMapping.h
ConfinementPotential.o: ConfinementPotential.h Basis.h D3vector.h UnitCell.h
Constraint.o: Constraint.h
ConstraintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
ConstraintCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
ConstraintCmd.o: ExtForceSet.h Wavefunction.h Control.h
ConstraintSet.o: ConstraintSet.h PositionConstraint.h Constraint.h
ConstraintSet.o: DistanceConstraint.h AngleConstraint.h D3vector.h
ConstraintSet.o: TorsionConstraint.h Atom.h AtomSet.h Context.h blacs.h
ConstraintSet.o: UnitCell.h D3tensor.h blas.h
Context.o: Context.h blacs.h
Context.o: blacs.h
Control.o: D3vector.h
D3tensor.o: blas.h D3vector.h
Debug.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Debug.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Debug.o: Control.h
DistanceCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
DistanceCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
DistanceCmd.o: ExtForceSet.h Wavefunction.h Control.h
DistanceConstraint.o: DistanceConstraint.h Constraint.h AtomSet.h Context.h
DistanceConstraint.o: blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
DistanceConstraint.o: Species.h
DistanceConstraint.o: Constraint.h
Dspin.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Dspin.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Dspin.o: Control.h
Dt.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Dt.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Dt.o: Control.h
Ecut.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Ecut.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Ecut.o: Control.h
Ecutprec.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Ecutprec.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Ecutprec.o: Control.h
Ecuts.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Ecuts.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Ecuts.o: Control.h
Efield.o: D3vector.h Sample.h AtomSet.h Context.h blacs.h Atom.h UnitCell.h
Efield.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Efield.o: Control.h
ElectricEnthalpy.o: Timer.h Context.h blacs.h Matrix.h Sample.h AtomSet.h
ElectricEnthalpy.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
ElectricEnthalpy.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ElectricEnthalpy.o: ElectricEnthalpy.h SlaterDet.h Basis.h MLWFTransform.h
ElectricEnthalpy.o: BasisMapping.h ChargeDensity.h FourierTransform.h
ElectricEnthalpy.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
ElectricEnthalpy.o: UnitCell.h SlaterDet.h Basis.h Timer.h Sample.h AtomSet.h
ElectricEnthalpy.o: Atom.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
ElectricEnthalpy.o: Control.h
Emass.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Emass.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Emass.o: Control.h
EnergyFunctional.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h
EnergyFunctional.o: Matrix.h Context.h blacs.h D3vector.h Wavefunction.h
EnergyFunctional.o: UnitCell.h SlaterDet.h Basis.h Timer.h Sample.h AtomSet.h
EnergyFunctional.o: Atom.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
EnergyFunctional.o: Control.h Species.h ChargeDensity.h FourierTransform.h
EnergyFunctional.o: BasisMapping.h XCOperator.h NonLocalPotential.h
EnergyFunctional.o: ConfinementPotential.h ExternalPotential.h
EnergyFunctional.o: StructureFactor.h ElectricEnthalpy.h Matrix.h Context.h
EnergyFunctional.o: blacs.h D3vector.h Wavefunction.h UnitCell.h SlaterDet.h
EnergyFunctional.o: Basis.h Timer.h Sample.h AtomSet.h Atom.h D3tensor.h
EnergyFunctional.o: blas.h ConstraintSet.h ExtForceSet.h Control.h
ExchangeOperator.o: VectorLess.h ExchangeOperator.h Sample.h AtomSet.h
ExchangeOperator.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
ExchangeOperator.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ExchangeOperator.o: Control.h SlaterDet.h Basis.h Matrix.h Timer.h
ExchangeOperator.o: FourierTransform.h BasisMapping.h Bisection.h
ExchangeOperator.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
ExchangeOperator.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h
ExchangeOperator.o: ExtForceSet.h Wavefunction.h Control.h SlaterDet.h
ExchangeOperator.o: Basis.h Matrix.h Timer.h FourierTransform.h
ExchangeOperator.o: BasisMapping.h
ExponentialIntegral.o: ExponentialIntegral.h
ExtForce.o: ExtForce.h D3vector.h
ExtForce.o: D3vector.h
ExtForceCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
ExtForceCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
ExtForceCmd.o: ExtForceSet.h Wavefunction.h Control.h
ExtForceSet.o: ExtForceSet.h AtomicExtForce.h ExtForce.h D3vector.h
ExtForceSet.o: PairExtForce.h GlobalExtForce.h Atom.h AtomSet.h Context.h
ExtForceSet.o: blacs.h UnitCell.h D3tensor.h blas.h
ExtStress.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
ExtStress.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
ExtStress.o: Wavefunction.h Control.h
ExternalPotential.o: Basis.h D3vector.h UnitCell.h ExternalPotential.h
ExternalPotential.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3tensor.h
ExternalPotential.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ExternalPotential.o: Control.h ChargeDensity.h Timer.h FourierTransform.h
ExternalPotential.o: BasisMapping.h Function3d.h Base64Transcoder.h
ExternalPotential.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
ExternalPotential.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h
ExternalPotential.o: ExtForceSet.h Wavefunction.h Control.h ChargeDensity.h
ExternalPotential.o: Timer.h
FermiTemp.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
FermiTemp.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
FermiTemp.o: Wavefunction.h Control.h
FoldInWsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
FoldInWsCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
FoldInWsCmd.o: ExtForceSet.h Wavefunction.h Control.h
ForceTol.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
ForceTol.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ForceTol.o: Control.h
FourierTransform.o: FourierTransform.h Timer.h BasisMapping.h Basis.h
FourierTransform.o: D3vector.h UnitCell.h blas.h
FourierTransform.o: Timer.h BasisMapping.h
Function3d.o: Function3d.h D3vector.h Function3dHandler.h StrX.h
Function3d.o: Base64Transcoder.h qbox_xmlns.h Timer.h
Function3d.o: D3vector.h
Function3dHandler.o: Function3d.h D3vector.h Function3dHandler.h StrX.h
Function3dHandler.o: Base64Transcoder.h Timer.h
Function3dHandler.o: Function3d.h D3vector.h StrX.h
GlobalExtForce.o: GlobalExtForce.h ExtForce.h D3vector.h AtomSet.h Context.h
GlobalExtForce.o: blacs.h Atom.h UnitCell.h D3tensor.h blas.h Species.h
GlobalExtForce.o: ExtForce.h D3vector.h
HSEFunctional.o: HSEFunctional.h XCFunctional.h ExponentialIntegral.h
HSEFunctional.o: XCFunctional.h
HelpCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
HelpCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
HelpCmd.o: ExtForceSet.h Wavefunction.h Control.h
IonicStepper.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
IonicStepper.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
IonicStepper.o: Wavefunction.h Control.h Species.h
IterCmd.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
IterCmd.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
IterCmd.o: Control.h
IterCmdPeriod.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
IterCmdPeriod.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
IterCmdPeriod.o: Wavefunction.h Control.h
JDWavefunctionStepper.o: JDWavefunctionStepper.h WavefunctionStepper.h
JDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
JDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
JDWavefunctionStepper.o: EnergyFunctional.h StructureFactor.h
JDWavefunctionStepper.o: ElectricEnthalpy.h Sample.h AtomSet.h Atom.h
JDWavefunctionStepper.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
JDWavefunctionStepper.o: Control.h Preconditioner.h
JDWavefunctionStepper.o: WavefunctionStepper.h Timer.h Wavefunction.h
JDWavefunctionStepper.o: D3vector.h UnitCell.h
KpointCmd.o: UserInterface.h D3vector.h Sample.h AtomSet.h Context.h blacs.h
KpointCmd.o: Atom.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
KpointCmd.o: ExtForceSet.h Wavefunction.h Control.h
LDAFunctional.o: LDAFunctional.h XCFunctional.h
LDAFunctional.o: XCFunctional.h
LineMinimizer.o: LineMinimizer.h
ListAtomsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
ListAtomsCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
ListAtomsCmd.o: ExtForceSet.h Wavefunction.h Control.h
ListSpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
ListSpeciesCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
ListSpeciesCmd.o: ExtForceSet.h Wavefunction.h Control.h
LoadCmd.o: LoadCmd.h UserInterface.h Sample.h AtomSet.h Context.h blacs.h
LoadCmd.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
LoadCmd.o: ExtForceSet.h Wavefunction.h Control.h SampleReader.h
LoadCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
LoadCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
LoadCmd.o: ExtForceSet.h Wavefunction.h Control.h
MDIonicStepper.o: MDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
MDIonicStepper.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
MDIonicStepper.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
MDIonicStepper.o: Control.h Species.h sampling.h
MDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h blacs.h Atom.h
MDIonicStepper.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
MDIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
MDWavefunctionStepper.o: MDWavefunctionStepper.h WavefunctionStepper.h
MDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
MDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
MDWavefunctionStepper.o: Sample.h AtomSet.h Atom.h D3tensor.h blas.h
MDWavefunctionStepper.o: ConstraintSet.h ExtForceSet.h Control.h
MDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
MLWFTransform.o: MLWFTransform.h D3vector.h BasisMapping.h Basis.h UnitCell.h
MLWFTransform.o: SlaterDet.h Context.h blacs.h Matrix.h Timer.h jade.h blas.h
MLWFTransform.o: D3vector.h BasisMapping.h
Matrix.o: Context.h blacs.h Matrix.h
Matrix.o: Context.h blacs.h
MoveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
MoveCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
MoveCmd.o: ExtForceSet.h Wavefunction.h Control.h
MuRSH.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
MuRSH.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
MuRSH.o: Control.h
MultidistanceConstraint.o: MultidistanceConstraint.h Constraint.h AtomSet.h
MultidistanceConstraint.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h
MultidistanceConstraint.o: D3tensor.h blas.h Species.h
MultidistanceConstraint.o: Constraint.h
Nempty.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Nempty.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Nempty.o: Control.h
NetCharge.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
NetCharge.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
NetCharge.o: Wavefunction.h Control.h
NonLocalPotential.o: NonLocalPotential.h AtomSet.h Context.h blacs.h Atom.h
NonLocalPotential.o: D3vector.h UnitCell.h D3tensor.h blas.h Basis.h
NonLocalPotential.o: SlaterDet.h Matrix.h Timer.h Species.h
NonLocalPotential.o: AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
NonLocalPotential.o: D3tensor.h blas.h Basis.h SlaterDet.h Matrix.h Timer.h
Nrowmax.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Nrowmax.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Nrowmax.o: Control.h
Nspin.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Nspin.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Nspin.o: Control.h
Occ.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Occ.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Occ.o: Control.h
PBEFunctional.o: PBEFunctional.h XCFunctional.h
PBEFunctional.o: XCFunctional.h
PSDAWavefunctionStepper.o: PSDAWavefunctionStepper.h WavefunctionStepper.h
PSDAWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
PSDAWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
PSDAWavefunctionStepper.o: Preconditioner.h
PSDAWavefunctionStepper.o: WavefunctionStepper.h Timer.h Wavefunction.h
PSDAWavefunctionStepper.o: D3vector.h UnitCell.h
PSDWavefunctionStepper.o: PSDWavefunctionStepper.h WavefunctionStepper.h
PSDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
PSDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
PSDWavefunctionStepper.o: Preconditioner.h
PSDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
PairExtForce.o: PairExtForce.h ExtForce.h D3vector.h AtomSet.h Context.h
PairExtForce.o: blacs.h Atom.h UnitCell.h D3tensor.h blas.h Species.h
PairExtForce.o: ExtForce.h D3vector.h
PartialChargeCmd.o: PartialChargeCmd.h UserInterface.h Sample.h AtomSet.h
PartialChargeCmd.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
PartialChargeCmd.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
PartialChargeCmd.o: Control.h Basis.h ChargeDensity.h Timer.h
PartialChargeCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h
PartialChargeCmd.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
PartialChargeCmd.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
PlotCmd.o: PlotCmd.h UserInterface.h Sample.h AtomSet.h Context.h blacs.h
PlotCmd.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
PlotCmd.o: ExtForceSet.h Wavefunction.h Control.h isodate.h release.h
PlotCmd.o: SampleReader.h Basis.h EnergyFunctional.h StructureFactor.h
PlotCmd.o: ElectricEnthalpy.h Matrix.h SlaterDet.h Timer.h FourierTransform.h
PlotCmd.o: BasisMapping.h Species.h ChargeDensity.h
PlotCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
PlotCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
PlotCmd.o: ExtForceSet.h Wavefunction.h Control.h
Polarization.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
Polarization.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
Polarization.o: Wavefunction.h Control.h
PositionConstraint.o: PositionConstraint.h Constraint.h AtomSet.h Context.h
PositionConstraint.o: blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
PositionConstraint.o: Species.h
PositionConstraint.o: Constraint.h
Preconditioner.o: Preconditioner.h Basis.h D3vector.h UnitCell.h
Preconditioner.o: Wavefunction.h EnergyFunctional.h StructureFactor.h
Preconditioner.o: ElectricEnthalpy.h Matrix.h Context.h blacs.h SlaterDet.h
Preconditioner.o: Timer.h Sample.h AtomSet.h Atom.h D3tensor.h blas.h
Preconditioner.o: ConstraintSet.h ExtForceSet.h Control.h
Preconditioner.o: ConfinementPotential.h
PrintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
PrintCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
PrintCmd.o: ExtForceSet.h Wavefunction.h Control.h
QuitCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
QuitCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
QuitCmd.o: ExtForceSet.h Wavefunction.h Control.h
RSHFunctional.o: RSHFunctional.h XCFunctional.h ExponentialIntegral.h
RSHFunctional.o: XCFunctional.h
RandomizeRCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
RandomizeRCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
RandomizeRCmd.o: ExtForceSet.h Wavefunction.h Control.h
RandomizeVCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
RandomizeVCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
RandomizeVCmd.o: ExtForceSet.h Wavefunction.h Control.h
RandomizeWfCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
RandomizeWfCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
RandomizeWfCmd.o: ExtForceSet.h Wavefunction.h Control.h
RefCell.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
RefCell.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
RefCell.o: Control.h
RescaleVCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
RescaleVCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
RescaleVCmd.o: ExtForceSet.h Wavefunction.h Control.h
ResetRotationCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h
ResetRotationCmd.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
ResetRotationCmd.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ResetVcmCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
ResetVcmCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
ResetVcmCmd.o: ExtForceSet.h Wavefunction.h Control.h
ResponseCmd.o: BOSampleStepper.h SampleStepper.h Timer.h EnergyFunctional.h
ResponseCmd.o: StructureFactor.h ElectricEnthalpy.h Matrix.h Context.h
ResponseCmd.o: blacs.h D3vector.h Wavefunction.h UnitCell.h SlaterDet.h
ResponseCmd.o: Basis.h Sample.h AtomSet.h Atom.h D3tensor.h blas.h
ResponseCmd.o: ConstraintSet.h ExtForceSet.h Control.h ChargeDensity.h
ResponseCmd.o: ExternalPotential.h FourierTransform.h BasisMapping.h
ResponseCmd.o: ResponseCmd.h UserInterface.h release.h isodate.h Species.h
ResponseCmd.o: Function3d.h
ResponseCmd.o: UserInterface.h
RseedCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
RseedCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
RseedCmd.o: ExtForceSet.h Wavefunction.h Control.h
RunCmd.o: RunCmd.h UserInterface.h BOSampleStepper.h SampleStepper.h Timer.h
RunCmd.o: EnergyFunctional.h StructureFactor.h ElectricEnthalpy.h Matrix.h
RunCmd.o: Context.h blacs.h D3vector.h Wavefunction.h UnitCell.h SlaterDet.h
RunCmd.o: Basis.h Sample.h AtomSet.h Atom.h D3tensor.h blas.h ConstraintSet.h
RunCmd.o: ExtForceSet.h Control.h ChargeDensity.h CPSampleStepper.h
RunCmd.o: UserInterface.h
SCANFunctional.o: SCANFunctional.h XCFunctional.h
SCANFunctional.o: XCFunctional.h
SDAIonicStepper.o: SDAIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDAIonicStepper.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
SDAIonicStepper.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
SDAIonicStepper.o: Control.h Species.h LineMinimizer.h
SDAIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h blacs.h Atom.h
SDAIonicStepper.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SDAIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
SDAIonicStepper.o: LineMinimizer.h
SDCellStepper.o: SDCellStepper.h CellStepper.h Sample.h AtomSet.h Context.h
SDCellStepper.o: blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
SDCellStepper.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
SDCellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h blacs.h Atom.h
SDCellStepper.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SDCellStepper.o: ExtForceSet.h Wavefunction.h Control.h
SDIonicStepper.o: SDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDIonicStepper.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
SDIonicStepper.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
SDIonicStepper.o: Control.h Species.h
SDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h blacs.h Atom.h
SDIonicStepper.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SDIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
SDWavefunctionStepper.o: SDWavefunctionStepper.h WavefunctionStepper.h
SDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
SDWavefunctionStepper.o: SlaterDet.h Context.h blacs.h Basis.h Matrix.h
SDWavefunctionStepper.o: Sample.h AtomSet.h Atom.h D3tensor.h blas.h
SDWavefunctionStepper.o: ConstraintSet.h ExtForceSet.h Control.h
SDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
Sample.o: AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
Sample.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
SampleHandler.o: SampleHandler.h StructureHandler.h Sample.h AtomSet.h
SampleHandler.o: Context.h blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h
SampleHandler.o: blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
SampleHandler.o: Control.h AtomSetHandler.h WavefunctionHandler.h SlaterDet.h
SampleHandler.o: Basis.h Matrix.h Timer.h StrX.h
SampleHandler.o: StructureHandler.h
SampleReader.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
SampleReader.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
SampleReader.o: Wavefunction.h Control.h SampleReader.h XMLGFPreprocessor.h
SampleReader.o: Matrix.h Timer.h SampleHandler.h StructureHandler.h
SampleReader.o: StructuredDocumentHandler.h StrX.h
SampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h StructureFactor.h
SampleStepper.o: ElectricEnthalpy.h Matrix.h Context.h blacs.h D3vector.h
SampleStepper.o: Wavefunction.h UnitCell.h SlaterDet.h Basis.h Sample.h
SampleStepper.o: AtomSet.h Atom.h D3tensor.h blas.h ConstraintSet.h
SampleStepper.o: ExtForceSet.h Control.h ChargeDensity.h Species.h
SampleStepper.o: Timer.h EnergyFunctional.h StructureFactor.h
SampleStepper.o: ElectricEnthalpy.h Matrix.h Context.h blacs.h D3vector.h
SampleStepper.o: Wavefunction.h UnitCell.h SlaterDet.h Basis.h Sample.h
SampleStepper.o: AtomSet.h Atom.h D3tensor.h blas.h ConstraintSet.h
SampleStepper.o: ExtForceSet.h Control.h ChargeDensity.h
SampleWriter.o: SampleWriter.h Context.h blacs.h Sample.h AtomSet.h Atom.h
SampleWriter.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SampleWriter.o: ExtForceSet.h Wavefunction.h Control.h qbox_xmlns.h Timer.h
SampleWriter.o: SharedFilePtr.h
SampleWriter.o: Context.h blacs.h
SaveCmd.o: SaveCmd.h UserInterface.h Sample.h AtomSet.h Context.h blacs.h
SaveCmd.o: Atom.h D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SaveCmd.o: ExtForceSet.h Wavefunction.h Control.h SampleWriter.h isodate.h
SaveCmd.o: release.h
SaveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
SaveCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SaveCmd.o: ExtForceSet.h Wavefunction.h Control.h
ScfTol.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
ScfTol.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ScfTol.o: Control.h
SetCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
SetCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SetCmd.o: ExtForceSet.h Wavefunction.h Control.h
SetVelocityCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
SetVelocityCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SetVelocityCmd.o: ExtForceSet.h Wavefunction.h Control.h
SlaterDet.o: SlaterDet.h Context.h blacs.h Basis.h D3vector.h UnitCell.h
SlaterDet.o: Matrix.h Timer.h FourierTransform.h BasisMapping.h blas.h
SlaterDet.o: Base64Transcoder.h SharedFilePtr.h
SlaterDet.o: Context.h blacs.h Basis.h D3vector.h UnitCell.h Matrix.h Timer.h
Species.o: Species.h spline.h sinft.h
SpeciesCmd.o: SpeciesCmd.h UserInterface.h Sample.h AtomSet.h Context.h
SpeciesCmd.o: blacs.h Atom.h D3vector.h UnitCell.h D3tensor.h blas.h
SpeciesCmd.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
SpeciesCmd.o: SpeciesReader.h Species.h
SpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
SpeciesCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
SpeciesCmd.o: ExtForceSet.h Wavefunction.h Control.h
SpeciesHandler.o: SpeciesHandler.h StructureHandler.h Species.h StrX.h
SpeciesHandler.o: StructureHandler.h
SpeciesReader.o: Species.h SpeciesReader.h StructuredDocumentHandler.h StrX.h
SpeciesReader.o: StructureHandler.h SpeciesHandler.h
SpectrumCmd.o: SpectrumCmd.h UserInterface.h Context.h blacs.h
SpectrumCmd.o: ChargeDensity.h Timer.h D3vector.h EnergyFunctional.h
SpectrumCmd.o: StructureFactor.h ElectricEnthalpy.h Matrix.h Wavefunction.h
SpectrumCmd.o: UnitCell.h SlaterDet.h Basis.h Sample.h AtomSet.h Atom.h
SpectrumCmd.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Control.h
SpectrumCmd.o: MLWFTransform.h BasisMapping.h
SpectrumCmd.o: UserInterface.h
StatusCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
StatusCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
StatusCmd.o: ExtForceSet.h Wavefunction.h Control.h ChargeDensity.h Timer.h
StatusCmd.o: FourierTransform.h BasisMapping.h
StrainCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
StrainCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
StrainCmd.o: ExtForceSet.h Wavefunction.h Control.h
Stress.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Stress.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Stress.o: Control.h
StressTol.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
StressTol.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
StressTol.o: Wavefunction.h Control.h
StructureFactor.o: StructureFactor.h Basis.h D3vector.h UnitCell.h
StructuredDocumentHandler.o: StructuredDocumentHandler.h StrX.h
StructuredDocumentHandler.o: StructureHandler.h
StructuredDocumentHandler.o: StrX.h StructureHandler.h
ThTemp.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
ThTemp.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ThTemp.o: Control.h
ThTime.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
ThTime.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ThTime.o: Control.h
ThWidth.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
ThWidth.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ThWidth.o: Control.h
Thermostat.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h
Thermostat.o: UnitCell.h D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
Thermostat.o: Wavefunction.h Control.h
TorsionCmd.o: UserInterface.h Sample.h AtomSet.h Context.h blacs.h Atom.h
TorsionCmd.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
TorsionCmd.o: ExtForceSet.h Wavefunction.h Control.h
TorsionConstraint.o: TorsionConstraint.h Constraint.h D3vector.h AtomSet.h
TorsionConstraint.o: Context.h blacs.h Atom.h UnitCell.h D3tensor.h blas.h
TorsionConstraint.o: Species.h
TorsionConstraint.o: Constraint.h D3vector.h
UnitCell.o: UnitCell.h D3vector.h
UnitCell.o: D3vector.h
UserInterface.o: UserInterface.h qbox_xmlns.h
VWNFunctional.o: VWNFunctional.h XCFunctional.h
VWNFunctional.o: XCFunctional.h
Vext.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Vext.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Vext.o: Control.h ExternalPotential.h ChargeDensity.h Timer.h
Wavefunction.o: Wavefunction.h D3vector.h UnitCell.h SlaterDet.h Context.h
Wavefunction.o: blacs.h Basis.h Matrix.h Timer.h FourierTransform.h
Wavefunction.o: BasisMapping.h jacobi.h SharedFilePtr.h
Wavefunction.o: D3vector.h UnitCell.h
WavefunctionHandler.o: WavefunctionHandler.h StructureHandler.h UnitCell.h
WavefunctionHandler.o: D3vector.h Wavefunction.h SlaterDet.h Context.h
WavefunctionHandler.o: blacs.h Basis.h Matrix.h Timer.h FourierTransform.h
WavefunctionHandler.o: BasisMapping.h SampleReader.h StrX.h
WavefunctionHandler.o: StructureHandler.h UnitCell.h D3vector.h
WavefunctionHandler.o: Wavefunction.h SlaterDet.h Context.h blacs.h Basis.h
WavefunctionHandler.o: Matrix.h Timer.h
WavefunctionStepper.o: Timer.h
WfDiag.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
WfDiag.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
WfDiag.o: Control.h
WfDyn.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
WfDyn.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
WfDyn.o: Control.h SlaterDet.h Basis.h Matrix.h Timer.h
XCOperator.o: XCOperator.h Sample.h AtomSet.h Context.h blacs.h Atom.h
XCOperator.o: D3vector.h UnitCell.h D3tensor.h blas.h ConstraintSet.h
XCOperator.o: ExtForceSet.h Wavefunction.h Control.h ChargeDensity.h Timer.h
XCOperator.o: XCPotential.h ExchangeOperator.h SlaterDet.h Basis.h Matrix.h
XCOperator.o: FourierTransform.h BasisMapping.h HSEFunctional.h
XCOperator.o: XCFunctional.h RSHFunctional.h
XCPotential.o: XCPotential.h ChargeDensity.h Timer.h Context.h blacs.h
XCPotential.o: D3vector.h LDAFunctional.h XCFunctional.h VWNFunctional.h
XCPotential.o: PBEFunctional.h BLYPFunctional.h HSEFunctional.h
XCPotential.o: RSHFunctional.h B3LYPFunctional.h BHandHLYPFunctional.h
XCPotential.o: SCANFunctional.h Sample.h AtomSet.h Atom.h UnitCell.h
XCPotential.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
XCPotential.o: Control.h SlaterDet.h Basis.h Matrix.h FourierTransform.h
XCPotential.o: BasisMapping.h
XCPotential.o: ChargeDensity.h Timer.h Context.h blacs.h D3vector.h
XMLGFPreprocessor.o: Timer.h Context.h blacs.h Base64Transcoder.h Matrix.h
XMLGFPreprocessor.o: XMLGFPreprocessor.h
XMLGFPreprocessor.o: Matrix.h Context.h blacs.h
Xc.o: Sample.h AtomSet.h Context.h blacs.h Atom.h D3vector.h UnitCell.h
Xc.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
Xc.o: Control.h
isodate.o: isodate.h
jacobi.o: blacs.h Context.h Matrix.h blas.h
jade.o: blacs.h Context.h Matrix.h blas.h Timer.h
qb.o: isodate.h release.h qbox_xmlns.h uuid_str.h Context.h blacs.h
qb.o: UserInterface.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
qb.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h Wavefunction.h
qb.o: Control.h Timer.h AngleCmd.h AtomCmd.h ComputeMLWFCmd.h MLWFTransform.h
qb.o: BasisMapping.h ConstraintCmd.h DistanceCmd.h ExtForceCmd.h
qb.o: FoldInWsCmd.h HelpCmd.h KpointCmd.h ListAtomsCmd.h ListSpeciesCmd.h
qb.o: LoadCmd.h MoveCmd.h PartialChargeCmd.h PlotCmd.h PrintCmd.h QuitCmd.h
qb.o: RandomizeRCmd.h RandomizeVCmd.h RandomizeWfCmd.h ResetRotationCmd.h
qb.o: ResetVcmCmd.h RescaleVCmd.h ResponseCmd.h RseedCmd.h RunCmd.h SaveCmd.h
qb.o: SetCmd.h SetVelocityCmd.h SpeciesCmd.h SpectrumCmd.h StatusCmd.h
qb.o: ChargeDensity.h FourierTransform.h StrainCmd.h TorsionCmd.h
qb.o: BisectionCmd.h Bisection.h SlaterDet.h Basis.h Matrix.h AlphaPBE0.h
qb.o: AlphaRSH.h AtomsDyn.h BetaRSH.h BlHF.h BtHF.h Cell.h CellDyn.h
qb.o: CellLock.h CellMass.h ChargeMixCoeff.h ChargeMixNdim.h ChargeMixRcut.h
qb.o: Debug.h Dspin.h Ecut.h Ecutprec.h Ecuts.h Efield.h ForceTol.h
qb.o: Polarization.h Emass.h ExtStress.h FermiTemp.h IterCmd.h
qb.o: IterCmdPeriod.h Dt.h MuRSH.h Nempty.h NetCharge.h Nrowmax.h Nspin.h
qb.o: Occ.h RefCell.h ScfTol.h Stress.h StressTol.h Thermostat.h ThTemp.h
qb.o: ThTime.h ThWidth.h Vext.h ExternalPotential.h WfDiag.h WfDyn.h Xc.h
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
testEnergyFunctional.o: Context.h blacs.h Sample.h AtomSet.h Atom.h
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
testMatrix.o: Timer.h Context.h blacs.h Matrix.h
testSample.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h UnitCell.h
testSample.o: Matrix.h Timer.h Sample.h AtomSet.h Atom.h D3tensor.h blas.h
testSample.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
testSampleReader.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h
testSampleReader.o: UnitCell.h Matrix.h Timer.h Sample.h AtomSet.h Atom.h
testSampleReader.o: D3tensor.h blas.h ConstraintSet.h ExtForceSet.h
testSampleReader.o: Wavefunction.h Control.h SampleReader.h
testSlaterDet.o: Context.h blacs.h SlaterDet.h Basis.h D3vector.h UnitCell.h
testSlaterDet.o: Matrix.h Timer.h FourierTransform.h BasisMapping.h
testSpecies.o: Species.h SpeciesReader.h
testUnitCell.o: UnitCell.h D3vector.h
testVWN.o: VWNFunctional.h XCFunctional.h LDAFunctional.h
testWavefunction.o: Context.h blacs.h Wavefunction.h D3vector.h UnitCell.h
testWavefunction.o: SlaterDet.h Basis.h Matrix.h Timer.h
testXCFunctional.o: LDAFunctional.h XCFunctional.h PBEFunctional.h Timer.h
testXMLGFPreprocessor.o: Context.h blacs.h Matrix.h XMLGFPreprocessor.h
test_fftw.o: Timer.h readTSC.h
test_vext.o: Function3d.h D3vector.h
testjacobi.o: Timer.h Context.h blacs.h Matrix.h jacobi.h
testjade.o: Timer.h Context.h blacs.h Matrix.h jade.h
uuid_str.o: uuid_str.h
