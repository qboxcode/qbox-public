#-------------------------------------------------------------------------------
#
# Copyright (c) 2008 The Regents of the University of California
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
# $Id: Makefile,v 1.70 2010-05-10 20:51:00 fgygi Exp $
#------------------------------------------------------------------------------
#
include $(TARGET).mk
#------------------------------------------------------------------------------
#
EXEC=qb
TESTEXECS=testMatrix testFTGrid testGridFunction testBasis \
testBlacsContext testSlaterDet testEnergyFunctional testSample \
testChargeDensity  testFourierTransform testSpecies testContext

OBJECTS=qb.o AtomSet.o Atom.o Species.o \
        Wavefunction.o SlaterDet.o \
        EnergyFunctional.o SampleStepper.o \
        Basis.o FourierTransform.o Matrix.o Context.o \
        sinft.o spline.o UnitCell.o StructureFactor.o \
        ChargeDensity.o UserInterface.o RunCmd.o \
        LoadCmd.o SaveCmd.o \
        SpeciesCmd.o SpeciesReader.o SpeciesHandler.o \
        XCOperator.o ExchangeOperator.o Bisection.o KPGridConnectivity.o \
	XCPotential.o LDAFunctional.o \
        PBEFunctional.o BLYPFunctional.o B3LYPFunctional.o  \
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
        SampleWriter.o ComputeMLWFCmd.o BasisMapping.o MLWFTransform.o \
        jade.o PlotCmd.o \
        ExtForceSet.o ExtForce.o AtomicExtForce.o PairExtForce.o \
        GlobalExtForce.o \
        uuid_str.o sampling.o CGOptimizer.o LineMinimizer.o \
        $(PLTOBJECTS)
CXXFLAGS += -DTARGET='"$(TARGET)"'
 $(EXEC):    $(OBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 lib:   $(OBJECTS)
	ar cr libqb.a $^
 kpgen:	kpgen.C
	$(CXX) -o $@ $^
 testLineMinimizer: testLineMinimizer.o LineMinimizer.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testCGOptimizer: testCGOptimizer.o CGOptimizer.o LineMinimizer.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testAndersonMixer: testAndersonMixer.o AndersonMixer.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSample: testSample.o AtomSet.o Atom.o Species.o \
        Wavefunction.o SlaterDet.o \
        Basis.o FourierTransform.o Matrix.o Context.o \
        sinft.o spline.o UnitCell.o \
        Base64Transcoder.o Constraint.o ConstraintSet.o DistanceConstraint.o \
        AngleConstraint.o TorsionConstraint.o PositionConstraint.o \
        ExtForceSet.o ExtForce.o PairExtForce.o AtomicExtForce.o \
        GlobalExtForce.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testreadSample: testreadSample.o \
         Context.o AtomSet.o Atom.o Basis.o SlaterDet.o \
         Constraint.o ConstraintSet.o Wavefunction.o UnitCell.o Matrix.o \
         SampleReader.o StructuredDocumentHandler.o FourierTransform.o \
         Species.o SpeciesReader.o sinft.o spline.o SpeciesHandler.o \
         SampleHandler.o WavefunctionHandler.o AtomSetHandler.o \
         Base64Transcoder.o ConstraintSet.o DistanceConstraint.o \
         AngleConstraint.o TorsionConstraint.o XMLGFPreprocessor.o \
         PositionConstraint.o \
        ExtForceSet.o ExtForce.o PairExtForce.o AtomicExtForce.o \
        GlobalExtForce.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 truncateMLWF: truncateMLWF.o \
         MLWFTransform.o BasisMapping.o jade.o \
         Context.o AtomSet.o Atom.o Basis.o SlaterDet.o \
         Constraint.o ConstraintSet.o Wavefunction.o UnitCell.o Matrix.o \
         SampleReader.o StructuredDocumentHandler.o FourierTransform.o \
         Species.o SpeciesReader.o sinft.o spline.o SpeciesHandler.o \
         SampleHandler.o WavefunctionHandler.o AtomSetHandler.o \
         Base64Transcoder.o ConstraintSet.o DistanceConstraint.o \
         AngleConstraint.o TorsionConstraint.o XMLGFPreprocessor.o \
         SampleWriter.o qbox_xmlns.o isodate.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 bisection: bisection.o jade.o \
         Context.o AtomSet.o Atom.o Basis.o SlaterDet.o \
         Constraint.o ConstraintSet.o Wavefunction.o UnitCell.o Matrix.o \
         SampleReader.o StructuredDocumentHandler.o FourierTransform.o \
         Species.o SpeciesReader.o sinft.o spline.o SpeciesHandler.o \
         SampleHandler.o WavefunctionHandler.o AtomSetHandler.o \
         Base64Transcoder.o ConstraintSet.o DistanceConstraint.o \
         AngleConstraint.o TorsionConstraint.o PositionConstraint.o \
         XMLGFPreprocessor.o \
         SampleWriter.o qbox_xmlns.o isodate.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 plotSample: plotSample.o \
         Context.o AtomSet.o Atom.o Basis.o SlaterDet.o \
         Constraint.o ConstraintSet.o Wavefunction.o UnitCell.o Matrix.o \
         SampleReader.o StructuredDocumentHandler.o FourierTransform.o \
         Species.o SpeciesReader.o sinft.o spline.o SpeciesHandler.o \
         SampleHandler.o WavefunctionHandler.o AtomSetHandler.o \
         Base64Transcoder.o ConstraintSet.o DistanceConstraint.o \
         AngleConstraint.o TorsionConstraint.o XMLGFPreprocessor.o \
         ChargeDensity.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testChargeDensity: testChargeDensity.o ChargeDensity.o \
        Wavefunction.o SlaterDet.o \
        Basis.o FourierTransform.o Matrix.o UnitCell.o Context.o \
        Base64Transcoder.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testWavefunction: testWavefunction.o Wavefunction.o SlaterDet.o \
        Basis.o FourierTransform.o Matrix.o UnitCell.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testEnergyFunctional: testEnergyFunctional.o EnergyFunctional.o Basis.o \
	SlaterDet.o Matrix.o UnitCell.o Context.o FourierTransform.o \
        Wavefunction.o Species.o Atom.o AtomSet.o StructureFactor.o \
        ChargeDensity.o \
        sinft.o spline.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSlaterDet: testSlaterDet.o SlaterDet.o FourierTransform.o \
        Basis.o UnitCell.o Matrix.o Context.o Base64Transcoder.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testBasisMapping: testBasisMapping.o BasisMapping.o Basis.o \
        Context.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSpecies: testSpecies.o Species.o sinft.o spline.o \
	SpeciesReader.o StructuredDocumentHandler.o SpeciesHandler.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testXCFunctional: testXCFunctional.o LDAFunctional.o PBEFunctional.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
#testMatrix: testMatrix.o Matrix.o Context.o pdgetri.o
 testMatrix: testMatrix.o Matrix.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testjacobi: testjacobi.o jacobi.o Matrix.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testjade: testjade.o jade.o Matrix.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testFourierTransform: testFourierTransform.o FourierTransform.o \
        Basis.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testFTGrid: testFTGrid.o FTGrid.o BlacsContext.o Basis.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testGridFunction: testGridFunction.o GridFunction.o BlacsContext.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testUnitCell: testUnitCell.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testBasis: testBasis.o Basis.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testMapping: testMapping.o Basis.o UnitCell.o Context.o ShiftMapping.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testBlacsContext: testBlacsContext.o BlacsContext.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testContext: testContext.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testHandler: testHandler.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testPParse: testPParse.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testMemParse: testMemParse.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 upf2qso: upf2qso.o PeriodicTable.o spline.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 xmlSpecies: xmlSpecies.o qbox_xmlns.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 xmlget: xmlget.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 xmlextract: xmlextract.o
	$(LD) -o xmlextract xmlextract.C -I$(XERCESCDIR)/include \
        -L$(XERCESCDIR)/lib -lxerces-c
 testXMLGFPreprocessor: testXMLGFPreprocessor.o XMLGFPreprocessor.o Context.o \
        Base64Transcoder.o Matrix.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testwrite: testwrite.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testStepper: testStepper.o AtomSet.o Context.o MDIonicStepper.o \
Wavefunction.o Species.o SlaterDet.o spline.o sinft.o Basis.o Matrix.o \
FourierTransform.o UnitCell.o Atom.o BOSampleStepper.o \
SampleStepper.o EnergyFunctional.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 test_fftw: test_fftw.o $(PLTOBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 test_fftw3: test_fftw3.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 test_sinft: test_sinft.o sinft.o $(PLTOBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 test_spline: test_spline.o spline.o $(PLTOBJECTS)
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

AlphaPBE0.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
AlphaPBE0.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
AndersonMixer.o: AndersonMixer.h blas.h
AngleCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
AngleCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
AngleConstraint.o: AngleConstraint.h Constraint.h D3vector.h AtomSet.h
AngleConstraint.o: Context.h Atom.h UnitCell.h Species.h
AngleConstraint.o: Constraint.h D3vector.h
Atom.o: Atom.h D3vector.h
AtomCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
AtomCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Atom.o: D3vector.h
AtomicExtForce.o: AtomicExtForce.h ExtForce.h D3vector.h AtomSet.h Context.h
AtomicExtForce.o: Atom.h UnitCell.h Species.h
AtomicExtForce.o: ExtForce.h D3vector.h
AtomsDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
AtomsDyn.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
AtomSet.o: AtomSet.h Context.h Atom.h D3vector.h UnitCell.h Species.h
AtomSet.o: NameOf.h sampling.h
AtomSet.o: Context.h Atom.h D3vector.h UnitCell.h
AtomSetHandler.o: AtomSetHandler.h StructureHandler.h D3vector.h AtomSet.h
AtomSetHandler.o: Context.h Atom.h UnitCell.h Species.h SpeciesHandler.h
AtomSetHandler.o: SpeciesReader.h StrX.h SampleReader.h
AtomSetHandler.o: StructureHandler.h D3vector.h
B3LYPFunctional.o: B3LYPFunctional.h XCFunctional.h
B3LYPFunctional.o: XCFunctional.h
Base64Transcoder.o: Base64Transcoder.h
Basis.o: Basis.h D3vector.h UnitCell.h
Basis.o: D3vector.h UnitCell.h
BasisMapping.o: Basis.h D3vector.h UnitCell.h Context.h BasisMapping.h
Bisection.o: Bisection.h Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h
Bisection.o: Matrix.h Timer.h jade.h FourierTransform.h
BisectionCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
BisectionCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
BisectionCmd.o: Wavefunction.h Control.h Bisection.h SlaterDet.h Basis.h
BisectionCmd.o: Matrix.h Timer.h
Bisection.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
Bisection.o: Timer.h
BlHF.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
BlHF.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
BLYPFunctional.o: BLYPFunctional.h XCFunctional.h
BLYPFunctional.o: XCFunctional.h
BMDIonicStepper.o: BMDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
BMDIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
BMDIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
BMDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
BMDIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
BMDIonicStepper.o: Wavefunction.h Control.h Species.h
BOSampleStepper.o: BOSampleStepper.h SampleStepper.h Timer.h
BOSampleStepper.o: EnergyFunctional.h StructureFactor.h Sample.h AtomSet.h
BOSampleStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
BOSampleStepper.o: ExtForceSet.h Wavefunction.h Control.h ChargeDensity.h
BOSampleStepper.o: SlaterDet.h Basis.h Matrix.h WavefunctionStepper.h
BOSampleStepper.o: SDWavefunctionStepper.h PSDWavefunctionStepper.h
BOSampleStepper.o: PSDAWavefunctionStepper.h JDWavefunctionStepper.h
BOSampleStepper.o: Preconditioner.h SDIonicStepper.h IonicStepper.h Species.h
BOSampleStepper.o: SDAIonicStepper.h LineMinimizer.h CGIonicStepper.h
BOSampleStepper.o: CGOptimizer.h MDIonicStepper.h BMDIonicStepper.h
BOSampleStepper.o: SDCellStepper.h CellStepper.h CGCellStepper.h
BOSampleStepper.o: AndersonMixer.h MLWFTransform.h BasisMapping.h
BOSampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h
BOSampleStepper.o: StructureFactor.h Sample.h AtomSet.h Context.h Atom.h
BOSampleStepper.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
BOSampleStepper.o: Wavefunction.h Control.h ChargeDensity.h
BtHF.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
BtHF.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
CellDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
CellDyn.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h SlaterDet.h
CellDyn.o: Basis.h Matrix.h Timer.h
Cell.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Cell.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
CellLock.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
CellLock.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
CellMass.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
CellMass.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
CellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
CellStepper.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
CellStepper.o: Control.h
CellStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
CellStepper.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
CGCellStepper.o: CGCellStepper.h CellStepper.h Sample.h AtomSet.h Context.h
CGCellStepper.o: Atom.h D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
CGCellStepper.o: Wavefunction.h Control.h CGOptimizer.h LineMinimizer.h
CGCellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
CGCellStepper.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
CGCellStepper.o: Control.h CGOptimizer.h LineMinimizer.h
CGIonicStepper.o: CGIonicStepper.h IonicStepper.h Sample.h AtomSet.h
CGIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
CGIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
CGIonicStepper.o: CGOptimizer.h LineMinimizer.h
CGIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
CGIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
CGIonicStepper.o: Wavefunction.h Control.h Species.h CGOptimizer.h
CGIonicStepper.o: LineMinimizer.h
CGOptimizer.o: CGOptimizer.h LineMinimizer.h blas.h
CGOptimizer.o: LineMinimizer.h
ChargeDensity.o: ChargeDensity.h Timer.h Context.h Basis.h D3vector.h
ChargeDensity.o: UnitCell.h Wavefunction.h FourierTransform.h SlaterDet.h
ChargeDensity.o: Matrix.h
ChargeDensity.o: Timer.h Context.h
ChargeMixCoeff.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ChargeMixCoeff.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ChargeMixNdim.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ChargeMixNdim.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ChargeMixRcut.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ChargeMixRcut.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ComputeMLWFCmd.o: ComputeMLWFCmd.h UserInterface.h Sample.h AtomSet.h
ComputeMLWFCmd.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
ComputeMLWFCmd.o: ExtForceSet.h Wavefunction.h Control.h MLWFTransform.h
ComputeMLWFCmd.o: BasisMapping.h SlaterDet.h Basis.h Matrix.h Timer.h
ComputeMLWFCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ComputeMLWFCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
ComputeMLWFCmd.o: Wavefunction.h Control.h MLWFTransform.h BasisMapping.h
ConfinementPotential.o: ConfinementPotential.h Basis.h D3vector.h UnitCell.h
Constraint.o: Constraint.h
ConstraintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ConstraintCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
ConstraintCmd.o: Wavefunction.h Control.h
ConstraintSet.o: ConstraintSet.h PositionConstraint.h Constraint.h
ConstraintSet.o: DistanceConstraint.h AngleConstraint.h D3vector.h
ConstraintSet.o: TorsionConstraint.h Atom.h AtomSet.h Context.h UnitCell.h
Context.o: Context.h blacs.h
CPSampleStepper.o: CPSampleStepper.h SampleStepper.h Timer.h
CPSampleStepper.o: EnergyFunctional.h StructureFactor.h ChargeDensity.h
CPSampleStepper.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
CPSampleStepper.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
CPSampleStepper.o: SlaterDet.h Basis.h Matrix.h MDWavefunctionStepper.h
CPSampleStepper.o: WavefunctionStepper.h MDIonicStepper.h IonicStepper.h
CPSampleStepper.o: Species.h SDCellStepper.h CellStepper.h CGCellStepper.h
CPSampleStepper.o: CGOptimizer.h LineMinimizer.h
CPSampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h
CPSampleStepper.o: StructureFactor.h ChargeDensity.h Context.h Sample.h
CPSampleStepper.o: AtomSet.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
CPSampleStepper.o: ExtForceSet.h Wavefunction.h Control.h
Debug.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Debug.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
DistanceCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
DistanceCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
DistanceCmd.o: Control.h
DistanceConstraint.o: DistanceConstraint.h Constraint.h AtomSet.h Context.h
DistanceConstraint.o: Atom.h D3vector.h UnitCell.h Species.h
DistanceConstraint.o: Constraint.h
Dspin.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Dspin.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Dt.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Dt.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Ecut.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Ecut.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Ecutprec.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Ecutprec.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Ecuts.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Ecuts.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Emass.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Emass.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
EnergyFunctional.o: EnergyFunctional.h StructureFactor.h Timer.h Sample.h
EnergyFunctional.o: AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
EnergyFunctional.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
EnergyFunctional.o: Species.h ChargeDensity.h SlaterDet.h Basis.h Matrix.h
EnergyFunctional.o: FourierTransform.h XCOperator.h NonLocalPotential.h
EnergyFunctional.o: ConfinementPotential.h blas.h
EnergyFunctional.o: StructureFactor.h Timer.h
ExchangeOperator.o: VectorLess.h ExchangeOperator.h Sample.h AtomSet.h
ExchangeOperator.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
ExchangeOperator.o: ExtForceSet.h Wavefunction.h Control.h SlaterDet.h
ExchangeOperator.o: Basis.h Matrix.h Timer.h FourierTransform.h
ExchangeOperator.o: KPGridConnectivity.h Bisection.h
ExchangeOperator.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ExchangeOperator.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ExchangeOperator.o: SlaterDet.h Basis.h Matrix.h Timer.h FourierTransform.h
ExchangeOperator.o: KPGridConnectivity.h
ExtForce.o: ExtForce.h D3vector.h
ExtForceCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
ExtForceCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ExtForceCmd.o: Control.h
ExtForce.o: D3vector.h
ExtForceSet.o: ExtForceSet.h AtomicExtForce.h ExtForce.h D3vector.h
ExtForceSet.o: PairExtForce.h GlobalExtForce.h Atom.h AtomSet.h Context.h
ExtForceSet.o: UnitCell.h
ExtStress.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ExtStress.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
FermiTemp.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
FermiTemp.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
FoldInWsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
FoldInWsCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
FoldInWsCmd.o: Control.h
FourierTransform.o: FourierTransform.h Timer.h Basis.h D3vector.h UnitCell.h
FourierTransform.o: blas.h
FourierTransform.o: Timer.h
GlobalExtForce.o: GlobalExtForce.h ExtForce.h D3vector.h AtomSet.h Context.h
GlobalExtForce.o: Atom.h UnitCell.h Species.h
GlobalExtForce.o: ExtForce.h D3vector.h
HelpCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
HelpCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
IonicStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
IonicStepper.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
IonicStepper.o: Species.h
isodate.o: isodate.h
jacobi.o: blacs.h Context.h Matrix.h blas.h
jade.o: blacs.h Context.h Matrix.h blas.h Timer.h
JDWavefunctionStepper.o: JDWavefunctionStepper.h WavefunctionStepper.h
JDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
JDWavefunctionStepper.o: SlaterDet.h Context.h Basis.h Matrix.h
JDWavefunctionStepper.o: EnergyFunctional.h StructureFactor.h
JDWavefunctionStepper.o: Preconditioner.h
JDWavefunctionStepper.o: WavefunctionStepper.h Timer.h Wavefunction.h
JDWavefunctionStepper.o: D3vector.h UnitCell.h
kpgen.o: D3vector.h
KPGridConnectivity.o: KPGridConnectivity.h Sample.h AtomSet.h Context.h
KPGridConnectivity.o: Atom.h D3vector.h UnitCell.h ConstraintSet.h
KPGridConnectivity.o: ExtForceSet.h Wavefunction.h Control.h SlaterDet.h
KPGridConnectivity.o: Basis.h Matrix.h Timer.h
KPGridConnectivity.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
KPGridConnectivity.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
KPGridConnectivity.o: Control.h SlaterDet.h Basis.h Matrix.h Timer.h
KpointCmd.o: UserInterface.h D3vector.h Sample.h AtomSet.h Context.h Atom.h
KpointCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
KpointCmd.o: Control.h
LDAFunctional.o: LDAFunctional.h XCFunctional.h
LDAFunctional.o: XCFunctional.h
LineMinimizer.o: LineMinimizer.h
ListAtomsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListAtomsCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
ListAtomsCmd.o: Wavefunction.h Control.h
ListSpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListSpeciesCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
ListSpeciesCmd.o: Wavefunction.h Control.h
LoadCmd.o: LoadCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
LoadCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
LoadCmd.o: Control.h SampleReader.h
LoadCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
LoadCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Matrix.o: Context.h blacs.h Matrix.h
MDIonicStepper.o: MDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
MDIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
MDIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h sampling.h
MDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
MDIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
MDIonicStepper.o: Wavefunction.h Control.h Species.h
MDWavefunctionStepper.o: MDWavefunctionStepper.h WavefunctionStepper.h
MDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
MDWavefunctionStepper.o: SlaterDet.h Context.h Basis.h Matrix.h Sample.h
MDWavefunctionStepper.o: AtomSet.h Atom.h ConstraintSet.h ExtForceSet.h
MDWavefunctionStepper.o: Control.h
MDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
MLWFTransform.o: MLWFTransform.h D3vector.h BasisMapping.h Basis.h UnitCell.h
MLWFTransform.o: SlaterDet.h Context.h Matrix.h Timer.h jade.h blas.h
MLWFTransform.o: D3vector.h BasisMapping.h
MoveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
MoveCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Nempty.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Nempty.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
NetCharge.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
NetCharge.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
NonLocalPotential.o: NonLocalPotential.h AtomSet.h Context.h Atom.h
NonLocalPotential.o: D3vector.h UnitCell.h Basis.h SlaterDet.h Matrix.h
NonLocalPotential.o: Timer.h Species.h blas.h
NonLocalPotential.o: AtomSet.h Context.h Atom.h D3vector.h UnitCell.h Basis.h
NonLocalPotential.o: SlaterDet.h Matrix.h Timer.h
Nrowmax.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Nrowmax.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
Nspin.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Nspin.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
PairExtForce.o: PairExtForce.h ExtForce.h D3vector.h AtomSet.h Context.h
PairExtForce.o: Atom.h UnitCell.h Species.h
PairExtForce.o: ExtForce.h D3vector.h
PBEFunctional.o: PBEFunctional.h XCFunctional.h
PBEFunctional.o: XCFunctional.h
PlotCmd.o: PlotCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
PlotCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
PlotCmd.o: Control.h isodate.h release.h SampleReader.h Basis.h
PlotCmd.o: FourierTransform.h Timer.h SlaterDet.h Matrix.h Species.h
PlotCmd.o: ChargeDensity.h
PlotCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
PlotCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
PositionConstraint.o: PositionConstraint.h Constraint.h AtomSet.h Context.h
PositionConstraint.o: Atom.h D3vector.h UnitCell.h Species.h
PositionConstraint.o: Constraint.h
Preconditioner.o: Preconditioner.h Basis.h D3vector.h UnitCell.h
Preconditioner.o: Wavefunction.h EnergyFunctional.h StructureFactor.h Timer.h
Preconditioner.o: ConfinementPotential.h SlaterDet.h Context.h Matrix.h
PrintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
PrintCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
PSDAWavefunctionStepper.o: PSDAWavefunctionStepper.h WavefunctionStepper.h
PSDAWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
PSDAWavefunctionStepper.o: SlaterDet.h Context.h Basis.h Matrix.h
PSDAWavefunctionStepper.o: Preconditioner.h
PSDAWavefunctionStepper.o: WavefunctionStepper.h Timer.h Wavefunction.h
PSDAWavefunctionStepper.o: D3vector.h UnitCell.h
PSDWavefunctionStepper.o: PSDWavefunctionStepper.h WavefunctionStepper.h
PSDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
PSDWavefunctionStepper.o: SlaterDet.h Context.h Basis.h Matrix.h
PSDWavefunctionStepper.o: Preconditioner.h
PSDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
qb.o: isodate.h release.h qbox_xmlns.h Context.h UserInterface.h Sample.h
qb.o: AtomSet.h Atom.h D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
qb.o: Wavefunction.h Control.h Timer.h AngleCmd.h AtomCmd.h ComputeMLWFCmd.h
qb.o: MLWFTransform.h BasisMapping.h ConstraintCmd.h DistanceCmd.h
qb.o: ExtForceCmd.h FoldInWsCmd.h HelpCmd.h KpointCmd.h ListAtomsCmd.h
qb.o: ListSpeciesCmd.h LoadCmd.h MoveCmd.h PlotCmd.h PrintCmd.h QuitCmd.h
qb.o: RandomizeRCmd.h RandomizeVCmd.h RandomizeWfCmd.h ResetVcmCmd.h
qb.o: RescaleVCmd.h RseedCmd.h RunCmd.h SaveCmd.h SetCmd.h SpeciesCmd.h
qb.o: StatusCmd.h StrainCmd.h TorsionCmd.h BisectionCmd.h Bisection.h
qb.o: SlaterDet.h Basis.h Matrix.h AlphaPBE0.h AtomsDyn.h BlHF.h BtHF.h
qb.o: Cell.h CellDyn.h CellLock.h CellMass.h ChargeMixCoeff.h ChargeMixNdim.h
qb.o: ChargeMixRcut.h Debug.h Dspin.h Ecut.h Ecutprec.h Ecuts.h Emass.h
qb.o: ExtStress.h FermiTemp.h Dt.h Nempty.h NetCharge.h Nrowmax.h Nspin.h
qb.o: RefCell.h ScfTol.h Stress.h Thermostat.h ThTemp.h ThTime.h ThWidth.h
qb.o: WfDiag.h WfDyn.h Xc.h
qbox_xmlns.o: qbox_xmlns.h
QuitCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
QuitCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
RandomizeRCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
RandomizeRCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
RandomizeRCmd.o: Wavefunction.h Control.h
RandomizeVCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
RandomizeVCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
RandomizeVCmd.o: Wavefunction.h Control.h
RandomizeWfCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
RandomizeWfCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
RandomizeWfCmd.o: Wavefunction.h Control.h
RefCell.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
RefCell.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
release.o: release.h
RescaleVCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
RescaleVCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
RescaleVCmd.o: Control.h
ResetVcmCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
ResetVcmCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
ResetVcmCmd.o: Control.h
RseedCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
RseedCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
RunCmd.o: RunCmd.h UserInterface.h BOSampleStepper.h SampleStepper.h Timer.h
RunCmd.o: EnergyFunctional.h StructureFactor.h Sample.h AtomSet.h Context.h
RunCmd.o: Atom.h D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
RunCmd.o: Wavefunction.h Control.h ChargeDensity.h CPSampleStepper.h
RunCmd.o: UserInterface.h
Sample.o: AtomSet.h Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
Sample.o: ExtForceSet.h Wavefunction.h Control.h
SampleHandler.o: SampleHandler.h StructureHandler.h Sample.h AtomSet.h
SampleHandler.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
SampleHandler.o: ExtForceSet.h Wavefunction.h Control.h AtomSetHandler.h
SampleHandler.o: WavefunctionHandler.h SlaterDet.h Basis.h Matrix.h Timer.h
SampleHandler.o: StrX.h
SampleHandler.o: StructureHandler.h
SampleReader.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
SampleReader.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
SampleReader.o: SampleReader.h XMLGFPreprocessor.h Matrix.h Timer.h
SampleReader.o: SampleHandler.h StructureHandler.h
SampleReader.o: StructuredDocumentHandler.h StrX.h
SampleStepper.o: SampleStepper.h Timer.h Sample.h AtomSet.h Context.h Atom.h
SampleStepper.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
SampleStepper.o: Wavefunction.h Control.h Species.h
SampleStepper.o: Timer.h
SampleWriter.o: SampleWriter.h Context.h Sample.h AtomSet.h Atom.h D3vector.h
SampleWriter.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
SampleWriter.o: Control.h qbox_xmlns.h Timer.h SharedFilePtr.h
SampleWriter.o: Context.h
SaveCmd.o: SaveCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
SaveCmd.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
SaveCmd.o: Control.h SampleWriter.h isodate.h release.h
SaveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SaveCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ScfTol.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ScfTol.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
SDAIonicStepper.o: SDAIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDAIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
SDAIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
SDAIonicStepper.o: LineMinimizer.h
SDAIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
SDAIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
SDAIonicStepper.o: Wavefunction.h Control.h Species.h LineMinimizer.h
SDCellStepper.o: SDCellStepper.h CellStepper.h Sample.h AtomSet.h Context.h
SDCellStepper.o: Atom.h D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
SDCellStepper.o: Wavefunction.h Control.h
SDCellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SDCellStepper.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
SDCellStepper.o: Control.h
SDIonicStepper.o: SDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
SDIonicStepper.o: ExtForceSet.h Wavefunction.h Control.h Species.h
SDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
SDIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
SDIonicStepper.o: Wavefunction.h Control.h Species.h
SDWavefunctionStepper.o: SDWavefunctionStepper.h WavefunctionStepper.h
SDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
SDWavefunctionStepper.o: SlaterDet.h Context.h Basis.h Matrix.h Sample.h
SDWavefunctionStepper.o: AtomSet.h Atom.h ConstraintSet.h ExtForceSet.h
SDWavefunctionStepper.o: Control.h
SDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
SetCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SetCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
sinft.o: sinft.h
SlaterDet.o: SlaterDet.h Context.h Basis.h D3vector.h UnitCell.h Matrix.h
SlaterDet.o: Timer.h FourierTransform.h blas.h Base64Transcoder.h
SlaterDet.o: SharedFilePtr.h
SlaterDet.o: Context.h Basis.h D3vector.h UnitCell.h Matrix.h Timer.h
Species.o: Species.h spline.h sinft.h
SpeciesCmd.o: SpeciesCmd.h UserInterface.h Sample.h AtomSet.h Context.h
SpeciesCmd.o: Atom.h D3vector.h UnitCell.h ConstraintSet.h ExtForceSet.h
SpeciesCmd.o: Wavefunction.h Control.h SpeciesReader.h Species.h
SpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SpeciesCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
SpeciesCmd.o: Control.h
SpeciesHandler.o: SpeciesHandler.h StructureHandler.h Species.h StrX.h
SpeciesHandler.o: StructureHandler.h
SpeciesReader.o: Species.h SpeciesReader.h StructuredDocumentHandler.h StrX.h
SpeciesReader.o: StructureHandler.h SpeciesHandler.h
spline.o: spline.h
StatusCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
StatusCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
StatusCmd.o: Control.h
StrainCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
StrainCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
StrainCmd.o: Control.h
Stress.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Stress.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
StructuredDocumentHandler.o: StructuredDocumentHandler.h StrX.h
StructuredDocumentHandler.o: StructureHandler.h
StructuredDocumentHandler.o: StrX.h StructureHandler.h
StructureFactor.o: StructureFactor.h Basis.h D3vector.h UnitCell.h
testAndersonMixer.o: Context.h AndersonMixer.h
testBase64Transcoder.o: Base64Transcoder.h
testBasis.o: Basis.h D3vector.h UnitCell.h
testBasisMapping.o: Context.h Basis.h D3vector.h UnitCell.h BasisMapping.h
testBasisMapping.o: Timer.h
testCGOptimizer.o: CGOptimizer.h LineMinimizer.h
testChargeDensity.o: Context.h Wavefunction.h D3vector.h UnitCell.h
testChargeDensity.o: ChargeDensity.h Timer.h SlaterDet.h Basis.h Matrix.h
testChargeDensity.o: FourierTransform.h
testContext.o: Context.h
testEnergyFunctional.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h
testEnergyFunctional.o: UnitCell.h ConstraintSet.h ExtForceSet.h
testEnergyFunctional.o: Wavefunction.h Control.h EnergyFunctional.h
testEnergyFunctional.o: StructureFactor.h Timer.h
test_fftw.o: Timer.h readTSC.h
testFourierTransform.o: Basis.h D3vector.h UnitCell.h FourierTransform.h
testFourierTransform.o: Timer.h
testjacobi.o: Timer.h Context.h Matrix.h jacobi.h
testjade.o: Timer.h Context.h Matrix.h jade.h
testLDAFunctional.o: LDAFunctional.h XCFunctional.h
testLineMinimizer.o: LineMinimizer.h
testMatrix.o: Timer.h Context.h Matrix.h
testSample.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
testSample.o: Timer.h Sample.h AtomSet.h Atom.h ConstraintSet.h ExtForceSet.h
testSample.o: Wavefunction.h Control.h
testSlaterDet.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
testSlaterDet.o: Timer.h FourierTransform.h
testSpecies.o: Species.h SpeciesReader.h
testUnitCell.o: UnitCell.h D3vector.h
testWavefunction.o: Context.h Wavefunction.h D3vector.h UnitCell.h
testWavefunction.o: SlaterDet.h Basis.h Matrix.h Timer.h
testXCFunctional.o: LDAFunctional.h XCFunctional.h PBEFunctional.h Timer.h
testXMLGFPreprocessor.o: Context.h Matrix.h XMLGFPreprocessor.h
Thermostat.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Thermostat.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ThTemp.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ThTemp.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ThTime.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ThTime.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
ThWidth.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ThWidth.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
TorsionCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
TorsionCmd.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
TorsionCmd.o: Control.h
TorsionConstraint.o: TorsionConstraint.h Constraint.h D3vector.h AtomSet.h
TorsionConstraint.o: Context.h Atom.h UnitCell.h Species.h
TorsionConstraint.o: Constraint.h D3vector.h
UnitCell.o: UnitCell.h D3vector.h
UnitCell.o: D3vector.h
UserInterface.o: UserInterface.h qbox_xmlns.h
Wavefunction.o: Wavefunction.h D3vector.h UnitCell.h SlaterDet.h Context.h
Wavefunction.o: Basis.h Matrix.h Timer.h jacobi.h SharedFilePtr.h
Wavefunction.o: D3vector.h UnitCell.h
WavefunctionHandler.o: WavefunctionHandler.h StructureHandler.h UnitCell.h
WavefunctionHandler.o: D3vector.h Wavefunction.h SlaterDet.h Context.h
WavefunctionHandler.o: Basis.h Matrix.h Timer.h FourierTransform.h
WavefunctionHandler.o: SampleReader.h StrX.h
WavefunctionHandler.o: StructureHandler.h UnitCell.h D3vector.h
WavefunctionHandler.o: Wavefunction.h SlaterDet.h Context.h Basis.h Matrix.h
WavefunctionHandler.o: Timer.h
WavefunctionStepper.o: Timer.h
WfDiag.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
WfDiag.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
WfDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
WfDyn.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h SlaterDet.h
WfDyn.o: Basis.h Matrix.h Timer.h
Xc.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Xc.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
XCOperator.o: XCOperator.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
XCOperator.o: UnitCell.h ConstraintSet.h ExtForceSet.h Wavefunction.h
XCOperator.o: Control.h ChargeDensity.h Timer.h XCPotential.h
XCOperator.o: ExchangeOperator.h SlaterDet.h Basis.h Matrix.h
XCOperator.o: FourierTransform.h KPGridConnectivity.h
XCOperator.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
XCOperator.o: ConstraintSet.h ExtForceSet.h Wavefunction.h Control.h
XCPotential.o: XCPotential.h ChargeDensity.h Timer.h Context.h
XCPotential.o: LDAFunctional.h XCFunctional.h PBEFunctional.h
XCPotential.o: BLYPFunctional.h B3LYPFunctional.h Basis.h D3vector.h
XCPotential.o: UnitCell.h FourierTransform.h blas.h
XCPotential.o: ChargeDensity.h Timer.h Context.h
XMLGFPreprocessor.o: Timer.h Context.h Base64Transcoder.h Matrix.h
XMLGFPreprocessor.o: XMLGFPreprocessor.h
XMLGFPreprocessor.o: Matrix.h
