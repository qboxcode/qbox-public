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
# $Id: Makefile,v 1.60 2009-03-08 01:16:33 fgygi Exp $
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
        ChargeDensity.o UserInterface.o RunCmd.o LDAFunctional.o \
        LoadCmd.o SaveCmd.o \
        SpeciesCmd.o SpeciesReader.o SpeciesHandler.o \
        XCPotential.o LDAFunctional.o PBEFunctional.o \
        BLYPFunctional.o NonLocalPotential.o \
	SampleReader.o StructuredDocumentHandler.o \
        SampleHandler.o AtomSetHandler.o WavefunctionHandler.o \
        XMLGFPreprocessor.o Base64Transcoder.o \
        CPSampleStepper.o BOSampleStepper.o \
        SDWavefunctionStepper.o \
        MDWavefunctionStepper.o SDIonicStepper.o MDIonicStepper.o \
        PSDWavefunctionStepper.o PSDAWavefunctionStepper.o \
        SDCellStepper.o ConfinementPotential.o Preconditioner.o \
        release.o qbox_xmlns.o isodate.o \
        AndersonMixer.o SDAIonicStepper.o CGIonicStepper.o \
        ConstraintSet.o Constraint.o DistanceConstraint.o \
        AngleConstraint.o TorsionConstraint.o jacobi.o \
        SampleWriter.o ComputeMLWFCmd.o BasisMapping.o MLWFTransform.o \
        jade.o \
        $(PLTOBJECTS)
CXXFLAGS += -DTARGET='"$(TARGET)"'
 $(EXEC):    $(OBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 lib:   $(OBJECTS)
	ar cr libqb.a $^
 SamplePrint: SamplePrint.o SamplePrintHandlers.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testMatrixColumnRotation : testMatrixColumnRotation.o Matrix.o Context.o \
        MatrixColumnRotation.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testAndersonMixer: testAndersonMixer.o AndersonMixer.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSample: testSample.o AtomSet.o Atom.o Species.o \
        Wavefunction.o SlaterDet.o \
        Basis.o FourierTransform.o Matrix.o Context.o \
        sinft.o spline.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testreadSample: testreadSample.o \
         Context.o AtomSet.o Atom.o Basis.o SlaterDet.o \
         Constraint.o ConstraintSet.o Wavefunction.o UnitCell.o Matrix.o \
         SampleReader.o StructuredDocumentHandler.o FourierTransform.o \
         Species.o SpeciesReader.o sinft.o spline.o SpeciesHandler.o \
         SampleHandler.o WavefunctionHandler.o AtomSetHandler.o \
         Base64Transcoder.o ConstraintSet.o DistanceConstraint.o \
         AngleConstraint.o TorsionConstraint.o XMLGFPreprocessor.o
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
         AngleConstraint.o TorsionConstraint.o XMLGFPreprocessor.o \
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
        Basis.o FourierTransform.o Matrix.o UnitCell.o Context.o
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
 testSpecies: testSpecies.o Species.o sinft.o spline.o Context.o \
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
        Context.o Basis.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testFTGrid: testFTGrid.o FTGrid.o BlacsContext.o Basis.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testGridFunction: testGridFunction.o GridFunction.o BlacsContext.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testUnitCell: testUnitCell.o UnitCell.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testBasis: testBasis.o Basis.o UnitCell.o Context.o
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
	makedepend -Y -D$(PLT) *.[cCh]
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
	etags -o tags *.[Ch]
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

AndersonMixer.o: AndersonMixer.h Context.h blas.h
AndersonMixer.o: Context.h
AngleCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
AngleCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
AngleCmd.o: Control.h
AngleConstraint.o: AngleConstraint.h Constraint.h D3vector.h AtomSet.h
AngleConstraint.o: Context.h Atom.h UnitCell.h Species.h
AngleConstraint.o: Constraint.h D3vector.h
Atom.o: Atom.h D3vector.h
AtomCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
AtomCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
AtomCmd.o: Control.h
Atom.o: D3vector.h
AtomsDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
AtomsDyn.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
AtomSet.o: AtomSet.h Context.h Atom.h D3vector.h UnitCell.h Species.h
AtomSet.o: NameOf.h
AtomSet.o: Context.h Atom.h D3vector.h UnitCell.h
AtomSetHandler.o: StructureHandler.h D3vector.h
Base64Transcoder.o: Base64Transcoder.h
Basis.o: Basis.h D3vector.h UnitCell.h Context.h
Basis.o: D3vector.h UnitCell.h Context.h
BasisMapping.o: Basis.h D3vector.h UnitCell.h Context.h BasisMapping.h
bisection1.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
bisection1.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
bisection1.o: SampleReader.h Basis.h FourierTransform.h Timer.h SlaterDet.h
bisection1.o: Matrix.h
bisection.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
bisection.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
bisection.o: SampleReader.h SampleWriter.h Basis.h FourierTransform.h Timer.h
bisection.o: SlaterDet.h Matrix.h isodate.h jade.h
bisection_level4.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
bisection_level4.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
bisection_level4.o: SampleReader.h SampleWriter.h Basis.h FourierTransform.h
bisection_level4.o: Timer.h SlaterDet.h Matrix.h isodate.h jade.h
bisection_save.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
bisection_save.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
bisection_save.o: SampleReader.h SampleWriter.h Basis.h FourierTransform.h
bisection_save.o: Timer.h SlaterDet.h Matrix.h isodate.h jade.h
BLYPFunctional.o: BLYPFunctional.h XCFunctional.h
BLYPFunctional.o: XCFunctional.h
BOSampleStepper.o: BOSampleStepper.h SampleStepper.h Timer.h
BOSampleStepper.o: EnergyFunctional.h ChargeDensity.h Context.h
BOSampleStepper.o: StructureFactor.h Sample.h AtomSet.h Atom.h D3vector.h
BOSampleStepper.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
BOSampleStepper.o: Control.h SlaterDet.h Basis.h Matrix.h
BOSampleStepper.o: WavefunctionStepper.h SDWavefunctionStepper.h
BOSampleStepper.o: PSDWavefunctionStepper.h PSDAWavefunctionStepper.h
BOSampleStepper.o: SDIonicStepper.h IonicStepper.h Species.h
BOSampleStepper.o: SDAIonicStepper.h LineMinimizer.h CGIonicStepper.h
BOSampleStepper.o: MDIonicStepper.h SDCellStepper.h CellStepper.h
BOSampleStepper.o: Preconditioner.h AndersonMixer.h MLWFTransform.h
BOSampleStepper.o: BasisMapping.h
BOSampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h ChargeDensity.h
BOSampleStepper.o: Context.h StructureFactor.h Sample.h AtomSet.h Atom.h
BOSampleStepper.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
BOSampleStepper.o: SharedFilePtr.h Control.h
CellDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
CellDyn.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
CellDyn.o: SlaterDet.h Basis.h Matrix.h Timer.h
Cell.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Cell.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
CellLock.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
CellLock.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
CellMass.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
CellMass.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
CellStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
CellStepper.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
CGIonicStepper.o: CGIonicStepper.h IonicStepper.h Sample.h AtomSet.h
CGIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
CGIonicStepper.o: Wavefunction.h SharedFilePtr.h Control.h Species.h
CGIonicStepper.o: LineMinimizer.h
CGIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
CGIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
CGIonicStepper.o: SharedFilePtr.h Control.h Species.h LineMinimizer.h
ChargeDensity.o: ChargeDensity.h Timer.h Context.h Basis.h D3vector.h
ChargeDensity.o: UnitCell.h Wavefunction.h SharedFilePtr.h FourierTransform.h
ChargeDensity.o: SlaterDet.h Matrix.h
ChargeDensity.o: Timer.h Context.h
ChargeMixCoeff.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ChargeMixCoeff.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
ChargeMixNdim.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ChargeMixNdim.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
ChargeMixRcut.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ChargeMixRcut.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
ComputeMLWFCmd.o: ComputeMLWFCmd.h UserInterface.h Sample.h AtomSet.h
ComputeMLWFCmd.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
ComputeMLWFCmd.o: Wavefunction.h SharedFilePtr.h Control.h MLWFTransform.h
ComputeMLWFCmd.o: BasisMapping.h SlaterDet.h Basis.h Matrix.h Timer.h
ComputeMLWFCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ComputeMLWFCmd.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
ComputeMLWFCmd.o: SharedFilePtr.h Control.h MLWFTransform.h BasisMapping.h
ConfinementPotential.o: ConfinementPotential.h Basis.h D3vector.h UnitCell.h
ConfinementPotential.o: Context.h
Constraint.o: Constraint.h
ConstraintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ConstraintCmd.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
ConstraintCmd.o: SharedFilePtr.h Control.h
ConstraintSet.o: ConstraintSet.h DistanceConstraint.h Constraint.h
ConstraintSet.o: AngleConstraint.h D3vector.h TorsionConstraint.h Atom.h
ConstraintSet.o: AtomSet.h Context.h UnitCell.h
Context.o: Context.h
CPSampleStepper.o: CPSampleStepper.h SampleStepper.h Timer.h
CPSampleStepper.o: EnergyFunctional.h ChargeDensity.h Context.h
CPSampleStepper.o: StructureFactor.h Sample.h AtomSet.h Atom.h D3vector.h
CPSampleStepper.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
CPSampleStepper.o: Control.h SlaterDet.h Basis.h Matrix.h
CPSampleStepper.o: MDWavefunctionStepper.h WavefunctionStepper.h
CPSampleStepper.o: MDIonicStepper.h IonicStepper.h Species.h SDCellStepper.h
CPSampleStepper.o: CellStepper.h
CPSampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h ChargeDensity.h
CPSampleStepper.o: Context.h StructureFactor.h Sample.h AtomSet.h Atom.h
CPSampleStepper.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
CPSampleStepper.o: SharedFilePtr.h Control.h
Debug.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Debug.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
DistanceCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
DistanceCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
DistanceCmd.o: Control.h
DistanceConstraint.o: DistanceConstraint.h Constraint.h AtomSet.h Context.h
DistanceConstraint.o: Atom.h D3vector.h UnitCell.h Species.h
DistanceConstraint.o: Constraint.h
Dt.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Dt.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
Ecut.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Ecut.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
Ecutprec.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Ecutprec.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
Ecuts.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Ecuts.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
Emass.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Emass.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
EnergyFunctional.o: EnergyFunctional.h ChargeDensity.h Timer.h Context.h
EnergyFunctional.o: StructureFactor.h Sample.h AtomSet.h Atom.h D3vector.h
EnergyFunctional.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
EnergyFunctional.o: Control.h Species.h SlaterDet.h Basis.h Matrix.h
EnergyFunctional.o: FourierTransform.h XCPotential.h LDAFunctional.h
EnergyFunctional.o: XCFunctional.h PBEFunctional.h BLYPFunctional.h
EnergyFunctional.o: NonLocalPotential.h ConfinementPotential.h blas.h
EnergyFunctional.o: ChargeDensity.h Timer.h Context.h StructureFactor.h
ExchangeOperator.o: ExchangeOperator.h Sample.h AtomSet.h Context.h Atom.h
ExchangeOperator.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
ExchangeOperator.o: SharedFilePtr.h Control.h SlaterDet.h Basis.h Matrix.h
ExchangeOperator.o: Timer.h FourierTransform.h
ExtStress.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ExtStress.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
FermiTemp.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
FermiTemp.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
FoldInWsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
FoldInWsCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
FoldInWsCmd.o: Control.h
FourierTransform.o: FourierTransform.h Timer.h Basis.h D3vector.h UnitCell.h
FourierTransform.o: Context.h
FourierTransform.o: Timer.h
HelpCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
HelpCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
HelpCmd.o: Control.h
IonicStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
IonicStepper.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
IonicStepper.o: Species.h
isodate.o: isodate.h
jacobi.o: Context.h Matrix.h blas.h
jade.o: Context.h Matrix.h blas.h Timer.h
kpgen.o: D3vector.h
KpointCmd.o: UserInterface.h D3vector.h Sample.h AtomSet.h Context.h Atom.h
KpointCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
KpointCmd.o: Control.h
LDAFunctional.o: LDAFunctional.h XCFunctional.h
LDAFunctional.o: XCFunctional.h
ListAtomsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListAtomsCmd.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
ListAtomsCmd.o: SharedFilePtr.h Control.h
ListConstraintsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListConstraintsCmd.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
ListConstraintsCmd.o: SharedFilePtr.h Control.h
ListSpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListSpeciesCmd.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
ListSpeciesCmd.o: SharedFilePtr.h Control.h
LoadCmd.o: LoadCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
LoadCmd.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
LoadCmd.o: SharedFilePtr.h Control.h SampleReader.h
LoadCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
LoadCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
LoadCmd.o: Control.h
Matrix.o: Context.h Matrix.h
MatrixColumnRotation.o: MatrixColumnRotation.h Context.h Matrix.h Timer.h
MatrixColumnRotation.o: Context.h Matrix.h
MDIonicStepper.o: MDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
MDIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
MDIonicStepper.o: Wavefunction.h SharedFilePtr.h Control.h Species.h
MDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
MDIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
MDIonicStepper.o: SharedFilePtr.h Control.h Species.h
MDWavefunctionStepper.o: MDWavefunctionStepper.h WavefunctionStepper.h
MDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
MDWavefunctionStepper.o: SharedFilePtr.h SlaterDet.h Context.h Basis.h
MDWavefunctionStepper.o: Matrix.h Sample.h AtomSet.h Atom.h ConstraintSet.h
MDWavefunctionStepper.o: Control.h
MDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
MLWFTransform.o: MLWFTransform.h D3vector.h BasisMapping.h Basis.h UnitCell.h
MLWFTransform.o: Context.h SlaterDet.h Matrix.h Timer.h jade.h blas.h
MLWFTransform.o: D3vector.h BasisMapping.h
MoveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
MoveCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
MoveCmd.o: Control.h
Nempty.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Nempty.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
NetCharge.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
NetCharge.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
NonLocalPotential.o: NonLocalPotential.h AtomSet.h Context.h Atom.h
NonLocalPotential.o: D3vector.h UnitCell.h Basis.h SlaterDet.h Matrix.h
NonLocalPotential.o: Timer.h Species.h blas.h
NonLocalPotential.o: AtomSet.h Context.h Atom.h D3vector.h UnitCell.h Basis.h
NonLocalPotential.o: SlaterDet.h Matrix.h Timer.h
Nrowmax.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Nrowmax.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
PBEFunctional.o: PBEFunctional.h XCFunctional.h
PBEFunctional.o: XCFunctional.h
plotSample.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
plotSample.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
plotSample.o: SampleReader.h Basis.h FourierTransform.h Timer.h SlaterDet.h
plotSample.o: Matrix.h Species.h ChargeDensity.h atomic_properties.h
Preconditioner.o: Preconditioner.h EnergyFunctional.h ChargeDensity.h Timer.h
Preconditioner.o: Context.h StructureFactor.h Sample.h AtomSet.h Atom.h
Preconditioner.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
Preconditioner.o: SharedFilePtr.h Control.h Basis.h SlaterDet.h Matrix.h
Preconditioner.o: ConfinementPotential.h
PrintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
PrintCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
PrintCmd.o: Control.h
PSDAWavefunctionStepper.o: PSDAWavefunctionStepper.h WavefunctionStepper.h
PSDAWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
PSDAWavefunctionStepper.o: SharedFilePtr.h SlaterDet.h Context.h Basis.h
PSDAWavefunctionStepper.o: Matrix.h Preconditioner.h
PSDAWavefunctionStepper.o: WavefunctionStepper.h Timer.h Wavefunction.h
PSDAWavefunctionStepper.o: D3vector.h UnitCell.h SharedFilePtr.h
PSDWavefunctionStepper.o: PSDWavefunctionStepper.h WavefunctionStepper.h
PSDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
PSDWavefunctionStepper.o: SharedFilePtr.h SlaterDet.h Context.h Basis.h
PSDWavefunctionStepper.o: Matrix.h Preconditioner.h
PSDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
qb.o: isodate.h release.h qbox_xmlns.h Context.h UserInterface.h Sample.h
qb.o: AtomSet.h Atom.h D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
qb.o: SharedFilePtr.h Control.h Timer.h AngleCmd.h AtomCmd.h ComputeMLWFCmd.h
qb.o: MLWFTransform.h BasisMapping.h ConstraintCmd.h DistanceCmd.h
qb.o: FoldInWsCmd.h HelpCmd.h KpointCmd.h ListAtomsCmd.h ListSpeciesCmd.h
qb.o: LoadCmd.h MoveCmd.h PrintCmd.h QuitCmd.h RandomizeWfCmd.h ResetVcmCmd.h
qb.o: RunCmd.h SaveCmd.h SetCmd.h SpeciesCmd.h StatusCmd.h StrainCmd.h
qb.o: TorsionCmd.h AtomsDyn.h Cell.h CellDyn.h SlaterDet.h Basis.h Matrix.h
qb.o: CellLock.h CellMass.h ChargeMixCoeff.h ChargeMixNdim.h ChargeMixRcut.h
qb.o: Debug.h Ecut.h Ecutprec.h Ecuts.h Emass.h ExtStress.h FermiTemp.h Dt.h
qb.o: Nempty.h NetCharge.h Nrowmax.h RefCell.h Stress.h Thermostat.h ThTemp.h
qb.o: ThTime.h ThWidth.h WfDiag.h WfDyn.h Xc.h
qbox_xmlns.o: qbox_xmlns.h
QuitCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
QuitCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
QuitCmd.o: Control.h
RandomizeWfCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
RandomizeWfCmd.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
RandomizeWfCmd.o: SharedFilePtr.h Control.h
RefCell.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
RefCell.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
release.o: release.h
ResetVcmCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
ResetVcmCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
ResetVcmCmd.o: Control.h
RunCmd.o: RunCmd.h UserInterface.h BOSampleStepper.h SampleStepper.h Timer.h
RunCmd.o: EnergyFunctional.h ChargeDensity.h Context.h StructureFactor.h
RunCmd.o: Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
RunCmd.o: Wavefunction.h SharedFilePtr.h Control.h CPSampleStepper.h
RunCmd.o: UserInterface.h
Sample.o: AtomSet.h Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
Sample.o: Wavefunction.h SharedFilePtr.h Control.h
SampleHandler.o: StructureHandler.h
SampleReader.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
SampleReader.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
SampleReader.o: SampleReader.h SpeciesReader.h Species.h Basis.h
SampleReader.o: FourierTransform.h Timer.h SlaterDet.h Matrix.h
SampleReader.o: XMLGFPreprocessor.h
SampleStepper.o: SampleStepper.h Timer.h Sample.h AtomSet.h Context.h Atom.h
SampleStepper.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
SampleStepper.o: SharedFilePtr.h Control.h Species.h
SampleStepper.o: Timer.h
SampleWriter.o: SampleWriter.h Context.h Sample.h AtomSet.h Atom.h D3vector.h
SampleWriter.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
SampleWriter.o: Control.h qbox_xmlns.h Timer.h
SampleWriter.o: Context.h
SaveCmd.o: SaveCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
SaveCmd.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
SaveCmd.o: SharedFilePtr.h Control.h SampleWriter.h isodate.h release.h
SaveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SaveCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
SaveCmd.o: Control.h
SDAIonicStepper.o: SDAIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDAIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
SDAIonicStepper.o: Wavefunction.h SharedFilePtr.h Control.h Species.h
SDAIonicStepper.o: LineMinimizer.h
SDAIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
SDAIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
SDAIonicStepper.o: SharedFilePtr.h Control.h Species.h LineMinimizer.h
SDCellStepper.o: SDCellStepper.h CellStepper.h Sample.h AtomSet.h Context.h
SDCellStepper.o: Atom.h D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
SDCellStepper.o: SharedFilePtr.h Control.h
SDCellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SDCellStepper.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
SDCellStepper.o: Control.h
SDIonicStepper.o: SDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDIonicStepper.o: Context.h Atom.h D3vector.h UnitCell.h ConstraintSet.h
SDIonicStepper.o: Wavefunction.h SharedFilePtr.h Control.h Species.h
SDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
SDIonicStepper.o: D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
SDIonicStepper.o: SharedFilePtr.h Control.h Species.h
SDWavefunctionStepper.o: SDWavefunctionStepper.h WavefunctionStepper.h
SDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
SDWavefunctionStepper.o: SharedFilePtr.h SlaterDet.h Context.h Basis.h
SDWavefunctionStepper.o: Matrix.h Sample.h AtomSet.h Atom.h ConstraintSet.h
SDWavefunctionStepper.o: Control.h
SDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
SetCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SetCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
sinft.o: sinft.h
SlaterDet.o: SlaterDet.h Context.h Basis.h D3vector.h UnitCell.h Matrix.h
SlaterDet.o: Timer.h FourierTransform.h blas.h Base64Transcoder.h
SlaterDet.o: SharedFilePtr.h
SlaterDet.o: Context.h Basis.h D3vector.h UnitCell.h Matrix.h Timer.h
Species.o: Species.h Context.h spline.h sinft.h
SpeciesCmd.o: SpeciesCmd.h UserInterface.h Sample.h AtomSet.h Context.h
SpeciesCmd.o: Atom.h D3vector.h UnitCell.h ConstraintSet.h Wavefunction.h
SpeciesCmd.o: SharedFilePtr.h Control.h SpeciesReader.h Species.h
SpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SpeciesCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
SpeciesCmd.o: Control.h
Species.o: Context.h
SpeciesHandler.o: StructureHandler.h
SpeciesReader.o: Species.h Context.h SpeciesReader.h
SpeciesReader.o: Context.h
spline.o: spline.h
StatusCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
StatusCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
StatusCmd.o: Control.h
StrainCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
StrainCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
StrainCmd.o: Control.h
Stress.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Stress.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
StructuredDocumentHandler.o: StrX.h StructureHandler.h
StructureFactor.o: StructureFactor.h Basis.h D3vector.h UnitCell.h Context.h
testAndersonMixer.o: Context.h AndersonMixer.h
testBase64Transcoder.o: Base64Transcoder.h
testBasis.o: Basis.h D3vector.h UnitCell.h Context.h
testBasisMapping.o: Context.h Basis.h D3vector.h UnitCell.h BasisMapping.h
testBasisMapping.o: Timer.h
testChargeDensity.o: Context.h Wavefunction.h D3vector.h UnitCell.h
testChargeDensity.o: SharedFilePtr.h ChargeDensity.h Timer.h SlaterDet.h
testChargeDensity.o: Basis.h Matrix.h FourierTransform.h
TestCmd.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
TestCmd.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
TestCmd.o: ExchangeOperator.h
testContext.o: Context.h
testEnergyFunctional.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h
testEnergyFunctional.o: UnitCell.h ConstraintSet.h Wavefunction.h
testEnergyFunctional.o: SharedFilePtr.h Control.h EnergyFunctional.h
testEnergyFunctional.o: ChargeDensity.h Timer.h StructureFactor.h
test_fftw3.o: Timer.h
test_fftw.o: Timer.h
testFourierTransform.o: Context.h Basis.h D3vector.h UnitCell.h
testFourierTransform.o: FourierTransform.h Timer.h
testHandler.o: StrX.h
testjacobi.o: Timer.h Context.h Matrix.h jacobi.h
testjade.o: Timer.h Context.h Matrix.h jade.h
testLDAFunctional.o: LDAFunctional.h XCFunctional.h
testLineMinimizer.o: LineMinimizer.h
testMapping.o: Basis.h D3vector.h UnitCell.h Context.h
testMatrix.o: Timer.h Context.h Matrix.h
testMatrixColumnRotation.o: MatrixColumnRotation.h Context.h Matrix.h Timer.h
testMemParse.o: StrX.h
testofstream.o: Context.h
testPParse.o: StrX.h
testreadSample.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
testreadSample.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
testreadSample.o: SampleReader.h
testSample.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
testSample.o: Timer.h Sample.h AtomSet.h Atom.h ConstraintSet.h
testSample.o: Wavefunction.h SharedFilePtr.h Control.h
test_sinft.o: sinft.h readTSC.h
testSlaterDet.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
testSlaterDet.o: Timer.h FourierTransform.h
testSpecies.o: Species.h Context.h SpeciesReader.h
test_spline.o: spline.h
testStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
testStepper.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
testStepper.o: BOSampleStepper.h SampleStepper.h Timer.h EnergyFunctional.h
testStepper.o: ChargeDensity.h StructureFactor.h
testUnitCell.o: UnitCell.h D3vector.h
testWavefunction.o: Context.h Wavefunction.h D3vector.h UnitCell.h
testWavefunction.o: SharedFilePtr.h SlaterDet.h Basis.h Matrix.h Timer.h
testwrite.o: Context.h
testXCFunctional.o: LDAFunctional.h XCFunctional.h PBEFunctional.h Timer.h
testXMLGFPreprocessor.o: Context.h Matrix.h XMLGFPreprocessor.h
Thermostat.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Thermostat.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
ThTemp.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ThTemp.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
ThTime.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ThTime.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
ThWidth.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
ThWidth.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
TorsionCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
TorsionCmd.o: UnitCell.h ConstraintSet.h Wavefunction.h SharedFilePtr.h
TorsionCmd.o: Control.h
TorsionConstraint.o: TorsionConstraint.h Constraint.h D3vector.h AtomSet.h
TorsionConstraint.o: Context.h Atom.h UnitCell.h Species.h
TorsionConstraint.o: Constraint.h D3vector.h
truncateMLWF.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h UnitCell.h
truncateMLWF.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
truncateMLWF.o: SampleReader.h SampleWriter.h Basis.h FourierTransform.h
truncateMLWF.o: Timer.h SlaterDet.h Matrix.h isodate.h jade.h MLWFTransform.h
truncateMLWF.o: BasisMapping.h ComputeMLWFCmd.h UserInterface.h
UnitCell.o: UnitCell.h D3vector.h
UnitCell.o: D3vector.h
UserInterface.o: UserInterface.h qbox_xmlns.h
Wavefunction.o: Wavefunction.h D3vector.h UnitCell.h SharedFilePtr.h
Wavefunction.o: SlaterDet.h Context.h Basis.h Matrix.h Timer.h jacobi.h
Wavefunction.o: D3vector.h UnitCell.h SharedFilePtr.h
WavefunctionHandler.o: StructureHandler.h UnitCell.h D3vector.h
WavefunctionHandler.o: Wavefunction.h SharedFilePtr.h SlaterDet.h Context.h
WavefunctionHandler.o: Basis.h Matrix.h Timer.h
WavefunctionStepper.o: Timer.h
WfDiag.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
WfDiag.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
WfDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
WfDyn.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h SlaterDet.h
WfDyn.o: Basis.h Matrix.h Timer.h
Xc.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h UnitCell.h
Xc.o: ConstraintSet.h Wavefunction.h SharedFilePtr.h Control.h
XCPotential.o: XCPotential.h ChargeDensity.h Timer.h Context.h
XCPotential.o: LDAFunctional.h XCFunctional.h PBEFunctional.h
XCPotential.o: BLYPFunctional.h Basis.h D3vector.h UnitCell.h
XCPotential.o: FourierTransform.h blas.h
XCPotential.o: ChargeDensity.h Timer.h Context.h LDAFunctional.h
XCPotential.o: XCFunctional.h PBEFunctional.h BLYPFunctional.h
xmlextract.o: StrX.h
XMLGFPreprocessor.o: Timer.h Context.h Base64Transcoder.h Matrix.h
XMLGFPreprocessor.o: XMLGFPreprocessor.h
XMLGFPreprocessor.o: Matrix.h
xmlSpecies.o: qbox_xmlns.h
