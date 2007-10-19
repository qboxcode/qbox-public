#-------------------------------------------------------------------------------
# $Id: Makefile,v 1.46 2007-10-19 16:24:36 fgygi Exp $
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
        AndersonMixer.o SDAIonicStepper.o \
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
	$(LD) -o xmlget xmlget.C -I$(XERCESCDIR)/include \
        -L$(XERCESCDIR)/lib -lxerces-c
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
AngleCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
AngleConstraint.o: AngleConstraint.h Constraint.h D3vector.h AtomSet.h
AngleConstraint.o: Context.h Atom.h Species.h
AngleConstraint.o: Constraint.h D3vector.h
Atom.o: Atom.h D3vector.h
AtomCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
AtomCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
Atom.o: D3vector.h
AtomsDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
AtomsDyn.o: Wavefunction.h UnitCell.h Control.h
AtomSet.o: AtomSet.h Context.h Atom.h D3vector.h Species.h NameOf.h
AtomSet.o: Context.h Atom.h D3vector.h
AtomSetHandler.o: StructureHandler.h D3vector.h
Base64Transcoder.o: Base64Transcoder.h
Basis.o: Basis.h D3vector.h UnitCell.h Context.h
Basis.o: D3vector.h UnitCell.h
BasisMapping.o: Basis.h D3vector.h UnitCell.h Context.h BasisMapping.h
bisection.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h ConstraintSet.h
bisection.o: Wavefunction.h UnitCell.h Control.h SampleReader.h
bisection.o: SampleWriter.h Basis.h FourierTransform.h Timer.h SlaterDet.h
bisection.o: Matrix.h isodate.h jade.h
BLACSContext.o: BLACSContext.h blacs.h
BLACSContext.o: blacs.h
BLYPFunctional.o: BLYPFunctional.h XCFunctional.h
BLYPFunctional.o: XCFunctional.h
BOSampleStepper.o: BOSampleStepper.h SampleStepper.h Timer.h
BOSampleStepper.o: EnergyFunctional.h ChargeDensity.h Context.h
BOSampleStepper.o: StructureFactor.h Sample.h AtomSet.h Atom.h D3vector.h
BOSampleStepper.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
BOSampleStepper.o: SlaterDet.h Basis.h Matrix.h WavefunctionStepper.h
BOSampleStepper.o: SDWavefunctionStepper.h PSDWavefunctionStepper.h
BOSampleStepper.o: PSDAWavefunctionStepper.h SDIonicStepper.h IonicStepper.h
BOSampleStepper.o: Species.h SDAIonicStepper.h AndersonMixer.h
BOSampleStepper.o: MDIonicStepper.h SDCellStepper.h CellStepper.h
BOSampleStepper.o: Preconditioner.h
BOSampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h ChargeDensity.h
BOSampleStepper.o: Context.h StructureFactor.h Sample.h AtomSet.h Atom.h
BOSampleStepper.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
BOSampleStepper.o: Control.h
CellDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
CellDyn.o: Wavefunction.h UnitCell.h Control.h SlaterDet.h Basis.h Matrix.h
CellDyn.o: Timer.h
Cell.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Cell.o: Wavefunction.h UnitCell.h Control.h
CellLock.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
CellLock.o: Wavefunction.h UnitCell.h Control.h
CellMass.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
CellMass.o: Wavefunction.h UnitCell.h Control.h
CellStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
CellStepper.o: Wavefunction.h UnitCell.h Control.h
ChargeDensity.o: ChargeDensity.h Timer.h Context.h Basis.h D3vector.h
ChargeDensity.o: UnitCell.h Wavefunction.h FourierTransform.h SlaterDet.h
ChargeDensity.o: Matrix.h
ChargeDensity.o: Timer.h Context.h
ChargeMixCoeff.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
ChargeMixCoeff.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
ChargeMixRcut.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
ChargeMixRcut.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
ComputeMLWFCmd.o: ComputeMLWFCmd.h UserInterface.h Sample.h AtomSet.h
ComputeMLWFCmd.o: Context.h Atom.h D3vector.h ConstraintSet.h Wavefunction.h
ComputeMLWFCmd.o: UnitCell.h Control.h MLWFTransform.h BasisMapping.h
ComputeMLWFCmd.o: SlaterDet.h Basis.h Matrix.h Timer.h
ComputeMLWFCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ComputeMLWFCmd.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
ComputeMLWFCmd.o: Control.h MLWFTransform.h BasisMapping.h
ConfinementPotential.o: ConfinementPotential.h Basis.h D3vector.h UnitCell.h
Constraint.o: Constraint.h
ConstraintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ConstraintCmd.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
ConstraintCmd.o: Control.h
ConstraintSet.o: ConstraintSet.h DistanceConstraint.h Constraint.h
ConstraintSet.o: AngleConstraint.h D3vector.h TorsionConstraint.h Atom.h
ConstraintSet.o: AtomSet.h Context.h
Context.o: Context.h
CPSampleStepper.o: CPSampleStepper.h SampleStepper.h Timer.h
CPSampleStepper.o: EnergyFunctional.h ChargeDensity.h Context.h
CPSampleStepper.o: StructureFactor.h Sample.h AtomSet.h Atom.h D3vector.h
CPSampleStepper.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
CPSampleStepper.o: SlaterDet.h Basis.h Matrix.h MDWavefunctionStepper.h
CPSampleStepper.o: WavefunctionStepper.h MDIonicStepper.h IonicStepper.h
CPSampleStepper.o: Species.h SDCellStepper.h CellStepper.h
CPSampleStepper.o: SampleStepper.h Timer.h EnergyFunctional.h ChargeDensity.h
CPSampleStepper.o: Context.h StructureFactor.h Sample.h AtomSet.h Atom.h
CPSampleStepper.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
CPSampleStepper.o: Control.h
Debug.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Debug.o: Wavefunction.h UnitCell.h Control.h
DiagCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
DiagCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
DistanceCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
DistanceCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
DistanceConstraint.o: DistanceConstraint.h Constraint.h AtomSet.h Context.h
DistanceConstraint.o: Atom.h D3vector.h Species.h
DistanceConstraint.o: Constraint.h
Dt.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Dt.o: Wavefunction.h UnitCell.h Control.h
Ecut.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Ecut.o: Wavefunction.h UnitCell.h Control.h
Ecutprec.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Ecutprec.o: Wavefunction.h UnitCell.h Control.h
Ecuts.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Ecuts.o: Wavefunction.h UnitCell.h Control.h
Emass.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Emass.o: Wavefunction.h UnitCell.h Control.h
EnergyFunctional.o: EnergyFunctional.h ChargeDensity.h Timer.h Context.h
EnergyFunctional.o: StructureFactor.h Sample.h AtomSet.h Atom.h D3vector.h
EnergyFunctional.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
EnergyFunctional.o: Species.h SlaterDet.h Basis.h Matrix.h FourierTransform.h
EnergyFunctional.o: XCPotential.h LDAFunctional.h XCFunctional.h
EnergyFunctional.o: PBEFunctional.h BLYPFunctional.h NonLocalPotential.h
EnergyFunctional.o: ConfinementPotential.h blas.h
EnergyFunctional.o: ChargeDensity.h Timer.h Context.h StructureFactor.h
ExtStress.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
ExtStress.o: Wavefunction.h UnitCell.h Control.h
FermiTemp.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
FermiTemp.o: Wavefunction.h UnitCell.h Control.h
FourierTransform.o: FourierTransform.h Timer.h Basis.h D3vector.h UnitCell.h
FourierTransform.o: Context.h
FourierTransform.o: Timer.h
HelpCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
HelpCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
IonicStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
IonicStepper.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h Species.h
isodate.o: isodate.h
jacobi.o: Context.h Matrix.h blas.h
jade.o: Context.h Matrix.h blas.h Timer.h
LDAFunctional.o: LDAFunctional.h XCFunctional.h
LDAFunctional.o: XCFunctional.h
ListAtomsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListAtomsCmd.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
ListAtomsCmd.o: Control.h
ListConstraintsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListConstraintsCmd.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
ListConstraintsCmd.o: Control.h
ListSpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListSpeciesCmd.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
ListSpeciesCmd.o: Control.h
LoadCmd.o: LoadCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
LoadCmd.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h Control.h
LoadCmd.o: SampleReader.h
LoadCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
LoadCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
Matrix.o: Context.h Matrix.h
MatrixColumnRotation.o: MatrixColumnRotation.h Context.h Matrix.h Timer.h
MatrixColumnRotation.o: Context.h Matrix.h
MDIonicStepper.o: MDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
MDIonicStepper.o: Context.h Atom.h D3vector.h ConstraintSet.h Wavefunction.h
MDIonicStepper.o: UnitCell.h Control.h Species.h
MDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
MDIonicStepper.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
MDIonicStepper.o: Control.h Species.h
MDWavefunctionStepper.o: MDWavefunctionStepper.h WavefunctionStepper.h
MDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
MDWavefunctionStepper.o: SlaterDet.h Context.h Basis.h Matrix.h Sample.h
MDWavefunctionStepper.o: AtomSet.h Atom.h ConstraintSet.h Control.h
MDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
MLWFTransform.o: MLWFTransform.h D3vector.h BasisMapping.h Basis.h UnitCell.h
MLWFTransform.o: SlaterDet.h Context.h Matrix.h Timer.h jade.h blas.h
MLWFTransform.o: D3vector.h BasisMapping.h
MoveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
MoveCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
MultiDistanceConstraint.o: MultiDistanceConstraint.h Constraint.h Atom.h
MultiDistanceConstraint.o: D3vector.h
MultiDistanceConstraint.o: Constraint.h
Nempty.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Nempty.o: Wavefunction.h UnitCell.h Control.h
NetCharge.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
NetCharge.o: Wavefunction.h UnitCell.h Control.h
NonLocalPotential.o: NonLocalPotential.h AtomSet.h Context.h Atom.h
NonLocalPotential.o: D3vector.h Basis.h UnitCell.h SlaterDet.h Matrix.h
NonLocalPotential.o: Timer.h Species.h blas.h
NonLocalPotential.o: AtomSet.h Context.h Atom.h D3vector.h Basis.h UnitCell.h
NonLocalPotential.o: SlaterDet.h Matrix.h Timer.h
Nrowmax.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Nrowmax.o: Wavefunction.h UnitCell.h Control.h
PBEFunctional.o: PBEFunctional.h XCFunctional.h
PBEFunctional.o: XCFunctional.h
plotSample.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h ConstraintSet.h
plotSample.o: Wavefunction.h UnitCell.h Control.h SampleReader.h Basis.h
plotSample.o: FourierTransform.h Timer.h SlaterDet.h Matrix.h Species.h
plotSample.o: ChargeDensity.h atomic_properties.h
Preconditioner.o: Preconditioner.h EnergyFunctional.h ChargeDensity.h Timer.h
Preconditioner.o: Context.h StructureFactor.h Sample.h AtomSet.h Atom.h
Preconditioner.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
Preconditioner.o: Control.h Basis.h SlaterDet.h Matrix.h
Preconditioner.o: ConfinementPotential.h
PrintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
PrintCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
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
qb.o: AtomSet.h Atom.h D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
qb.o: Control.h Timer.h AngleCmd.h AtomCmd.h ComputeMLWFCmd.h MLWFTransform.h
qb.o: BasisMapping.h ConstraintCmd.h DistanceCmd.h HelpCmd.h ListAtomsCmd.h
qb.o: ListSpeciesCmd.h LoadCmd.h MoveCmd.h PrintCmd.h QuitCmd.h
qb.o: RandomizeWfCmd.h ResetVcmCmd.h RunCmd.h SaveCmd.h SetCmd.h SpeciesCmd.h
qb.o: StatusCmd.h TorsionCmd.h AtomsDyn.h Cell.h CellDyn.h SlaterDet.h
qb.o: Basis.h Matrix.h CellLock.h CellMass.h ChargeMixCoeff.h ChargeMixRcut.h
qb.o: Debug.h Ecut.h Ecutprec.h Ecuts.h Emass.h ExtStress.h FermiTemp.h Dt.h
qb.o: Nempty.h NetCharge.h Nrowmax.h RefCell.h Stress.h Thermostat.h ThTemp.h
qb.o: ThTime.h ThWidth.h WfDiag.h WfDyn.h Xc.h
qbox_xmlns.o: qbox_xmlns.h
QuitCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
QuitCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
RandomizeWfCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
RandomizeWfCmd.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
RandomizeWfCmd.o: Control.h
RefCell.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
RefCell.o: Wavefunction.h UnitCell.h Control.h
release.o: release.h
ResetVcmCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
ResetVcmCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
RunCmd.o: RunCmd.h UserInterface.h BOSampleStepper.h SampleStepper.h Timer.h
RunCmd.o: EnergyFunctional.h ChargeDensity.h Context.h StructureFactor.h
RunCmd.o: Sample.h AtomSet.h Atom.h D3vector.h ConstraintSet.h Wavefunction.h
RunCmd.o: UnitCell.h Control.h CPSampleStepper.h
RunCmd.o: UserInterface.h
Sample.o: AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Sample.o: Wavefunction.h UnitCell.h Control.h
SampleHandler.o: StructureHandler.h
SampleReader.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
SampleReader.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
SampleReader.o: SampleReader.h SpeciesReader.h Species.h Basis.h
SampleReader.o: FourierTransform.h Timer.h SlaterDet.h Matrix.h
SampleReader.o: XMLGFPreprocessor.h
SampleStepper.o: SampleStepper.h Timer.h Sample.h AtomSet.h Context.h Atom.h
SampleStepper.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
SampleStepper.o: Control.h Species.h
SampleStepper.o: Timer.h
SampleViewer.o: SampleViewer.h D3vector.h Sample.h AtomSet.h Context.h Atom.h
SampleViewer.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
SampleViewer.o: D3vector.h
SampleWriter.o: SampleWriter.h Context.h Sample.h AtomSet.h Atom.h D3vector.h
SampleWriter.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
SampleWriter.o: qbox_xmlns.h
SampleWriter.o: Context.h
SaveCmd.o: SaveCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
SaveCmd.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h Control.h
SaveCmd.o: SampleWriter.h isodate.h release.h
SaveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SaveCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
scansample.o: Timer.h
SDAIonicStepper.o: SDAIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDAIonicStepper.o: Context.h Atom.h D3vector.h ConstraintSet.h Wavefunction.h
SDAIonicStepper.o: UnitCell.h Control.h Species.h AndersonMixer.h
SDAIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
SDAIonicStepper.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
SDAIonicStepper.o: Control.h Species.h AndersonMixer.h
SDCellStepper.o: SDCellStepper.h CellStepper.h Sample.h AtomSet.h Context.h
SDCellStepper.o: Atom.h D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
SDCellStepper.o: Control.h
SDCellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SDCellStepper.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
SDIonicStepper.o: SDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDIonicStepper.o: Context.h Atom.h D3vector.h ConstraintSet.h Wavefunction.h
SDIonicStepper.o: UnitCell.h Control.h Species.h
SDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
SDIonicStepper.o: D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
SDIonicStepper.o: Control.h Species.h
SDWavefunctionStepper.o: SDWavefunctionStepper.h WavefunctionStepper.h
SDWavefunctionStepper.o: Timer.h Wavefunction.h D3vector.h UnitCell.h
SDWavefunctionStepper.o: SlaterDet.h Context.h Basis.h Matrix.h Sample.h
SDWavefunctionStepper.o: AtomSet.h Atom.h ConstraintSet.h Control.h
SDWavefunctionStepper.o: WavefunctionStepper.h Timer.h
SetCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SetCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
sinft.o: sinft.h
SlaterDet.o: SlaterDet.h Context.h Basis.h D3vector.h UnitCell.h Matrix.h
SlaterDet.o: Timer.h FourierTransform.h blas.h Base64Transcoder.h
SlaterDet.o: Context.h Basis.h D3vector.h UnitCell.h Matrix.h Timer.h
Species.o: Species.h Context.h spline.h sinft.h
SpeciesCmd.o: SpeciesCmd.h UserInterface.h Sample.h AtomSet.h Context.h
SpeciesCmd.o: Atom.h D3vector.h ConstraintSet.h Wavefunction.h UnitCell.h
SpeciesCmd.o: Control.h SpeciesReader.h Species.h
SpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SpeciesCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
Species.o: Context.h
SpeciesHandler.o: StructureHandler.h
SpeciesReader.o: Species.h Context.h SpeciesReader.h
SpeciesReader.o: Context.h
spline.o: spline.h
StatusCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
StatusCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
Stress.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Stress.o: Wavefunction.h UnitCell.h Control.h
StructuredDocumentHandler.o: StrX.h StructureHandler.h
StructureFactor.o: StructureFactor.h Basis.h D3vector.h UnitCell.h
testAndersonMixer.o: Context.h AndersonMixer.h
testBase64Transcoder.o: Base64Transcoder.h
testBasis.o: Basis.h D3vector.h UnitCell.h Context.h
testBasisMapping.o: Context.h Basis.h D3vector.h UnitCell.h BasisMapping.h
testBasisMapping.o: Timer.h
testChargeDensity.o: Context.h Wavefunction.h D3vector.h UnitCell.h
testChargeDensity.o: ChargeDensity.h Timer.h SlaterDet.h Basis.h Matrix.h
testChargeDensity.o: FourierTransform.h
testContext.o: Context.h
testEnergyFunctional.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h
testEnergyFunctional.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
testEnergyFunctional.o: EnergyFunctional.h ChargeDensity.h Timer.h
testEnergyFunctional.o: StructureFactor.h
test_fftw3.o: Timer.h
test_fftw.o: Timer.h
testFourierTransform.o: Context.h Basis.h D3vector.h UnitCell.h
testFourierTransform.o: FourierTransform.h Timer.h
testHandler.o: StrX.h
TestHandler.o: StrX.h
testjacobi.o: Timer.h Context.h Matrix.h jacobi.h
testjade.o: Timer.h Context.h Matrix.h jade.h
testLDAFunctional.o: LDAFunctional.h XCFunctional.h
testMatrix.o: Timer.h Context.h Matrix.h
testMatrixColumnRotation.o: MatrixColumnRotation.h Context.h Matrix.h Timer.h
testMemParse.o: StrX.h
testofstream.o: Context.h
testPParse.o: StrX.h
testreadSample.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h
testreadSample.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
testreadSample.o: SampleReader.h
testSample.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
testSample.o: Timer.h Sample.h AtomSet.h Atom.h ConstraintSet.h
testSample.o: Wavefunction.h Control.h
test_sinft.o: sinft.h readTSC.h
testSlaterDet.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
testSlaterDet.o: Timer.h FourierTransform.h
testSpecies.o: Species.h Context.h SpeciesReader.h
test_spline.o: spline.h
testStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
testStepper.o: Wavefunction.h UnitCell.h Control.h BOSampleStepper.h
testStepper.o: SampleStepper.h Timer.h EnergyFunctional.h ChargeDensity.h
testStepper.o: StructureFactor.h
testUnitCell.o: UnitCell.h D3vector.h
testWavefunction.o: Context.h Wavefunction.h D3vector.h UnitCell.h
testWavefunction.o: SlaterDet.h Basis.h Matrix.h Timer.h
testwrite.o: Context.h
testXCFunctional.o: LDAFunctional.h XCFunctional.h PBEFunctional.h Timer.h
testXMLGFPreprocessor.o: Context.h Matrix.h XMLGFPreprocessor.h
Thermostat.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Thermostat.o: Wavefunction.h UnitCell.h Control.h
ThTemp.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
ThTemp.o: Wavefunction.h UnitCell.h Control.h
ThTime.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
ThTime.o: Wavefunction.h UnitCell.h Control.h
ThWidth.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
ThWidth.o: Wavefunction.h UnitCell.h Control.h
TorsionCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
TorsionCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
TorsionConstraint.o: TorsionConstraint.h Constraint.h D3vector.h AtomSet.h
TorsionConstraint.o: Context.h Atom.h Species.h
TorsionConstraint.o: Constraint.h D3vector.h
UnitCell.o: UnitCell.h D3vector.h
UnitCell.o: D3vector.h
UserInterface.o: UserInterface.h
ViewCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
ViewCmd.o: ConstraintSet.h Wavefunction.h UnitCell.h Control.h
Wavefunction.o: Wavefunction.h D3vector.h UnitCell.h SlaterDet.h Context.h
Wavefunction.o: Basis.h Matrix.h Timer.h jacobi.h
Wavefunction.o: D3vector.h UnitCell.h
WavefunctionHandler.o: StructureHandler.h UnitCell.h D3vector.h
WavefunctionHandler.o: Wavefunction.h SlaterDet.h Context.h Basis.h Matrix.h
WavefunctionHandler.o: Timer.h
WavefunctionStepper.o: Timer.h
WfDiag.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
WfDiag.o: Wavefunction.h UnitCell.h Control.h
WfDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
WfDyn.o: Wavefunction.h UnitCell.h Control.h SlaterDet.h Basis.h Matrix.h
WfDyn.o: Timer.h
Xc.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h ConstraintSet.h
Xc.o: Wavefunction.h UnitCell.h Control.h
XCPotential.o: XCPotential.h ChargeDensity.h Timer.h Context.h
XCPotential.o: LDAFunctional.h XCFunctional.h PBEFunctional.h
XCPotential.o: BLYPFunctional.h Basis.h D3vector.h UnitCell.h
XCPotential.o: FourierTransform.h blas.h
XCPotential.o: ChargeDensity.h Timer.h Context.h LDAFunctional.h
XCPotential.o: XCFunctional.h PBEFunctional.h BLYPFunctional.h
xml2xy.o: StrX.h
xmlextract.o: StrX.h
xmlget.o: StrX.h
XMLGFPreprocessor.o: Timer.h Context.h Base64Transcoder.h Matrix.h
XMLGFPreprocessor.o: XMLGFPreprocessor.h
XMLGFPreprocessor.o: Matrix.h
xmlSpecies.o: qbox_xmlns.h
