#-------------------------------------------------------------------------------
# $Id: Makefile,v 1.27 2004-09-14 22:24:11 fgygi Exp $
#------------------------------------------------------------------------------
#
include $(TARGET).mk
#------------------------------------------------------------------------------
#
EXECS=qb testMatrix testFTGrid testGridFunction testBasis \
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
        WavefunctionStepper.o SDWavefunctionStepper.o \
        MDWavefunctionStepper.o SDIonicStepper.o MDIonicStepper.o \
        PSDWavefunctionStepper.o PSDAWavefunctionStepper.o \
        SDCellStepper.o ConfinementPotential.o Preconditioner.o \
        release.o qbox_xmlns.o isodate.o $(PLTOBJECTS)
CXXFLAGS += -DTARGET='"$(TARGET)"'
 qb:    $(OBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 lib:   $(OBJECTS)
	ar cq libqb.a $^
 SamplePrint: SamplePrint.o SamplePrintHandlers.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSample: testSample.o AtomSet.o Atom.o Species.o \
        Wavefunction.o SlaterDet.o \
        Basis.o FourierTransform.o Matrix.o Context.o \
        sinft.o spline.o UnitCell.o 
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
        Basis.o UnitCell.o Matrix.o Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testSpecies: testSpecies.o Species.o sinft.o spline.o Context.o \
	SpeciesReader.o StructuredDocumentHandler.o SpeciesHandler.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testXCFunctional: testXCFunctional.o LDAFunctional.o PBEFunctional.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testMatrix: testMatrix.o Matrix.o Context.o
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
 test_fftw: test_fftw.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 test_fftw3: test_fftw3.o
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
	rm -f $(EXECS)
	rm -rf ti_files
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
	tar czf qbsrc.tgz *.[Ch] Makefile *.mk
#------------------------------------------------------------------------------
# DO NOT DELETE

Atom.o: Atom.h D3vector.h
AtomCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
AtomCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
Atom.o: D3vector.h
AtomsDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
AtomsDyn.o: Wavefunction.h UnitCell.h Control.h
AtomSet.o: AtomSet.h Context.h Atom.h D3vector.h Species.h NameOf.h
AtomSet.o: Context.h Atom.h D3vector.h Species.h
AtomSetHandler.o: AtomSetHandler.h StructureHandler.h D3vector.h AtomSet.h
AtomSetHandler.o: Context.h Atom.h Species.h SpeciesHandler.h SpeciesReader.h
AtomSetHandler.o: StrX.h
AtomSetHandler.o: StructureHandler.h D3vector.h
Base64Transcoder.o: Base64Transcoder.h
Basis.o: Basis.h D3vector.h UnitCell.h Context.h
Basis.o: D3vector.h UnitCell.h
BLYPFunctional.o: BLYPFunctional.h XCFunctional.h
BLYPFunctional.o: XCFunctional.h
BOSampleStepper.o: BOSampleStepper.h SampleStepper.h Sample.h AtomSet.h
BOSampleStepper.o: Context.h Atom.h D3vector.h Species.h Wavefunction.h
BOSampleStepper.o: UnitCell.h Control.h Timer.h EnergyFunctional.h
BOSampleStepper.o: ChargeDensity.h StructureFactor.h SlaterDet.h Basis.h
BOSampleStepper.o: Matrix.h WavefunctionStepper.h SDWavefunctionStepper.h
BOSampleStepper.o: PSDWavefunctionStepper.h PSDAWavefunctionStepper.h
BOSampleStepper.o: SDIonicStepper.h IonicStepper.h MDIonicStepper.h
BOSampleStepper.o: SDCellStepper.h CellStepper.h Preconditioner.h
BOSampleStepper.o: SampleStepper.h Sample.h AtomSet.h Context.h Atom.h
BOSampleStepper.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
BOSampleStepper.o: Timer.h EnergyFunctional.h ChargeDensity.h
BOSampleStepper.o: StructureFactor.h
CellDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
CellDyn.o: Wavefunction.h UnitCell.h Control.h SlaterDet.h Basis.h Matrix.h
CellDyn.o: Timer.h
Cell.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Cell.o: Wavefunction.h UnitCell.h Control.h
CellLock.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
CellLock.o: Wavefunction.h UnitCell.h Control.h
CellMass.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
CellMass.o: Wavefunction.h UnitCell.h Control.h
CellStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
CellStepper.o: Wavefunction.h UnitCell.h Control.h
ChargeDensity.o: ChargeDensity.h Context.h Basis.h D3vector.h UnitCell.h
ChargeDensity.o: Wavefunction.h FourierTransform.h SlaterDet.h Matrix.h
ChargeDensity.o: Timer.h
ChargeDensity.o: Context.h
ConfinementPotential.o: ConfinementPotential.h Basis.h D3vector.h UnitCell.h
Context.o: Context.h
CPSampleStepper.o: CPSampleStepper.h SampleStepper.h Sample.h AtomSet.h
CPSampleStepper.o: Context.h Atom.h D3vector.h Species.h Wavefunction.h
CPSampleStepper.o: UnitCell.h Control.h Timer.h EnergyFunctional.h
CPSampleStepper.o: ChargeDensity.h StructureFactor.h SlaterDet.h Basis.h
CPSampleStepper.o: Matrix.h MDWavefunctionStepper.h WavefunctionStepper.h
CPSampleStepper.o: MDIonicStepper.h IonicStepper.h SDCellStepper.h
CPSampleStepper.o: CellStepper.h
CPSampleStepper.o: SampleStepper.h Sample.h AtomSet.h Context.h Atom.h
CPSampleStepper.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
CPSampleStepper.o: Timer.h EnergyFunctional.h ChargeDensity.h
CPSampleStepper.o: StructureFactor.h
Debug.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Debug.o: Wavefunction.h UnitCell.h Control.h
Dt.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h Wavefunction.h
Dt.o: UnitCell.h Control.h
Ecut.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Ecut.o: Wavefunction.h UnitCell.h Control.h
Ecutprec.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Ecutprec.o: Wavefunction.h UnitCell.h Control.h
Ecuts.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Ecuts.o: Wavefunction.h UnitCell.h Control.h
Emass.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Emass.o: Wavefunction.h UnitCell.h Control.h
EnergyFunctional.o: EnergyFunctional.h ChargeDensity.h Context.h
EnergyFunctional.o: StructureFactor.h Timer.h Sample.h AtomSet.h Atom.h
EnergyFunctional.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
EnergyFunctional.o: SlaterDet.h Basis.h Matrix.h FourierTransform.h
EnergyFunctional.o: XCPotential.h LDAFunctional.h XCFunctional.h
EnergyFunctional.o: PBEFunctional.h BLYPFunctional.h NonLocalPotential.h
EnergyFunctional.o: ConfinementPotential.h blas.h
EnergyFunctional.o: ChargeDensity.h Context.h StructureFactor.h Timer.h
ExtStress.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
ExtStress.o: Wavefunction.h UnitCell.h Control.h
FermiTemp.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
FermiTemp.o: Wavefunction.h UnitCell.h Control.h
FourierTransform.o: FourierTransform.h Basis.h D3vector.h UnitCell.h
FourierTransform.o: Context.h Timer.h
HelpCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
HelpCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
IonicStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
IonicStepper.o: Wavefunction.h UnitCell.h Control.h
isodate.o: isodate.h
LDAFunctional.o: LDAFunctional.h XCFunctional.h
LDAFunctional.o: XCFunctional.h
ListAtomsCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListAtomsCmd.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
ListSpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
ListSpeciesCmd.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
LoadCmd.o: LoadCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
LoadCmd.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
LoadCmd.o: SampleReader.h
LoadCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
LoadCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
Matrix.o: Context.h Matrix.h
MDIonicStepper.o: MDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
MDIonicStepper.o: Context.h Atom.h D3vector.h Species.h Wavefunction.h
MDIonicStepper.o: UnitCell.h Control.h
MDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
MDIonicStepper.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
MDWavefunctionStepper.o: MDWavefunctionStepper.h WavefunctionStepper.h
MDWavefunctionStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
MDWavefunctionStepper.o: Species.h Wavefunction.h UnitCell.h Control.h
MDWavefunctionStepper.o: Timer.h SlaterDet.h Basis.h Matrix.h
MDWavefunctionStepper.o: WavefunctionStepper.h Sample.h AtomSet.h Context.h
MDWavefunctionStepper.o: Atom.h D3vector.h Species.h Wavefunction.h
MDWavefunctionStepper.o: UnitCell.h Control.h Timer.h
Nempty.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Nempty.o: Wavefunction.h UnitCell.h Control.h
NonLocalPotential.o: NonLocalPotential.h AtomSet.h Context.h Atom.h
NonLocalPotential.o: D3vector.h Species.h Basis.h UnitCell.h SlaterDet.h
NonLocalPotential.o: Matrix.h Timer.h blas.h
NonLocalPotential.o: AtomSet.h Context.h Atom.h D3vector.h Species.h Basis.h
NonLocalPotential.o: UnitCell.h SlaterDet.h Matrix.h Timer.h
Nrowmax.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Nrowmax.o: Wavefunction.h UnitCell.h Control.h
PBEFunctional.o: PBEFunctional.h XCFunctional.h
PBEFunctional.o: XCFunctional.h
Preconditioner.o: Preconditioner.h EnergyFunctional.h ChargeDensity.h
Preconditioner.o: Context.h StructureFactor.h Timer.h Sample.h AtomSet.h
Preconditioner.o: Atom.h D3vector.h Species.h Wavefunction.h UnitCell.h
Preconditioner.o: Control.h Basis.h SlaterDet.h Matrix.h
Preconditioner.o: ConfinementPotential.h
PrintCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
PrintCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
PSDAWavefunctionStepper.o: PSDAWavefunctionStepper.h WavefunctionStepper.h
PSDAWavefunctionStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
PSDAWavefunctionStepper.o: Species.h Wavefunction.h UnitCell.h Control.h
PSDAWavefunctionStepper.o: Timer.h SlaterDet.h Basis.h Matrix.h
PSDAWavefunctionStepper.o: Preconditioner.h
PSDAWavefunctionStepper.o: WavefunctionStepper.h Sample.h AtomSet.h Context.h
PSDAWavefunctionStepper.o: Atom.h D3vector.h Species.h Wavefunction.h
PSDAWavefunctionStepper.o: UnitCell.h Control.h Timer.h
PSDWavefunctionStepper.o: PSDWavefunctionStepper.h WavefunctionStepper.h
PSDWavefunctionStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
PSDWavefunctionStepper.o: Species.h Wavefunction.h UnitCell.h Control.h
PSDWavefunctionStepper.o: Timer.h SlaterDet.h Basis.h Matrix.h
PSDWavefunctionStepper.o: Preconditioner.h
PSDWavefunctionStepper.o: WavefunctionStepper.h Sample.h AtomSet.h Context.h
PSDWavefunctionStepper.o: Atom.h D3vector.h Species.h Wavefunction.h
PSDWavefunctionStepper.o: UnitCell.h Control.h Timer.h
qb.o: isodate.h release.h qbox_xmlns.h Context.h UserInterface.h Sample.h
qb.o: AtomSet.h Atom.h D3vector.h Species.h Wavefunction.h UnitCell.h
qb.o: Control.h Timer.h AtomCmd.h HelpCmd.h ListAtomsCmd.h ListSpeciesCmd.h
qb.o: LoadCmd.h PrintCmd.h QuitCmd.h RandomizeWfCmd.h RunCmd.h SaveCmd.h
qb.o: SetCmd.h SpeciesCmd.h StatusCmd.h AtomsDyn.h Cell.h CellDyn.h
qb.o: SlaterDet.h Basis.h Matrix.h CellLock.h CellMass.h Debug.h Ecut.h
qb.o: Ecutprec.h Ecuts.h Emass.h ExtStress.h FermiTemp.h Dt.h Nempty.h
qb.o: Nrowmax.h RefCell.h Stress.h Thermostat.h ThTemp.h ThTime.h ThWidth.h
qb.o: WfDiag.h WfDyn.h Xc.h
qbox_xmlns.o: qbox_xmlns.h
QuitCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
QuitCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
RandomizeWfCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h
RandomizeWfCmd.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
RefCell.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
RefCell.o: Wavefunction.h UnitCell.h Control.h
release.o: release.h
RunCmd.o: RunCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
RunCmd.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
RunCmd.o: BOSampleStepper.h SampleStepper.h Timer.h EnergyFunctional.h
RunCmd.o: ChargeDensity.h StructureFactor.h CPSampleStepper.h
RunCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
RunCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
Sample.o: AtomSet.h Context.h Atom.h D3vector.h Species.h Wavefunction.h
Sample.o: UnitCell.h Control.h
SampleHandler.o: SampleHandler.h StructureHandler.h Sample.h AtomSet.h
SampleHandler.o: Context.h Atom.h D3vector.h Species.h Wavefunction.h
SampleHandler.o: UnitCell.h Control.h AtomSetHandler.h WavefunctionHandler.h
SampleHandler.o: SlaterDet.h Basis.h Matrix.h Timer.h StrX.h
SampleHandler.o: StructureHandler.h
SampleReader.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
SampleReader.o: Wavefunction.h UnitCell.h Control.h SampleReader.h
SampleReader.o: SpeciesReader.h StructuredDocumentHandler.h StrX.h
SampleReader.o: StructureHandler.h SampleHandler.h Basis.h FourierTransform.h
SampleReader.o: SlaterDet.h Matrix.h Timer.h XMLGFPreprocessor.h
SampleReader.o: Context.h
SampleStepper.o: SampleStepper.h Sample.h AtomSet.h Context.h Atom.h
SampleStepper.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
SampleStepper.o: Timer.h
SampleStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
SampleStepper.o: Wavefunction.h UnitCell.h Control.h Timer.h
SaveCmd.o: SaveCmd.h UserInterface.h Sample.h AtomSet.h Context.h Atom.h
SaveCmd.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h isodate.h
SaveCmd.o: release.h qbox_xmlns.h
SaveCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SaveCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
SDCellStepper.o: SDCellStepper.h CellStepper.h Sample.h AtomSet.h Context.h
SDCellStepper.o: Atom.h D3vector.h Species.h Wavefunction.h UnitCell.h
SDCellStepper.o: Control.h
SDCellStepper.o: CellStepper.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SDCellStepper.o: Species.h Wavefunction.h UnitCell.h Control.h
SDIonicStepper.o: SDIonicStepper.h IonicStepper.h Sample.h AtomSet.h
SDIonicStepper.o: Context.h Atom.h D3vector.h Species.h Wavefunction.h
SDIonicStepper.o: UnitCell.h Control.h
SDIonicStepper.o: IonicStepper.h Sample.h AtomSet.h Context.h Atom.h
SDIonicStepper.o: D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
SDWavefunctionStepper.o: SDWavefunctionStepper.h WavefunctionStepper.h
SDWavefunctionStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
SDWavefunctionStepper.o: Species.h Wavefunction.h UnitCell.h Control.h
SDWavefunctionStepper.o: Timer.h SlaterDet.h Basis.h Matrix.h
SDWavefunctionStepper.o: WavefunctionStepper.h Sample.h AtomSet.h Context.h
SDWavefunctionStepper.o: Atom.h D3vector.h Species.h Wavefunction.h
SDWavefunctionStepper.o: UnitCell.h Control.h Timer.h
SetCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SetCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
sinft.o: sinft.h
SlaterDet.o: SlaterDet.h Context.h Basis.h D3vector.h UnitCell.h Matrix.h
SlaterDet.o: Timer.h FourierTransform.h blas.h
SlaterDet.o: Context.h Basis.h D3vector.h UnitCell.h Matrix.h Timer.h
Species.o: Species.h Context.h spline.h sinft.h
SpeciesCmd.o: SpeciesCmd.h UserInterface.h Sample.h AtomSet.h Context.h
SpeciesCmd.o: Atom.h D3vector.h Species.h Wavefunction.h UnitCell.h Control.h
SpeciesCmd.o: SpeciesReader.h
SpeciesCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
SpeciesCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
Species.o: Context.h
SpeciesHandler.o: SpeciesHandler.h StructureHandler.h Species.h Context.h
SpeciesHandler.o: StrX.h
SpeciesHandler.o: StructureHandler.h
SpeciesReader.o: Species.h Context.h SpeciesReader.h
SpeciesReader.o: StructuredDocumentHandler.h StrX.h StructureHandler.h
SpeciesReader.o: SpeciesHandler.h
SpeciesReader.o: Context.h
SpeciesReaderHandler.o: SpeciesReaderHandler.h Species.h Context.h StrX.h
SpeciesReaderHandler.o: Species.h Context.h StrX.h
spline.o: spline.h
StatusCmd.o: UserInterface.h Sample.h AtomSet.h Context.h Atom.h D3vector.h
StatusCmd.o: Species.h Wavefunction.h UnitCell.h Control.h
Stress.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Stress.o: Wavefunction.h UnitCell.h Control.h
StructuredDocumentHandler.o: StructuredDocumentHandler.h StrX.h
StructuredDocumentHandler.o: StructureHandler.h
StructuredDocumentHandler.o: StrX.h StructureHandler.h
StructureFactor.o: StructureFactor.h Basis.h D3vector.h UnitCell.h
testBasis.o: Basis.h D3vector.h UnitCell.h Context.h
testChargeDensity.o: Context.h Wavefunction.h D3vector.h UnitCell.h
testChargeDensity.o: ChargeDensity.h SlaterDet.h Basis.h Matrix.h Timer.h
testChargeDensity.o: FourierTransform.h
testContext.o: Context.h
testEnergyFunctional.o: Context.h Sample.h AtomSet.h Atom.h D3vector.h
testEnergyFunctional.o: Species.h Wavefunction.h UnitCell.h Control.h
testEnergyFunctional.o: EnergyFunctional.h ChargeDensity.h StructureFactor.h
testEnergyFunctional.o: Timer.h
test_fftw3.o: Timer.h
test_fftw.o: Timer.h
testFourierTransform.o: Context.h Basis.h D3vector.h UnitCell.h
testFourierTransform.o: FourierTransform.h Timer.h
testHandler.o: StrX.h
TestHandler.o: StrX.h
testLDAFunctional.o: LDAFunctional.h XCFunctional.h
testMatrix.o: Timer.h Context.h Matrix.h
testMemParse.o: StrX.h
testofstream.o: Context.h
testPParse.o: StrX.h
testSample.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
testSample.o: Timer.h Sample.h AtomSet.h Atom.h Species.h Wavefunction.h
testSample.o: Control.h
testSlaterDet.o: Context.h SlaterDet.h Basis.h D3vector.h UnitCell.h Matrix.h
testSlaterDet.o: Timer.h FourierTransform.h
testSpecies.o: Species.h Context.h SpeciesReader.h
testStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
testStepper.o: Wavefunction.h UnitCell.h Control.h BOSampleStepper.h
testStepper.o: SampleStepper.h Timer.h EnergyFunctional.h ChargeDensity.h
testStepper.o: StructureFactor.h
testUnitCell.o: UnitCell.h D3vector.h
testWavefunction.o: Context.h Wavefunction.h D3vector.h UnitCell.h
testWavefunction.o: SlaterDet.h Basis.h Matrix.h Timer.h
testwrite.o: Context.h
testXCFunctional.o: LDAFunctional.h XCFunctional.h PBEFunctional.h Timer.h
testXMLGFPreprocessor.o: Context.h Matrix.h XMLGFPreprocessor.h
Thermostat.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
Thermostat.o: Wavefunction.h UnitCell.h Control.h
ThTemp.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
ThTemp.o: Wavefunction.h UnitCell.h Control.h
ThTime.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
ThTime.o: Wavefunction.h UnitCell.h Control.h
ThWidth.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
ThWidth.o: Wavefunction.h UnitCell.h Control.h
UnitCell.o: UnitCell.h D3vector.h
UnitCell.o: D3vector.h
UserInterface.o: UserInterface.h
Wavefunction.o: Wavefunction.h D3vector.h UnitCell.h SlaterDet.h Context.h
Wavefunction.o: Basis.h Matrix.h Timer.h
Wavefunction.o: D3vector.h UnitCell.h
WavefunctionHandler.o: WavefunctionHandler.h StructureHandler.h UnitCell.h
WavefunctionHandler.o: D3vector.h Wavefunction.h SlaterDet.h Context.h
WavefunctionHandler.o: Basis.h Matrix.h Timer.h FourierTransform.h StrX.h
WavefunctionHandler.o: StructureHandler.h UnitCell.h D3vector.h
WavefunctionHandler.o: Wavefunction.h SlaterDet.h Context.h Basis.h Matrix.h
WavefunctionHandler.o: Timer.h
WavefunctionStepper.o: WavefunctionStepper.h Sample.h AtomSet.h Context.h
WavefunctionStepper.o: Atom.h D3vector.h Species.h Wavefunction.h UnitCell.h
WavefunctionStepper.o: Control.h Timer.h SlaterDet.h Basis.h Matrix.h
WavefunctionStepper.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h
WavefunctionStepper.o: Species.h Wavefunction.h UnitCell.h Control.h Timer.h
WfDiag.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
WfDiag.o: Wavefunction.h UnitCell.h Control.h
WfDyn.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h
WfDyn.o: Wavefunction.h UnitCell.h Control.h SlaterDet.h Basis.h Matrix.h
WfDyn.o: Timer.h
Xc.o: Sample.h AtomSet.h Context.h Atom.h D3vector.h Species.h Wavefunction.h
Xc.o: UnitCell.h Control.h
XCPotential.o: XCPotential.h ChargeDensity.h Context.h LDAFunctional.h
XCPotential.o: XCFunctional.h PBEFunctional.h BLYPFunctional.h Basis.h
XCPotential.o: D3vector.h UnitCell.h FourierTransform.h blas.h
XCPotential.o: ChargeDensity.h Context.h LDAFunctional.h XCFunctional.h
XCPotential.o: PBEFunctional.h BLYPFunctional.h
xmlextract.o: StrX.h
xmlget.o: StrX.h
XMLGFPreprocessor.o: Timer.h Context.h Base64Transcoder.h Matrix.h
XMLGFPreprocessor.o: XMLGFPreprocessor.h
XMLGFPreprocessor.o: Matrix.h
xmlSpecies.o: qbox_xmlns.h
