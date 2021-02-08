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
        ANDIonicStepper.o \
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
        GlobalExtForce.o sampling.o MPIdata.o cout0.o
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
        XMLGFPreprocessor.o Base64Transcoder.o MPIdata.o cout0.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testChargeDensity: testChargeDensity.o ChargeDensity.o \
        Wavefunction.o SlaterDet.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o UnitCell.o \
        Context.o Base64Transcoder.o MPIdata.o cout0.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testWavefunction: testWavefunction.o Wavefunction.o SlaterDet.o \
        Basis.o BasisMapping.o FourierTransform.o Matrix.o UnitCell.o \
        Context.o Base64Transcoder.o MPIdata.o cout0.o
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
        jade.o MPIdata.o cout0.o
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
