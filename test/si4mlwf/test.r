<?xml version="1.0" encoding="UTF-8"?>
<fpmd:simulation xmlns:fpmd="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0">

                   ============================
                   I qbox wrk                 I
                   I                          I
                   I                          I
                   I                          I
                   I                          I
                   I                          I
                   I                          I
                   I                          I
                   I                          I
                   I                          I
                   I                          I
                   I                          I
                   I http://eslab.ucdavis.edu I
                   ============================


<release> wrk pavane </release>
<svn_version> 1174 </svn_version>
<user> fgygi </user>
<sysname> Linux </sysname>
<nodename> pavane.ucdavis.edu </nodename>
<start_time> 2012-04-25T17:38:09Z </start_time>
<mpi_processes count="4">
<process id="0"> pavane.ucdavis.edu </process>
<process id="1"> pavane.ucdavis.edu </process>
<process id="2"> pavane.ucdavis.edu </process>
<process id="3"> pavane.ucdavis.edu </process>
</mpi_processes>
[qbox] <cmd># Si4 MLWF</cmd>
[qbox] <cmd>load ../si4gs/test.xml</cmd>
 LoadCmd: loading from ../si4gs/test.xml
 SampleReader: reading from file ../si4gs/test.xml size: 239854
 XMLGFPreprocessor: read time: 0.0001121
 XMLGFPreprocessor: local read rate: 510.3 MB/s  aggregate read rate: 2041 MB/s
 XMLGFPreprocessor: tag fixing time: 0.0003521
 XMLGFPreprocessor: segment definition time: 0.004858
 XMLGFPreprocessor: boundary adjustment time: 0.001927
 XMLGFPreprocessor: transcoding time: 0.0001631
 XMLGFPreprocessor: data redistribution time: 0.001187
 XMLGFPreprocessor: XML compacting time: 0.0007529
 XMLGFPreprocessor: total time: 0.009787
 xmlcontent.size(): 2406
 Starting XML parsing
 SpeciesHandler: found href in species definition
 name=silicon href=Si_VBC_LDA-1.0.xml

 species silicon:
<species name="silicon" href="Si_VBC_LDA-1.0.xml">
 <description>
Translated from UPF format by upf2qso
Generated using unknown code
Author: Von Barth-Car ( 1984)
Info: automatically converted from PWSCF format
    0        The Pseudo was generated with a Non-Relativistic Calculation
  0.00000000000E+00    Local Potential cutoff radius
nl pn  l   occ               Rcut            Rcut US             E pseu
3S  0  0  2.00      0.00000000000      0.00000000000      0.00000000000
3P  0  1  2.00      0.00000000000      0.00000000000      0.00000000000
SLA PZ NOGX NOGC
</description>
 <symbol>Si</symbol>
 <atomic_number>14</atomic_number>
 <mass>28.09</mass>
 <norm_conserving_pseudopotential>
 <valence_charge>4</valence_charge>
 <lmax>2</lmax>
 <llocal>2</llocal>
 <nquad>0</nquad>
 <rquad>0</rquad>
 <mesh_spacing>0.01</mesh_spacing>
 </norm_conserving_pseudopotential>
</species>
 Kleinman-Bylander potential
 rcps_ =   1.500000
 WavefunctionHandler::startElement: wavefunction nspin=1 nel=16 nempty=0
 kpoint=0.000000 0.000000 0.000000 weight=1.000000 size=8
 XML parsing done
 wavefunction was read
 SampleReader::readSample: grid nx,ny,nz=14 14 14
 SampleReader: read time: 0.058581 s
[qbox] <cmd>compute_mlwf</cmd>
 <mlwf_set size="8">
   <mlwf center="   -1.642349   -1.956538    0.000125 " spread=" 2.272888 "/>
   <mlwf center="   -1.642344    1.956555    0.000119 " spread=" 2.272892 "/>
   <mlwf center="   -4.661109   -0.000001   -0.000015 " spread=" 2.248696 "/>
   <mlwf center="    0.000132   -0.000006    0.000225 " spread=" 2.752357 "/>
   <mlwf center="    4.661105    0.000001   -0.000017 " spread=" 2.248691 "/>
   <mlwf center="    1.642252   -1.956737    0.000173 " spread=" 2.272917 "/>
   <mlwf center="    1.642258    1.956729    0.000158 " spread=" 2.272918 "/>
   <mlwf center="    0.000003    0.000000   -0.001036 " spread=" 3.100348 "/>
 </mlwf_set>
 <electronic_dipole> 0.000104 -0.000005 0.000535 </electronic_dipole>
 <ionic_dipole> 0.000000 0.000000 0.000000 </ionic_dipole>
 <total_dipole> 0.000104 -0.000005 0.000535 </total_dipole>
 <total_dipole_length> 0.000545 </total_dipole_length>
[qbox]  End of command stream 
<real_time> 0.808292 </real_time>
<end_time> 2012-04-25T17:38:09Z </end_time>
</fpmd:simulation>
