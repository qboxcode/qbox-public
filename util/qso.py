#!/usr/bin/env python3
import xml.sax
import numpy as np
# quantum-simulation.org (QSO) definitions
class Atom:
  def __init__(self,name,species,pos,vel):
    self.name = name
    self.species = species
    self.position = pos
    self.velocity = vel

class Species:
  def __init__(self,name,href,symbol,atomic_number,mass):
    self.name = name
    self.href = href
    self.symbol = symbol
    self.atomic_number = atomic_number
    self.mass = mass

class UnitCell:
  def __init__(self,ax=0.0,ay=0.0,az=0.0,bx=0.0,by=0.0,bz=0.0,cx=0.0,cy=0.0,cz=0.0):
    self.a = [ax,ay,az]
    self.b = [bx,by,bz]
    self.c = [cx,cy,cz]
    self.update()

  def update(self):
    self.av = np.array(self.a)
    self.bv = np.array(self.b)
    self.cv = np.array(self.c)
    self.an = [np.zeros(3) for i in range(13)]
    self.an[0] = self.av
    self.an[1] = self.bv
    self.an[2] = self.cv
    self.an[3] = self.av + self.bv
    self.an[4] = self.av - self.bv
    self.an[5] = self.bv + self.cv
    self.an[6] = self.bv - self.cv
    self.an[7] = self.cv + self.av
    self.an[8] = self.cv - self.av
    self.an[9] = self.av + self.bv + self.cv
    self.an[10] = self.av - self.bv - self.cv
    self.an[11] = self.av + self.bv - self.cv
    self.an[12] = self.av - self.bv + self.cv
    self.a_norm = np.array([np.linalg.norm(self.an[i]) for i in range(13)])
    self.an2h = np.array([0.5*self.a_norm[i]*self.a_norm[i] for i in range(13)])

  def fold_in_ws(self,x,y,z):
    v = np.array([x,y,z])
    epsilon=1.e-10
    done = False
    maxiter = 10
    iter = 0
    while ( not done and iter < maxiter ):
      done = True
      for i in range(13):
        sp = np.dot(v,self.an[i])
        if ( sp > self.an2h[i] + epsilon ):
          done = False
          v -= self.an[i]
          while ( np.dot(v,self.an[i]) > self.an2h[i] + epsilon ):
            v -= self.an[i]
        elif ( sp < -self.an2h[i] - epsilon ):
          done = False
          v += self.an[i]
          while ( np.dot(v,self.an[i]) < -self.an2h[i] - epsilon ):
            v += self.an[i]
      iter += 1
    assert( iter < maxiter)
    return (v[0],v[1],v[2])

class AtomSet:
  def __init__(self):
    self.atom_list = []
    self.species_list = []
    self.cell = UnitCell()

  def add_atom(self,atom):
    self.atom_list.append(atom)

  def add_species(self,species):
    # add species to the species list, only if name is not a duplicate
    if not self.find_species(species.name):
      self.species_list.append(species)

  def find_species(self,species_name):
    found = False
    for sp in range(len(self.species_list)):
      found |= (self.species_list[sp].name == species_name)
    return found

class Sample:
  def __init__(self):
    self.atoms = AtomSet()

# The following handler processes the <atomset> element of
# an XML document and updates the AtomSet data of the Sample
# If multiple instances of <atomset> are found, the
# handler overwrites the AtomSet data
class QSOAtomSetHandler(xml.sax.handler.ContentHandler):
  def __init__(self,sample):
    self.s = sample
    self.inAtomSet = False
    self.inAtom = False
    self.inPosition = False
    self.inVelocity = False
    self.inSpecies = False
    self.inSymbol = False
    self.inAtomicNumber = False
    self.inMass = False
    self.buffer = ""
    # flag to signal that the first <atomset> has been processed
    self.done_first = False

  def startElement(self, name, attributes):
    if name == "atomset":
      self.s.atoms.atom_list = []
      self.inAtomSet = True
    elif (name == "unit_cell") and self.inAtomSet:
      self.s.atoms.cell.a = [float(s) for s in attributes["a"].split()]
      self.s.atoms.cell.b = [float(s) for s in attributes["b"].split()]
      self.s.atoms.cell.c = [float(s) for s in attributes["c"].split()]
      self.cell.update()
    elif (name == "species"):
      self.inSpecies = True
      self.species_name = "species_name"
      if "name" in attributes:
        self.species_name = attributes["name"]
      self.species_href = self.species_name+"_href"
      if "href" in attributes:
        self.species_href = attributes["href"]
      self.species_symbol = self.species_name+"_symbol"
      self.species_atomic_number = self.species_name+"_atomic_number"
      self.species_mass = self.species_name+"_mass"
    elif (name == "atom") and self.inAtomSet:
      self.inAtom = True
      self.atom_name = attributes["name"]
      self.atom_species = attributes["species"]
      sp = Species(self.species_name,self.species_href,self.species_symbol,
           self.species_atomic_number,self.species_mass)
      self.s.atoms.add_species(sp)
      self.atom_position = []
      self.atom_velocity = []
    elif (name == "position") and self.inAtom:
      self.buffer = ""
      self.inPosition = True
    elif (name == "velocity") and self.inAtom:
      self.buffer = ""
      self.inVelocity = True
    elif (name == "symbol") and self.inSpecies:
      self.buffer = ""
      self.inSymbol = True
    elif (name == "atomic_number") and self.inSpecies:
      self.buffer = ""
      self.inAtomicNumber = True
    elif (name == "mass") and self.inSpecies:
      self.buffer = ""
      self.inMass = True

  def characters(self, data):
    if self.inPosition or self.inVelocity or self.inSymbol or self.inAtomicNumber or self.inMass:
      self.buffer += data

  def endElement(self, name):
    if (name == "atom") and self.inAtomSet:
      a = Atom(self.atom_name,self.atom_species,self.atom_position,
      self.atom_velocity)
      self.s.atoms.add_atom(a)
      self.inAtom = False
    if (name == "position") and self.inAtom:
      pos = self.buffer.split()
      x = float(pos[0])
      y = float(pos[1])
      z = float(pos[2])
      self.atom_position = [x,y,z]
      self.inPosition = False
    if (name == "velocity") and self.inAtom:
      vel = self.buffer.split()
      vx = float(vel[0])
      vy = float(vel[1])
      vz = float(vel[2])
      self.atom_velocity = [vx,vy,vz]
      self.inVelocity = False
    elif name == "atomset":
      self.done_first = True
      self.inAtomSet = False
    elif name == "species":
      sp = Species(self.species_name,self.species_href,self.species_symbol,
           self.species_atomic_number,self.species_mass)
      self.s.atoms.add_species(sp)
      self.inSpecies = False
    elif (name == "symbol") and self.inSpecies:
      self.species_symbol = self.buffer
      self.inSymbol = False
    elif name == ("atomic_number") and self.inSpecies:
      self.species_atomic_number = int(self.buffer)
      self.inAtomicNumber = False
    elif (name == "mass") and self.inSpecies:
      self.species_mass = float(self.buffer)
      self.inMass = False

