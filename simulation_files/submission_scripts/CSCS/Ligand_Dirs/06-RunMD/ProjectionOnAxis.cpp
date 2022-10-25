/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Angle.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR PROJECTION_ON_AXIS
/*
Calculate the projection on an axis.

\par Examples

*/
//+ENDPLUMEDOC

class ProjectionOnAxis : public Colvar {
  bool pbc;

public:
  explicit ProjectionOnAxis(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(ProjectionOnAxis,"PROJECTION_ON_AXIS")

void ProjectionOnAxis::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","AXIS_ATOMS","the atoms that define the direction of the axis of interest");
  keys.add("atoms","ATOM","the atom whose position we want to project on the axis of interest");
}

ProjectionOnAxis::ProjectionOnAxis(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true)
{
  vector<AtomNumber> axis_atoms;
  parseAtomList("AXIS_ATOMS",axis_atoms);
  if( axis_atoms.size()!=2 ) error("should only be two atoms specified to AXIS_ATOMS keyword");
  vector<AtomNumber> atom; 
  parseAtomList("ATOM",atom);
  if( atom.size()!=1 ) error("should only be one atom specified to ATOM keyword");
  log.printf("  calculating projection of vector connecting atom %d and atom %d on vector connecting atom %d and atom %d \n",
                axis_atoms[0].serial(), atom[0].serial(), axis_atoms[0].serial(), axis_atoms[1].serial() ); 
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  not using periodic boundary conditions\n");

  // Add values to store data
  addComponentWithDerivatives("proj"); componentIsNotPeriodic("proj"); 
  addComponentWithDerivatives("ext"); componentIsNotPeriodic("ext");
  // Get all the atom positions 
  axis_atoms.push_back( atom[0] ); 
  requestAtoms(axis_atoms);
  checkRead();
}

// calculator
void ProjectionOnAxis::calculate() {

  Vector rik, rjk;
  if( pbc ) {
      rik = pbcDistance( getPosition(2), getPosition(0) );
      rjk = pbcDistance( getPosition(2), getPosition(1) );
  } else {
      rik = delta( getPosition(2), getPosition(0) );
      rjk = delta( getPosition(2), getPosition(1) );
  }
  Vector rij = delta( rik, rjk ); double dij = rij.modulo(); 
  Vector nij = (1.0/dij)*rij; Tensor dij_a1;
  // Derivative of director connecting atom1 - atom2 wrt the position of atom 1
  dij_a1(0,0) = ( -(nij[1]*nij[1]+nij[2]*nij[2])/dij );   // dx/dx
  dij_a1(0,1) = (  nij[0]*nij[1]/dij );                   // dx/dy
  dij_a1(0,2) = (  nij[0]*nij[2]/dij );                   // dx/dz
  dij_a1(1,0) = (  nij[1]*nij[0]/dij );                   // dy/dx
  dij_a1(1,1) = ( -(nij[0]*nij[0]+nij[2]*nij[2])/dij );   // dy/dy
  dij_a1(1,2) = (  nij[1]*nij[2]/dij );
  dij_a1(2,0) = (  nij[2]*nij[0]/dij );
  dij_a1(2,1) = (  nij[2]*nij[1]/dij );
  dij_a1(2,2) = ( -(nij[1]*nij[1]+nij[0]*nij[0])/dij );

  // Calculate dot product and derivatives
  double d = dotProduct( -rik, nij );
  Vector dd1 = matmul(-rik, dij_a1) - nij;
  Vector dd2 = matmul(rik, dij_a1);
  Vector dd3 = nij;
  Value* pval=getPntrToComponent("proj"); pval->set( d );
  setAtomsDerivatives( pval, 0, dd1 );
  setAtomsDerivatives( pval, 1, dd2 );
  setAtomsDerivatives( pval, 2, dd3 );
  setBoxDerivatives( pval, -Tensor( rik, dd1 ) - Tensor( rjk, dd2 ) );
  // Calculate derivatives of perpendicular distance from axis
  double c = sqrt( rik.modulo2() - d*d ); double invc = (1.0/c);
  // Calculate derivatives of the other thing
  Vector der1 = invc*(rik - d*dd1);
  Vector der2 = invc*(-d*dd2);
  Vector der3 = invc*(-rik - d*dd3);

  Value* cval=getPntrToComponent("ext"); cval->set( c ); 
  setAtomsDerivatives( cval, 0, der1 );
  setAtomsDerivatives( cval, 1, der2 );
  setAtomsDerivatives( cval, 2, der3 );
  setBoxDerivatives( cval, -Tensor( rik, der1 ) - Tensor( rjk, der2 ) );
}

}
}



