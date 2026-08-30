/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Johanna Moedl
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef ADVECTION_DIFFUSSION_REACTION_COUPLING_POST_PROCESSOR_3D_H
#define ADVECTION_DIFFUSSION_REACTION_COUPLING_POST_PROCESSOR_3D_H

#include "core/blockStructure.h"
#include "core/postProcessing.h"
#include "core/util.h"
#include "descriptor/descriptor.h"
#include "utilities/omath.h"


namespace olb {


/**
* Coupling of ADlattice[0] with the other AD lattices (tpartners)
*/

//======================================================================
// ======== AD coupling with Concentration 3D ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class ConcentrationAdvectionDiffusionCouplingPostProcessor3D
 : public LocalPostProcessor3D<T,DESCRIPTOR> {

public:
  ConcentrationAdvectionDiffusionCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    const std::vector<T>& stoichiometricCoeff_,
    const std::vector<T> latticeReactionCoeff_,
    const std::vector<T>& react_order_,
    std::vector<BlockStructureD<3>* > partners_)
   : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), stoichiometricCoeff(stoichiometricCoeff_), latticeReactionCoeff(latticeReactionCoeff_),
     react_order(react_order_), partners(partners_) {
    this->getName() = "ConcentrationAdvectionDiffusionCouplingPostProcessor3D";
    reaction_number = static_cast<int>(latticeReactionCoeff.size());
    component_number = static_cast<int>(partners.size())+1;
    for (int i = 0; i<component_number;i++){
      tpartners.emplace_back(
        static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[i]));
    }
  }

  int extent() const override {
    return 0;
  }

  int extent(int whichDirection) const override {
    return 0;
  }

  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override {
    processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
  }

  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)  override {

    int newX0, newX1, newY0, newY1, newZ0, newZ1;
    if ( util::intersect (
           x0, x1, y0, y1, z0, z1,
           x0_, x1_, y0_, y1_, z0_, z1_,
           newX0, newX1, newY0, newY1, newZ0, newZ1) ) {

      for (int iX=newX0; iX<=newX1; ++iX) {
        for (int iY=newY0; iY<=newY1; ++iY) {
          for (int iZ=newZ0; iZ<=newZ1; ++iZ) {

            std::vector<T> conc;
            conc.emplace_back(blockLattice.get(iX,iY, iZ).computeRho());
            for (int iter_component = 0; iter_component<component_number-1; ++ iter_component){
              conc.emplace_back(tpartners[iter_component]->get(iX,iY,iZ).computeRho());
              }

            T lambda[reaction_number];
            T reaction_rate;
            for (int iter_reaction = 0; iter_reaction<reaction_number; ++ iter_reaction){
              lambda[iter_reaction] = 0;
              reaction_rate = 1;
              for(int iter_component = 0; iter_component <component_number; ++ iter_component){

                reaction_rate = reaction_rate*(util::pow(conc[iter_component],react_order[iter_reaction*component_number+iter_component]));
              }

              lambda[iter_reaction] = reaction_rate*latticeReactionCoeff[iter_reaction];
           }
            T temp_source;
            for (int iter_component = 0; iter_component<component_number; ++ iter_component){
              temp_source = 0;
              for (int iter_reaction = 0; iter_reaction<reaction_number; ++ iter_reaction){
                temp_source = temp_source + stoichiometricCoeff[iter_reaction*component_number+iter_component]*lambda[iter_reaction];
              }
              if (iter_component == 0){
                blockLattice.get(iX,iY,iZ).template setField<descriptors::SOURCE>(
                temp_source);
              }
              else {
                tpartners[iter_component-1]->get(iX,iY,iZ).template setField<descriptors::SOURCE>(
                temp_source);
                }
            }
         }
         }
      }
    }
  }

private:
  int x0, x1, y0, y1, z0, z1;
  int reaction_number;
  int component_number;
  const std::vector<T>& stoichiometricCoeff;
  const std::vector<T> latticeReactionCoeff;
  const std::vector<T>& react_order;
  std::vector<BlockLattice<T,DESCRIPTOR>*> tpartners;
  std::vector<BlockStructureD<3>* > partners;
};

template<typename T, typename DESCRIPTOR>
class ConcentrationAdvectionDiffusionCouplingGenerator3D
 : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {

public:
  ConcentrationAdvectionDiffusionCouplingGenerator3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    const std::vector<T>& stoichiometricCoeff_,
    const std::vector<T> latticeReactionCoeff_,
    const std::vector<T>& react_order_)
   : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
     stoichiometricCoeff(stoichiometricCoeff_), latticeReactionCoeff(latticeReactionCoeff_), react_order(react_order_)
  { }

  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<BlockStructureD<3>* > partners) const override {
    return new ConcentrationAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>(
      this->x0,this->x1,this->y0,this->y1, this->z0, this->z1, stoichiometricCoeff, latticeReactionCoeff, react_order, partners);
  }

  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override {
    return new ConcentrationAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>(*this);
  }

private:
  const std::vector<T>& stoichiometricCoeff;
  const std::vector<T> latticeReactionCoeff;
  const std::vector<T>& react_order;
};




}

#endif
