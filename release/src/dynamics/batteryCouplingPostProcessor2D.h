/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Paul Neugebauer, Lukas Richter, Kevin Schuelein
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

#ifndef BATTERY_COUPLING_POST_PROCESSOR_2D_H
#define BATTERY_COUPLING_POST_PROCESSOR_2D_H

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
// ======== AD coupling with Concentration 2D ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class batteryCouplingPostProcessor2D
 : public LocalPostProcessor2D<T,DESCRIPTOR> {
    private:
  mutable OstreamManager            clout {std::cout, "Reaction2dSolver"};

public:
  batteryCouplingPostProcessor2D(
    int x0_, int x1_, int y0_, int y1_,
    const T conversionDiffusivity_,
    const T diffusion_S_,
    const T kappa_S_,
    const T factor1_,
    const T factor2_,
    const T fac_sep_,
    std::vector<BlockStructureD<2>* > partners_)
   : x0(x0_), x1(x1_), y0(y0_), y1(y1_),
   conversionDiffusivity(conversionDiffusivity_), diffusion_S(diffusion_S_),
   kappa_S(kappa_S_), factor1(factor1_), factor2(factor2_), fac_sep(fac_sep_),
   partners(partners_) {
    this->getName() = "batteryCouplingPostProcessor2D";
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
    processSubDomain(blockLattice, x0, x1, y0, y1);
  }

  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_) override {

    int newX0, newX1, newY0, newY1;
    if ( util::intersect (
           x0, x1, y0, y1,
           x0_, x1_, y0_, y1_,
           newX0, newX1, newY0, newY1 ) ) {

      for (int iX=newX0; iX<=newX1; ++iX) {
        for (int iY=newY0; iY<=newY1; ++iY) {
          T materialNumber = blockLattice.get(iX,iY).template getField<descriptors::CELL_TYPE>();


      //if ((int) materialNumber == 6 ||(int) materialNumber == 7 || (int) materialNumber == 1 || (int) materialNumber == 3)
          //{

            T physDiffusion_c;
            T tmp_source_c;
            T physDiffusion_phi;
            T tmp_source_phi;

            std::vector<T> conc;
            std::vector<T> conc_plusX;
            std::vector<T> conc_minusX;
            std::vector<T> conc_plusY;
            std::vector<T> conc_minusY;

            conc.emplace_back(blockLattice.get(iX,iY).computeRho());
            for (int iter_component = 0; iter_component<component_number-1; ++ iter_component){
              conc.emplace_back(tpartners[iter_component]->get(iX,iY).computeRho());
              }

            T tmp_Diffusion_E = 1.2e-21*util::pow(conc[0],4)- 6.5e-18*util::pow(conc[0],3)+ 1.14e-14*util::pow(conc[0],2)-8.06e-12*conc[0]+2.24e-9;
            T kappa_E = -2.39e-11*util::pow(conc[0],4)+1.21e-7*util::pow(conc[0],3)-2.89e-4*util::pow(conc[0],2)+0.32*conc[0]-2.789;

          if ((int) materialNumber == 6 || (int) materialNumber == 7 || (int) materialNumber == 1 || (int) materialNumber == 3 || (int) materialNumber == 4 || (int) materialNumber == 5){

            if (iX== newX0 ){
              conc_minusX.emplace_back(0.0);
              conc_minusX.emplace_back(0.0);
            }
            else {
              conc_minusX.emplace_back(blockLattice.get(iX-1,iY).computeRho());
              conc_minusX.emplace_back(tpartners[0]->get(iX-1,iY).computeRho());
            }

            if (iX== newX1){
              conc_plusX.emplace_back(0.0);
              conc_plusX.emplace_back(0.0);
            }
            else {
              conc_plusX.emplace_back(blockLattice.get(iX+1,iY).computeRho());
              conc_plusX.emplace_back(tpartners[0]->get(iX+1,iY).computeRho());
            }

            if (iY== newY0){
              conc_minusY.emplace_back(blockLattice.get(iX, newY1).computeRho());
              conc_minusY.emplace_back(tpartners[0]->get(iX, newY1).computeRho());
            }
            else {
              conc_minusY.emplace_back(blockLattice.get(iX,iY-1).computeRho());
              conc_minusY.emplace_back(tpartners[0]->get(iX,iY-1).computeRho());
            }

            if (iY == newY1){
              conc_plusY.emplace_back(blockLattice.get(iX,newY0).computeRho());
              conc_plusY.emplace_back(tpartners[0]->get(iX,newY0).computeRho());
            }
            else {
              conc_plusY.emplace_back(blockLattice.get(iX,iY+1).computeRho());
              conc_plusY.emplace_back(tpartners[0]->get(iX,iY+1).computeRho());
            }

            T kappa_E_minusX = -2.39e-11*util::pow(conc_minusX[0],4)+1.21e-7*util::pow(conc_minusX[0],3)-2.89e-4*util::pow(conc_minusX[0],2)+0.32*conc_minusX[0]-2.789;
            T kappa_E_plusX  = -2.39e-11*util::pow(conc_plusX[0],4)+1.21e-7*util::pow(conc_plusX[0],3)-2.89e-4*util::pow(conc_plusX[0],2)+0.32*conc_plusX[0]-2.789;
            T kappa_E_minusY = -2.39e-11*util::pow(conc_minusY[0],4)+1.21e-7*util::pow(conc_minusY[0],3)-2.89e-4*util::pow(conc_minusY[0],2)+0.32*conc_minusY[0]-2.789;
            T kappa_E_plusY  = -2.39e-11*util::pow(conc_plusY[0],4)+1.21e-7*util::pow(conc_plusY[0],3)-2.89e-4*util::pow(conc_plusY[0],2)+0.32*conc_plusY[0]-2.789;

            if ((int) materialNumber == 6 || (int)materialNumber == 5){ //Solid
                physDiffusion_c = diffusion_S;
                physDiffusion_phi = kappa_S;
                tmp_source_c = 0.0;
                tmp_source_phi= 0.0;
            }
            else if ((int) materialNumber == 7 || (int) materialNumber == 3){ //Separator
                physDiffusion_c = fac_sep*(tmp_Diffusion_E- factor1*factor2*(kappa_E/conc[0]));
                physDiffusion_phi = kappa_E;
                tmp_source_c = fac_sep*(factor1* (1./4.)*((kappa_E_plusX-kappa_E_minusX)*(conc_plusX[1]-conc_minusX[1])+
                                + (kappa_E_plusY-kappa_E_minusY)*(conc_plusY[1]-conc_minusY[1]))
                                -kappa_E*(conc_minusX[1]+conc_plusX[1]+ conc_minusY[1]+conc_plusY[1]-4*conc[1]));
                tmp_source_phi = factor2*((1./4.)*((kappa_E_plusX/conc_plusX[0]-kappa_E_minusX/conc_minusX[0])*(conc_plusX[0]-conc_minusX[0])
                                  + (kappa_E_plusY/conc_plusY[0]-kappa_E_minusY/conc_minusY[0])*(conc_plusY[0]-conc_minusY[0]))
                                  -(kappa_E/conc[0])*(conc_minusX[0]+conc_plusX[0]+conc_minusY[0]+conc_plusY[0]-4*conc[0]));
            }
            else if((int) materialNumber == 1 || (int) materialNumber == 4) { //Electrolyte
              physDiffusion_c = tmp_Diffusion_E- factor1*factor2*(kappa_E/conc[0]);
              physDiffusion_phi = kappa_E;
              tmp_source_c = (factor1* (1./4.)*((kappa_E_plusX-kappa_E_minusX)*(conc_plusX[1]-conc_minusX[1])+
                                + (kappa_E_plusY-kappa_E_minusY)*(conc_plusY[1]-conc_minusY[1]))
                                -kappa_E*(conc_minusX[1]+conc_plusX[1]+ conc_minusY[1]+conc_plusY[1]-4*conc[1]));
              tmp_source_phi = factor2*((1./4.)*((kappa_E_plusX/conc_plusX[0]-kappa_E_minusX/conc_minusX[0])*(conc_plusX[0]-conc_minusX[0])
                                  + (kappa_E_plusY/conc_plusY[0]-kappa_E_minusY/conc_minusY[0])*(conc_plusY[0]-conc_minusY[0]))
                                  -(kappa_E/conc[0])*(conc_minusX[0]+conc_plusX[0]+conc_minusY[0]+conc_plusY[0]-4*conc[0]));

            }


            T tau_c = (physDiffusion_c/ conversionDiffusivity * descriptors::invCs2<T,DESCRIPTOR>())+0.5;
            T tmp_omega_c = 1.0/tau_c;

            T tau_phi = (physDiffusion_phi/ conversionDiffusivity * descriptors::invCs2<T,DESCRIPTOR>())+0.5;
            T tmp_omega_phi = 1.0/tau_phi;


            blockLattice.get(iX,iY).template setField<descriptors::SOURCE>(tmp_source_c*2.5e-13);
            blockLattice.get(iX,iY).template setField<descriptors::OMEGA>(tmp_omega_c);
            blockLattice.get(iX,iY).template setField<collision::TRT::MAGIC>(0.1*physDiffusion_c);

            tpartners[0]->get(iX,iY).template setField<descriptors::SOURCE>(tmp_source_phi*2.5e-13);
            tpartners[0]->get(iX,iY).template setField<descriptors::OMEGA>(tmp_omega_phi);
            tpartners[0]->get(iX,iY).template setField<collision::TRT::MAGIC>(0.1*physDiffusion_phi);

          }

            //Berechnung von BOUNDARY für MN 3 für BlockLattices
            else{
            T boundaryInflow = factor1/tmp_Diffusion_E*(kappa_E*(tpartners[0]->get(iX+1,iY).computeRho()-conc[1])
                              +kappa_S/conc[0]*(blockLattice.get(iX+1,iY).computeRho()-conc[0]));
            blockLattice.get(iX,iY).template setField<descriptors::BOUNDARY> (boundaryInflow);
            }

         //}
        }
      }
    }
  }

private:
  int x0, x1, y0, y1;
  int component_number;
  T conversionDiffusivity;
  T diffusion_S;
  T kappa_S;
  T factor1;
  T factor2;
  T fac_sep;
  std::vector<BlockLattice<T,DESCRIPTOR>*> tpartners;
  std::vector<BlockStructureD<2>* > partners;
};

template<typename T, typename DESCRIPTOR>
class batteryCouplingGenerator2D
 : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {

public:
  batteryCouplingGenerator2D(
    int x0_, int x1_, int y0_, int y1_,
    T conversionDiffusivity_, T diffusion_S_, T kappa_S_, T factor1_, T factor2_, T fac_sep_)
   : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
     conversionDiffusivity(conversionDiffusivity_), diffusion_S(diffusion_S_),
     kappa_S(kappa_S_), factor1(factor1_), factor2(factor2_), fac_sep(fac_sep_)
  { }

  PostProcessor2D<T,DESCRIPTOR>* generate(
    std::vector<BlockStructureD<2>* > partners) const override {
    return new batteryCouplingPostProcessor2D<T,DESCRIPTOR>(
      this->x0,this->x1,this->y0,this->y1, conversionDiffusivity, diffusion_S,
      kappa_S, factor1, factor2, fac_sep, partners);
  }

  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override {
    return new batteryCouplingGenerator2D<T,DESCRIPTOR>(*this);
  }

private:
  T conversionDiffusivity;
  T diffusion_S;
  T kappa_S;
  T factor1;
  T factor2;
  T fac_sep;
};




}

#endif
