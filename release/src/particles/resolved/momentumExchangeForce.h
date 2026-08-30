/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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



#ifndef MOMENTUM_EXCHANGE_FORCE_H
#define MOMENTUM_EXCHANGE_FORCE_H
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

namespace particles {

namespace resolved {


template<typename T, unsigned D, bool useLadd=false>
struct population_momentum_exchange;
template<typename T>
//Momentum Exchange Wen2014
struct population_momentum_exchange<T,2,false> {
  static constexpr T calculate( T f1, T f2, T c, T pVel, T deltaR )
  {
    return (f1*(c-pVel) + f2*(c+pVel))*(1./deltaR) ;
  }
};
template<typename T>
struct population_momentum_exchange<T,3,false> {
  static constexpr T calculate( T f1, T f2, T c, T pVel, T deltaR )
  {
    return (f1*(c-pVel) + f2*(c+pVel));
  }
};
//Momentum Exchange Ladd1994
template<typename T>
struct population_momentum_exchange<T,2,true> {
  static constexpr T calculate( T f1, T f2, T c, T pVel, T deltaR )
  {
    return ((f1 + f2)*c)*(1./deltaR) ;
  }
};
template<typename T>
struct population_momentum_exchange<T,3,true> {
  static constexpr T calculate( T f1, T f2, T c, T pVel, T deltaR )
  {
    return ((f1 + f2)*c);
  }
};



//Evaluate local momentum exchange contribution and return whether particle and bulk found
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool momentumExchangeAtSurfaceLocation(
  T momentumExchange[],
  PhysR<T,DESCRIPTOR::d>& lever,
  const LatticeR<DESCRIPTOR::d>& latticeRinner,
  const BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  Particle<T,PARTICLETYPE>& particle,
  const int bulkMaterial = 1
)
{

  using namespace descriptors;
  constexpr unsigned D = DESCRIPTOR::d;
  auto position = particles::access::getPosition(particle);

  //Retrieve grid spacing
  const T deltaR = blockGeometry.getDeltaR();
  //Get inner phys position
  T physRinner[D] = { };
  Vector<T,D> physRV;
  blockGeometry.getPhysR(physRV, latticeRinner);
  for (unsigned iD = 0; iD < D; ++iD) {
    physRinner[iD] = physRV[iD];
  }
  //Retrieve inner porosity
  T porosityInner[1] = { };
  particles::resolved::evalSolidVolumeFraction(porosityInner, physRinner, particle);
  //Check whether particle and bulk is existent at location
  if ( !util::nearZero(porosityInner[0]) && blockGeometry.get(latticeRinner)==bulkMaterial ) {
    //Loop over distribution functions
    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      //Calculate outer lattice position (get next cell located in current direction)
      const Vector<int,D> c = descriptors::c<DESCRIPTOR>(iPop);
      const LatticeR<D> latticeRouter = latticeRinner + c;
      //Retrieve outer phys position
      T physRouter[D] = {0.};
      Vector<T,D> physRouterV;
      blockGeometry.getPhysR(physRouterV, latticeRouter);
      for (unsigned iD = 0; iD < D; ++iD) {
        physRouter[iD] = physRouterV[iD];
      }
      //Retrieve outer porosity
      T porosityOuter[1] = {0.};
      particles::resolved::evalSolidVolumeFraction(porosityOuter, physRouter, particle);
      //if not both cells are in the full solid domain calculate force
      if ( !(porosityInner[0]==1 && porosityOuter[0]==1) ) {
        //Momentum Exchange Wen
        const T f1 = blockLattice.get( latticeRouter )[iPop];
        const T f2 = blockLattice.get( latticeRinner )[descriptors::opposite<DESCRIPTOR>(iPop)];
        T pVel[D] = {0.};
        for (unsigned iDim=0; iDim<D; ++iDim) {
          pVel[iDim] = converter.getLatticeVelocity(particles::dynamics::calculateLocalVelocity(particle, PhysR<T,D>(physRinner))[iDim]);
          momentumExchange[iDim] -= converter.getPhysForce(
                                      population_momentum_exchange<T,D>::calculate( f1, f2, c[iDim], pVel[iDim], deltaR ) );
        }

        // Optional correction of the lubrication force
        if constexpr (particles::access::providesIsInContact<PARTICLETYPE>()) {
          auto contact = particle.template getField<NUMERICPROPERTIES,IS_IN_CONTACT>();

          if (contact == 1) {
            T vn = 0.0;
            for (unsigned iDim=0; iDim<D; ++iDim) {
                //std::cout << momentumExchange[iDim];
                vn = vn + util::pow(momentumExchange[iDim],2);
            }

            vn = util::sqrt(vn);

            auto cell = blockLattice.get(latticeRinner);
            Vector<T,DESCRIPTOR::d> localLatticeFluidVelocity;
            Vector<T,DESCRIPTOR::d> localPhysicalFluidVelocity;
            cell.computeU(localLatticeFluidVelocity.data());

            for(unsigned iD=0; iD<D; ++iD) {
                localPhysicalFluidVelocity[iD] = converter.getPhysVelocity(localLatticeFluidVelocity[iD]);
            }

            T U = 0.0;
            for (unsigned iDim=0; iDim<D; ++iDim) {
                U = U + util::pow(pVel[iDim]-localLatticeFluidVelocity[iDim],2.0);
            }
            U = util::sqrt(U);


            T deltaT = converter.getPhysDeltaT();
            T deltaTLattice = converter.getLatticeTime(deltaT);
            T rho_fluid = converter.getPhysDensity();

            T factor = 1.0;

            T maxForce = converter.getPhysForce(U * converter.getLatticeDensity(rho_fluid)/deltaTLattice)*factor;

            if (util::fabs(vn) == 0){
                vn = 0.5 * maxForce;
            }
            if (util::fabs(vn) > util::fabs(maxForce)) {
                T mp = maxForce/vn;
                for (unsigned iDim=0; iDim<D; ++iDim) {
                    //Set the force to the maximal force
                    momentumExchange[iDim] *= mp;
                }
            }
          }
        }
      }

    }
    //Calculate lever (necessary for torque calculation)
    lever = PhysR<T,D>(physRinner)-position;
    return true;
  }
  else {
    return false;
  }
}






} //namespace resolved

} //namespace particles

} //namespace olb

#endif
