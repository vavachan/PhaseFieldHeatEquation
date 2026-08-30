/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle
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

#ifndef ADSORPTION_COUPLING_POST_PROCESSOR_3D_H
#define ADSORPTION_COUPLING_POST_PROCESSOR_3D_H

#include "core/blockStructure.h"
#include "core/postProcessing.h"
#include "core/util.h"
#include "descriptor/descriptor.h"
#include "utilities/omath.h"

#include "adsorptionReaction.h"

namespace olb {

/// Coupling for adsorption reaction
template <typename ADSORPTION_REACTION, typename FORCE>
struct AdsorptionFullCoupling3D {
  using reaction_t = ADSORPTION_REACTION;
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  // Combine the parameters of ADSORPTION_REACTION e FORCE
  using parameters = meta::merge<typename ADSORPTION_REACTION::parameters,
                                  typename FORCE::parameters>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::Concentration0>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Concentration0>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellParticle = cells.template get<names::Concentration0>();
    auto& cellSolute = cells.template get<names::Concentration1>();
    auto& cellLoading = cells.template get<names::Concentration2>();

    // computation of particle velocity for particle and load lattices
    V velocity[3], velocityRight[3], velocityLeft[3], velocityUp[3], velocityDown[3], velocityFront[3], velocityBack[3] {V(0)};
    cellNS.computeU(velocity);
    // x-velocities
    cellNS.neighbor({1,0,0}).computeU(velocityRight);
    cellNS.neighbor({-1,0,0}).computeU(velocityLeft);
    // y-velocities
    cellNS.neighbor({0,1,0}).computeU(velocityUp);
    cellNS.neighbor({0,-1,0}).computeU(velocityDown);
    // z-velocities
    cellNS.neighbor({0,0,1}).computeU(velocityFront);
    cellNS.neighbor({0,0,-1}).computeU(velocityBack);
    // velocity gradients
    V velocityGrad[3] = {0., 0., 0.};
    velocityGrad[0] = 0.5*( velocity[0]*(velocityRight[0] - velocityLeft[0])
                          + velocity[1]*(velocityUp[0] - velocityDown[0])
                          + velocity[2]*(velocityFront[0] - velocityBack[0]) );
    velocityGrad[1] = 0.5*( velocity[0]*(velocityRight[1] - velocityLeft[1])
                          + velocity[1]*(velocityUp[1] - velocityDown[1])
                          + velocity[2]*(velocityFront[1] - velocityBack[1]) );
    velocityGrad[2] = 0.5*( velocity[0]*(velocityRight[2] - velocityLeft[2])
                          + velocity[1]*(velocityUp[2] - velocityDown[2])
                          + velocity[2]*(velocityFront[2] - velocityBack[2]) );

    V forceValue[3] = {0.,0.,0.};

    FORCE().applyForce(forceValue,cells,velocity,parameters);

    // compute new particle velocity under action of forces
    V newVel[3];
    for (int i=0; i < DESCRIPTOR::d; i++) {
      newVel[i] = velocity[i] + forceValue[i] - velocityGrad[i];
    }
    cellParticle.template setField<descriptors::VELOCITY>(newVel);
    cellLoading.template setField<descriptors::VELOCITY>(newVel);

    //setting of the solute velocity to the carrier fluid velocity
    cellSolute.template setField<descriptors::VELOCITY>(velocity);

    //computation of the adsorption source terms for solute and load lattices
    V soluteConcentration = cellSolute.computeRho();
    V particleConcentration = cellParticle.computeRho();
    V particleLoading = cellLoading.computeRho();

    Vector<V, 2> reactionRates = ADSORPTION_REACTION().getReactionRate(soluteConcentration, particleLoading, particleConcentration,parameters);
    cellLoading.template setField<descriptors::SOURCE>(reactionRates[1]);
    cellSolute.template setField<descriptors::SOURCE>(-reactionRates[0]);

   }
};

}

#endif //OLB_ADSORPTIONCOUPLING3D_H
