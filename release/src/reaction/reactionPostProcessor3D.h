/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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

/** \file
 * 3D postprocessor to locally perform a generic chemical reactions
 *  -- header file
 */
#ifndef REACTION_POST_PROCESSOR_3D_H
#define REACTION_POST_PROCESSOR_3D_H

#include "rate.h"
#include "reactingSpecies3D.h"

namespace olb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Postprocessor to perform a generic chemical reactions
 */
template<typename T, typename DESCRIPTOR, typename REACTIONS>
class ReactionPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ReactionPostProcessor3D ( int x0, int x1, int y0, int y1, int z0, int z1,
                            std::vector<std::shared_ptr<Rate<T>>> rate, std::shared_ptr<REACTIONS> reactions,
                            std::vector<BlockStructureD<3>*> partners );
  ReactionPostProcessor3D ( std::vector<std::shared_ptr<Rate<T>>> rate, std::shared_ptr<REACTIONS> reactions,
                            std::vector<BlockStructureD<3>*> partners );
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice<T,DESCRIPTOR>& blockLattice,
                          int x0, int x1, int y0, int y1, int z0, int z1 ) override;
private:
  template <typename VECT_TYPE, typename F>
  void functOverReactions(std::vector<VECT_TYPE>& vect, F&& f);
  int _x0, _x1, _y0, _y1, _z0, _z1;
  std::vector<size_t> _sizes;
  std::vector<std::shared_ptr<Rate<T>>> _rate;
  std::shared_ptr<REACTIONS> _reactions;
  std::vector<BlockStructureD<3>*> _partners;
};

template<typename T, typename DESCRIPTOR, typename REACTIONS>
class ReactionGenerator3D final : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  ReactionGenerator3D ( int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                        std::vector<std::shared_ptr<Rate<T>>> rate, REACTIONS&& reactions );
  ReactionGenerator3D ( std::vector<std::shared_ptr<Rate<T>>> rate, REACTIONS&& reactions );
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<3>*> partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  std::vector<std::shared_ptr<Rate<T>>> _rate;
  std::shared_ptr<REACTIONS> _reactions;
};


}  // namespace olb

#endif
