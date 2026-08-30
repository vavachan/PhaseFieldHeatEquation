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
 * 2D postprocessor to locally perform a generic chemical reaction
 *  -- generic implementation
 */
#ifndef REACTION_POST_PROCESSOR_2D_HH
#define REACTION_POST_PROCESSOR_2D_HH

namespace olb {


///////////////////////////////////// class ReactionPostProcessor2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename REACTIONS>
ReactionPostProcessor2D<T,DESCRIPTOR,REACTIONS>::ReactionPostProcessor2D ( int x0, int x1, int y0, int y1,
  std::vector<std::shared_ptr<Rate<T>>> rate, std::shared_ptr<REACTIONS> reactions,
  std::vector<BlockStructureD<2>*> partners )
  : _x0(x0), _x1(x1), _y0(y0), _y1(y1),
    _rate(rate), _reactions(reactions), _partners(partners)
{
  std::size_t iReaction = 0;
  _sizes = std::vector<std::size_t>(std::tuple_size_v<REACTIONS>);
  meta::tuple_for_each(*_reactions.get(), [&](auto& line){
    using line_type = typename std::remove_reference_t<decltype(line)>;
    _sizes[iReaction++] = std::tuple_size_v<line_type>;
  });
  if (std::accumulate(_sizes.begin(), _sizes.end(), std::size_t(0)) != partners.size()) {
    throw std::invalid_argument("The number of species must equate the number of input lattices.");
  }
}

template <typename T, typename DESCRIPTOR, typename REACTIONS>
ReactionPostProcessor2D<T,DESCRIPTOR,REACTIONS>::ReactionPostProcessor2D (
  std::vector<std::shared_ptr<Rate<T>>> rate, std::shared_ptr<REACTIONS> reactions,
  std::vector<BlockStructureD<2>*> partners )
  : ReactionPostProcessor2D<T,DESCRIPTOR,REACTIONS>(0,0,0,0,rate,reactions,partners)
{}

template <typename T, typename DESCRIPTOR, typename REACTIONS>
template <typename VECT_TYPE, typename F>
void ReactionPostProcessor2D<T,DESCRIPTOR,REACTIONS>::functOverReactions(std::vector<VECT_TYPE>& vect, F&& f)
{
  std::size_t iReaction = 0;
  std::vector<size_t> sizes(std::tuple_size_v<REACTIONS>);
  meta::tuple_for_each(*_reactions.get(), [&](auto& line){
    std::size_t iReagent = 0;
    meta::tuple_for_each(line, [&](auto& elem){
      f(elem, vect[std::accumulate(_sizes.begin(), _sizes.begin()+iReaction, 0) + (iReagent++)]);
    });
    iReaction++;
  });
}

template <typename T, typename DESCRIPTOR, typename REACTIONS>
void ReactionPostProcessor2D<T,DESCRIPTOR,REACTIONS>::process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, _x0, _x1, _y0, _y1);
}

template <typename T, typename DESCRIPTOR, typename REACTIONS>
void ReactionPostProcessor2D<T,DESCRIPTOR,REACTIONS>::processSubDomain ( BlockLattice<T,DESCRIPTOR>& blockLattice,
    int x0, int x1, int y0, int y1 )
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( _x0, _x1, _y0, _y1,
                         x0, x1, y0, y1,
                         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {

        // storing the old values of the local field from the partner lattices & resetting the sources
        std::vector<T> fields_asSingleVector;
        functOverReactions(_partners, [&](auto& elem, auto& component){
          fields_asSingleVector.push_back( elem.getField(component, iX, iY) );
          elem.resetSource(component, iX, iY);
        });

        // computing the rates from the old values of the local fields as per rate's own law
        std::size_t iPartner = 0;
        std::vector<std::tuple<T,BlockStructureD<2>*>> ratePartner;
        for (std::size_t iReaction=0; iReaction<_sizes.size(); ++iReaction) {
          auto reagents_begin = fields_asSingleVector.begin() + std::accumulate(_sizes.begin(), _sizes.begin()+iReaction, 0);
          std::vector<T> reagents = { reagents_begin, reagents_begin + _sizes[iReaction] };
          T rate_thisReaction = _rate[iReaction]->compute(reagents);
          for (std::size_t iReagent=0; iReagent<_sizes[iReaction]; ++iReagent) {
            ratePartner.push_back(std::make_tuple(rate_thisReaction, _partners[iPartner++]));
          }
        }

        // updating local field in the partner lattice: source = sum_i rate_i*coeff
        functOverReactions(ratePartner, [&](auto& elem, auto& component){
          elem.incrementSource(std::get<1>(component), std::get<0>(component) * elem.getStoichioCoeff(), iX, iY);
        });
      }
    }
  }
}


///////////////////////////////////// class ReactionGenerator2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename REACTIONS>
ReactionGenerator2D<T,DESCRIPTOR,REACTIONS>::ReactionGenerator2D ( int x0_, int x1_, int y0_, int y1_,
  std::vector<std::shared_ptr<Rate<T>>> rate, REACTIONS&& reactions)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    _rate(rate), _reactions(std::make_shared<REACTIONS>(reactions))
{}

template <typename T, typename DESCRIPTOR, typename REACTIONS>
ReactionGenerator2D<T,DESCRIPTOR,REACTIONS>::ReactionGenerator2D (
  std::vector<std::shared_ptr<Rate<T>>> rate, REACTIONS&& reactions)
  : ReactionGenerator2D<T,DESCRIPTOR,REACTIONS>(0,0,0,0,rate,std::move(reactions))
{}

template<typename T, typename DESCRIPTOR, typename REACTIONS>
PostProcessor2D<T,DESCRIPTOR>* ReactionGenerator2D<T,DESCRIPTOR,REACTIONS>::generate(std::vector<BlockStructureD<2>*> partners) const
{
  return new ReactionPostProcessor2D<T,DESCRIPTOR,REACTIONS>(
           this->x0,this->x1,this->y0,this->y1, _rate, _reactions, partners);
}

template<typename T, typename DESCRIPTOR, typename REACTIONS>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* ReactionGenerator2D<T,DESCRIPTOR,REACTIONS>::clone() const
{
  return new ReactionGenerator2D<T,DESCRIPTOR,REACTIONS>(*this);
}

}  // namespace olb

#endif
