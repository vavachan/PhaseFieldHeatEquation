/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef REFINEMENT_COARSE_CELL_H
#define REFINEMENT_COARSE_CELL_H

namespace olb {

namespace refinement {

template <typename CELL>
class CoarseCell {
private:
  CELL _cell;
  LatticeR<CELL::descriptor_t::d> _fineLatticeR;

public:
  using value_t = typename CELL::value_t;
  using descriptor_t = typename CELL::descriptor_t;

  CoarseCell(CELL&& cell, const LatticeR<CELL::descriptor_t::d>& fineLatticeR) any_platform
  : _cell{std::move(cell)}
  , _fineLatticeR{fineLatticeR}
  {
    _cell.setLatticeR(_fineLatticeR / 2);
  }

  CoarseCell(const CELL& cell, const LatticeR<CELL::descriptor_t::d>& fineLatticeR) any_platform
  : _cell{cell}
  , _fineLatticeR{fineLatticeR}
  {
    _cell.setLatticeR(_fineLatticeR / 2);
  }

  CoarseCell(const CoarseCell& rhs) any_platform
  : _cell(rhs._cell)
  , _fineLatticeR(rhs._fineLatticeR)
  { }

  /// Returns true iff present fine location is co-incident with a coarse cell
  operator bool() const any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    bool coIncident = true;
    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      coIncident &= !(_fineLatticeR[iD] & 1);
    }
    return coIncident;
  }

  /// Return reference to coarse cell (only well defined if isCoIncident)
  CELL& operator*() any_platform {
    return _cell;
  }

  /// Return neighbor along fine offset
  CoarseCell neighbor(LatticeR<descriptor_t::d> offset) any_platform {
    return CoarseCell{_cell, _fineLatticeR + offset};
  }

};

template <typename DATA>
class ContextData {
private:
  DATA _data;

public:
  using value_t = typename DATA::value_t;
  using descriptor_t = typename DATA::descriptor_t;

  ContextData(DATA&& data) any_platform
  : _data{std::move(data)}
  { }

  ContextData(const DATA& data) any_platform
  : _data{data}
  { }

  ContextData(const ContextData& rhs) any_platform
  : _data(rhs._data)
  { }

  DATA& operator*() any_platform {
    return _data;
  }

  DATA* operator->() any_platform {
    return &_data;
  }

};

template <typename DATA>
class ContextDataWithNeighbors {
private:
  DATA _data;

public:
  using value_t = typename DATA::value_t;
  using descriptor_t = typename DATA::descriptor_t;

  ContextDataWithNeighbors(DATA&& data) any_platform
  : _data{std::move(data)}
  { }

  ContextDataWithNeighbors(const DATA& data) any_platform
  : _data{data}
  { }

  ContextDataWithNeighbors(const ContextDataWithNeighbors& rhs) any_platform
  : _data(rhs._data)
  { }

  operator bool() const any_platform {
    return _data.getCellId() != std::numeric_limits<CellID>::max();
  }

  DATA& operator*() any_platform {
    return _data;
  }

  DATA* operator->() any_platform {
    return &_data;
  }

  ContextDataWithNeighbors neighbor(unsigned iN) any_platform {
    auto data = _data;
    const CellID index = data.template getFieldComponent<fields::refinement::CONTEXT_NEIGHBORS>(iN);
    data.setCellId(index);
    return ContextDataWithNeighbors{data};
  }

};

}

}

#endif
