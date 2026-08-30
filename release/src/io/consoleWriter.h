/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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



#ifndef CONSOLE_WRITER_H
#define CONSOLE_WRITER_H

namespace olb {


template<typename T, unsigned D, typename DATA>
class ConsoleWriter{
  using loc = int;
  using Pos = Vector<loc,D>; //Has dimension of source data
private:
  DATA& _data;
  Pos _extent = Pos(0);
  std::ostream& os;
  const char* _charEmpty=" ";
  Vector<unsigned,2> _spacing = {2,0};

  //Make sequence by using multiplier n
  std::string sequence(std::string a, unsigned n){
    std::string seq="";
    for (unsigned i=0; i<n; ++i){
      seq+=a;
    }
    return seq;
  }

public:

  //Construct from DATA and ostream
  ConsoleWriter(std::ostream& str, DATA& data)
    : _data(data), os(str)
  {
    _extent = data.getExtent();
  }

  //Print Writer info
  void print(){
    os << "DataType: " << typeid(DATA).name() << std::endl;
    os << "Extent=" << _extent << std::endl;
  }

  void adjustSettings(Vector<unsigned,2> spacing)
  {
    _spacing = spacing;
  }

  //Write method with function f
  template<typename F>
  void write(F f, unsigned sliceAxis=2, int slicePosOnAxis=0){
    for (loc iY=0; iY < _extent[1]; ++iY) {
      for (loc iX=0; iX < _extent[0]; ++iX) {
        if constexpr (D == 3) {
          int iZ = slicePosOnAxis;
          loc iXeval;
          loc iYeval;
          loc iZeval;
          if (sliceAxis==0){            //YZ-Slice
            iXeval = iZ;
            iYeval = iX;
            iZeval = _extent[1]-iY-1;
          } else if (sliceAxis==1){     //XZ-Slice
            iXeval = iX;
            iYeval = iZ;
            iZeval = _extent[1]-iY-1;
          } else {                      //XY-Slice
            iXeval = iX;
            iYeval = _extent[1]-iY-1;
            iZeval = iZ;
          }
          //Evaluate function f
          if constexpr (std::is_invocable_v<F, DATA&, LatticeR<D>>) {
            os << f(_data,{iXeval,iYeval,iZeval}); //Version B: call with braced enclosed dimensions
          } else {
            os << f(_data,iXeval,iYeval,iZeval);   //Version A: call with dimensions separat
          }
        } else {
          //Project coordinates
          loc iXeval = iX;
          loc iYeval = _extent[1]-iY-1;
          //Evaluate function f
          if constexpr (std::is_invocable_v<F, DATA&, LatticeR<D>>) {
            os << f(_data,{iXeval,iYeval}); //Version B: call with braced enclosed dimensions
          } else {
            os << f(_data,iXeval,iYeval);   //Version A: call with dimensions separat
          }
        }
        os << (iX==_extent[0]-1? "" : sequence(_charEmpty,_spacing[0]));
      }
      os << std::endl;
      //Add _spacing (if needed)
      if (iY<_extent[1]-1){
        for (unsigned i=0; i<_spacing[1]; ++i){
          for (loc iX=0; iX < _extent[0]; ++iX) {
            os << (iX==_extent[0]-1? _charEmpty : sequence(_charEmpty,_spacing[0]+1));
          }
          os << std::endl;
        }
      }
    }
  }

  /// Non-generic common functions

  //Write method for separate GETABLE (e.g. SuperGeometry)
  template<typename GETABLE>
  void writeGetable(GETABLE& getable, unsigned sliceAxis=2, int slicePosOnAxis=0){
    write([&](auto& data, LatticeR<D> pos){
      auto res = getable.get(pos);
      return res;
    });
  }

  //Write method for separate indicator (e.g. BlockIndicator)
  template<typename INDICATOR>
  void writeIndicator(INDICATOR& indicator, unsigned sliceAxis=2, int slicePosOnAxis=0){
    write([&](auto& data, LatticeR<D> pos){
      bool out[1];
      indicator( out, pos.data() );
      return out[0];
    });
  }

  //Write method for separate phys indicator (e.g. IndicatorCircle)
  template<typename INDICATOR>
  void writeIndicator(BlockGeometry<T,D>& blockGeometry, INDICATOR& indicator, unsigned sliceAxis=2, int slicePosOnAxis=0){
    write([&](auto& data, LatticeR<D> pos){
      bool out[1];
      Vector<T,D> physPos;
      blockGeometry.getPhysR( physPos.data(), pos );
      indicator( out, physPos.data() );
      return out[0];
    });
  }

  //Write method for separate functor (e.g. BlockLatticeVelocity)
  template<unsigned DIM, typename FUNCTOR, typename S=T>
  void writeFunctor(FUNCTOR& functor, unsigned sliceAxis=2, int slicePosOnAxis=0){
    write([&](auto& data, LatticeR<D> pos){
      Vector<S,DIM> out;
      functor( out.data(), pos.data() );
      if constexpr (DIM==1){
        return out[0];
      } else {
        return out;
      }
    });
  }

};


} // namespace olb

#endif
