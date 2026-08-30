/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender, Mathias J. Krause
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

#ifndef BLOCK_COMMUNICATOR_H
#define BLOCK_COMMUNICATOR_H

#include <map>
#include <set>
#include <typeindex>

#include "core/platform/platform.h"

#include "communication/loadBalancer.h"
#include "communication/blockCommunicationNeighborhood.h"

namespace olb {

template <typename T>
class SuperCommunicationTagCoordinator;

/// Generic communicator for the overlap neighborhood of a block
/**
 * Managed by SuperCommunicator
 *
 * Reconstructed during SuperCommunicator::exchangeRequests
 **/
struct BlockCommunicator {
  virtual ~BlockCommunicator() { };

#ifdef PARALLEL_MODE_MPI
  virtual void receive() = 0;
  virtual void send()    = 0;
  virtual void unpack()  = 0;
  virtual void wait()    = 0;
#else
  virtual void copy()    = 0;
#endif

};

template <typename BLOCK>
class ConcreteBlockCommunicator final : public BlockCommunicator {
private:
  const int _iC;
#ifdef PARALLEL_MODE_MPI
  MPI_Comm _mpiCommunicator;
#endif

#ifdef PARALLEL_MODE_MPI
  class SendTask;
  class RecvTask;

  std::vector<SendTask> _sendTasks;
  std::vector<RecvTask> _recvTasks;
#else
  class CopyTask;

  std::vector<CopyTask> _copyTasks;
#endif

public:
  template <typename T, typename SUPER>
  ConcreteBlockCommunicator(SUPER& super,
                            LoadBalancer<T>& loadBalancer,
#ifdef PARALLEL_MODE_MPI
                            SuperCommunicationTagCoordinator<T>& tagCoordinator,
                            MPI_Comm comm,
#endif
                            int iC,
                            const BlockCommunicationNeighborhood<T,SUPER::d>& neighborhood);

#ifdef PARALLEL_MODE_MPI
  void receive() override;
  void send() override;
  void unpack() override;
  void wait() override;
#else
  void copy() override;
#endif

};


}

#endif
