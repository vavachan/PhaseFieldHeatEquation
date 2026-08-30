/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause
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
 * Wrapper functions that simplify the use of MPI, generic template code
 */

#ifndef MPI_MANAGER_AD_HH
#define MPI_MANAGER_AD_HH

#ifdef PARALLEL_MODE_MPI



//#include <adolc/adouble.h>
#include "utilities/aDiff.h"
#include "communication/mpiManager.h"


//using namespace adtl;


namespace olb {

namespace singleton {


template <typename T, unsigned DIM>
void MpiManager::send(util::ADf<T,DIM> *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), (sizeof(util::ADf<T,DIM>)/8)*count, MPI_DOUBLE, dest, tag, comm);
}

template <typename T, unsigned DIM>
void MpiManager::sendInit(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Send_init(static_cast<void*>(buf), (sizeof(util::ADf<T,DIM>)/8)*count, MPI_DOUBLE, dest, tag, comm, request);
  }
}

template <typename T, unsigned DIM>
void MpiManager::receive(util::ADf<T,DIM> *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), (sizeof(util::ADf<T,DIM>)/8)*count, MPI_DOUBLE, source, tag, comm, &status);
}

template <typename T, unsigned DIM>
void MpiManager::recvInit(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Recv_init(static_cast<void*>(buf), (sizeof(util::ADf<T,DIM>)/8)*count, MPI_DOUBLE, dest, tag, comm, request);
  }
}


template <typename T, unsigned DIM>
void MpiManager::iSend
(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), (sizeof(util::ADf<T,DIM>)/8)*count, MPI_DOUBLE, dest, tag, comm, request);
  }
}





template <typename T, unsigned DIM>
void MpiManager::iRecv(util::ADf<T,DIM> *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), (sizeof(util::ADf<T,DIM>)/8)*count, MPI_DOUBLE, source, tag, comm, request);
  }
}

template <typename T, unsigned DIM>
void MpiManager::sendRecv(util::ADf<T,DIM> *sendBuf, util::ADf<T,DIM> *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               (sizeof(util::ADf<T,DIM>)/8)*count,
               MPI_DOUBLE, dest, tag,
               static_cast<void*>(recvBuf),
               (sizeof(util::ADf<T,DIM>)/8)*count,
               MPI_DOUBLE, source, tag, comm, &status);
}




template <typename T, unsigned DIM>
void MpiManager::bCast(util::ADf<T,DIM>* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            (sizeof(util::ADf<T,DIM>)/8)*sendCount, MPI_DOUBLE, root, comm);
}

template <typename T, unsigned DIM>
void MpiManager::bCast(
  BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& sendData,
  int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
//  MPI_Bcast(static_cast<void*>(sendData.getRawData()),
//            sendData.getDataSize(), MPI_DOUBLE, root, comm);
  for (unsigned iD=0; iD < sendData.getSize(); ++iD) {
    MPI_Bcast(static_cast<void*>(sendData.getColumn(iD).data()),
              sendData.getNcells(), MPI_DOUBLE, root, comm);
  }
}

template <typename T, unsigned DIM>
void MpiManager::reduce(util::ADf<T,DIM>& sendVal, util::ADf<T,DIM>& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }

  int sizeADouble = sizeof(util::ADf<T,DIM>)/8-1;

  MPI_Reduce(static_cast<void*>(&sendVal.v()),
             static_cast<void*>(&recvVal.v()), 1, MPI_DOUBLE, op, root, comm);

  for (int i=0; i<sizeADouble; i++) {
    MPI_Reduce(static_cast<void*>(&sendVal.d(i)),
               static_cast<void*>(&recvVal.d(i)), 1, MPI_DOUBLE, op, root, comm);
  }
}

template <typename T, unsigned DIM>
void MpiManager::reduce(
  BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& sendVal,
  BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& recvVal,
  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
//  MPI_Reduce(static_cast<void*>(sendVal.getRawData()),
//             static_cast<void*>(recvVal.getRawData()),
//             sendVal.getDataSize(), MPI_DOUBLE, op, root, comm);
  for (unsigned iD=0; iD < sendVal.getSize(); ++iD) {
    MPI_Reduce(static_cast<void*>(sendVal.getColumn(iD).data()),
               static_cast<void*>(recvVal.getColumn(iD).data()),
               sendVal.getNcells(), MPI_DOUBLE, op, root, comm);
  }
}


template <typename T, unsigned DIM>
void MpiManager::bCastThroughMaster(util::ADf<T,DIM>* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <typename T, unsigned DIM>
void MpiManager::reduceAndBcast(util::ADf<T,DIM>& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  util::ADf<T,DIM> recvVal;
  reduce(reductVal, recvVal, op, root, comm);

  //MPI_Reduce(&reductVal, &recvVal, 1, MPI_DOUBLE, op, root, comm);
  reductVal = recvVal;
  bCast(&reductVal, 1, root, comm);

  //MPI_Bcast(&reductVal, 1, MPI_DOUBLE, root, comm);

}


/*

template <typename T>
void MpiManager::wait(MPI_Request* request, MPI_Status* status, void* ptr, T* writeBack, int writeBackSize)
{
    if (!ok) return;
    MPI_Wait(request, status);

    if (ptr!=NULL) {
        delete [] ptr;
        ptr = NULL;
    }
}
*/




}   // namespace singleton


}  // namespace olb



#endif  // PARALLEL_MODE_MPI

#endif  // MPI_MANAGER_AD_HH
