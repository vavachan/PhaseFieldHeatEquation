/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 The OpenLB project
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
 * Wrapper functions that simplify the use of MPI
 */

#ifndef MPI_MANAGER_H
#define MPI_MANAGER_H

#include "core/platform/platform.h"

#ifdef PARALLEL_MODE_MPI
#include "mpi.h"

#include <vector>
#include <memory>
#endif

#include <string>
#include "io/ostreamManager.h"
#include "utilities/aDiff.h"

namespace olb {

template <unsigned D, typename T, typename U>
class BlockData;

namespace singleton {

#ifdef PARALLEL_MODE_MPI

/// Helper class for non blocking MPI communication
class MpiNonBlockingHelper {
private:
  /// Size of the vector _mpiRequest/_mpiStatus
  unsigned _size;
  /// vector of MPI_Request
  std::unique_ptr<MPI_Request[]> _mpiRequest;
  /// vector of MPI_Status
  std::unique_ptr<MPI_Status[]> _mpiStatus;
public:
  MpiNonBlockingHelper();
  ~MpiNonBlockingHelper() = default;

  MpiNonBlockingHelper(MpiNonBlockingHelper&& rhs) = default;
  MpiNonBlockingHelper(const MpiNonBlockingHelper&) = delete;
  MpiNonBlockingHelper& operator=(const MpiNonBlockingHelper&) = delete;

  /// Allocates memory
  void allocate(unsigned i);
  /// Reset
  void free();

  /// Returns the size of the vector _mpiRequest/_mpiStatus
  unsigned get_size() const;

  /// Get the specified request object
  MPI_Request* get_mpiRequest(int i=0) const;
  /// Get the specified status object
  MPI_Status* get_mpiStatus(int i=0) const;

  void start(int i);
  void wait(int i);
  bool isDone(int i);

  /// Swap method
  void swap(MpiNonBlockingHelper& rhs);
};

/// Wrapper functions that simplify the use of MPI

class MpiManager {
public:
  /// Initializes the mpi manager
  void init(int *argc, char ***argv, bool verbose=true);
  /// Returns the number of processes
  int getSize() const;
  /// Returns the process ID
  int getRank() const;
  /// Returns process ID of main processor
  int bossId() const;
  /// Tells whether current processor is main processor
  bool isMainProcessor() const;
  /// Returns universal MPI-time in seconds
  double getTime() const;

  /// Synchronizes the processes
  void barrier(MPI_Comm comm = MPI_COMM_WORLD);

  /// Synchronizes the processes and wait to ensure correct cout order
  void synchronizeIO(unsigned tDelay = 100, MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data at *buf, blocking
  template <typename T>
  void send(T *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void send(util::ADf<T,DIM> *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template<typename... args>
  void send(std::vector<args...>& vec, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD){
    send( vec.data(), vec.size(), dest, tag, comm );
  }
  template<class T, std::size_t N>
  void send(std::array<T,N>& array, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD){
    send( array.data(), array.size(), dest, tag, comm );
  }

  /// Initialize persistent non-blocking send
  template <typename T>
  void sendInit(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void sendInit(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data at *buf, non blocking
  template <typename T>
  void iSend(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void iSend(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data at *buf, non blocking and buffered
  template <typename T>
  void ibSend(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void ibSend(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Probe size of incoming message
  std::size_t probeReceiveSize(int source, MPI_Datatype type, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  /// Probe size of incoming message with TYPE
  template <typename TYPE>
  std::size_t probeReceiveSize(int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Receives data at *buf, blocking
  template <typename T>
  void receive(T *buf, int count, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void receive(util::ADf<T,DIM> *buf, int count, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template<typename... args>
  void receive(std::vector<args...>& vec, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD){
    receive( vec.data(), vec.size(), source, tag, comm );
  }
  template<class T, std::size_t N>
  void receive(std::array<T,N>& array, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD){
    receive( array.data(), array.size(), source, tag, comm );
  }

  /// Initialize persistent non-blocking receive
  template <typename T>
  void recvInit(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void recvInit(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Receives data at *buf, non blocking
  template <typename T>
  void iRecv(T *buf, int count, int source, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void iRecv(util::ADf<T,DIM> *buf, int count, int source, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Send and receive data between two partners
  template <typename T>
  void sendRecv(T *sendBuf, T *recvBuf, int count, int dest, int source, int tag = 0,
                MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void sendRecv(util::ADf<T,DIM> *sendBuf, util::ADf<T,DIM> *recvBuf, int count,
                int dest, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data to master processor
  template <typename T>
  void sendToMaster(T* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);

  /// Scatter data from one processor over multiple processors
  template <typename T>
  void scatterv(T *sendBuf, int* sendCounts, int* displs,
                T* recvBuf, int recvCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Gather data from multiple processors to one processor
  template <typename T>
  void gather(T* sendBuf, int sendCount, T* recvBuf, int recvCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Gather data from multiple processors to every processor
  template <typename T>
  void allGather(T* sendBuf, int sendCount, T* recvBuf, int recvCount, MPI_Comm comm = MPI_COMM_WORLD);

  /// Gather data from multiple processors to one processor
  template <typename T>
  void gatherv(T* sendBuf, int sendCount, T* recvBuf, int* recvCounts, int* displs, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Gather data from multiple processors to every processor
  template <typename T>
  void allGatherv(T* sendBuf, int sendCount, T* recvBuf, int* recvCounts, int* displs, MPI_Comm comm = MPI_COMM_WORLD);

  /// Broadcast data from one processor to multiple processors
  template <typename T>
  void bCast(T* sendBuf, int sendCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void bCast(util::ADf<T,DIM>* sendBuf, int sendCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void bCast(BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& sendData, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T>
  void bCast(T& sendVal, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void bCast(Vector<T,DIM>& sendData, int root = 0, MPI_Comm comm = MPI_COMM_WORLD) {
    bCast(sendData.data(), DIM, root, comm);
  }

  /// Broadcast data when root is unknown to other processors
  template <typename T>
  void bCastThroughMaster(T* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void bCastThroughMaster(util::ADf<T,DIM>* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);

  /// Special case for broadcasting strings. Memory handling is automatic.
  void bCast(std::string& message, int root = 0);
  /// Special case for broadcasting BlockData2D
  void bCast(BlockData<2,double,double>& sendData, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  /// Special case for broadcasting BlockData2D
  void bCast(BlockData<2,float,float>& sendData, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Reduction operation toward one processor
  template <typename T>
  void reduce(T& sendVal, T& recvVal, MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);
  template <typename T>
  void reduce(T* sendVal, T* recvVal, int count, MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void reduce(util::ADf<T,DIM>& sendVal, util::ADf<T,DIM>& recvVal,
              MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void reduce(BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& sendVal,
              BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& recvVal,
              MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Element-per-element reduction of a vector of data
  template <typename T>
  void reduceVect(std::vector<T>& sendVal, std::vector<T>& recvVal,
                  MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Reduction operation, followed by a broadcast
  template <typename T>
  void reduceAndBcast(T& reductVal, MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void reduceAndBcast(util::ADf<T,DIM>& reductVal, MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  template <typename T>
  void allreduce(const T* in, T* out, int count, MPI_Op op, MPI_Comm comm = MPI_COMM_WORLD);

  /// All reduction operation of a vector data
  template <typename T>
  void allReduce(T& reductVal, MPI_Op op, MPI_Comm comm = MPI_COMM_WORLD);

  /// Element-per-element all-reduction of a vector of data (in-place)
  template <typename T>
  void allReduceVect(std::vector<T>& reductVal, MPI_Op op, MPI_Comm comm = MPI_COMM_WORLD);

  /// Complete a non-blocking MPI operation
  void wait(MPI_Request* request, MPI_Status* status);

  /// Complete a series of non-blocking MPI operations
  void waitAll(MpiNonBlockingHelper& mpiNbHelper);

private:
  MpiManager();
  ~MpiManager();
private:
  int numTasks, taskId;
  bool ok;
  mutable OstreamManager clout;

  friend MpiManager& mpi();
};

#else

class MpiManager {
public:
  /// Initializes the mpi manager
  void init(int *argc, char ***argv, bool verbose=false) { }
  /// Returns the number of processes
  int getSize() const
  {
    return 1;
  }
  /// Returns the process ID
  int getRank() const
  {
    return 0;
  }
  /// Returns process ID of main processor
  int bossId() const
  {
    return 0;
  }
  /// Tells whether current processor is main processor
  bool isMainProcessor() const
  {
    return true;
  }

  /// Synchronizes the processes
  void barrier() const {};

  friend MpiManager& mpi();
};

#endif

MpiManager& mpi();

}  // namespace singleton

}  // namespace olb


#endif  // MPI_MANAGER_H
