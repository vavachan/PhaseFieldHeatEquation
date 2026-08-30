/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 František Prinz
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


/**
 *  \file eulerRotation.h
 *  \brief $(Helper functions for subgrid scale spheroid rotation computation in moving reference frame)
 *  This file contains helper functions for computation of movement of subgrid scale rotational elipsoids (spheroid).
 *  This method is also called Euler-Lagrange-Euler-Rotation.
 *  Definitions and conventions for the euler angles and quaternions used here are based on
 *  Lin Tian,
 *  Computational modeling of fiber transport in human respiratory airways—A review
 *  Matrices are treated as 9 dimensional vectors with following component ordering:
 *   / 0 1 2 \
 *   | 3 4 5 |
 *   \ 6 7 8 /
 *  By matrix vector multiplication the column vectors are used!
 *  This method is appliable just for the 3dimensional case (due to quaternion definition)
 *  DOI https://doi.org/10.1007/s42757-020-0061 -7
 *  \author František Prinz
 */

/* This file contains functions used for the calculation of particle related dynamics.
 *
*/

#ifndef EULER_ROTATION_H
#define EULER_ROTATION_H


#include "core/vector.h"


namespace olb {

  namespace eler{ //Euler-Lagrange-Euler-Rotation

//BE CAREFULL ABOUT THE QUATERNION NOTATION xyzomega

///normalize quaternion (size==4) to magnitude 1
template <typename T>
void qnormalize(Vector <T,4> & vec)
{
  T norm = 0.;
  for(int i= 0; i < 4;i++)
   norm += vec[i]*vec[i];

  if (norm != 0) {
    norm = olb::util::sqrt(norm);
    for(int i=0; i< 4;i++)
      vec[i] /= norm;
    }
  }

///compute transformation matrix from eul angles using x-convention Goldstein(1980), in degrees
template <typename T, unsigned D>
  Vector<T,D*D> computeTransformMatrixFromEulAng(olb::Vector<T,D> eulAng)//x-convention Goldstein(1980), in degrees
  {
    T phi = olb::util::degreeToRadian(eulAng[0]);
    T theta = olb::util::degreeToRadian(eulAng[1]);
    T varphi = olb::util::degreeToRadian(eulAng[2]);

   /* / 0 1 2 \
    | 3 4 5 |
     \ 6 7 8 /*/
using namespace olb::util;
Vector<T,9> trafoMatrix;
trafoMatrix[0] = cos(varphi)*cos(phi) - cos(theta)*sin(phi)*sin(varphi);
trafoMatrix[3] = -1.*sin(varphi)*cos(phi) - cos(theta)*sin(phi)*cos(varphi);
trafoMatrix[6] = sin(theta)*sin(phi);

trafoMatrix[1] = cos(varphi)*sin(phi)+cos(theta)*cos(phi)*sin(varphi);
trafoMatrix[4] = -1.*sin(varphi)*sin(phi)+cos(theta)*cos(phi)*cos(varphi);
trafoMatrix[7] =  -1.*sin(theta)*cos(phi);

trafoMatrix[2] =    sin(varphi)*sin(theta);
trafoMatrix[5] =    cos(varphi)*sin(theta);
trafoMatrix[8] =     cos(theta);

     return trafoMatrix;
}

///compute quaternion from eul angles using x-convention Goldstein(1980), in degrees
template <typename T>
Vector<T,4> computeEulQuatFromEulAng(olb::Vector<T,3> eulAng)//x-convention Goldstein(1980)
{
Vector<T,4> quat (0.,0.,0.,0.);
Vector<T,9> tMat = computeTransformMatrixFromEulAng(eulAng);
 if((1+tMat[0] + tMat[4] + tMat[8]) > 1.0e-8)
 {
    quat[3] = 1./2.*(olb::util::sqrt(1+tMat[0] + tMat[4] + tMat[8]));
    quat[0] = 1./(4*quat[3])*(tMat[5]-tMat[7]);
    quat[1] = 1./(4*quat[3])*(tMat[6]-tMat[2]);
    quat[2] = 1./(4*quat[3])*(tMat[1]-tMat[3]);
 }
 else
    {
    quat[3] = 0;
    quat[0] = util::sqrt((1+tMat[0])/2.);
    quat[1] = tMat[1]/(2*quat[0]);
    quat[2] = tMat[5]/(2*quat[1]);
    }
 return quat;
}

///compute transformation matrix from eul angles using x-convention Goldstein(1980), in degrees
template <typename T>
Vector<T,9> computeTransformMatrixFromEulQuat(olb::Vector<T,4> eulQuat)//x-convention Goldstein(1980)
  {
   /* / 0 1 2 \
    | 3 4 5 |
     \ 6 7 8 /*/
using namespace olb::util;
Vector<T,9> trafoMatrix;
trafoMatrix[0] = 1 - 2*(eulQuat[1]*eulQuat[1] + eulQuat[2]*eulQuat[2]);
trafoMatrix[3] = 2*(eulQuat[1]*eulQuat[0]-eulQuat[2]*eulQuat[3]);
trafoMatrix[6] = 2*(eulQuat[2]*eulQuat[0]+eulQuat[1]*eulQuat[3]);
trafoMatrix[1] = 2*(eulQuat[1]*eulQuat[0]+eulQuat[2]*eulQuat[3]);
trafoMatrix[4] = 1 - 2*(eulQuat[2]*eulQuat[2] + eulQuat[0]*eulQuat[0]);
trafoMatrix[7] =  2*(eulQuat[1]*eulQuat[2]-eulQuat[0]*eulQuat[3]);
trafoMatrix[2] =   2*(eulQuat[2]*eulQuat[0]-eulQuat[1]*eulQuat[3]);
trafoMatrix[5] =   2*(eulQuat[1]*eulQuat[2]+eulQuat[0]*eulQuat[3]);
trafoMatrix[8] = 1 - 2*(eulQuat[0]*eulQuat[0] + eulQuat[1]*eulQuat[1]);

return trafoMatrix;
 }
    template <typename T>
  Vector<T,3> rotate (Vector<T,9> tMat, Vector<T,3> pos)
  {
    Vector<T,3> out;
    out[0] = tMat[0]*pos[0] + tMat[1]*pos[1] + tMat[2]*pos[2];
    out[1] = tMat[3]*pos[0] + tMat[4]*pos[1] + tMat[5]*pos[2];
    out[2] = tMat[6]*pos[0] + tMat[7]*pos[1] + tMat[8]*pos[2];

return out;
  }

///compute inverse transformation matrix (working just for orthogonal trafo matrices)
template <typename T>
Vector<T,9> inverseTransformMatrix(Vector<T,9>  tMat)
{
/* / 0 1 2 \
   | 3 4 5 |
  \ 6 7 8 /*/
Vector <T,9> invTMat;
invTMat [0]= tMat[0];
invTMat [4]= tMat[4];
invTMat [8]= tMat[8];
invTMat[1] = tMat[3];
invTMat[2] = tMat[6];
invTMat[3] = tMat[1];
invTMat[5] = tMat[7];
invTMat[6] = tMat[2];
invTMat[7] = tMat[5];
return invTMat;
}

///compute rotational dynamics of a object from the time step and vector of rotational velocity (rad/s) by explicit Euler
   template <typename T>
   void quaternionRotate ( Vector <T,4>  & quat, T time_step, Vector<T,3>  revo)
   {
     Vector<T,4> temp (0.,0.,0.,0.);
     temp[0] = 0.5*(quat[3]*revo[0]-quat[2]*revo[1] + quat[1]*revo[2])*time_step;
     temp[1] = 0.5*(quat[2]*revo[0]+quat[3]*revo[1] - quat[0]*revo[2])*time_step;
     temp[2] = 0.5*(-1.*quat[1]*revo[0]+quat[0]*revo[1] + quat[3]*revo[2])*time_step;
     temp[3] = 0.5*(-1.*quat[0]*revo[0]-quat[1]*revo[1] - quat[2]*revo[2])*time_step;
     quat =quat + temp;
     qnormalize(quat);
     return;
   }

   ///compute rotational dynamics of a object from the time step and vector of rotational velocity (rad/s) by Leap Frog
   template <typename T>
   void quaternionRotateLeapFrog ( Vector <T,4>  & quat, T time_step, Vector<T,3>  revo)
   {
     Vector<T,4> temp (0.,0.,0.,0.);
     //temp  = quat;
     Vector <T,4> quat_temp (0.,0.,0.,0.);
     quat_temp = quat;
     for(int i=0; i<5;i++)
     {
     temp[0] = 0.5*0.5*(quat_temp[3]*revo[0]-quat_temp[2]*revo[1] + quat_temp[1]*revo[2])*time_step;
     temp[1] = 0.5*0.5*(quat_temp[2]*revo[0]+quat_temp[3]*revo[1] - quat_temp[0]*revo[2])*time_step;
     temp[2] = 0.5*0.5*(-1.*quat_temp[1]*revo[0]+quat_temp[0]*revo[1] + quat_temp[3]*revo[2])*time_step;
     temp[3] = 0.5*0.5*(-1.*quat_temp[0]*revo[0]-quat_temp[1]*revo[1] - quat_temp[2]*revo[2])*time_step;
     quat_temp = quat + temp;
     }
     quat =quat + temp;
     qnormalize(quat);
     return;
   }

   //https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
    template <typename T>
   void quaternionRotateRungeKutta ( Vector <T,4>  & quat, T time_step, Vector<T,3>  revo)
   {
     Vector<T,4> K1,K2,K3,K4;
     K1[0] = 0.5*(quat[3]*revo[0]-quat[2]*revo[1] + quat[1]*revo[2]);
     K1[1] = 0.5*(quat[2]*revo[0]+quat[3]*revo[1] - quat[0]*revo[2]);
     K1[2] = 0.5*(-1.*quat[1]*revo[0]+quat[0]*revo[1] + quat[3]*revo[2]);
     K1[3] = 0.5*(-1.*quat[0]*revo[0]-quat[1]*revo[1] - quat[2]*revo[2]);
     K2[0] = 0.5*((quat[3]+time_step*K1[3]*0.5)*revo[0]-(quat[2]+time_step*K1[2]*0.5)*revo[1] + (quat[1]+time_step*K1[1]*0.5)*revo[2]);
     K2[1] = 0.5*((quat[2]+time_step*K1[2]*0.5)*revo[0]+(quat[3]+time_step*K1[3]*0.5)*revo[1] - (quat[0]+time_step*K1[0]*0.5)*revo[2]);
     K2[2] = 0.5*(-1.*(quat[1]+time_step*K1[1]*0.5)*revo[0]+(quat[0]+time_step*K1[0]*0.5)*revo[1] + (quat[3]+time_step*K1[3]*0.5)*revo[2]);
     K2[3] = 0.5*(-1.*(quat[0]+time_step*K1[0]*0.5)*revo[0]-(quat[1]+time_step*K1[1]*0.5)*revo[1] - (quat[2]+time_step*K1[2]*0.5)*revo[2]);
     K3[0] = 0.5*((quat[3]+time_step*K2[3]*0.5)*revo[0]-(quat[2]+time_step*K2[2]*0.5)*revo[1] + (quat[1]+time_step*K2[1]*0.5)*revo[2]);
     K3[1] = 0.5*((quat[2]+time_step*K2[2]*0.5)*revo[0]+(quat[3]+time_step*K2[3]*0.5)*revo[1] - (quat[0]+time_step*K2[0]*0.5)*revo[2]);
     K3[2] = 0.5*(-1.*(quat[1]+time_step*K2[1]*0.5)*revo[0]+(quat[0]+time_step*K2[0]*0.5)*revo[1] + (quat[3]+time_step*K2[3]*0.5)*revo[2]);
     K3[3] = 0.5*(-1.*(quat[0]+time_step*K2[0]*0.5)*revo[0]-(quat[1]+time_step*K2[1]*0.5)*revo[1] - (quat[2]+time_step*K2[2]*0.5)*revo[2]);
     K4[0] = 0.5*((quat[3]+time_step*K3[3])*revo[0]-(quat[2]+time_step*K3[2])*revo[1] + (quat[1]+time_step*K3[1])*revo[2]);
     K4[1] = 0.5*((quat[2]+time_step*K3[2])*revo[0]+(quat[3]+time_step*K3[3])*revo[1] - (quat[0]+time_step*K3[0])*revo[2]);
     K4[2] = 0.5*(-1.*(quat[1]+time_step*K3[1])*revo[0]+(quat[0]+time_step*K3[0])*revo[1] + (quat[3]+time_step*K3[3])*revo[2]);
     K4[3] = 0.5*(-1.*(quat[0]+time_step*K3[0])*revo[0]-(quat[1]+time_step*K3[1])*revo[1] - (quat[2]+time_step*K3[2])*revo[2]);
     quat[0] = quat[0]+time_step/6.*(K1[0]+2*K2[0]+2*K3[0]+K4[0]);
     quat[1] = quat[1]+time_step/6.*(K1[1]+2*K2[1]+2*K3[1]+K4[1]);
     quat[2] = quat[2]+time_step/6.*(K1[2]+2*K2[2]+2*K3[2]+K4[2]);
     quat[3] = quat[3]+time_step/6.*(K1[3]+2*K2[3]+2*K3[3]+K4[3]);
     return;


   }



   ///3-dimensional matrix multiplication
   template <typename T>
    Vector<T,9> matrixMultiply (Vector <T,9>  m1 , Vector <T,9>  m2)
    {
       /*/ 0 1 2 \ / 0 1 2 \
       | 3 4 5 |  | 3 4 5 |
        \ 6 7 8/ \ 6 7 8 /*/
      Vector<T,9> res (0.);
      res[0] = m1[0]*m2[0] + m1[1]*m2[3] + m1[2]*m2[6];
      res[1] = m1[0]*m2[1] + m1[1]*m2[4] + m1[2]*m2[7];
      res[2] = m1[0]*m2[2] + m1[1]*m2[5] + m1[2]*m2[8];
      res[3] = m1[3]*m2[0] + m1[4]*m2[3] + m1[5]*m2[6];
      res[4] = m1[3]*m2[1] + m1[4]*m2[4] + m1[5]*m2[7];
      res[5] = m1[3]*m2[2] + m1[4]*m2[5] + m1[5]*m2[8];
      res[6] = m1[6]*m2[0] + m1[7]*m2[3] + m1[8]*m2[6];
      res[7] = m1[6]*m2[1] + m1[7]*m2[4] + m1[8]*m2[7];
      res[8] = m1[6]*m2[2] + m1[7]*m2[5] + m1[8]*m2[8];
      return res;
    }
    ///3-dimensional matrix columnvector multiplication matrix*vector
    template <typename T>
    Vector<T,3> matrixVectorMultiply (Vector <T,9>  mat, Vector <T,3>  vec)
    {
      Vector<T,3> res (0.);
      res[0] = mat[0]*vec[0] + mat[1]*vec[1] + mat[2]*vec[2];
      res[1] = mat[3]*vec[0] + mat[4]*vec[1] + mat[5]*vec[2];
      res[2] = mat[6]*vec[0] + mat[7]*vec[1] + mat[8]*vec[2];
      return res;
    }

    ///update angular velocity ang_vel from torque , mom_inertia and time-step using explicit Euler
    template<typename T>
    void computeAngularVelocity(Vector<T,3> & ang_vel,Vector <T,3>   torque, Vector<T,3>   mom_inertia, T delta_t )
    {
      Vector<T,3> res (0.);
      //std::cout << "delta_t " << delta_t << std::endl;
      res[0] = (torque[0]/mom_inertia[0] + ang_vel[1]*ang_vel[2]*(mom_inertia[1]-mom_inertia[2])/mom_inertia[0])*delta_t+ ang_vel[0];
      res[1] = (torque[1]/mom_inertia[1] + ang_vel[2]*ang_vel[0]*(mom_inertia[2]-mom_inertia[0])/mom_inertia[1])*delta_t+ ang_vel[1];
      res[2] = (torque[2]/mom_inertia[2] + ang_vel[0]*ang_vel[1]*(mom_inertia[0]-mom_inertia[1])/mom_inertia[2])*delta_t+ ang_vel[2];
      ang_vel = res;
    }

    ///update angular velocity ang_vel from torque , mom_inertia and time-step using explicit Euler Leap Frog
    template<typename T>
    void computeAngularVelocityLeapFrog(Vector<T,3> & ang_vel,Vector <T,3>   torque, Vector<T,3>   mom_inertia, T delta_t )
    {
      Vector<T,3>ang_vel_new (ang_vel);
      Vector<T,3> ang_temp (0.);


      for (int i=0;i<5;i++)///hardcoded 5 inner iterations!
      {
      ang_temp[0] = (torque[0]/mom_inertia[0] + ang_vel_new[1]*ang_vel_new[2]*(mom_inertia[1]-mom_inertia[2])/mom_inertia[0])*0.5*delta_t+ ang_vel[0];
      ang_temp[1] = (torque[1]/mom_inertia[1] + ang_vel_new[2]*ang_vel_new[0]*(mom_inertia[2]-mom_inertia[0])/mom_inertia[1])*0.5*delta_t+ ang_vel[1];
      ang_temp[2] = (torque[2]/mom_inertia[2] + ang_vel_new[0]*ang_vel_new[1]*(mom_inertia[0]-mom_inertia[1])/mom_inertia[2])*0.5*delta_t+ ang_vel[2];
      ang_vel_new  = ang_temp;
      //std::cout << "ang temp iteration " << i << " " << ang_temp[0] << " " << ang_temp[1] << " " <<ang_temp[2] << std::endl;
      }
      ang_temp[0] = (torque[0]/mom_inertia[0] + ang_vel_new[1]*ang_vel_new[2]*(mom_inertia[1]-mom_inertia[2])/mom_inertia[0])*delta_t+ ang_vel[0];
      ang_temp[1] = (torque[1]/mom_inertia[1] + ang_vel_new[2]*ang_vel_new[0]*(mom_inertia[2]-mom_inertia[0])/mom_inertia[1])*delta_t+ ang_vel[1];
      ang_temp[2] = (torque[2]/mom_inertia[2] + ang_vel_new[0]*ang_vel_new[1]*(mom_inertia[0]-mom_inertia[1])/mom_inertia[2])*delta_t+ ang_vel[2];
      ang_vel = ang_temp;
     // std::cout << ang_vel[0] << " " << ang_vel[1] << " " <<ang_vel[2] << std::endl;
    }

    template<typename T>
    void computeAngularVelocityRungeKutta(Vector<T,3> & ang_vel,Vector <T,3>   torque, Vector<T,3>   mom_inertia, T time_step )
    {
    Vector<T,3> K1,K2,K3,K4;
    K1[0] = (torque[0]/mom_inertia[0] + ang_vel[1]*ang_vel[2]*(mom_inertia[1]-mom_inertia[2])/mom_inertia[0]);
    K1[1] = (torque[1]/mom_inertia[1] + ang_vel[2]*ang_vel[0]*(mom_inertia[2]-mom_inertia[0])/mom_inertia[1]);
    K1[2] = (torque[2]/mom_inertia[2] + ang_vel[0]*ang_vel[1]*(mom_inertia[0]-mom_inertia[1])/mom_inertia[2]);
    K2[0] = (torque[0]/mom_inertia[0] + (ang_vel[1]+time_step*K1[1]*0.5)*(ang_vel[2]*time_step*K1[2]*0.5)*(mom_inertia[1]-mom_inertia[2])/mom_inertia[0]);
    K2[1] = (torque[1]/mom_inertia[1] + (ang_vel[2]+time_step*K1[2]*0.5)*(ang_vel[0]+time_step*K1[0]*0.5)*(mom_inertia[2]-mom_inertia[0])/mom_inertia[1]);
    K2[2] = (torque[2]/mom_inertia[2] + (ang_vel[0]+time_step*K1[0]*0.5)*(ang_vel[1]+time_step*K1[1]*0.5)*(mom_inertia[0]-mom_inertia[1])/mom_inertia[2]);
    K3[0] = (torque[0]/mom_inertia[0] + (ang_vel[1]+time_step*K2[1]*0.5)*(ang_vel[2]*time_step*K2[2]*0.5)*(mom_inertia[1]-mom_inertia[2])/mom_inertia[0]);
    K3[1] = (torque[1]/mom_inertia[1] + (ang_vel[2]+time_step*K2[2]*0.5)*(ang_vel[0]+time_step*K2[0]*0.5)*(mom_inertia[2]-mom_inertia[0])/mom_inertia[1]);
    K3[2] = (torque[2]/mom_inertia[2] + (ang_vel[0]+time_step*K2[0]*0.5)*(ang_vel[1]+time_step*K2[1]*0.5)*(mom_inertia[0]-mom_inertia[1])/mom_inertia[2]);
    K4[0] = (torque[0]/mom_inertia[0] + (ang_vel[1]+time_step*K3[1])*(ang_vel[2]*time_step*K3[2])*(mom_inertia[1]-mom_inertia[2])/mom_inertia[0]);
    K4[1] = (torque[1]/mom_inertia[1] + (ang_vel[2]+time_step*K3[2])*(ang_vel[0]+time_step*K3[0])*(mom_inertia[2]-mom_inertia[0])/mom_inertia[1]);
    K4[2] = (torque[2]/mom_inertia[2] + (ang_vel[0]+time_step*K3[0])*(ang_vel[1]+time_step*K3[1])*(mom_inertia[0]-mom_inertia[1])/mom_inertia[2]);
    ang_vel[0] = ang_vel[0] + time_step/6.*(K1[0]+2*K2[0]+2*K3[0]+K4[0]);
    ang_vel[1] = ang_vel[1] + time_step/6.*(K1[1]+2*K2[1]+2*K3[1]+K4[1]);
    ang_vel[2] = ang_vel[2] + time_step/6.*(K1[2]+2*K2[2]+2*K3[2]+K4[2]);
    return;
    }





    ///computes velocity gradient in given direction in positionArray by central finite difference using points in the distance of deltaX
    template<typename T, typename DESCRIPTOR>
    void veloGradientInterpFD ( Vector<T,3> &res, BlockLatticeInterpPhysVelocity<T,DESCRIPTOR> & blockInterpPhysVelF , Vector<T,3> & direction, Vector<T,3> & positionArray, T deltaX)
    {

      T fluidVelArray2 [3];
      T fluidVelArray1 [3];
      auto neighbour1 = positionArray+(deltaX*direction);
      blockInterpPhysVelF(fluidVelArray2,      neighbour1.data());
      auto neighbour2 = positionArray-(deltaX*direction);
       blockInterpPhysVelF(fluidVelArray1, neighbour2.data());
      res[0] = (fluidVelArray2[0] - fluidVelArray1[0])/(2*deltaX);
      res[1] = (fluidVelArray2[1] - fluidVelArray1[1])/(2*deltaX);
      res[2] = (fluidVelArray2[2] - fluidVelArray1[2])/(2*deltaX);

    }

    template<typename T>
    Vector<T,9> veloGradientAnalyticalPoiuseuille3D (Vector<T,3>& orientation, Vector<T,3>& pos, T radius, T velo_max, Vector<T,2> midp_pos_proj = (0.)  )
    {
      Vector <T,9> res (0.);
      if (orientation[0] > 0.9999 && orientation[1] < 0.0001 && orientation[2] < 0.0001)
      {
        res[1] = 8*velo_max*(pos[1]-midp_pos_proj[0])/radius^2;
        res[2] = 8*velo_max*(pos[2]-midp_pos_proj[1])/radius^2;
      }
      else if (orientation[0] < 0.0001 && orientation[1] > 0.9999 && orientation[2] < 0.0001)
      {
        res[3] = 8*velo_max*(pos[0]-midp_pos_proj[0])/radius^2;
        res[5] = 8*velo_max*(pos[2]-midp_pos_proj[1])/radius^2;
      }
      else if (orientation[0] < 0.0001 && orientation[1] < 0.0001 && orientation[2] > 0.9999)
      {
        res[6] = 8*velo_max*(pos[0]-midp_pos_proj[0])/radius^2;
        res[7] = 8*velo_max*(pos[1]-midp_pos_proj[1])/radius^2;
      }
      else
      { std::cout << "WARNING! VeloGradientAnalytialPoiseuille3D: No discrete normal!" << std::endl;}

    }

    template<typename T, typename DESCRIPTOR>
    void veloGradientInterpFDTurbulenceAnalytical( Vector<T,3> &res,CirclePoiseuille3D<T> & pois,  Vector<T,3> & direction, Vector<T,3> & positionArray, T deltaX, Vector<T,3> & std_dev)
    {
      bool verbose =false;
      std::random_device rd;
      if (direction[0] == 1) {
        if(verbose) { std::cout << "X" << std::endl; }
      }
      else if (direction[1] == 1) {
        if (verbose) { std::cout << "Y" << std::endl; }
      }
      else if(direction[2] == 1) {
        if(verbose) { std::cout << "Z" << std::endl; }
      }
      std::mt19937 generator (rd());
      T fluidVelArray2 [3];
      T fluidVelArray1 [3];
      if(verbose) std::cout.precision(15);
      auto neighbour1 = positionArray+(deltaX*direction);
     // std::cout << "data neighbour1 " << neighbour1 << std::endl;
     // std::cout << neighbour1.data()[0] << " " << neighbour1.data()[1] << " " << neighbour1.data()[2] <<std::endl;
      T pos[3] = {neighbour1[0], neighbour1[1], neighbour1[2]};
      pois(fluidVelArray1, neighbour1.data());
      if(verbose) std::cout  << fluidVelArray1[0] << " " <<   fluidVelArray1[1] << " " <<  fluidVelArray1[2] << std::endl;
      //blockInterpPhysVelF(fluidVelArray2,neighbour1.data());
      auto neighbour2 = positionArray-(deltaX*direction);
      //std::cout << "data neighbour2 " << neighbour2 << std::endl;
      pois(fluidVelArray2,neighbour2.data());
      if(verbose) std::cout << fluidVelArray2[0] << " " <<   fluidVelArray2[1] << " " <<  fluidVelArray2[2] << std::endl;
      std::normal_distribution<T> distx ((fluidVelArray2[0] - fluidVelArray1[0])/(2*deltaX), util::sqrt(2.)*std_dev[0]/(2*deltaX));
      std::normal_distribution<T> disty ((fluidVelArray2[1] - fluidVelArray1[1])/(2*deltaX), util::sqrt(2.)*std_dev[1]/(2*deltaX));
      std::normal_distribution<T> distz ((fluidVelArray2[2] - fluidVelArray1[2])/(2*deltaX), util::sqrt(2.)*std_dev[2]/(2*deltaX));
      auto DeltaVec = 2*deltaX*direction;
      if(verbose) std::cout << "thoery " << (fluidVelArray1[0] - fluidVelArray2[0])/(2*deltaX) << "," <<(fluidVelArray1[1] - fluidVelArray2[1])/(2*deltaX)<< "," <<(fluidVelArray1[2] - fluidVelArray2[2])/(2*deltaX) << std::endl;
      res[0] = distx(generator);
      res[0] = util::min(util::max(res[0],(fluidVelArray2[0] - fluidVelArray1[0])/(2*deltaX)-2.*util::sqrt(2.)*std_dev[0]/(2*deltaX)),(fluidVelArray2[0] - fluidVelArray1[0])/(2*deltaX)+2.*util::sqrt(2.)*std_dev[0]/(2*deltaX));
      res[1] = util::min(util::max(res[1],(fluidVelArray2[1] - fluidVelArray1[1])/(2*deltaX)-2.*util::sqrt(2.)*std_dev[1]/(2*deltaX)),(fluidVelArray2[1] - fluidVelArray1[1])/(2*deltaX)+2.*util::sqrt(2.)*std_dev[1]/(2*deltaX));
      res[2] = util::min(util::max(res[2],(fluidVelArray2[2] - fluidVelArray1[2])/(2*deltaX)-2.*util::sqrt(2.)*std_dev[2]/(2*deltaX)),(fluidVelArray2[2] - fluidVelArray1[2])/(2*deltaX)+2.*util::sqrt(2.)*std_dev[2]/(2*deltaX));
      if(verbose) std::cout << res << std::endl;
    }





    ///computes velocity gradient in given direction in positionArray by central finite difference using points in the distance of deltaX
    template<typename T, typename DESCRIPTOR>
    void veloGradientInterpFDTurbulence ( Vector<T,3> &res, BlockLatticeInterpPhysVelocity<T,DESCRIPTOR> & blockInterpPhysVelF , Vector<T,3> & direction, Vector<T,3> & positionArray, T deltaX, Vector<T,3> & std_dev)
    {
      std::random_device rd;
      std::mt19937 generator (rd());
      T fluidVelArray2 [3];
      T fluidVelArray1 [3];
      auto neighbour1 = positionArray+(deltaX*direction);
      blockInterpPhysVelF(fluidVelArray2,neighbour1.data());
      auto neighbour2 = positionArray-(deltaX*direction);
      blockInterpPhysVelF(fluidVelArray1, neighbour2.data());
      std::normal_distribution<T> distx ((fluidVelArray2[0] - fluidVelArray1[0])/(2*deltaX), util::sqrt(2*std_dev[0]/(2*deltaX)));
      std::normal_distribution<T> disty ((fluidVelArray2[1] - fluidVelArray1[1])/(2*deltaX), util::sqrt(2*std_dev[1]/(2*deltaX)));
      std::normal_distribution<T> distz ((fluidVelArray2[2] - fluidVelArray1[2])/(2*deltaX), util::sqrt(2*std_dev[2]/(2*deltaX)));

      res[0] = distx(generator);
      res[1] = disty(generator);
      res[2] = distz(generator);
    }

    ///computes velocity gradient in given direction in positionArray by central finite difference using points in the distance of deltaX
    template<typename T, typename DESCRIPTOR>
    void veloGradientInterpFDTurbulenceAnalytical ( Vector<T,3> & res , Vector<T,3> & direction, Vector<T,3> & positionArray, T deltaX, Vector<T,3> & std_dev)
    {
      std::random_device rd;
      std::mt19937 generator (rd());
      T fluidVelArray2 [3];
      T fluidVelArray1 [3];

      std::normal_distribution<T> distx ((fluidVelArray2[0] - fluidVelArray1[0])/(2*deltaX), util::sqrt(2*std_dev[0]/(2*deltaX)));
      std::normal_distribution<T> disty ((fluidVelArray2[1] - fluidVelArray1[1])/(2*deltaX), util::sqrt(2*std_dev[1]/(2*deltaX)));
      std::normal_distribution<T> distz ((fluidVelArray2[2] - fluidVelArray1[2])/(2*deltaX), util::sqrt(2*std_dev[2]/(2*deltaX)));


      res[0] = distx(generator);
      res[1] = disty(generator);
      res[2] = distz(generator);


    }




    ///return quaternion composition of two rotations defined by quaternions
      template<typename T>
    Vector<T,4> multiplyQuaternion(Vector<T,4> & a, Vector<T,4> & b)
    {
      Vector <T,4> res (0.);
     Vector<T,3> v (a[0], a[1], a[2]);
     Vector<T,3> w (b[0], b[1], b[2]);
     T s = a[3];
     T t = b[3];
     res[3] = s*t - v*w;
     Vector<T,3> res_vec = s*w + t*v+ crossProduct3D(v,w);
      res[0] = res_vec[0];
      res[1] = res_vec[1];
      res[2] = res_vec[2];
      return res;
    }

   ///this function uses the different original orientation of the object (1,0,0) instead of (0,0,1)
   ///and also another notation
   ///used to correcty visualize spheroid orientation in Paraview due to different convention there
    template<typename T>
    Vector<T,4> getParaquaternion(Vector<T,4>  quat )
    {
      Vector<T,4> quat_to_z (0.,util::sqrt(2.)/2.,0.,util::sqrt(2.)/2.);
      Vector<T,4> res2 (0.);
      res2 = multiplyQuaternion(quat, quat_to_z);
      Vector<T,4> res (res2[3],res2[0],res2[1], res2[2]);
      res = normalize(res);
      return res;
    }

    ///get uniformly distributed random Orientation Vector defined in the Euler Angles in degrees
    template<typename T>
    Vector<T,3> getRandomOrientation()
    {
      //std::srand((unsigned) time(NULL));
     T n1 =  static_cast<double>(std::rand()) / RAND_MAX;
     T n2 = static_cast<double>(std::rand()) / RAND_MAX;
     T n3 = static_cast<double>(std::rand()) / RAND_MAX;
     Vector<T,3> result;
     result[0] = 360*n1;
     result[1] = util::acos((2*n2-1))*180/M_PI;
     result[2] = 360*n3;
     return result;
    }

  ///this is used for detection of the exact position of particle when located near to the wall
  ///DOI 10.1007/s42757-020-0061-7 Tian and Ahmadi 2021, Chapter 2.5.2, (3) scennario, modification with respect to the exact location of the spheroid
   template<typename T>
    bool elipsoidDeposition(T semi_minor, T aspect_ratio, T distance, Vector<T,3>& orientation, Vector<T,3>& normal_to_surface)
    {
      T a = semi_minor;
      T b = semi_minor*aspect_ratio;
      using namespace util;
      orientation = normalize(orientation);
      normal_to_surface = normalize(normal_to_surface);
      T alpha = std::acos(util::dotProduct3D(orientation, normal_to_surface));
      if(alpha>3.141592653/2.)
        alpha = 3.141592653-alpha;
      alpha = 3.141592653/2.-alpha;
      T zeta = util::sqrt((util::pow<T,int>((b*b-a*a),2)*(util::pow<T,int>(cos(alpha)/a,2)+util::pow<T,int>((sin(alpha)/b),2))*util::pow<T,int>(sin(alpha),2)*util::pow<T,int>(cos(alpha),2))/
      (1.+a*a*b*b*util::pow<T,int>((1./(a*a)-1./(b*b)),2)*util::pow<T,int>(sin(alpha),2)*util::pow<T,int>(cos(alpha),2)));
      T theta = (-1.*zeta*(1./(a*a)-1./(b*b))*sin(alpha)*cos(alpha)-
      sqrt(util::pow<T,int>(cos(alpha)/a,2)+util::pow<T,int>(sin(alpha)/b,2)-util::pow<T,int>(zeta/(a*b),2)))/
      (util::pow<T,int>(cos(alpha)/a,2)+util::pow<T,int>(sin(alpha)/b,2));
      if(util::abs(theta)< distance)
        return false;
      else
      {
        std::cout << "DEPOSITED DUE TO ORIENTATION! theta " << util::abs(theta) << " distance " << distance << std::endl;
        return true;
      }

    }

    // Function to generate a random unit quaternion
    template<typename T>
    Vector<T, 4> randomUnitQuaternion() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(0.0, 1.0);

    T u[3] = {dis(gen), dis(gen), dis(gen)};
    T q[4];
    q[0] = sqrt(1 - u[0]) * std::sin(2 * M_PI * u[1]);
    q[1] = sqrt(1 - u[0]) * std::cos(2 * M_PI * u[1]);
    q[2] = sqrt(u[0]) * std::sin(2 * M_PI * u[2]);
    q[3] = sqrt(u[0]) * std::cos(2 * M_PI * u[2]);

    // Normalize quaternion
    T norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    Vector<T, 4> result (0.);
    result[0] = q[0] / norm;
    result[1] = q[1] / norm;
    result[2] = q[2] / norm;
    result[3] = q[3] / norm;
    return result;
}

  } //eler
}


#endif
