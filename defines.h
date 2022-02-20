/*
 Copyright © 2021  TokiNoBug
This file is part of UniSimulator.

    UniSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    UniSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with UniSimulator.  If not, see <https://www.gnu.org/licenses/>.

    Contact with me:
    github:https://github.com/ToKiNoBug
    bilibili:https://space.bilibili.com/351429231
*/

#ifndef DEFINES_H
#define DEFINES_H


//#define EIGEN_NO_DEBUG

#include <D:/CppLibs/eigen-3.4.0/Eigen/Dense>
#include <D:/CppLibs/eigen-3.4.0/unsupported/Eigen/CXX11/Tensor>

#define DIM_COUNT 2
#define BODY_COUNT 3

using Position=Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>>;

using Velocity=Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>>;

using Acceleration=Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>>;

using Interaction=Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT,BODY_COUNT>>;

using Time=double;

typedef Eigen::Array<double,BODY_COUNT,1> BodyVector ;

typedef Eigen::Array<double,DIM_COUNT,1> DimVector ;

typedef Eigen::Array<double,BODY_COUNT,BODY_COUNT> DistanceMat ;

typedef std::pair<Position,Velocity> Statue ;

typedef  std::pair<Velocity,Acceleration> Derivative ;

typedef std::pair<Time,Time> TimeSpan ;

typedef  std::pair<Time,Statue> Point ;


#if (BODY_COUNT==2) && (DIM_COUNT==2)
#define TEST_BODY2DIM2
#endif

#if (BODY_COUNT==2) && (DIM_COUNT==3)
#define TEST_BODY2DIM3
#endif

#if (BODY_COUNT==3) && (DIM_COUNT==3)
#define BODY3_DIM3
#define TEST_BODY3DIM3
#endif

#if DIM_COUNT <=0
#error DIM_COUNT should be a positive integer
#endif

#if BODY_COUNT <2
#error BODY_COUNT should be a positive integer not less than 2
#endif

const bool DoMotionAlign=false;

#endif // DEFINES_H
