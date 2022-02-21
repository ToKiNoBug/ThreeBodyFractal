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

#ifndef SIMULATOR_H
#define SIMULATOR_H

#define _USE_MATH_DEFINES
//#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include <tuple>
#include <fstream>
#include <string>

#include "defines.h"

extern const double G; //gravity constant
extern const double Ms; //standard mass (solar mass)
extern const double omega_s; //standard angle speed
extern const double rho; //solar density
extern const double year; //seconds of a year
extern const double rs; //standard distance (1AU)
extern const double vs; //standard speed
extern const double as; //standard accelerate


class Simulator
{
public:
    Simulator();

    enum Algorithm {
        Euler,
        RK4Fixed,
        RK4Var1,
        //ODE45
    };
    static const std::string paraSuffix;
    static const std::string dataSuffix;
public:
    void simulateEuler(const double,
                       TimeSpan,
                       Statue,
                       bool * noCollide=nullptr);

    void simulateRK4Fixed(const double,
                          TimeSpan,
                          Statue,
                          bool * noCollide=nullptr);
    
    template<bool record=true>
    void simulateRK4Var1(double step,
                        TimeSpan tSpan, 
                        Statue y, 
                        bool *noCollide,
                        double precision=1e-8) {
    if(tSpan.second==tSpan.first) {
        mexErrMsgTxt("the begging time shouldn't equal to end time");
        return;
    }

    if(tSpan.second<tSpan.first) {
        std::swap(tSpan.first,tSpan.second);
    }

    if(tSpan.second-tSpan.first<step) {
        mexErrMsgTxt("time span should be greater than time step");
        return;
    }

    clear();

    positonAlign(y.first);
    motionAlign(mass,y.second);

    DistanceMat safeDistance;
    calculateSafeDistance(mass,safeDistance);

    Interaction GM;
    calculateGM(mass,GM);

    Acceleration acc;
    Time curTime=tSpan.first;

    Statue y_h,y_h_2;

    static const double searchRatio=0.5;
    static const int rank=4;
    static const double ratio=std::pow(2,rank)-1;

    while(true) {
        if(record)
            sol.emplace_back(std::make_pair(curTime,y));

        if(curTime>tSpan.second) {
            break;
        }

        double minStep=16*std::nextafter(curTime,curTime*2)-16*curTime;

        bool isOk=true;

        isOk=RK4(step,y,GM,safeDistance,y_h);
        if(!isOk) {
            if(noCollide!=nullptr) {
                *noCollide=isOk;
            }
            //std::cerr<<"stars will collide (Line"<<__LINE__<<")\n";
            break;
        }

        RK4(step/2,y,GM,safeDistance,y_h_2);

        RK4(step/2,y_h_2,GM,safeDistance,y_h_2);

        if(isErrorTolerantable(y_h,y_h_2)) {
            //error is tolerantable, scale up until next value is untolerantable
            while(true) {
                step/=searchRatio;
                //std::cerr<<"step increasing to "<<step<<std::endl;
                RK4(step,y,GM,safeDistance,y_h);
                RK4(step/2,y,GM,safeDistance,y_h_2);
                RK4(step/2,y_h_2,GM,safeDistance,y_h_2);
                if(!isErrorTolerantable(y_h,y_h_2,precision)) {
                    step*=searchRatio;
                    RK4(step,y,GM,safeDistance,y_h);
                    RK4(step/2,y,GM,safeDistance,y_h_2);
                    RK4(step/2,y_h_2,GM,safeDistance,y_h_2);
                    break;
                }
            }
        } else {
            while(true) {
                step*=searchRatio;
                //std::cerr<<"step shrinking to "<<step<<std::endl;
                if(step<=minStep) {
                    if(noCollide!=nullptr) {
                        *noCollide=false;
                    }
                    //std::cerr<<"stars will collide\n";
                    break;
                }
                RK4(step,y,GM,safeDistance,y_h);
                RK4(step/2,y,GM,safeDistance,y_h_2);
                RK4(step/2,y_h_2,GM,safeDistance,y_h_2);
                if(isErrorTolerantable(y_h,y_h_2,precision)) {
                    break;
                }
            }
        }
        //std::cerr<<"step accepted ="<<step<<std::endl;
        Statue yh2_error;
        yh2_error.first=(y_h_2.first-y_h.first)/ratio;
        yh2_error.second=(y_h_2.second-y_h.second)/ratio;

        y.first=y_h_2.first+yh2_error.first;
        y.second=y_h_2.second+yh2_error.second;
        curTime+=step;

    }

    if(!record)
        sol.emplace_back(std::make_pair(curTime,y));

}

    void clear();

    void setMass(const BodyVector &);

    const BodyVector & getMass() const;

    const std::list<Point> & getResult() const;

    double calculateKinetic(const Statue & it) const;

    double calculatePotential(const Statue & it) const;

    double calculateEnergy(const Statue & it) const;

    void calculateTotalMotion(const Statue & it,
                                                DimVector &) const;

    void savePath(const char * fileName);

public:
    static void calculateSafeDistance(const BodyVector &,
                                                            DistanceMat & dest);

    static void calculateGM(const BodyVector &,
                                                            Interaction &);

    static bool calculateDiff(const Position & y,
                             const Interaction &,
                             const DistanceMat &,
                             Acceleration & dy);
    //To avoid useless deep copying, dy.first will not be used.

    static bool RK4(const Time h,
                    const Statue & y,
                    const Interaction &,
                    const DistanceMat& ,
                    Statue & y_n1);

    static bool isErrorTolerantable(const Statue & y_h,
                                    const Statue & y_h_2,
                                    double errorRatio=1e-8);

    static void deval(const Simulator * source,
            Simulator * dest,
                          const Eigen::ArrayXd & timeQueried);
    static void motionAlign(const BodyVector & mass, Velocity & velocity);
    static void positonAlign(Position &);
    static void saveParameters(const char * fileName,
                               const BodyVector & mass,
                               const Statue & y0,
                               const TimeSpan ts,
                               const double step);
    static bool loadParameters(const char * fileName,
                               BodyVector & mass,
                               Statue & y0,
                               TimeSpan & ts,
                               double & step);
    void saveAsData(const char * fileName) const;
    bool loadFromData(const char * fileName);

#ifdef BODY3_DIM3
    static void positionShrink(Position &);
#endif
private:
    std::list<Point> sol;

    BodyVector mass;


};

#endif // SIMULATOR_H
