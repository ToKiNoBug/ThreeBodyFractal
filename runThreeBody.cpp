#include <mex.h>
#include "Simulator.h"



inline bool checkSize(const mxArray * mxA,int rows,int cols) {
    if (mxGetM(mxA)==rows)
        return true;
    if (mxGetN(mxA)==cols)
        return true;
    return false;
}

inline bool checkSize(const mxArray * mxA,int Size) {
    if (mxGetNumberOfElements(mxA)==Size)
        return true;
    return false;
}

static const size_t Rows=DIM_COUNT*BODY_COUNT;

void mexFunction(
    int           outC,           /* number of expected outputs */
    mxArray       *outV[],        /* array of pointers to output arguments */
    int           inC,           /* number of inputs */
    const mxArray *inV[]         /* array of pointers to input arguments */
    ) 
{
    if(inC<4) {
        //mexPrintf("%s","Fatal error : Too few inputs!\n");
        mexErrMsgTxt("Too few inputs!");
        return;
    }

    if(inC>5) {
        mexErrMsgTxt("Too much inputs!");
    }

    if(!checkSize(inV[0],BODY_COUNT))
        mexErrMsgTxt("The first input should be a 3 dim vector!");
    Eigen::Map<const BodyVector> mass(mxGetPr(inV[0]));

    if(!checkSize(inV[1],Rows))
        mexErrMsgTxt("The second input should have 6 elements!");
    Eigen::TensorMap<const Position> begPos(mxGetPr(inV[1]),DIM_COUNT,BODY_COUNT);

    if(!checkSize(inV[2],Rows))
        mexErrMsgTxt("The third input should have 6 elements!");
    Eigen::TensorMap<const Velocity> begVelocity(mxGetPr(inV[2]),DIM_COUNT,BODY_COUNT);

    if(!checkSize(inV[3],2))
        mexErrMsgTxt("The forth input should have 2 elements!");
    TimeSpan ts(year**mxGetPr(inV[3]),year**(mxGetPr(inV[3])+1));

    Simulator s;

    s.setMass(mass*Ms);
    
    bool noCollide=true;
    s.simulateRK4Var1(1e-4*year,ts,Statue(begPos*rs,begVelocity*vs),&noCollide);

    mexPrintf("%s","Finished\n");
    
/*
Format:
[Pos,Velocity,TimeQ,noCollide]=runThreeBody(MassVec,BegPos,BegVelocity,tSpan)
[Pos,Velocity,noCollide]=runThreeBody(MassVec,BegPos,BegVelocity,tSpan,TimeQ)
*/

    if(inC==4) {    
        //[Pos,Velocity,TimeQ,noCollide]=runThreeBody(MassVec,BegPos,BegVelocity,tSpan)
        const size_t Cols=s.getResult().size();
        outV[0]=mxCreateDoubleMatrix(Rows,Cols,mxREAL);
        Eigen::Map<Eigen::Array<double,Rows,Eigen::Dynamic>> 
            posDst(mxGetPr(outV[0]),Rows,Cols);
        outV[1]=mxCreateDoubleMatrix(Rows,Cols,mxREAL);
        Eigen::Map<Eigen::Array<double,Rows,Eigen::Dynamic>> 
            velocityDst(mxGetPr(outV[1]),Rows,Cols);
        outV[2]=mxCreateDoubleMatrix(1,Cols,mxREAL);
        double * timeDst=mxGetPr(outV[2]);
        size_t c=0;
        for(const auto & i : s.getResult()) {
            timeDst[c]=i.first/year;
            for(size_t r=0;r<Rows;r++) {
                posDst(r,c)=i.second.first(r)/rs;
                velocityDst(r,c)=i.second.second(r)/vs;
            }
            c++;
        }

        outV[3]=mxCreateLogicalScalar(noCollide);

        return;
    }
    else {

    }
    
}