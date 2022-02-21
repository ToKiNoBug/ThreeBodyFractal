#include <mex.h>
#include "defines.h"
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

const size_t Rows=DIM_COUNT*BODY_COUNT;

void mexFunction(
    int           outC,           /* number of expected outputs */
    mxArray       *outV[],        /* array of pointers to output arguments */
    int           inC,           /* number of inputs */
    const mxArray *inV[]         /* array of pointers to input arguments */
    )
{
    if(inC<=0) {
        mexPrintf("%s\n",
"Format:\n \
[noCollide,lastTime]=threeBodyFast(MassVec,BegPos,BegVelocity,tSpan);\n"
        );
        return;
    }

    if(outC<=1) {
        return;
    }

    if(inC<4) {
        mexErrMsgTxt("Too few inputs!");
        return;
    }

    if(inC>5) {
        mexErrMsgTxt("Too much inputs!");
        return;
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
    s.simulateRK4Var1<false>(1e-4*year,ts,Statue(begPos*rs,begVelocity*vs),&noCollide,1e-4);

    //mexPrintf("%s%d\n","size of result list = ",s.getResult().size());

    if(inC>=1) {
        outV[0]=mxCreateLogicalScalar(noCollide);
    }

    if(inC>=2) {
        outV[1]=mxCreateDoubleScalar(std::min(s.getResult().back().first,ts.second)/year);
    }
}