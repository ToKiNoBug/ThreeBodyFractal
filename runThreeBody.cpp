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

const size_t Rows=DIM_COUNT*BODY_COUNT;

void storeResult(const Simulator & src,
    const bool noCollide,
    mxArray ** pos,
    mxArray ** velocity,
    mxArray ** timeQ,
    mxArray ** noCollide_dst) {


    const size_t Cols=src.getResult().size();
    *pos=mxCreateDoubleMatrix(Rows,Cols,mxREAL);
    *velocity=mxCreateDoubleMatrix(Rows,Cols,mxREAL);
    *timeQ=mxCreateDoubleMatrix(1,Cols,mxREAL);
        
    Eigen::Map<Eigen::Array<double,Rows,Eigen::Dynamic>> 
        posDst(mxGetPr(*pos),Rows,Cols);
    Eigen::Map<Eigen::Array<double,Rows,Eigen::Dynamic>> 
        velocityDst(mxGetPr(*velocity),Rows,Cols);
    double * timeDst=mxGetPr(*timeQ);
    size_t c=0;
    for(const auto & i : src.getResult()) {
        timeDst[c]=i.first/year;
        for(size_t r=0;r<Rows;r++) {
            posDst(r,c)=i.second.first(r)/rs;
            velocityDst(r,c)=i.second.second(r)/vs;
        }
        c++;
    }

    *noCollide_dst=mxCreateLogicalScalar(noCollide);
}

void storeResult(const Simulator & src,
    const bool noCollide,
    mxArray ** pos,
    mxArray ** velocity,
    mxArray ** noCollide_dst) {


    const size_t Cols=src.getResult().size();
    *pos=mxCreateDoubleMatrix(Rows,Cols,mxREAL);
    *velocity=mxCreateDoubleMatrix(Rows,Cols,mxREAL);
        
    Eigen::Map<Eigen::Array<double,Rows,Eigen::Dynamic>> 
        posDst(mxGetPr(*pos),Rows,Cols);
    Eigen::Map<Eigen::Array<double,Rows,Eigen::Dynamic>> 
        velocityDst(mxGetPr(*velocity),Rows,Cols);
    size_t c=0;
    for(const auto & i : src.getResult()) {
        for(size_t r=0;r<Rows;r++) {
            posDst(r,c)=i.second.first(r)/rs;
            velocityDst(r,c)=i.second.second(r)/vs;
        }
        c++;
    }

    *noCollide_dst=mxCreateLogicalScalar(noCollide);
}

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
[Pos,Velocity,TimeQ,noCollide]=runThreeBody(MassVec,BegPos,BegVelocity,tSpan);\n \
[Pos,Velocity,noCollide]=runThreeBody(MassVec,BegPos,BegVelocity,tSpan,TimeQ);\n \
[noCollide]=runThreeBody(___);\n \
[noCollide,lastTime]=runThreeBody(___);\n"

        );
        return;
    }

    if(outC<=0) {
        return;
    }

    if(inC<4) {
        //mexPrintf("%s","Fatal error : Too few inputs!\n");
        mexErrMsgTxt("Too few inputs!");
        return;
    }

    if(inC>5) {
        mexErrMsgTxt("Too much inputs!");
        return;
    }
    
    std::clock_t c=std::clock();
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
    s.simulateRK4Var1<true>(1e-4*year,ts,Statue(begPos*rs,begVelocity*vs),&noCollide,1e-8);

    if(std::clock()-c>=2*CLOCKS_PER_SEC)
        mexPrintf("%s","Finished\n");
    

    if(outC<=2) {
        if(outC==1) {
            outV[0]=mxCreateLogicalScalar(noCollide);
            return;
        }

        if(outC==2) {
            outV[0]=mxCreateLogicalScalar(noCollide);
            outV[1]=mxCreateDoubleScalar(std::min(s.getResult().back().first,ts.second)/year);
            return;
        }
    }

    if(inC==4) {    
        //[Pos,Velocity,TimeQ,noCollide]=runThreeBody(MassVec,BegPos,BegVelocity,tSpan)
        storeResult(s,noCollide,outV+0,outV+1,outV+2,outV+3);
        return;
    }
    else {
        //[Pos,Velocity,noCollide]=runThreeBody(MassVec,BegPos,BegVelocity,tSpan,TimeQ)
        const size_t Cols=mxGetNumberOfElements(inV[4]);
        
        mxAssert(mxGetM(inV[4])==1||mxGetN(inV[4])==1,"You should input a vector as timeQ");
        Simulator res;
        Eigen::Map<const Eigen::ArrayXd> tQ(mxGetPr(inV[4]),Cols);

        Simulator::deval(&s,&res,tQ*year);
        storeResult(res,noCollide,outV+0,outV+1,outV+2);
        return;
    }
    
}