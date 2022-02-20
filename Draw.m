Rows=100;
Cols=100;
load('BegStatue00.mat')

Ms1=1;
Ms2=logspace(-2,2,Rows)*MassVec(2);
Ms3=logspace(-2,2,Cols)*MassVec(3);

lastTimeMat=zeros(Rows,Cols);
noCollideMat=false(Rows,Cols);

for r=1:Rows
   for c=1:Cols
      [noCollideMat(r,c),lastTimeMat(r,c)]=runThreeBody([Ms1,Ms2(r),Ms3(c)],BegPos,BegVelocity,tSpan); 
   end
   disp(num2str(r))
end


imagesc(lastTimeMat);
