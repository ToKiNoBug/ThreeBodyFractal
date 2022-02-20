clear

Rows=100;
Cols=100;
load('BegStatue00.mat')
%global Ms1 Ms2 Ms3
Ms1=1;
Ms2=logspace(-4,2,Rows)*MassVec(2);

Ms3=logspace(-4,2,Cols)*MassVec(3);
lastTimeMat=zeros(Rows,Cols);
noCollideMat=false(Rows,Cols);


%poolObj=parpool("local",2);
tic


% parfor idx=1:Rows*Cols
%     [r,c]=ind2sub([Rows,Cols],idx);
%     [noCollideMat(idx),lastTimeMat(idx)]=runThreeBody([Ms1,Ms2(r),Ms3(c)],BegPos,BegVelocity,tSpan); 
% end

parfor r=1:Rows
   for c=1:Cols
       temp_Ms3=logspace(-4,2,Cols)*MassVec(3);
       [noCollideMat(r,c),lastTimeMat(r,c)]=runThreeBody([Ms1,Ms2(r),temp_Ms3(c)],BegPos,BegVelocity,tSpan); 
   end
   disp(num2str(r))
end


toc
imagesc(lastTimeMat);
