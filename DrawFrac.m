%clear

Rows=100;
Cols=100;
%load('BegStatue00.mat')
%global Ms1 Ms2 Ms3

%center=[-1,-1];
%rSpan=6;
%cSpan=6;

minPos=center-[rSpan,cSpan]/2;
maxPos=center+[rSpan,cSpan]/2;

Ms1=1;
Ms2=logspace(minPos(1),maxPos(2),Rows)*MassVec(2);

Ms3=logspace(minPos(1),maxPos(2),Cols)*MassVec(3);

clear minPos
clear maxPos

lastTimeMat=zeros(Rows,Cols);
noCollideMat=false(Rows,Cols);


%poolObj=parpool("local",2);
tic


% parfor idx=1:Rows*Cols
%     [r,c]=ind2sub([Rows,Cols],idx);
%     [noCollideMat(idx),lastTimeMat(idx)]=runThreeBody([Ms1,Ms2(r),Ms3(c)],BegPos,BegVelocity,tSpan); 
% end

parfor r=1:Rows
   temp_Ms3=Ms3;
   for c=1:Cols
       [noCollideMat(r,c),lastTimeMat(r,c)]=threeBodyFast([Ms1,Ms2(r),temp_Ms3(c)],BegPos,BegVelocity,tSpan,1e-6); 
   end
   %disp(num2str(r))
end

toc

clear temp_Ms3
clear Ms3
clear Ms2
clear Ms1
imagesc(lastTimeMat);
