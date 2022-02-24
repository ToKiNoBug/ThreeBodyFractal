function [center,rSpan,cSpan] = updateWindow(center,rSpan,cSpan,newR,newC,ratio)
%UPDATEWINDOW 此处显示有关此函数的摘要
%   此处显示详细说明
minPos=center-[rSpan,cSpan]/2;
maxPos=center+[rSpan,cSpan]/2;
center(1)=interp1([1,100],[minPos(1),maxPos(1)],newR);
center(2)=interp1([1,100],[minPos(2),maxPos(2)],newC);
rSpan=rSpan/ratio;
cSpan=cSpan/ratio;
end

