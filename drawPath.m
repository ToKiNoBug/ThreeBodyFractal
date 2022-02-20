function [] = drawPath(Pos,varargin)

BodyN=3;
%DimN=2;

if nargin==2
     TimeQ=varargin{1};
end
for b=1:BodyN
    if nargin~=2
        plot(Pos(b,:),Pos(b+BodyN,:))
    else
        plot3(TimeQ,Pos(b,:),Pos(b+BodyN,:))
    end
    hold on
end

if nargin~=2
    xlabel('x')
    ylabel('y')
else
    xlabel('t')
    ylabel('x')
    zlabel('y')
end

hold off
end

