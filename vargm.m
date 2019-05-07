function [ cov] = vargm( Coord1,Coord2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
theta=2*pi*45/360;
R=[cos(theta) sin(theta);-sin(theta) cos(theta)];
Coord1=R*Coord1';
Coord2=R*Coord2';

sill=1;
ax=14;
ay=56;

hx=Coord1(1)-Coord2(1);
hy=Coord1(2)-Coord2(2);
h=sqrt(((hx/ax)^2)+((hy/ay)^2));
if h<1
gamma=sill*(1.5*h-0.5*h^3);
else
    gamma=sill;
end

cov=sill-gamma;



% sill=1;
% a1=16;
% h1=sqrt((Coord1(1)-Coord2(1))^2+(Coord1(2)-Coord2(2))^2);
% if h1<=a1
% gamma=sill*(1.5*(h1/a1)-0.5*(h1/a1)^3);
% else
%     gamma=sill;
% end
% 
% cov=sill-gamma;

end

