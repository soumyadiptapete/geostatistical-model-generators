function [ cov ] = vargm_nst( Coord1,Coord2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

theta=2*pi*315/360;
R=[cos(theta) -sin(theta);sin(theta) cos(theta)];
Coord1=R*Coord1';
Coord2=R*Coord2';

sill1=0.6;
a1=16;
h1=sqrt((Coord1(1)-Coord2(1))^2+(Coord1(2)-Coord2(2))^2);
if h1<=a1
gamma1=sill1*(1.5*(h1/a1)-0.5*(h1/a1)^3);
else
    gamma1=sill1;
end

sill2=0.9;
a2=12;
h2=abs((Coord1(2)-Coord2(2)));
if h2<=a2
gamma2=(sill2-sill1)*(1.5*(h2/a2)-0.5*(h2/a2)^3);
else
    gamma2=(sill2-sill1);
end

gamma=gamma1+gamma2;
cov=sill2-gamma;
end





