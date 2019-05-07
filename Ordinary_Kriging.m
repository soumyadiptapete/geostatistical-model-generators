clc;close all;clear all;
%% ordinary kriging
NST;
load('Data');
N=length(Data.x);
Data=Data_nst;
sill=1;

%% formulate the left hand side matrix

for i=1:N
    for j=1:N
        Coord1=[Data.x(i) Data.y(i)];
        Coord2=[Data.x(j) Data.y(j)];
        
        % get covariance from function
        cov=vargm(Coord1,Coord2);
        A(i,j)=cov;
    end
end

% modify A for oridnary kriging

A(N+1,N+1)=0;
A(N+1,1:N)=1;
A(1:N,N+1)=-1;


%% calculate grid centers
grid_dimensions

% formulate the right hand matrix for all the grid centers

Ncenter=length(Grid.x);

for i=1:Ncenter
    
    for j=1:N
        Coord1=[Grid.x(i) Grid.y(i)];
        Coord2=[Data.x(j) Data.y(j)];
        
        %get covariance
        cov=vargm(Coord1,Coord2);
        
        B(j,1)=cov;
    end
    B(N+1)=1;
    
    lambda=A\B;
    
    mean=0;
   OK(i)= lambda(1:N)'*Data.lnperm;
   OK_var(i)=sill- lambda'*B -mean;
end
%back transform lnperm
pt=normcdf(OK);
OK=interp1(P,rearrdata,pt);

% transform logk to k
OK=exp(OK);

OK=reshape(OK,Nx,Ny);
OK=OK';
imagesc((OK));
set(gca,'YDir','Normal');
xlabel('East');
ylabel('North');
title('Ordinary Kriging');
colorbar;

figure;
OK_var=reshape(OK_var,Nx,Ny);
OK_var=OK_var';

imagesc((OK_var));
set(gca,'YDir','Normal');
xlabel('East');
ylabel('North');
title('Ordinary Kriging Variance');
colorbar;
   

