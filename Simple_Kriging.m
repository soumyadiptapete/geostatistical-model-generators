clear all; close all; clc;
%%% Simple kriging
NST;
load('Data');
N=length(Data.x);
sill=1;

Data=Data_nst;

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

%% cholesky factorization of A

L=chol(A,'lower');

%% calculate global mean
u= 0;

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
    
    
   %solve Lz=B;
   z=L\B;
   
   %solve L'x=z;
   lambda=L'\z;
   
   
   u=0;
   lambda0=u*(1-sum(lambda));
   SK(i)=lambda0+ lambda'*Data.lnperm;
   SK_var(i)=sill- lambda'*B;
   
end
%back transform lnperm
pt=normcdf(SK);
SK=interp1(P,rearrdata,pt);

% transform logk to k
SK=exp(SK);

SK=reshape(SK,Nx,Ny);
SK=SK';

imagesc((SK));
set(gca,'YDir','Normal');
xlabel('East');
ylabel('North');
title('Simple Kriging');
colorbar;

figure;
SK_var=reshape(SK_var,Nx,Ny);
SK_var=SK_var';

imagesc((SK_var));
set(gca,'YDir','Normal');
xlabel('East');
ylabel('North');
title('Simple Kriging Variance');
colorbar;
   
    
        

        
