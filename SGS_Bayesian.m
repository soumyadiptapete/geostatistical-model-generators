%% SGS with Bayesian updating
NST;%%normal score transform the data
nrlzn=5;%%% number of realizations

for t=1:nrlzn
Data=Data_nst;
Nx=40;
Ny=40;
dx=1;
dy=1;
L=[1:Nx*Ny];
Data_SGS=Data;%%initialize the data structure to which 
%new conditioning data will be added
Data_SGS.varperm=zeros(length(Data.x),1);
sill=1;
rho=0.8; %coefficient of correlation

iter=length(Data.x);

%% start SGS loop

while L %%while there are locations left in grid
    
%%choose 1st location
    pos=randi(length(L));
    index=L(pos);
    L(pos)=[];
       
%% calculate x and y index and co-ordinates of randomly selected grid centre
 if mod(index,Nx)==0
     xindex=Nx;
 else
     xindex=mod(index,Nx);
 end
 
 yindex=((index-xindex)/Ny) + 1;
 
 x_coord=(xindex-1)*dx+dx/2;
 y_coord=(yindex-1)*dy+dy/2;
 
 %% extract the k nearest neighbours
 N=40;
 List=[Data_SGS.x Data_SGS.y];
 Query=[x_coord y_coord];
 idx=knnsearch(List,Query,'k',N);
 idx=idx(5:end);
 N=N-4;
 %% extract a matrix of the known nearest points which will be used for kriging
 Dummy.x=Data_SGS.x(idx);
 Dummy.y=Data_SGS.y(idx);
 Dummy.lnperm=Data_SGS.lnperm(idx);
 Dummy.varperm=Data_SGS.varperm(idx);

%% formulate the left hand side matrix for kriging

for i=1:N
    for j=1:N
        Coord1=[Dummy.x(i) Dummy.y(i)];
        Coord2=[Dummy.x(j) Dummy.y(j)];
        
        % get covariance from function
        cov=vargm(Coord1,Coord2);
        A(i,j)=cov;
    end
end

%% cholesky decomposition of A
LA=chol(A,'lower');

%% calculate global mean
u= mean(Data_SGS.lnperm);

%% calculate grid centers
grid_dimensions

%% formulate the right hand matrix for all the grid centers
    
  for j=1:N
        Coord1=[x_coord y_coord];
        Coord2=[Dummy.x(j) Dummy.y(j)];
        
        %get covariance
        cov=vargm(Coord1,Coord2);
        
        B(j,1)=cov;
  end
  
    %solve Lz=B;
   z=LA\B;
   
   %solve L'x=z;
   lambda=LA'\z;
   
   lambda0=u*(1-sum(lambda));
   SK_est=lambda0+ lambda'*Dummy.lnperm;%%kriging mean
   SK_var=sill- lambda'*B;%% kriging variance
   
   %% bayesian update
   zi=Data_sec_nst(yindex,xindex);
   SK_est=(rho*zi*SK_var+ SK_est*(1-rho^2))/(rho^2*(SK_var-1)+1);%%updated mean
   SK_var=SK_var*(1-rho^2)/(rho^2*(SK_var-1)+1);%%updated variance
   
   %%Random selction of location data value
   p=rand(1); %% choose a random value b/w 0&1;
   Gauss_est=norminv(p,SK_est,SK_var);
   
   
   %% update dataset to include newly calculated data value
   iter=iter+1;
   Data_SGS.x(iter,1)=x_coord;
   Data_SGS.y(iter,1)=y_coord;
   Data_SGS.xindex(iter,1)=xindex;
   Data_SGS.yindex(iter,1)=yindex;
   Data_SGS.lnperm(iter,1)=Gauss_est;
   Data_SGS.varperm(iter,1)=SK_var;
   
   iter-64
end

RLZN=zeros(Ny,Nx);
  for i=65:iter
      %assign the dat values to the location on the grid for a particular
      %realization
      xindex=Data_SGS.xindex(i);
      yindex=Data_SGS.yindex(i);
      lnperm=Data_SGS.lnperm(i);
      
      %back transforming to original data space
      pt=normcdf(lnperm);
      lnperm=interp1(P,rearrdata,pt);
      
      RLZN(yindex,xindex)=lnperm;
  end
  Realization(t).RLZN=exp(RLZN);%% convert to permeability
  disp('realization complete');
end

%% plotting
for i=1:nrlzn
     subplot(3,2,i);
     imagesc((Realization(i).RLZN));
     if i==1
        cl = caxis; %# get color limits from the 1st image
    else
        caxis(cl) %# apply the same color limits to other images
    end
set(gca,'YDir','Normal');
xlabel('East');
ylabel('North');
s=strcat('Realization',num2str(i));
title(s);
colorbar;
 end
    
    

