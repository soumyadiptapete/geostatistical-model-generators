clear all; clc; close all;
%% SGS with bayesian updating
load('Data');
load('Data_sec');

%% normal score transform the primary and secondary data

[D,I]=sort(Data.lnperm);
N=length(Data.x);
for i=1:N
    P(i)=(i-0.5)/N;
end

D_nst_dummy=norminv(P,0,1);
D_nst=zeros(N,1);
for i=1:N
    D_nst(I(i))=D_nst_dummy(i);
end
Data_nst.x=Data.x;
Data_nst.y=Data.y;
Data_nst.lnperm=D_nst;

%%% back tranforming
pt=normcdf(D_nst);
rearrdata=Data.lnperm(I);
kh=interp1(P,rearrdata,pt);


[D,I]=sort(Data_sec);
N=length(Data_sec);
for i=1:N
    P1(i)=(i-0.5)/N;
end

D_nst_dummy=norminv(P1,0,1);
D_nst=zeros(N,1);

for i=1:N
    D_nst(I(i))=D_nst_dummy(i);
end

Data_sec_nst=D_nst;
Data_sec_nst=reshape(Data_sec_nst,40,40);
Data_sec_nst=Data_sec_nst';
