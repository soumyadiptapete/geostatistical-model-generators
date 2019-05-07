close all;
Mean=zeros(40,40);
for i=1:5
    Mean=Mean+Realization(i).RLZN;
end
Mean=Mean/5;
A=(Realization(1).RLZN-Mean);
B=(Realization(2).RLZN-Mean);
C=(Realization(3).RLZN-Mean);
D=(Realization(4).RLZN-Mean);
E=(Realization(5).RLZN-Mean);
A=A.^2;B=B.^2;C=C.^2;D=D.^2;E=E.^2;
VAR=(A+B+C+D+E)/5;

subplot(1,2,1)
imagesc((Mean));
set(gca,'YDir','Normal');
xlabel('East');
ylabel('North');
title('Mean of SGSCOSIM');
colorbar;

subplot(1,2,2);
imagesc((VAR));
set(gca,'YDir','Normal');
xlabel('East');
ylabel('North');
title('Variance of SGSCOSIM');
colorbar;


