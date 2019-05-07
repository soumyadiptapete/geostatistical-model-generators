%% grid is 40X40; each cell is 1x1
% calculate the center of each grid cell
Nx=40;
Ny=40;
dx=1;
dy=1;
a=0;
for i=1:Ny
    for j=1:Nx
        a=a+1;
        Grid.x(a)=(j-1)*dx+dx/2;
        Grid.y(a)=(i-1)*dy+dy/2;
    end
end