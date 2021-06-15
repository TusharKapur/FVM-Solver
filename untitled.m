format long g
Z = peaks(100); 
U = dlmread('data-u-tdma.txt', '', 0, 0);
V = dlmread('data-v-tdma.txt', '', 0, 0);
%a=max(max(data)); b=min(min(data));
%a(i,j);
%levels = prctile(data(:), linspace(0, 10, 10));
%contour(data, levels, 'k')
%imagesc(data)
%title('Streamlines','fontsize',30)
%colorbar
%axis ij
contour(U)
x = 1:50;
y = 1:50;
[X,Y] = meshgrid(x,y);
vx = U;
vy = V;
figure
%pcolor(X,Y,hypot(vx,vy))
shading interp
set ( gca, 'ydir', 'reverse' )