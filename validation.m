format long g
u_n= dlmread('u_n.txt', '', 0, 0);
u_ansys= dlmread('u_ansys.txt', '', 0, 0);
u_n(:,2)=flip(u_n(:,2));
plot(u_n(:,1),u_n(:,2))
hold on
scatter(u_ansys(:,1),u_ansys(:,2),'*')
%title('Validation for Horizontal Velocity','fontsize',15)
legend({'FORTRAN Code','ANSYS Fluent'},'Location','northwest')
xlabel({'x/l','(a)'}) 
ylabel('u/U')

%a(i,j);
%contour(data, levels, 'k')
%imagesc(data)
%title('Streamlines','fontsize',30)
%colorbar
%x = 1:50;
%y = 1:50;
%[X,Y] = meshgrid(x,y);
%vx = U;
%vy = V;
%figure
%pcolor(X,Y,hypot(vx,vy))
%shading interp
%N = 20;
%xstart = 0.1:0.1:1; 
%ystart = ones(size(startx));
%h=streamline(X,Y,vx,vy,xstart,ystart);
%set(h,'color','red')