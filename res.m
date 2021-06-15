clear all
res= dlmread('res.txt', '', 0, 0);
res(:,1) = flip(res(:,1));
%semilogx(res(:,1),res(:,3));
%semilogx(res(:,1),res(:,2));
plot(res(:,1),res(:,3))

hold on

title('Comparison Between LU Decomposition and TDMA','fontsize',15)
%legend({'y = sin(x)','y = cos(x)'},'Location','northeast')