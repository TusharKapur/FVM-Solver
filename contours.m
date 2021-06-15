data = dlmread('data-psi-tdma.txt', '', 0, 0);
figure
[C,h] = contour(data, 20, 'k');
set(gca,'XTick',[], 'YTick', [])
%clabel(C,h,'FontSize',30,'labelspacing', 700)
%imagesc(data)
title('Re = 10','fontsize',12)
%colorbar
axis ij
xlabel('(a)')