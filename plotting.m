data = dlmread('data-psi-tdma-30.txt', '', 0, 0);
a=max(max(data)); b=min(min(data))
c = 0.000000035738:100:0.000000070569;
x = 0.00000003:0.000000007;

contourf(data, x)