format long g
v_n= dlmread('v_n.txt', '', 0, 0);
v_ansys= dlmread('v_ansys100.txt', '', 0, 0);
v_n(:,1)=(v_n(:,1));
plot(v_n(:,2),-v_n(:,1))
hold on
scatter(v_ansys(:,1),v_ansys(:,2),'*')
%title('Validation for Vertical Velocity','fontsize',15)
legend({'FORTRAN Code','ANSYS Fluent'},'Location','northeast')
xlabel({'x/l', '(b)'}) 
ylabel('v/V')