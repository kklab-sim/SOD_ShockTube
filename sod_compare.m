%A demo script file to show the use of analytic_sod.m

time = 0.2;
dt = 0.001;
data = analytic_sod(time);
a = load('SOD_RESULTS.dat');
x_node_left = a(:,1);
x_node_right = a(:,2);
x_cell = a(:,3);
rho = a(:,4);
u = a(:,5);
p = a(:,6);
T = a(:,7);
e = a(:,8);
h = a(:,9);
Ncell = length(a(:,1));

figure,
str = strcat({'Sod results for TSIM= '},num2str(time),' s, Ncell= ',num2str(Ncell),'and DT=',num2str(dt),' s.');

subplot(2,2,1),
plot(data.x,data.rho,'-b',x_cell,rho,'-r','LineWidth',1);
xlabel('x (m)');
ylabel('Density (kg/m^3)');
title('Plot of Density vs Position');
grid on;

subplot(2,2,2),
plot(data.x,data.P,'-b',x_cell,p,'-r','LineWidth',1);
xlabel('x (m)');
ylabel('Pressure (Pa)');
title('Plot of Pressure vs Position');
grid on;

subplot(2,2,3),
plot(data.x,data.u,'-b',x_cell,u,'-r','LineWidth',1);
xlabel('x (m)');
ylabel('Velocity (m/s)');
title('Plot of Velocity vs Position');
grid on;

subplot(2,2,4),
plot(data.x,data.e,'-b',x_cell,p./(1.4-1)./rho,'-r','LineWidth',1);
xlabel('x (m)');
ylabel('Specific Internal Energy (J/kg)');
title('Plot of Internal Energy vs Position');
grid on;

suptitle(str);