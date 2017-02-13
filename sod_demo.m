%A demo script file to show the use of analytic_sod.m
time = 0.2;
data = analytic_sod(time);
load output.dat;
figure,

subplot(2,2,1),
plot(data.x,data.rho,'-b','LineWidth',2);
xlabel('x (m)');
ylabel('Density (kg/m^3)');
title('Plot of Density vs Position');
grid on;

subplot(2,2,2),
plot(data.x,data.P,'-b','LineWidth',2);
xlabel('x (m)');
ylabel('Pressure (Pa)');
title('Plot of Pressure vs Position');
grid on;

subplot(2,2,3),
plot(data.x,data.u,'-b','LineWidth',2);
xlabel('x (m)');
ylabel('Velocity (m/s)');
title('Plot of Velocity vs Position');
grid on;

subplot(2,2,4),
plot(data.x,data.e,'-b','LineWidth',2);
xlabel('x (m)');
ylabel('Specific Internal Energy (J/kg)');
title('Plot of Internal Energy vs Position');
grid on;
