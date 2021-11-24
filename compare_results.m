
clc;clear;
Newmark=load('Newmark.mat');
Runge=load('Runge_Kutta.mat');


subplot(3,1,1)
plot(Newmark.u(1,:)); hold on;
plot(Runge.u(1,:)); hold on;
ylabel('Displacement')
title('Newmark vs Runge Kutta')

subplot(3,1,2)
plot(Newmark.v(1,:)); hold on;
plot(Runge.v(1,:)); hold on;
ylabel('Velocity')

subplot(3,1,3)
plot(Newmark.a(1,:)); hold on;
plot(Runge.a(1,:)); hold on;
ylabel('Acceleration')
xlabel('Time');




