% in this file we compute analytical solution and compare it to 
% our DG-2D code results.
% clc; clear;

% load("data1.mat");

x = data(:,1);
y = data(:,2);

% DG resutls
p = data(:,3);

% analytical
p_a = exp(-x-y.^2);

figure(1);clf;
plot3(x,y,p, 'x');

figure(2);clf;
plot3(x,y,p_a,'x');

figure(3); clf;
[xsamp, ysamp] = meshgrid(0:0.01:1,0:0.01:1);
p_fine = exp(-xsamp-ysamp.^2);
plot3(xsamp,ysamp,p_fine,'x');

figure(4); clf;
plot(p, p_a, 'x');