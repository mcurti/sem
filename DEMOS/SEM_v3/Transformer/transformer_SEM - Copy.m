%% SEM Test Transformer Cartesian
clear variables;
close all;
clc;

%% SEM

% Load data and prepare for solving
S = SEM('Transformer.xml');

% Solve problem
S = S.solve;

%% Extract results

% Flux in the primary coil [Wb]
phi_p = S.PProcessing.flux_linkage(27);

% Flux in the secondary coil [Wb]
phi_s = S.PProcessing.flux_linkage(15);

%% Plot results

% Plot the vector potential
figure(1)
hold on;
% Plot points of geometry
S.Geometry.plot_geometry('markersize',12,'color','k','points');
% Plot lines of geometry
S.Geometry.plot_geometry('color','k','lines');
% Plot the vector potential
S.PProcessing.plot_contour();
hold off
c = colorbar;
caxis([min(S.PProcessing.PHI) max(S.PProcessing.PHI)])
title(c,'A_z')
xlabel('x [mm]')
ylabel('y [mm]')
axis equal