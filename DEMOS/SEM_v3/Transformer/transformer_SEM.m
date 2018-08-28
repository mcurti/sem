%% SEM Test Transformer Cartesian
clear variables;
% close all;
clc;

%% SEM

% Load data and prepare for solving
S = SEM('Transformer_V2.xml');

% Solve problem
S = S.solve;

%% Extract results

% Flux in the primary coil [Wb]
phi_p = S.PProcessing.flux_linkage(27);

% Flux in the secondary coil [Wb]
phi_s = S.PProcessing.flux_linkage(15);

% Current density in plate
N_el = 42;
f = 1;
sigma = 8.41e3;
J_e = cell(1,N_el);
for k = 1:N_el
    J_e{k} = real(S.PProcessing.Potential{k}*1j*2*pi*f*sigma);
end

%% Plot results

% Plot the vector potential
figure(1)
clf
hold on;
% Plot points of geometry
S.Geometry.plot_geometry('markersize',12,'color','k','points');
% Plot lines of geometry
S.Geometry.plot_geometry('color','k','lines');
% Plot the vector potential
% S.PProcessing.plot_contour(linspace(0,max(real(S.PProcessing.PHI))));
S.PProcessing.plot_surf_var(J_e);
hold off
c = colorbar;
% caxis([min(real(S.PProcessing.PHI))*0 max(real(S.PProcessing.PHI))])
title(c,'A_z')
xlabel('x [mm]')
ylabel('y [mm]')
daspect([1,1,1e4])