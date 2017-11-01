%==========================================================================
% SEM_v3.m
% Created: 29.10.2017 - 22:12:30
% By: M. Curti
%
% The third version of the SEM demo file where the global class SEM is used
%==========================================================================

clear; clc;

S = SEM('ISEF_corner.xml');

% Preparing the current density vector
J = linspace(.1,7,10)*20;
Flux = zeros(1,length(J)); remFlux = zeros(1,length(J));

for k = 1:length(J)
    new_J = sprintf('%.6f',J(k)*pi*4e-7);
    S.Physics = S.Physics.setParameter('J',new_J);
    S.Physics = S.Physics.load_current_density;
    
    S.Problem = S.Problem.updateSources(S.Physics.Sources);
    
    S = S.solve;
    
    %% Post processing
    [Flux(k), remFlux(k)] = S.PProcessing.flux_linkage([12, 14]);
    
end
%% Plot the results

% Plot the vector potential
figure(1)
clf
hold on
S.Geometry.plot_geometry('markersize',12,'color','k','points')
S.Geometry.plot_geometry('color','k','lines')
% G.plot_geometry('color','k','markersize',8,'ElementGrid')
S.PProcessing.plot_contour
c = colorbar;
% Pp.plot_surf
% Pp.plot_trisurf(dataFlux)
hold off
% caxis([min(PHI) max(PHI)]*pi*4e-7)
% figure_config(1,8.2,6.8,8)
title(c,'Az')
xlabel('x, mm')
ylabel('y, mm')
xlim([0,40])
daspect([10 10 2])

figure_config(2,10,10,8)

% Plot the convergence

figure(2)
clf
hold on
plot(J,Flux)
plot(J,Flux,'ok')
plot(J*0,remFlux,'*r')
hold off