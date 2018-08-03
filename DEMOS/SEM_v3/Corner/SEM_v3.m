%==========================================================================
% SEM_v3.m
% Created: 29.10.2017 - 22:12:30
% By: M. Curti
%
% The third version of the SEM demo file where the global class SEM is used
%==========================================================================

clear; clc;

S = SEM('ISEF_corner.xml');
NonlinearSolver.tol = 1e-4;
NonlinearSolver.max_iter = 30;
NonlinearSolver.display = 'on';
set(S,'NonlinearSolver',NonlinearSolver);
% Preparing the current density vector
% J = linspace(.1,8,40);
J  = [1e4 3e4 5e4 1e5 2e5 4e5 6e5 8e5 1e6]*1e-6;
I  = [0.5 1.5 2.5 5 10 20 30 40 50];
Flux = zeros(1,length(J)); remFlux = zeros(1,length(J));
Flux_inc = Flux; L_inc = Flux; L = Flux;
for k = 1:length(J)
    new_J = sprintf('%g',J(k)*pi*4e-7);
    S.Physics = S.Physics.setParameter('J',new_J);
    S.Physics = S.Physics.load_current_density;
    
    S.Problem = S.Problem.updateSources(S.Physics.Sources);
    
    S = S.solve;
    
    %% Post processing
    [Flux(k), remFlux(k)] = S.PProcessing.flux_linkage([12, 14]);
    
%     Flux_inc(k) = (Flux(k) - remFlux(k))/J(k);
    L(k) = Flux(k)*1e-3/I(k);
    L_inc(k) = (Flux(k) - remFlux(k))/I(k)*1e-3;
end
% Flux = Flux*1e-3; remFlux = remFlux*1e-3; L_inc = L_inc*1e-3;

%% Plot the results
SB  = spline(I,Flux*1e-3);


f   = @(x) ppval(SB,x);
df  = derivest(f,I);

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
plot(I,Flux)
plot(I,Flux,'ok')
% plot(J*0,remFlux,'*r')
% plot(I,L_inc,'o-k')
hold off
xlabel('Current, A')
ylabel('Flux linkage, Wb')

figure(3)
clf
hold on
% plot(I,Flux)
% plot(I,Flux,'ok')
% plot(J*0,remFlux,'*r')
plot(I,L_inc,'o-k')
plot(I,df,'o-g')
% plot(I,df,'o-g')
hold off
xlabel('Current, A')
ylabel('Inductance, Wb')