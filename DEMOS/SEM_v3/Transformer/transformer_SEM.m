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
N_el = 49;                                                                          % Number of elements
sigma_r = [6,12,18,24,30,36,42,43,44,45,46,47,48,49];                               % Elements having non-zero conductivity
f = 1;                                                                              % Frequency [Hz]
sigma_1 = 0;                                                                        % zero conductivity
sigma_2 = 8.41e3;                                                                   % Non-zero conductivity [S/mm]
J_e = cell(1,N_el);
for k = 1:N_el
    if any(k == sigma_r)
        J_e{k} = real(S.PProcessing.Potential{k}*1j*2*pi*f*sigma_2);
    else
        J_e{k} = real(S.PProcessing.Potential{k}*1j*2*pi*f*sigma_1);
    end
end

%% Plot results

% Plot the vector potential
figure(1)
hold on;
% Plot points of geometry
S.Geometry.plot_geometry('markersize',12,'color','k','points');
% Plot lines of geometry
S.Geometry.plot_geometry('color','k','lines');
% Plot the vector potential
S.PProcessing.plot_contour(linspace(0,max(real(S.PProcessing.PHI))));
hold off;
c = colorbar;
caxis([min(real(S.PProcessing.PHI))*0 max(real(S.PProcessing.PHI))])
title(c,'A_z')
xlabel('x [mm]')
ylabel('y [mm]')
axis equal

% Plot the current density
figure(2);
hold on;
% Plot points of geometry
S.Geometry.plot_geometry('markersize',12,'color','k','points');
% Plot lines of geometry
S.Geometry.plot_geometry('color','k','lines');
% Plot the current density in the geometry
S.PProcessing.plot_surf_var(J_e);
hold off;
c = colorbar;
title(c,'J_z')
xlabel('x [mm]')
ylabel('y [mm]')
axis equal


