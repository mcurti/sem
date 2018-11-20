%% SEM Test Transformer Cartesian
clear variables;
close all;
clc;

%% SEM
tic
% Load data and prepare for solving
S = SEM('Transformer_V2.xml');

% Solve problem
S = S.solve;

% Depth of the domain [m]
D = 10e-3;                                                                          

%% Extract results

% Flux in the primary coil [Wb]
phi_p = S.PProcessing.flux_linkage(27);

% Flux in the secondary coil [Wb]
phi_s = S.PProcessing.flux_linkage(15);

% Current density in plate
N_el = length(S.Geometry.elements.lines);                                           % Number of elements
sigma_r = [6,12,18,24,30,36,42];                                                    % Elements having non-zero conductivity
f = abs(str2num(cell2mat(S.Physics.parameters(3,2)))/(4e-7*pi*2*pi));               % Frequency [Hz]
sigma_1 = 0;                                                                        % zero conductivity
% Non-zero conductivity [S/mm]
str_1 = strsplit(cell2mat(S.Physics.parameters(8,2)),'*');
sigma_2 = str2num(cell2mat(str_1(1)));

% Skin depth
delta_s = sqrt(1/(pi*f*4e-7*pi*100*sigma_2));

J_e = cell(1,N_el);
J_int = cell(1,N_el);
J_int_2 = cell(1,N_el);
P_e = zeros(1,N_el);
for k = 1:N_el
    if any(k == sigma_r)
        J_e{k} = real(S.PProcessing.Potential{k}*1i*2*pi*f*sigma_2);
        J_int{k} = J_e{k}.*S.Geometry.metrics.W{k}.*S.Geometry.metrics.J{k};
        J_int_2{k} = abs(J_e{k}).*S.Geometry.metrics.W{k}.*S.Geometry.metrics.J{k};
        P_e(k) = sum(sum((1/sigma_2)*abs(S.PProcessing.Potential{k}*1i*2*pi*f*sigma_2./sqrt(2)).^2.*S.Geometry.metrics.W{k}.*S.Geometry.metrics.J{k}));
    else
        J_e{k} = real(S.PProcessing.Potential{k}*1i*2*pi*f*sigma_1);
        J_int{k} = real(S.PProcessing.Potential{k}*1i*2*pi*f*sigma_1);
        J_int_2{k} = real(S.PProcessing.Potential{k}*1i*2*pi*f*sigma_1);
    end
end
toc
P_j_p = sum(sum((cell2mat(J_int) + cell2mat(J_int_2))./2));
P_j_n = sum(sum((-cell2mat(J_int) + cell2mat(J_int_2))./2));
e_j = P_j_n./(P_j_n+P_j_p)*100;

%% Plot results

% Plot the vector potential
figure(1)
hold on;
% Plot points of geometry
S.Geometry.plot_geometry('markersize',12,'color','k','points');
% Plot lines of geometry
S.Geometry.plot_geometry('color','k','lines');
% Plot the vector potential
S.PProcessing.plot_contour(linspace(min(abs(S.PProcessing.PHI)),max(abs(S.PProcessing.PHI))));
hold off;
c = colorbar;
caxis([min(abs(S.PProcessing.PHI)) max(abs(S.PProcessing.PHI))])
title(c,'A_z')
xlabel('x [m]')
ylabel('y [m]')
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
J_c = cell2mat(J_e);
caxis([min(J_c(:)) max(J_c(:))]);
title(c,'J_z');
xlabel('x [m]');
ylabel('y [m]');
%axis equal
xlim([7.65e-3 7.65e-3+3*delta_s]);
ylim([0 18.5e-3]);

%% Plate losses

P_eddy = sum(P_e)*D
P_FEM = 20.5494270709388
e = abs(P_eddy - P_FEM)/P_FEM*100


