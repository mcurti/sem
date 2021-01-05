clear; clc;

G  = Geometry('ISEF_corner.xml');
Ph = Physics(G);

Ph = Ph.load_magnetic_materials;
Ph = Ph.load_current_density;

tic
Pr = Problem(G,Ph);
Pr = Pr.load_Y_sources;
Pr = Pr.building_Y_vector;
Pr = Pr.global_matrix;

PHI = sparse(Pr.Global_Matrix)\Pr.Y.vector';

Pp = PostProcessing(G, Pr, PHI);
Pp = Pp.compute_B;
toc
%

% Plot the results

mpy = @(s) s*1e3*pi*4e-7;
%%
% Converting to Wb/mm
A = cellfun(mpy,Pp.Potential,'UniformOutput',false);

figure(1)
clf
hold on
G.plot_geometry('markersize',12,'color','k','points')
G.plot_geometry('color','k','lines')
Pp.plot_contour_var(A)
c = colorbar;
hold off
title(c,'Az, Wb/mm')
xlabel('x, mm')
ylabel('y, mm')
xlim([0,40])
daspect([10 10 2])
%%
% Converting to tesla
B = cellfun(mpy,Pp.Flux.abs,'UniformOutput',false);

figure(2)
clf
hold on
Pp.plot_surf_var(B)
hold off
c = colorbar;
title(c,'B, T')
xlabel('x, mm')
ylabel('y, mm')
xlim([0,40])
daspect([10 10 1])
