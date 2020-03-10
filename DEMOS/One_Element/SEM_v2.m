%==========================================================================
% SEM_v2.m
% Created: 26.10.2017 - 18:09:30
% By: M. Curti
%
% The demo file which can be used to start your own SEM model
%==========================================================================
clear; clc;

G  = Geometry('quarter_motor.xml');
Ph = Physics(G);
% 
bh = @(H) H+4e-3/pi.*atan((pi*(1000-1)*H)./4e-3);
Ph = Ph.load_magnetic_materials;
Ph = Ph.load_current_density;

dataProblem.ProblemData = Ph.ProblemData;
dataProblem.Sources     = Ph.Sources;
dataProblem.metrics     = G.metrics;
dataProblem.xmlContent  = G.GeometryElement;
dataProblem.Materials   = Ph.Materials;
tic
Pr = Problem(dataProblem);




Pr = Pr.load_Y_sources;
Pr = Pr.building_Y_vector;
Pr = Pr.global_matrix;

PHI = Pr.Global_Matrix\Pr.Y.vector';
%% pause
dataPostProcessing.xmlContent   = G.GeometryElement;
dataPostProcessing.PHI          = PHI;
dataPostProcessing.ProblemData  = Ph.ProblemData;
dataPostProcessing.mappings     = G.mappings;
dataPostProcessing.metrics      = G.metrics;
dataPostProcessing.lines        = G.lines;
dataPostProcessing.points       = G.points;
dataPostProcessing.xi           = G.xi;
dataPostProcessing.elements     = G.elements;
dataPostProcessing.Permeability = Pr.Materials.Permeability;

%



Pp = PostProcessing(dataPostProcessing);
Pp = Pp.compute_B;
%     new_source.CurrentDensity = Ph.Sources.CurrentDensity;



% Plot the results
              
%%
figure(1)
clf
hold on
G.plot_geometry('markersize',12,'color','k','points')
G.plot_geometry('color','k','lines')
% G.plot_geometry('color','k','markersize',8,'ElementGrid')
Pp.plot_contour

Pp.plot_surf_var(Pp.Flux.abs)

% c = colorbar;
% Pp.plot_surf
% Pp.plot_trisurf(dataFlux)
hold off
% caxis([min(PHI) max(PHI)]*pi*4e-7)
% figure_config(1,8.2,6.8,8)
% title(c,'Az')
xlabel('x, mm')
ylabel('y, mm')

daspect([10 10 2])
figure_config(1,10,10,8)




figure(2)
clf
hold on
Pp.plot_surf_var(Pp.Flux.abs)
hold off
daspect([40 40 max(Pp.Flux.abs{1}(:))])
colorbar

figure_config(2,10,10,8)
