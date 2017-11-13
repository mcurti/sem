G  = Geometry('LPMSM.xml');


figure(1)
clf
hold on
G.plot_geometry('markersize',12,'color','k','points')
G.plot_geometry('color','k','lines')
G.plot_geometry('color','k','markersize',8,'ElementGrid')
axis equal

%%
clc; clear all;
S = SEM('LPMSM.xml');
S = S.solve;
figure(1)
clf
hold on
S.Geometry.plot_geometry('markersize',12,'color','k','points')
S.Geometry.plot_geometry('color','k','lines')
% S.Geometry.plot_geometry('color','k','markersize',8,'ElementGrid')
colorbar
axis equal

S.PProcessing.plot_contour


