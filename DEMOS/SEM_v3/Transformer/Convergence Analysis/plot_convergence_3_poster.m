%% Plot convergence 3
clear variables;
close all;
clc;

% Figure settings (presentation)
w_f = 19;                                                                          % Width [cm]
h_f = 10.0;                                                                          % Height [cm]
l_w = 2.0;                                                                          % Linewidth
l_o = 2.0;                                                                          % Linewidth optimum values
f_s = 14.0;                                                                          % Fontsize
m_sl = 10.0;                                                                         % Markersize for legend
m_sp = 4.0;                                                                         % Markersize
m_so = 10.0;                                                                         % Markersize for optimum values

% Pareto weights
w_n = 0.5;                                                                          % Weight for number of d.o.f. [-]
w_e = 0.5;                                                                          % Weight for error rate [-]

%% FEM data 1

% Import data
load('C:\Users\kbastiae\Documents\GitHub\RepositoryForPDEs\DEMOS\SEM_v3\Transformer\Convergence Analysis\2D_FEM_dof_pot_core_final_V5.mat');

% Total number of nodes and eddy current losses
N_f_1 = N_f_t(:);
P_e_1 = P_e_t(:);
t_c_1 = t_c_t(:);

% Reference values
[Y,I] = max(N_f_1);
P_e_FEM = P_e_1(I);

% Difference
e_FEM_1 = abs(P_e_1-P_e_FEM)./P_e_FEM*100;

%% FEM data 2

% Import data
load('C:\Users\kbastiae\Documents\GitHub\RepositoryForPDEs\DEMOS\SEM_v3\Transformer\Convergence Analysis\2D_FEM_dof_pot_core_final_V4.mat');

% Total number of nodes and eddy current losses (exclude first point)
N_f_t = N_f_t(2:end,:,:,:,:); 
P_e_t = P_e_t(2:end,:,:,:,:);
t_c_t = t_c_t(2:end,:,:,:,:);
N_f_2 = N_f_t(:);
P_e_2 = P_e_t(:);
t_c_2 = t_c_t(:);

% Difference
e_FEM_2 = abs(P_e_2-P_e_FEM)./P_e_FEM*100;

%% Pareto front FEM

% Combine datasets
N_f_F = [N_f_1; N_f_2];
e_FEM = [e_FEM_1; e_FEM_2];
P_e_F = [P_e_1; P_e_2];
t_c_F = [t_c_1; t_c_2];

% Pareto matrix
P = [-N_f_F -e_FEM];
% Calculate front
[p,idxs] = paretoFront(P);
% Return front
P_f_F = [N_f_F(idxs) e_FEM(idxs)];
% Sort matrix
P_f_F = sortrows(P_f_F,1);

% Pareto optimum data set 1
f_n1_F = (N_f_F-min(N_f_F))./(max(N_f_F)-min(N_f_F));                                 % Normalized number of d.o.f. for SEM
f_n2_F = (e_FEM-min(e_FEM))./(max(e_FEM)-min(e_FEM));                                 % Normalized relative error 
P_n = [-f_n1_F -f_n2_F];                                                                % Pareto matrix
[p_n,idxs_n] = paretoFront(P_n);                                                    % Pareto indeces
% Sort matrix
p_n_F = sortrows(abs(p_n),1);

F_FEM = f_n1_F(idxs_n)*w_n + f_n2_F(idxs_n)*w_e;                                        % Pareto front 
[Y,I] = min(F_FEM);                                                                 % Find pareto minium
I_FEM = idxs_n(I);                                                                  % Index of optimum
display(sprintf(['Pareto optimum values for FEM (weights ',num2str(w_n),' d.o.f. ',num2str(w_e),' discrepancy','): \n',num2str(N_f_F(I_FEM)),' d.o.f., ','eddy current losses: ',num2str(P_e_F(I_FEM)),' W ','relative discrepancy: ',num2str(e_FEM(I_FEM)),'%%\n','computation time: ',num2str(t_c_F(I_FEM)),' seconds\n']));

% Find optimum mesh discretization
[Y_1,I_1] = min(abs(e_FEM_1-e_FEM(I_FEM)));
[Y_2,I_2] = min(abs(e_FEM_2-e_FEM(I_FEM)));
if Y_1 == 0
    load('C:\Users\kbastiae\Documents\GitHub\RepositoryForPDEs\DEMOS\SEM_v3\Transformer\Convergence Analysis\2D_FEM_dof_pot_core_final_V5.mat');
else
    load('C:\Users\kbastiae\Documents\GitHub\RepositoryForPDEs\DEMOS\SEM_v3\Transformer\Convergence Analysis\2D_FEM_dof_pot_core_final_V4.mat');
end
I_FEM_2 = find(P_e_t(:) == P_e_F(I_FEM));
[I_1,I_2,I_3,I_4,I_5] = ind2sub(size(P_e_t),I_FEM_2);
display(sprintf(['Optimum mesh discretization for FEM:\nextra small: ',num2str(m_es_1(I_1)*1e3),' mm\n','Small 1: ',num2str(m_es_2(I_2)),' mm\n','Small 2: ',num2str(m_s(I_3)),' mm\n','Medium: ',num2str(m_m(I_4)),'mm\n','Large: ',num2str(m_l(I_5)),' mm\n']));

%% SEM data 

% Import data
load('C:\Users\kbastiae\Documents\GitHub\RepositoryForPDEs\DEMOS\SEM_v3\Transformer\Convergence Analysis\2D_SEM_convergence_V2_final.mat');

% Total number of nodes and losses
N_f_s = N_f_t(:);
P_e_s = P_e_t(:);
t_c_s = t_c_t(:);

% Difference
e_SEM = abs(P_e_s-P_e_FEM)./P_e_FEM*100;

%% Pareto front SEM

% Pareto front
P = [-N_f_s -e_SEM];
[p,idxs] = paretoFront(P);
P_f_S = [N_f_s(idxs) e_SEM(idxs)];
P_f_S = sortrows(P_f_S,1);

% Pareto optimum
f_n1_S = (N_f_s-min(N_f_s))./(max(N_f_s)-min(N_f_s));                                 % Normalized number of d.o.f. for SEM
f_n2_S = (e_SEM-min(e_SEM))./(max(e_SEM)-min(e_SEM));                                 % Normalized relative error 
P_n = [-f_n1_S -f_n2_S];                                                                % Pareto matrix
[p_n,idxs_n] = paretoFront(P_n);                                                    % Pareto indeces
% Sort matrix
p_n_S = sortrows(abs(p_n),1);
F_SEM = f_n1_S(idxs_n)*w_n + f_n2_S(idxs_n)*w_e;                                        % Pareto front 
[Y,I] = min(F_SEM);                                                                 % Find pareto minium
I_SEM = idxs_n(I);                                                                  % Index of optimum
display(sprintf(['Pareto optimum values for SEM (weights ',num2str(w_n),' d.o.f. ',num2str(w_e),' discrepancy','): \n',num2str(N_f_s(I_SEM)),' d.o.f., ','eddy current losses: ',num2str(P_e_s(I_SEM)),' W ','relative discrepancy: ',num2str(e_SEM(I_SEM)),'%%\n','computation time: ',num2str(t_c_s(I_SEM)),' seconds\n']));

% Optimum discretization
[I_1,I_2,I_3] = ind2sub(size(P_e_t),I_SEM);
display(sprintf(['Optimum polynomial degree discretization for SEM:\nTangential 1: ',num2str(N_x(I_1)),' -\n','Normal: ',num2str(N_y(I_2)),' -\n','Tangential 2: ',num2str(N_px(I_3)),' -\n']));

%% Plot results loglog

fig_1 = figure(1);
% FEM pareto front
h_0 = loglog(P_f_F(:,1),P_f_F(:,2),'k-x');
set(h_0,'LineWidth',l_w);
set(h_0,'MarkerSize',8);
hold on;
% FEM data point
h_1 = loglog(NaN,NaN,'.','Color',[192 192 192]./256);
set(h_1,'MarkerSize',m_sl);
% FEM optimum
h_2 = loglog(N_f_F(I_FEM),e_FEM(I_FEM),'o','Color',[34 139 34]./256);
set(h_2,'MarkerSize',m_so);
set(h_2,'LineWidth',l_o);
% SEM pareto front
h_3 = loglog(P_f_S(:,1),P_f_S(:,2),'Marker','x','LineStyle','-','Color',[200 25 25]./256);
set(h_3,'LineWidth',l_w);
set(h_3,'MarkerSize',8);
% SEM data point
h_4 = loglog(NaN,NaN,'.','Color',[105 105 105]./256);
set(h_4,'MarkerSize',m_sl);
% SEM optimum
h_5 = loglog(N_f_s(I_SEM),e_SEM(I_SEM),'o','Color',[16 16 115]./256);
set(h_5,'MarkerSize',m_so);
set(h_5,'LineWidth',l_o);

%% Figure settings

set(gca,'Ticklabelinterpreter','Latex','FontSize',f_s);
xlabel('d.o.f. [-]','Interpreter','LaTex','FontSize',f_s);
ylabel('$\epsilon$ [\%]','Interpreter','LaTex','FontSize',f_s);
ylim([10^(-5) 10^2])
xlim([10^1.5 10^4.5]);
set(gca,'Xtick',[10^2 10^3 10^4]);
set(gca,'Ytick',[10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^0 10^1 10^2]);
%set(gca,'Ytick',[10^(-5) 10^(-3) 10^(-1) 10^1]);
grid on;

%set desired output size
set(gcf,'color','w')
set(fig_1, 'Units','centimeters')
height = h_f;
width = w_f;

set(fig_1, 'Position',[25 10 width height],...
       'PaperSize',[width height],...
       'PaperPositionMode','auto',...
       'InvertHardcopy', 'off',...
       'Renderer','painters'...     %recommended if there are no alphamaps
   );

% Legend entries
S = cellstr({'Pareto front FEM','FEM data point','Optimum FEM','Pareto front SEM','SEM data point','Optimum SEM'});

% Legend
[legend_h,object_h,plot_h,text_strings] = columnlegend(1,S,'Location','NorthEastOutside','Fontsize',f_s,'Interpreter','LaTex','Padding',0.5,'boxoff');

%l_1 = legend('Pareto front FEM','FEM data point','Optimum FEM','Pareto front SEM','SEM data point','Optimum SEM','Location','NorthEastOutside');
%set(l_1,'Fontsize',f_s,'Interpreter','LaTex');

%% Plot remaining data

figure(1);
% FEM data points
h_6 = loglog(N_f_F,e_FEM,'.','Color',[192 192 192]./256);
set(h_6,'MarkerSize',m_sp);
% SEM data points
h_7 = loglog(N_f_s,e_SEM,'.','Color',[105 105 105]./256);
set(h_7,'MarkerSize',m_sp);

%% Plotting and positioning

h = get(gca,'Children');
set(gca,'Children',[h(3) h(5) h(4) h(6) h(8) h(7) h(1) h(2)]);

% set(legend_h,'Position',[0.0797 0.5465 0.8979 0.45])
% set(gca,'Position',[0.1523 0.1379 0.7527 0.625])
%saveas(fig_1,'FEM_SEM_convergence','pdf')

%% Plate thickess

l_w = 2.0;
f_s = 14;
h_f = 5.5;                                                                            % Height of the figure [cm]
w_f = 8.8 + 0.5;                                                                          % Width of the figure [cm]

% Import FEM data
load('pot_core_plate.mat');

% Calculate discrepancy [%]
e_p = abs(P_e_t-P_e_t(end-1))/P_e_t(end-1)*100;

fig_3 = figure(3);
h_8 = semilogy(d_p_loop,e_p,'LineStyle','-','Marker','x','Color',[200 25 25]./256);
set(h_8,'MarkerSize',8);
hold on;
set(h_8,'LineWidth',l_w);
grid on;
set(gca,'Ticklabelinterpreter','Latex','FontSize',f_s);
xlabel('$t_{p}$ [mm]','Interpreter','LaTex','FontSize',f_s);
ylabel('$\epsilon$ [\%]','Interpreter','LaTex','FontSize',f_s);
xlim([0 0.3]);
ylim([10^(-7) 10^2]);
h_9 = semilogy(d_p_loop(3),e_p(3),'o','Color',[16 16 115]./256);
set(h_9,'MarkerSize',m_so);
set(h_9,'LineWidth',l_o);
set(gca,'Ytick',[10^(-7) 10^(-5) 10^(-3) 10^(-1) 10^1]);
% l_1 = legend('FEM data','Selected value','Location','NorthEast'); 
% set(l_1,'Fontsize',f_s,'Interpreter','LaTex');

%set desired output size
set(gcf,'color','w')
set(fig_3, 'Units','centimeters')
height = h_f;
width = w_f;

set(fig_3, 'Position',[25 10 width height],...
       'PaperSize',[width height],...
       'PaperPositionMode','auto',...
       'InvertHardcopy', 'off',...
       'Renderer','painters'...     %recommended if there are no alphamaps
   );

% Legend entries
% S = cellstr({'FEM data','Selected value'});
% 
% % Legend
% [legend_h,object_h,plot_h,text_strings] = columnlegend(1,S,'Location','NorthOutside' ,'Fontsize',f_s,'Interpreter','LaTex');
% 
% % Adjust position
% set(legend_h,'Position',[0.3003 0.79 0.5129 0.2240])
% set(gca,'Position',[0.2086 0.2313 0.6964 0.55])

%saveas(fig_3,'FEM_plate_thickness','pdf');
