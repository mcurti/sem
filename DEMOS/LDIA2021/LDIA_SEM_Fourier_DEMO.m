clear; clc
%% Initial data
n           = struct;
n.Q         = 100;
n.Nbiron    = 4;
n.Nd_teethy = 3;
n.Ng        = [26 41];
load('Fx_FEM.mat','Fx')
[Fx_SEM, t, dof] = LM_SEM_fourier(n);
xj = linspace(0,24e-3,240); 

%% Figures
figure(1)
clf
plot(xj, Fx,'.k','markersize',10)
hold on
plot(xj, Fx_SEM,'.r','markersize',10)
hold off
ylabel('F_x [N]')
xlabel('position, [mm]')
grid on
legend('FEM','SEM')
figure_config(1,8,6,8)

%% Functions
function [Fx, t, dof] = LM_SEM_fourier(N)
% This function takes the refinement data for SEM elements and number of
% harmonics to build the hybride model

filename = 'LinearMachine_airgap.xml';


Q = N.Q;
%%
XMLnode = xmlread(filename);

% Writting the parameters
    XMLnode.getElementsByTagName('Parameters').item(0). ...
                  setAttribute('Nd',num2str(N.Ng(1)));
    XMLnode.getElementsByTagName('Parameters').item(0). ...
                  setAttribute('Nd_teeth',num2str(N.Ng(2)));
    %----------------------------------------------------------------------
    XMLnode.getElementsByTagName('Parameters').item(0). ...
                  setAttribute('Nd_teethy',num2str(N.Nd_teethy));
    XMLnode.getElementsByTagName('Parameters').item(0). ...
                  setAttribute('Nbiron',num2str(N.Nbiron));
xmlwrite(filename,XMLnode);
%-------------------------------------------------
% start the time
tic
%-------------------------------------------------
G  = Geometry(filename);
%% Changing the phase of the current density
param = linspace(0,24,240);
Fx    = zeros(1,length(param));
dof = zeros(1,length(param));
x_off = 0;
for ii = 1:length(param)
Ph = Physics(G);
% 
Ph = Ph.load_magnetic_materials;
Ph = Ph.load_current_density;
%% Fourier implementation

X1 = 0 + x_off; X2 = 48 + x_off;
h1 = 5; h2 = 9; h3 = 9.9; % sizes of Fourier domain
taup = 12;
alpham = 2/3;


Elements           = struct;
Elements.x_start   = X1;
Elements.tau       = X2-X1;
Elements.Harmonics = Q;
Elements.heights   = [h3, h2
                      h2, h1];
Nfel               = 2;
F    = fourierElements(Elements);

% Initiating Fourier Data
FourierData = struct; Gf = struct;
FourierData.connectedLines    = 58:66;
FourierData.connectedElements = 19:27;
Gf.lines             = G.lines;
Gf.xi                = G.xi;
Gf.metrics           = G.metrics;
Gf.Permeability      = Ph.Materials.Permeability;
% Magnetic sources in the fourier region
tauk = @(k) .5*taup*(k - (1 + (-1)^k)/2 + alpham*(-1)^k);
w_i  = F.w_n;
mu0  = pi*4e-7;
M    = 1.3e-3/mu0;

a_M  = zeros(Q,1); b_M  = zeros(Q,1);

for k = 2:2:8
    a_M = a_M + 2*M/(X2-X1)*(-1)^(k/2)./w_i.*(cos(w_i*tauk(k)) - cos(w_i*tauk(k-1)));
    b_M = b_M - 2*M/(X2-X1)*(-1)^(k/2)./w_i.*(sin(w_i*tauk(k)) - sin(w_i*tauk(k-1)));
end

abs_M = a_M;


%% Building the problem
Pr = Problem(G, Ph);
Pr = Pr.global_matrix;
Pr = Pr.load_Y_sources;
Pr = Pr.building_Y_vector;
% Adding Fourier terms in the matrix
E = Pr.Global_Matrix;
    

[Espace, ~] = F.fourier_space_matrix(Pr.ProblemData,Gf,...
                 FourierData);
Efrequency = F.fourier_frequency_matrix(Pr.ProblemData,Gf,...
                 FourierData);
%     
% Eglobal = [E Espace; Efrequency];

    x_off = param(ii);
    F = F.update_x_start(x_off);
    
    
    a_M = abs_M.*cos(-2*pi*(1:Q)'/Elements.tau*x_off);
    b_M = abs_M.*sin(-2*pi*(1:Q)'/Elements.tau*x_off);

    new_phi = sprintf('%.4f',-pi/(taup)*x_off - 50*pi/180);
    Ph = Ph.setParameter('phi',new_phi);
    Ph = Ph.load_current_density;

    Pr = Pr.updateSources(Ph.Sources);
    Pr = Pr.load_Y_sources;
    Pr = Pr.building_Y_vector;
    Y = [Pr.Y.vector zeros(1,Q) -b_M'./w_i' zeros(1,Q) a_M'./w_i' zeros(1,4*Q) 0];
    
    % Solving the linear system
    Eglobal = [E Espace; Efrequency];

    S = sparse(Eglobal);

    PHI = S\Y'*mu0;

    
    X = PHI(length(Pr.Y.vector)+1:end);
    
    c1 = zeros(1,Q*Nfel); c2 = zeros(1,Q*Nfel); c3 = zeros(1,Q*Nfel);
    c4 = zeros(1,Q*Nfel); %Bx0 = zeros(1,1);
    for k = 1:Nfel
        c1((1:Q) + (k-1)*Q) = X((1:Q) + 4*Q*(k-1));  
        c2((1:Q) + (k-1)*Q) = X((1:Q) + Q*(1 + 4*(k-1)));
        c3((1:Q) + (k-1)*Q) = X((1:Q) + Q*(2 + 4*(k-1)));
        c4((1:Q) + (k-1)*Q) = X((1:Q) + Q*(3 + 4*(k-1)));

    end
    
    Bx0 = X(end);
    
    Az0 = X(end-1)*0;
    
    s  = [zeros(1,2*Q) a_M'*mu0 b_M'*mu0];
    F = F.update_coefficients(s, c1, c2, c3, c4, Az0, Bx0);
    
    % Computing forces
    Fx(ii) = F.MaxwellForce(1);
    dof(ii) = length(Y);
    %% Plot the solution
    Pp = PostProcessing(G, Pr, PHI);
    
    % solution for Fourier regions
    [dx1, dy1, F1] = F.fourier2space(100,10,1);
    [dx2, dy2, F2] = F.fourier2space(100,50,2);
    figure(11)
    clf
    hold on
    G.plot_geometry('markersize',12,'color','k','points')
    G.plot_geometry('color','k','lines')
    Pp.plot_contour_var(Pp.Potential)
    contour(dx1,dy1,F1,linspace(min(PHI(1:length(Pr.Y.vector))), max(PHI(1:length(Pr.Y.vector))), 20))
    contour(dx2,dy2,F2,linspace(min(PHI(1:length(Pr.Y.vector))), max(PHI(1:length(Pr.Y.vector))), 20))
    c = colorbar;
    hold off
    title(c,'Az, Wb/mm')
    xlabel('x, mm')
    ylabel('y, mm')
    axis equal
    pause(1)
end
t = toc; dof = dof(end);
end

