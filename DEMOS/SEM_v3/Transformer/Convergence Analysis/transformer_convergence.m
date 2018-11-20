%% SEM Convergence analysis
clear variables;
close all;
clc;

%% Parameters

% Degree
N_x = 9;%1:1:10;
N_y = 8;%1:1:30;
N_px = 8;%1:1:30;
N_y1 = 8;%1:1:10;

% Geometry
D = 10e-3;                                                                          % Depth of the domain [m]
sigma_r = [6,12,18,24,30,36,42];                                                    % Elements having non-zero conductivity
N_el = 42;                                                                          % Number of elements [-]
f = 1e6;                                                                            % Frequency [Hz]
sigma = 8.41e6;                                                                     % Plate conductivity [S/m]
mu_r = 100;                                                                         % Relative permeability plate 
mu_0 = 4e-7*pi;                                                                     % Permeability of free space
delta_s = sqrt(1/(pi*f*mu_0*mu_r*sigma));                                           % Skin depth [m]

% Memory allocation
P_e_t = NaN(length(N_x),length(N_y),length(N_px),length(N_y1));
N_f_t = NaN(length(N_x),length(N_y),length(N_px),length(N_y1));
t_c_t = NaN(length(N_x),length(N_y),length(N_px),length(N_y1));
J_e = cell(1,N_el);
P_e = zeros(1,N_el);

% Total number of points and estimated time
x_tot = length(N_x)*length(N_y)*length(N_px)*length(N_y1);
t_est = x_tot*0.6586/3600/24

% Iteration variable
i = 1;
j = 1;
k = 1;
l = 1;
x = 1;

% Create save file
%save('2D_SEM_convergence_V5.mat','P_e_t','N_f_t','t_c_t','N_x','N_y','N_px','N_y1','i','j','k','l','x');

% Degree sweep 
while l <= length(N_y1);
    
    % Import xml file
    A = xmlread('Transformer_N.xml');
    
    % Set degree
    A.getElementsByTagName('Parameters').item(0).setAttribute('N_x',num2str(N_x(i)));
    A.getElementsByTagName('Parameters').item(0).setAttribute('N_y',num2str(N_y(j)));
    A.getElementsByTagName('Parameters').item(0).setAttribute('N_px',num2str(N_px(k)));
    A.getElementsByTagName('Parameters').item(0).setAttribute('N_y1',num2str(N_y1(l)));
    %A.getElementsByTagName('Parameters').item(0).setAttribute('N_y1',num2str(N_y(j)));
    
    % Write to xml file
    xmlwrite('Transformer_N.xml',A);
    
    % Start timer
    tic;
    
    % Load data and prepare for solving
    S = SEM('Transformer_N.xml');
    
    % Solve problem
    S = S.solve;
    
    % Calculate eddy current losses
    for m = 1:N_el
        if any(m == sigma_r)
            J_e{m} = real(S.PProcessing.Potential{m}*1i*2*pi*f*sigma);
            P_e(m) = sum(sum((1/sigma)*abs(S.PProcessing.Potential{m}*1i*2*pi*f*sigma./sqrt(2)).^2.*S.Geometry.metrics.W{m}.*S.Geometry.metrics.J{m}));
        else
            J_e{m} = 0;
        end
    end
    
    % End timer
    t = toc;
    
    % Store data
    P_e_t(i,j,k,l) = sum(P_e)*D;
    N_f_t(i,j,k,l) = length(S.Problem.Y.vector);
    t_c_t(i,j,k,l) = t;
    
    % Next iteration
    i = i + 1;
    if i > length(N_x);
        i = 1;
        j = j + 1;
        if j > length(N_y);
            j = 1;
            k = k + 1;
            if k > length(N_px);
                k = 1;
                l = l + 1;
            end
        end
    end
    clc;
    
    % Append to save file
    %save('2D_SEM_convergence_V5.mat','P_e_t','N_f_t','t_c_t','i','j','k','l','x','-append');
    
    % Print progress
    x_pro = [x x_tot]
    x = x + 1;
end

%% Save data
%save('2D_SEM_convergence_V5_final.mat','P_e_t','N_f_t','t_c_t','N_x','N_y','N_px','N_y1');
    
    
    
    