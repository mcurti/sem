%% CEFC transformer current sources
clear variables;
close all;
clc;

%% Parameters

% Input parameters
V_p = 48/sqrt(2);                                                                   % Primary voltage [Vrms]
R_L = 23.5;                                                                         % Load resistance [ohm]
f = 1e6;                                                                            % Frequency [Hz]
w = 2*pi*f;                                                                         % Frequency [rad/s]

% Geometrical parameters
R_2 = 3.0e-3;
R_3 = 5.8e-3;
h_w = 2.8e-3;

% Resistances (neglected)
R_p = 1e-9;
R_s = 1e-9;

% Inductances (per turn)
L_M = 1.4374160783334447E-8;                                                       % Magnetizing inductance [H]   
L_LKP = 2.052993273544207E-8 - L_M;                                                % Primary leakage inductance [H]
L_LKS = 2.0531848656223115E-8 - L_M;                                               % Seconadry leakage inductance [H]

% Number of turns
N_p_loop = 4;%:1:10;                                                                  % Primary number of turns [-]
N_s_loop = 7;%:1:10;                                                                  % Secondary number of turns [-]

% Memory allocation
I_p_t = NaN(length(N_p_loop),length(N_s_loop));
I_s_t = NaN(length(N_p_loop),length(N_s_loop));
P_o_t = NaN(length(N_p_loop),length(N_s_loop));

% Iteration variable
i = 1;
j = 1;

while j <= length(N_s_loop);
    
    % Set number of turns
    N_p = N_p_loop(i);
    N_s = N_s_loop(j);
    a_t = N_p/N_s;                                                                      
    
    % Total inductances
    L_m = L_M*N_p^2;                                                                   
    L_lkp = L_LKP*N_p^2;                                                              
    L_lks = L_LKS*N_s^2;                                                               
    
    % Impedances (no resonance)
    Z_1 = R_p + 1i*w*L_lkp;
    Z_2 = (a_t^2)*(R_L + R_s + 1i*w*L_lks);
    Z_3 = 1i*w*L_m;
    
    % Current calculation
    i_p = ((V_p*(Z_2 + Z_3))/(Z_1*(Z_2 + Z_3) + Z_2*Z_3));
    i_s = (a_t*(V_p*-Z_3)/(Z_1*(Z_2 + Z_3) + Z_2*Z_3));
    
    % Output power
    P_o_t(i,j) = abs(i_s)^2*R_L;
    
    % Currents
    I_p_t(i,j) = abs(i_p);
    I_s_t(i,j) = abs(i_s);
    
    % Next iteration
    i = i + 1;
    if i > length(N_p_loop);
        i = 1;
        j = j + 1;
    end
end

% Primary current 
J_p = (abs(i_p)*N_p)/((R_3-R_2)*h_w)
ph_p_rad = angle(i_p)
ph_p_deg = angle(i_p)*180/pi

% Secondary current
J_s = (abs(i_s)*N_s)/((R_3-R_2)*h_w)
ph_s_rad = angle(i_s)
ph_s_deg = angle(i_s)*180/pi
    
    
    