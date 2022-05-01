clc
clear all
close all

% Date: 1.5.2022
% Author: Radim Zedka, radim.zedka@vut.cz (Brno University of Technology)
% Description: In this script I attempt to simulate formula (6) from [1]
% and approximate it with Gamma probability-density function (PDF). 
% The script allows to set an arbitrary antenna count "M".
% All details can be found in:
% documentation/Full_rate_STLC_for_Four_Receive_Antennas___Notes.pdf



% References:

% [1] S. -chan Lim and J. Joung, “Full-Rate Space–Time Line Code for Four Receive Antennas”, IEEE wireless communications letters, pp. 1-1, 2021. 

N_sym = 1e6; % Num of samples of the random variable vector

M = 128; % Num. of TX antennas


%% 1) General variables (for all systems)

gamma_acc_vec = zeros(1, N_sym); % contains (A0^2 + A1^2  + A2^2  + A3^2 + B0^2 + B1^2  + B2^2  + B3^2)

Re_acc_vec = zeros(1, N_sym); % contains (A0*A3 + B0*B3 - A1*A2 - B1*B2)

for i = 1:M
    C_mat = sqrt(1/2)*randn(8, N_sym); % (8 x N_sym)
    
    gamma_temp_vec = sum(abs(C_mat).^2,1);
    Re_temp_vec = C_mat(1,:).*C_mat(4,:) + C_mat(5,:).*C_mat(8,:) ...
                - C_mat(2,:).*C_mat(3,:) - C_mat(6,:).*C_mat(7,:); % (1 x N_sym)
    
    gamma_acc_vec = gamma_acc_vec + gamma_temp_vec;
    Re_acc_vec   = Re_acc_vec   + Re_temp_vec;
    

    fprintf('progress: %d/%d \n',i,M)
end   

clear C_mat

%% formula (6) and (7) simulation (raw data):

% A) Original form:
Xi_A = 1/4*(gamma_acc_vec - 2*Re_acc_vec).^2./gamma_acc_vec; % (1 x N_sym) % (1 x N_sym)

% B) Approximated - after using Lemma 1:
Xi_B = 1/4*gamma_acc_vec; % (1 x N_sym)

% clear gamma_acc_vec 
% clear Re_acc_vec


%% PDF fitting parameters for Gamma distribution

switch M
    case 1
         xi_rng = [0 10];% (15) for M = 1,
         mpar_A = 1.530445; npar_A = 0.949676;
    case 2
         xi_rng = [0 14]; % (15) for M = 2,
         mpar_A = 1.574221; npar_A = 0.886598;

    case 4
         xi_rng = [0 20]; % (15) for M = 4, 
         mpar_A = 1.596783; npar_A = 0.847111;
        
    case 8
         xi_rng = [0 26]; % (15) for M = 8, 
         mpar_A = 1.600477; npar_A = 0.824861;
         
    case 16
         xi_rng = [0 40]; % (15) for M = 16, 
         mpar_A = 1.606456; npar_A = 0.815408;
         
    case 32
         xi_rng = [8 64]; % (15) for M = 32, 
         mpar_A = 1.610120; npar_A = 0.811129;
         
    case 64
         xi_rng = [30 110];  % (15)  for M = 64, 
         mpar_A = 1.608200; npar_A = 0.807196;
         
    case 128
         xi_rng = [70 190]; % (15)  for M = 128,
         mpar_A = 1.613510; npar_A = 0.808271;
         
    case 256
         xi_rng = [170 340]; % (15)  for M = 256,
         mpar_A = 1.612725; npar_A = 0.807231;
         
    case 512
        xi_rng = [400 640]; % (15)  for M = 512,
        mpar_A = 1.610609; npar_A = 0.805685;
        
    otherwise
        xi_rng = [0 4];% (15) for M = 1,
        
end

npar_B = 4;
mpar_B = 8;

%%  Histogram calculations

N_bins = 500;  % num. of histogram bins

xi_cent = linspace(xi_rng(1),xi_rng(2),N_bins); % histogram bin-center vector

d_xi = xi_cent(2) - xi_cent(1);

Xi_A_pdf = hist(Xi_A, xi_cent)/N_sym/d_xi; % (1 x N_bins)

Xi_B_pdf = hist(Xi_B, xi_cent)/N_sym/d_xi; % (1 x N_bins)

% Gamma distribution fit:
f_gamma_pdf_A = gampdf(xi_cent,npar_A*M,2/mpar_A);  % (1 x N_bins)

f_gamma_pdf_B = gampdf(xi_cent,npar_B*M,2/mpar_B);  % (1 x N_bins)

%% PDFs

idx = 10;
f1 = figure(1);

plot(xi_cent, Xi_A_pdf, 'b')
grid on
hold on
plot(xi_cent(1:8:end), f_gamma_pdf_A(1:8:end), 'b*')
plot(xi_cent, Xi_B_pdf, 'r')
plot(xi_cent(1:4:end), f_gamma_pdf_B(1:4:end), 'r*')

hold off
x1 = xlabel('$\xi$ [-]','fontsize',14);
y1 = ylabel('$f_{Xi}(\xi)$ [-]','fontsize',14);


lgd1 = legend('(6)',...
              'Gamma fit of (6)',...
              '(7)',...
              'Gamma fit of (7)',...
              'fontsize',12);

str1 = sprintf('M = %d',M);
title(str1,'fontsize',14)

set(f1, 'position',[700 200 600 500])
set(x1,'interpreter','latex')
set(y1,'interpreter','latex')
set(lgd1,'interpreter','latex')
set(lgd1,'location','best')
