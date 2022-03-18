clc
clear all
close all

% Date: 18.3.2022
% Author: Radim Zedka
% Description: In this script I attempt to simulate formula (6) from [1]
% and approximate it with Chi-squared probability-density function (PDF). 
% At high TX antenna counts (M > 64) the original Chi-squared formula does
% not work in MATLAB so I switch to the Normal distribution PDF which is
% equivalent to Chi-squared with lots of freedom degrees.



% References:

% [1] S. -chan Lim and J. Joung, “Full-Rate Space–Time Line Code for Four Receive Antennas”, IEEE wireless communications letters, pp. 1-1, 2021. 

N_smp = 1e6; % Num of samples of the random variable vector

M = 32; % Num. of TX antennas


%% 1) General variables (for all systems)

gamma_mat = zeros(M, N_smp); % contains (A0^2 + A1^2  + A2^2  + A3^2 + B0^2 + B1^2  + B2^2  + B3^2)

Re_mat = zeros(M, N_smp); % contains (A0*A3 + B0*B3 - A1*A2 - B1*B2)

for i = 1:M
    C_mat = sqrt(1/2)*randn(8, N_smp); % (8 x N_sym)
    
    gamma_mat(i,:) = sum(abs(C_mat).^2,1); % (1 x N_sym)
    
    Re_mat(i,:) = C_mat(1,:).*C_mat(3,:) + C_mat(5,:).*C_mat(7,:) ...
                + C_mat(2,:).*C_mat(4,:) + C_mat(6,:).*C_mat(8,:); % (1 x N_sym)
    fprintf('progress: %d/%d \n',i,M)
end   

clear C_mat

%% formula (6) and (7) simulation (raw data):

% A) Original form:
Xi_A = 1/4*(sum(gamma_mat,1) - 2*sum(Re_mat,1)).^2./sum(gamma_mat,1); % (1 x N_sym)

% B) Approximated - after using Lemma 1:
Xi_B = 1/4*sum(gamma_mat,1); % (1 x N_sym)


clear gamma_mat 
clear Re_mat

%%  Histogram calculations

N_bins = 400;  % num. of histogram bins
% x_max = 200;   % 
x_rng = [-M M]/1 + M;
x_cent = linspace(x_rng(1),x_rng(2),N_bins); % histogram bin-center vector
% x_cent = linspace(0,x_max,N_bins);
d_x = x_cent(2) - x_cent(1);

rho = 1; % TX SNR [-]

Xi_A_pdf = hist(rho*Xi_A, x_cent)/N_smp/d_x; % (1 x N_bins)

Xi_B_pdf = hist(rho*Xi_B, x_cent)/N_smp/d_x; % (1 x N_bins)

% sum(Xi_A_pdf)*d_x
% sum(Xi_B_pdf)*d_x

% Chi-squared distribution PDF for (c = rho/2, k = 2*M):
f_Xi_A2_theor_rho_chi = 1/gamma(M)/rho^M.*x_cent.^(M-1).*exp(-x_cent/rho); 

% Chi-squared distribution PDF for (c = rho/4, k = 4*M):
f_Xi_A4_theor_rho_chi = 1/(gamma(2*M)*(rho/2)^(2*M)).*x_cent.^(2*M-1).*exp(-x_cent*2/rho); 

% Chi-squared distribution PDF for (c = rho/8, k = 8*M):
f_Xi_A8_theor_rho_chi = 1/(gamma(4*M)*(rho/4)^(4*M)).*x_cent.^(4*M-1).*exp(-x_cent*4/rho); 

%% Normal approximations of Chi-squared distribution

var_A2 = M; % Equivalent of the scaled chi-squared distribution for (c = rho/2, k = 2*M)
f_Xi_A2_theor_rho_norm = 1/sqrt(2*pi*var_A2).*exp(-(x_cent-M).^2/2/var_A2);

% Normal approximation of chi-squared  distribution (valid for high M):
var_A4 = M/2;
f_Xi_A4_theor_rho_norm = 1/sqrt(2*pi*var_A4).*exp(-(x_cent-M).^2/2/var_A4); 

% Normal approximation of chi-squared distribution:
var_A8 = M/4;
f_Xi_A8_theor_rho_norm = 1/sqrt(2*pi*var_A8).*exp(-(x_cent-M).^2/2/var_A8); 

%% PDFs

idx = 10;
f1 = figure(1);

plot(x_cent, Xi_A_pdf, 'b')
grid on
hold on
% plot(x_cent, f_Xi_A_theor_rho, 'b--')
% plot(x_cent(1:4:end), f_Xi_A2_theor_rho_norm(1:4:end), 'bx')
plot(x_cent(1:4:end), f_Xi_A2_theor_rho_chi(1:4:end), 'bx')
plot(x_cent, f_Xi_A4_theor_rho_norm, 'r')
plot(x_cent(1:4:end), Xi_B_pdf(1:4:end), 'm')
plot(x_cent, f_Xi_A8_theor_rho_norm, 'mx')
% plot(x_cent, f_Xi_A8_theor_rho_chi, 'mx')

hold off
x1 = xlabel('$\xi$ [-]','fontsize',14);
y1 = ylabel('$f_{Xi}(\xi)$ [-]','fontsize',14);

sigmm = '\sigma';
mii = '\mi';

lgd1 = legend('(6)',...
              '(9) $k=2M,  c=\rho/2$',...
              '(9) $k=4M,  c=\rho/4$',...
              '(7)',...
              '(9) $k=8M, c=\rho/8$',...
              'fontsize',12);

str1 = sprintf('M = %d',M);
title(str1,'fontsize',14)
          
xlim([-M M]/2 + M)
% ylim([0 1.2])
set(f1, 'position',[700 200 600 500])
set(x1,'interpreter','latex')
set(y1,'interpreter','latex')
set(lgd1,'interpreter','latex')
set(lgd1,'location','best')