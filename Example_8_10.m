clc
close all
clearvars

m=0.01; zeta=0.03; f_n=10; w_n=2*pi*f_n; w_d=sqrt(1-zeta^2)*w_n;
x_std=2;

f_s=10*f_n;
N=250;
K=N;
[Delta_f,T,Delta_t]=samplingParameters_fs_N(f_s,N);
t_row=(0:K-1)*Delta_t;
f_row=(0:N-1)*Delta_f;

h_exact_row=exp(-zeta*w_n*t_row)/m.*sin(w_d*t_row)/w_d;

RecordLengthMultiplier=800; % 40 and 800
T_long=RecordLengthMultiplier*T;
K_long=RecordLengthMultiplier*K;
rng(0);
x_long_row=x_std*randn(1,K_long);    %std(x_long_row)=x_std

y_long_row=filter(h_exact_row,1,x_long_row)*Delta_t;

% Study Linear Correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_xx_lin_row=xcorr(x_long_row,x_long_row,K)*Delta_t;    %r_xx_lin_row is 2*K+1 length. The biased and unbiased scalings is equivalent to multiplying by Delta_t/T. This is the suitable option for periodic and random signals
r_yy_lin_row=xcorr(y_long_row,y_long_row,K)*Delta_t;    %r_yy_lin_row is 2*K+1 length
r_xy_lin_row=xcorr(y_long_row,x_long_row,K)*Delta_t;    %r_xy_lin_row is 2*K+1 length
tau_sym_row=(-K:K)*Delta_t;

figure
subplot(3,1,1)
plot(tau_sym_row,r_xx_lin_row/T_long)
title(['$K=',int2str(K_long),'$'],'interpreter','latex')
legend('$r_{xx}^{\mathrm{lin}} (\tau)/T$','interpreter','latex')
set(gca,'XTickLabel',[]);

subplot(3,1,2)
r_hh_lin_row=xcorr(h_exact_row,h_exact_row,K)*Delta_t;
plot(tau_sym_row,[r_yy_lin_row/T_long;r_hh_lin_row*x_std^2*Delta_t])
legend(["$r_{yy}^{\mathrm{lin}} (\tau)/T$","$r_{hh}(\tau)$"],'interpreter','latex')
set(gca,'XTickLabel',[]);

subplot(3,1,3)
plot(tau_sym_row,r_xy_lin_row/T_long,t_row,h_exact_row*x_std^2*Delta_t);
legend(["$r_{xy}^{\mathrm{lin}} (\tau)/T$","$h^{\mathrm{exact}}(t)$"],'interpreter','latex')
xlabel('$t$ \& $\tau$ (s)','interpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_row1=(0:2*K-2)*Delta_f/2;
H_exact_row=1./(w_n^2-(2*pi*f_row1).^2+2*1i*zeta*w_n*(2*pi*f_row1))/m;

% Periodogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following is equivalent to:
% H_raw=fft(y_long_row(1:K).')./fft(x_long_row(1:K).');
R_XX_circ_Corr=cpsd(x_long_row(1:K),x_long_row(1:K),rectwin(K),0,K,'twosided')*2*pi*K;
R_XY_circ_Corr=cpsd(y_long_row(1:K),x_long_row(1:K),rectwin(K),0,K,'twosided')*2*pi*K;
H_raw=R_XY_circ_Corr./R_XX_circ_Corr;

% Welch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_XX_Welch=cpsd(x_long_row,x_long_row,hanning(2*K-1),floor((2*K-1)/2),2*K-1,f_s,'twosided').';
R_XY_Welch=cpsd(y_long_row,x_long_row,hanning(2*K-1),floor((2*K-1)/2),2*K-1,f_s,'twosided').';
H_Welch=R_XY_Welch./R_XX_Welch;

figure
[ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(f_row1,H_exact_row,false);
hold(ax_mag_h,'on');
hold(ax_phase_h,'on');
plot_FRF_mag_phase(f_row,H_raw.',false,ax_mag_h,ax_phase_h);
plot_FRF_mag_phase(f_row1,H_Welch,false,ax_mag_h,ax_phase_h);
legend(ax_mag_h,["$H(f)^{\mathrm{exact}}$","$H_{1}^{\mathrm{periodogram}}$","$H_{1}^{\mathrm{Welch}}$"],'interpreter','latex')