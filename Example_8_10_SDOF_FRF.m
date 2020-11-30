clearvars
close all
clc


fs=100;
N=250;
K=N;
[T,Delta_t,Delta_f]=samplingParameters_fs_N(fs,N);
t_row=(0:N-1)*Delta_t;
f_row=(0:N-1)*Delta_f;

m=.01; zeta=0.03; f_n=10; w_n=2*pi*f_n; w_d=sqrt(1-zeta^2)*w_n;
h_exact_vec=exp(-zeta*w_n*t_row)/m.*sin(w_d*t_row)/w_d;

ScaleFactor=5; % 40 and 800
TT=ScaleFactor*T;
NN=ScaleFactor*N;
rng(0);
x_long_row=2*randn(1,NN); 

y_long_row=conv(h_exact_vec,x_long_row)*Delta_t;  % y_long_row=filter(h_exact_vec,1,x_long_row);
y_long_row=y_long_row(1:NN);

kappa_max=K;
r_xx_lin=xcorr(x_long_row,x_long_row,kappa_max,'unbiased');
r_yy_lin=xcorr(y_long_row,y_long_row,kappa_max,'unbiased');
r_xy_lin=xcorr(y_long_row,x_long_row,kappa_max,'unbiased');
tau_vec=(-kappa_max:kappa_max)*Delta_t;

figure
subplot(3,1,1)
plot(tau_vec,r_xx_lin)
xlabel('$\tau$','interpreter','latex'); ylabel('$r_{xx}(\tau)$','interpreter','latex')

subplot(3,1,2)
plot(tau_vec,r_yy_lin)
xlabel('$\tau$','interpreter','latex'); ylabel('$r_{yy}(\tau)$','interpreter','latex')

subplot(3,1,3)
plot(tau_vec,r_xy_lin)
xlabel('$\tau$','interpreter','latex'); ylabel('$r_{xy}(\tau)$','interpreter','latex')

N_z=2*kappa_max;
T_z=2*T;
D_f_z=Delta_f/2;
f_pad=(0:N_z-1)*D_f_z;

Sxx=fft(r_xx_lin(1:N_z));
Syy=fft(r_yy_lin(1:N_z));
Sxy=fft(r_xy_lin(1:N_z));

H1=Sxy./Sxx;
H_exact=fft(h_exact_vec,N_z)/fs;

figure
subplot(3,1,1:2)
plot(f_pad(1:N),20*log10(abs(H_exact(1:N))),f_pad(1:N),20*log10(abs(H1(1:N))),'r');
xlabel('$f$ (Hz)','interpreter','latex');
ylabel('$|H(f)|$ (dB)','interpreter','latex')
legend({'$H(f)$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

subplot(3,1,3)
plot(f_pad(1:N),unwrap(angle(H_exact(1:N))),f_pad(1:N),unwrap(angle(H1(1:N))),'r');
xlabel('$f$ (Hz)','interpreter','latex');
ylabel('$\angle H(f)$ (rad)','interpreter','latex')
legend({'$H(f)$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

%% Comments 1 (use T=2000) %%

% figure
% plot(t_row,4*h_exact_vec); hold on
% plot(tau_vec(kappa_max+1:end),r_xy_lin(kappa_max+1:end),'r:'); hold off
% xlabel('Time (s) and lag (\it\tau\rm)'); ylabel('Amplitude')
% 
% figure
% Rhh=xcorr(h_exact_vec,h_exact_vec,kappa_max);
% plot(tau_vec,4*Rhh); hold on
% plot(tau_vec,r_yy_lin,'r:'); hold off
% xlabel('Lag (\it\tau\rm)'); ylabel('Amplitude')

%% Comments 2 (use T=100) %%
Sxx_circ_Corr=cpsd(x_long_row,x_long_row,ones([1,NN]),0,N_z,fs,'twosided');
Sxy_circ_Corr=cpsd(y_long_row,x_long_row,ones([1,NN]),0,N_z,fs,'twosided');
H1_circ_Corr=Sxy_circ_Corr./Sxx_circ_Corr;

figure
subplot(3,1,1:2)
plot(f_pad(1:N),20*log10(abs(H1_circ_Corr(1:N))),f_pad(1:N),20*log10(abs(H1(1:N))),'r');
xlabel('$f$ (Hz)','interpreter','latex');
ylabel('$|H(f)|$ (dB)','interpreter','latex')
legend({'$H_{1}^{\mathrm{circ corr}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

subplot(3,1,3)
plot(f_pad(1:N),unwrap(angle(H1_circ_Corr(1:N))),f_pad(1:N),unwrap(angle(H1(1:N))),'r'); hold on;set(gca,'ColorOrderIndex',1);
xlabel('$f$ (Hz)','interpreter','latex');
ylabel('$\angle H(f)$ (rad)','interpreter','latex')
legend({'$H_{1}^{\mathrm{circ corr}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

Sxx_w=cpsd(x_long_row,x_long_row,hanning(N_z),N_z/2,N_z,fs,'twosided');
Sxy_w=cpsd(y_long_row,x_long_row,hanning(N_z),N_z/2,N_z,fs,'twosided');
H1_w=Sxy_w./Sxx_w;

figure
subplot(3,1,1:2)
plot(f_pad(1:N),20*log10(abs(H1_w(1:N))),f_pad(1:N),20*log10(abs(H1(1:N))),'r');
xlabel('$f$ (Hz)','interpreter','latex');
ylabel('$|H(f)|$ (dB)','interpreter','latex')
legend({'$H_{1}^{\mathrm{Welch}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

subplot(3,1,3)
plot(f_pad(1:N),unwrap(angle(H1_w(1:N))),f_pad(1:N),unwrap(angle(H1(1:N))),'r'); hold on;set(gca,'ColorOrderIndex',1);
xlabel('$f$ (Hz)','interpreter','latex');
ylabel('$\angle H(f)$ (rad)','interpreter','latex')
legend({'$H_{1}^{\mathrm{Welch}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')
