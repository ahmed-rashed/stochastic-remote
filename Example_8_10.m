% Although the rectangular window is included in Table 10.2, it is rarely used since its spectral
% window side lobes cause large leakage effects.However, in Chapter 8, the rectangular window
% function is applied to estimate the spectral density functions in MATLAB Examples 8.8-8.10
% (note that we have used a very long data length T).

clc
close all
clearvars

m=0.01; zeta=0.03; f_n=10; w_n=2*pi*f_n; w_d=sqrt(1-zeta^2)*w_n;
x_std=2;

f_s=10*f_n;
N=250;
K=N;
[T,Delta_t,Delta_f]=samplingParameters_fs_N(f_s,N);

t_row=(0:K-1)*Delta_t;
f_row=(0:N-1)*Delta_f;

h_exact_vec=exp(-zeta*w_n*t_row)/m.*sin(w_d*t_row)/w_d;

RecordLengthMultiplier=800; % 40 and 800
TT=RecordLengthMultiplier*T;
KK=RecordLengthMultiplier*K;
rng(0);
x_long_row=x_std*randn(1,KK);    %std(x_long_row)=x_std

y_long_row=conv(h_exact_vec,x_long_row)*Delta_t;    %Linear convolution. This yields y_long_row of length KK+K+1
% y_long_row=filter(h_exact_vec,1,x_long_row);
y_long_row=y_long_row(1:KK);  %Use same number of elements as x_long_row

kappa_max=K;
r_xx_lin=xcorr(x_long_row,x_long_row,kappa_max,'unbiased')*Delta_t;    %r_xx_lin is 2*kappa_max+1 length
r_yy_lin=xcorr(y_long_row,y_long_row,kappa_max,'unbiased')*Delta_t;    %r_yy_lin is 2*kappa_max+1 length
r_xy_lin=xcorr(y_long_row,x_long_row,kappa_max,'unbiased')*Delta_t;    %r_xy_lin is 2*kappa_max+1 length
tau_sym_row=(-kappa_max:kappa_max)*Delta_t;

figure
subplot(3,1,1)
plot(tau_sym_row,r_xx_lin)
ylabel('$r_{xx}^{\mathrm{lin}} (\tau)$','interpreter', 'latex')
set(gca,'XTickLabel',[]);

subplot(3,1,2)
plot(tau_sym_row,r_yy_lin)
ylabel('$r_{yy}^{\mathrm{lin}} (\tau)$','interpreter', 'latex')
set(gca,'XTickLabel',[]);

subplot(3,1,3)
plot(tau_sym_row,r_xy_lin)
ylabel('$r_{xy}^{\mathrm{lin}} (\tau)$','interpreter', 'latex')
xlabel('$\tau$','interpreter', 'latex');

K_z=2*K-1;
N_z=2*N-1;
T_z=K_z/K*T;
Delta_f_z=N/N_z*Delta_f;
f_z=(0:N_z-1)*Delta_f_z;

R_XX_lin=fft(r_xx_lin(1:end-1))/(2*K-1)*2;    %R_XX_lin is 2K length
R_YY_lin=fft(r_yy_lin(1:end-1))/(2*K-1)*2;    %R_YY_lin is 2K length
R_XY_lin=fft(r_xy_lin(1:end-1))/(2*K-1)*2;    %R_XY_lin is 2K length

H1=R_XY_lin./R_XX_lin;

H_exact=1/m./(w_n^2-(2*pi*f_z).^2+2*1i*zeta*w_n*(2*pi*f_z));

figure
plot_FRF_mag_phase(f_z(1:N), [H_exact(1:N);H1(1:N)],0);
legend({'$H(f)^{\mathrm{exact}}$','$H_{1}^{\mathrm{lin}}$'},'interpreter','latex')

%%Comments 1 (use T=2000) %%
figure
plot(t_row,h_exact_vec,tau_sym_row, r_xy_lin);
legend({'$h(t)$','$r_{xy}(\tau)$'},'interpreter','latex')
xlabel('$t$ \& $\tau$ (s)','interpreter','latex');
ylabel('Amplitude')

% figure(2)
% Rhh=xcorr(h,h,N);
% plot(tau_sym_row,4*Rhh); hold on
% plot(tau_sym_row, Ryy, 'r:'); hold off
% xlabel('Lag (\it\tau_sym_row\rm)'); ylabel('Amplitude')



% Periodogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [R_XX_circ_Corr,ff]=cpsd(x_long_row,x_long_row,rectwin(KK),0,KK, f_s, 'twosided');
% [R_XY_circ_Corr,ff]=cpsd(y_long_row,x_long_row,rectwin(KK),0,KK, f_s, 'twosided');
% H1_circ_Corr=R_XY_circ_Corr./R_XX_circ_Corr;
[R_XX_circ_Corr,ff]=cpsd(x_long_row(1:K),x_long_row(1:K),rectwin(K),0,K, f_s, 'twosided');
[R_XY_circ_Corr,ff]=cpsd(y_long_row(1:K),x_long_row(1:K),rectwin(K),0,K, f_s, 'twosided');
H1_circ_Corr=R_XY_circ_Corr./R_XX_circ_Corr;

% figure
% [ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(ff(1:floor(KK/2)), H1_circ_Corr(1:floor(KK/2)),0);
% legend(ax_mag_h,{'$H_{1}^{\mathrm{circ corr}}$','$H_{1}^{\mathrm{lin}}$'},'interpreter','latex')
% grid(ax_mag_h,'off')
% grid(ax_phase_h,'off')

% Welch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_XX_Welch=cpsd(x_long_row,x_long_row, hanning(N_z),floor(N_z/2), N_z, f_s, 'twosided').';
R_XY_Welch=cpsd(y_long_row,x_long_row, hanning(N_z),floor(N_z/2), N_z, f_s, 'twosided').';
H1_Welch=R_XY_Welch./R_XX_Welch;

figure
plot_FRF_mag_phase(f_z(1:N), [H1_circ_Corr.';H1_Welch(1:N);H1(1:N)],0);
legend({'$H_{1}^{\mathrm{periodogram}}$','$H_{1}^{\mathrm{Welch}}$','$H_{1}^{\mathrm{lin}}$'},'interpreter','latex')