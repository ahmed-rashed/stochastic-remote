clc
close all
clearvars

%DFT parameters
a=1; b=0.8; c=0.75;
D_T_1=1;
D_T_2=D_T_1+0.5;
%D_T_2=D_T_1+0.07;
SNR_y=10;

D_T=abs(D_T_2-D_T_1);
T_start=max(D_T_1,D_T_2);
f_s=50*(1/D_T);
K_Tot=50e3;
[T_Tot,D_t]=samplingParameters_fs_N(f_s,K_Tot);
T_s=T_Tot+T_start;
K_s=T_s/D_t;
K=K_Tot/10;
N=K;
[T,D_t,D_f]=samplingParameters_fs_N(f_s,K);

%s(t) parameters
rng(0);
s_vec=randn(1,K_s);   %white noise
f_c=f_s/10;

if f_c<(1/(D_T_2-D_T_1)),warning('s(t) is not wide band enough for detecting D_T=D_T_1 & D_T_2'),end

[b_filt,a_filt]=butter(9,f_c/(f_s/2));  %designs a 9th-order low-pass digital Butterworth filter (IIR),where b is a vector containing coefficients of a moving average part and a is a vector containing coefficients of an auto-regressive part of the transfer function (see Equation (6.12) of Shin's book).
s_vec=filtfilt(b_filt,a_filt,s_vec);  %Filter s(t) with (full) bandwidth approximately 20 Hz (- fc to fc).

s_vec=s_vec-mean(s_vec);    % Makes mean(s)=0
s_vec=s_vec/std(s_vec);     % Makes std(s)=1

k_start=T_start/D_t;
k_1=D_T_1/D_t;k_2=D_T_2/D_t;

x_vec =a*s_vec (k_start+1:K_s);
y_vec=b*s_vec((k_start+1:K_s)-k_1)+c*s_vec((k_start+1:K_s)-k_2);

rng(10);
y_vec=y_vec+randn(1,K_Tot)/SNR_y;

%kappa_max=round(1.5*max(k_1,k_2));
r_xy=xcorr(y_vec,x_vec,'unbiased');
tau=(-(K_Tot-1):K_Tot-1)*D_t;

plot(tau,r_xy)
set(gca,'XGrid','on')
xlabel('$\tau$ (sec.)','interpreter','latex')
ylabel('$r_{xy}(\tau)$','interpreter','latex')
xlim([0,2*max(D_T_1,D_T_2)])

%Welch with hanning window and 50% overlap
win_col=window(@hann,K,'periodic');
R_XX=cpsd(x_vec,x_vec,win_col,K/2,K,f_s);
R_YY=cpsd(y_vec,y_vec,win_col,K/2,K,f_s);
R_XY=cpsd(y_vec,x_vec,win_col,K/2,K,f_s);
Gamma_2_XY=abs(R_XY).^2./(R_XX.*R_YY);

f_vec=(0:N/2)*D_f;
figure
plot(f_vec,unwrap(angle(R_XY)))
xlabel('$f$ (Hz)','interpreter','latex')
ylabel('arg\itG_x_y\rm(\itf\rm) (rad)')

figure
plot(f_vec,Gamma_2_XY)
xlabel('$f$ (Hz)','interpreter','latex')
ylabel('Coherence function')
