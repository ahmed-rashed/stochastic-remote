clc
close all
clearvars

%DFT parameters
a=1; b=0.8; c=0.75;
% D_T_1=1.5;
% D_T_2=2;
D_T_1=1;
D_T_2=1.5;
%D_T_2=D_T_1+0.07;
SNR_y=10;


D_T=abs(D_T_2-D_T_1);
f_s=50*(1/D_T);
D_t=1/f_s;
T_start=max(D_T_1,D_T_2);
K_Tot=50*1e3;
T_Tot=samplingParameters_fs_N(f_s,K_Tot);
T_s_vec=T_Tot+T_start;
K_s_vec=T_s_vec/D_t;
K=K_Tot/10;
N=K;
[T,D_t,D_f]=samplingParameters_fs_N(f_s,K);

%s(t) parameters
rng(0);
s_vec=randn(1,K_s_vec);   %white noise
f_c=f_s/10;

if f_c<(1/(D_T_2-D_T_1)),warning('s(t) is not wide band enough for detecting D_T=D_T_1 & D_T_2'),end

[b_filt,a_filt]=butter(9,f_c/(f_s/2));s_vec=filtfilt(b_filt,a_filt,s_vec);  %Filter s(t) with (full) bandwidth approximately 20 Hz (- fc to fc).

s_vec=s_vec-mean(s_vec); s_vec=s_vec/std(s_vec); % Makes mean(s)=0 & std(s)=1;

k_start=T_start/D_t;
k_1=D_T_1/D_t;k_2=D_T_2/D_t;

x_vec =a*s_vec (k_start+1:K_s_vec);
y1_vec=b*s_vec((k_start+1:K_s_vec)-k_1); 
y2_vec=c*s_vec((k_start+1:K_s_vec)-k_2); 
y_vec=y1_vec+y2_vec;

rng(10);
y_vec=y_vec+randn(1,K_Tot)/SNR_y;

%kappa_max=round(1.5*max(k_1,k_2));
r_xy=xcorr(y_vec,x_vec, 'unbiased');
tau=(-(K_Tot-1):K_Tot-1)*D_t;

figure
plot(tau,r_xy)
set(gca,'XGrid','on')
xlabel('$\tau$ (sec.)', 'interpreter', 'latex')
ylabel('$r_{xy}(\tau)$', 'interpreter', 'latex')
xlim([0,2*max(D_T_1,D_T_2)])

R_XX=cpsd(x_vec,x_vec, hanning(K),K/2, K, f_s);  %Welch 50% overlap
R_YY=cpsd(y_vec,y_vec, hanning(K),K/2, K, f_s);
R_XY=cpsd(y_vec,x_vec, hanning(K),K/2, K, f_s);
Gamma_2_XY=abs(R_XY).^2./(R_XX.*R_YY);

f_vec=(0:N/2)*D_f;
figure
plot(f_vec,unwrap(angle(R_XY)))
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('arg\itG_x_y\rm(\itf\rm) (rad)')

figure
plot(f_vec,Gamma_2_XY)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('Coherence function')
