function Example_10_1()
clc
close all

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultLineLineWidth',1);

A_r=[-0.31837,-0.15916];
f_r=[5,15];
zeta_r=[0.02,0.01];

f_s=100/15*max(f_r);
T=10;   %T and T_r should be merged in one variable called T (T_r is the notation of shin, and T is my notation)
T_r=4;
%T_r=20;

[D_t,K,D_f]=samplingParameters_T_fs(T,f_s);
N=K;
t_col=(0:K-1).'*D_t;
f_col=(0:N-1).'*D_f;

K_Tot=200*K;     %T_Tot=200*T;
%K_Tot=1e3*K;    %T_Tot=1e3*T;

w_r=2*pi*f_r;
w_d_r=sqrt(1-zeta_r.^2).*w_r;
h_exact=(-2*A_r*(exp(-(zeta_r.*w_r).'*t_col.').*sin(w_d_r.'*t_col.'))).';
%H_exact=(-2*(w_d_r.*A_r)*(1./(w_r.'.^2*ones(1,N)-ones(2,1)*(2*pi*f_col.').^2+2*1i*(zeta_r.*w_r).'*2*pi*f_col.'))).';
H_DFT=fft(h_exact)/f_s;

rng(0); 
x_vec_long=randn(1,K_Tot);
y_vec_long=filter(h_exact,1,x_vec_long)/f_s; %scaled appropriately. 

%%Comments on the segment averaging method: (use T_Tot=2000 and T_r=20)
% x=[x 2*x 3*x 4*x 5*x]; x=x-mean(x); x=x/std(x);

N_r=T_r/D_t;
R_XX=cpsd(x_vec_long,x_vec_long, hanning(N_r),N_r/2, N_r, f_s, 'twosided')*f_s; %Welch with 50% overlap
R_YY=cpsd(y_vec_long,y_vec_long, hanning(N_r),N_r/2, N_r, f_s, 'twosided')*f_s; %Welch with 50% overlap
[R_XY,f]=cpsd(y_vec_long,x_vec_long, hanning(N_r),N_r/2, N_r, f_s, 'twosided'); %Welch with 50% overlap
R_XY=R_XY*f_s;

figure
semilogy(f,R_XX, f, ones(size(f)))
xlabel('$f$ (Hz)','interpreter','latex')
legend({'$R_{\hat{X}\hat{X}}(f)$','$R_{XX}(f)$'},'interpreter','latex')
xlim([0,f_s/2])

figure
semilogy(f,R_YY, f_col,abs(H_DFT.^2))
xlabel('$f$ (Hz)','interpreter','latex')
legend({'$R_{\hat{Y}\hat{Y}}(f)$','$R_{YY}(f)$'},'interpreter','latex')
xlim([0,f_s/2])

figure
[ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(f,R_XY,false);
hold(ax_mag_h, 'all');hold(ax_phase_h, 'all')
plot_FRF_mag_phase(f_col,H_DFT,false,ax_mag_h,ax_phase_h,'','R_{XY}(f)');
legend(ax_mag_h,{'$R_{\hat{X}\hat{Y}}(f)$','$R_{XY}(f)$'},'interpreter','latex')
xlim(ax_mag_h,[0,f_s/2]);xlim(ax_phase_h,[0,f_s/2])

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove');