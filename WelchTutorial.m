function WelchTutorial()
close all
clc

%%%
f_s=1000;
D_t=1/f_s;

N=200;
T=N*D_t;
D_f=f_s/N;
f_col=(0:N-1).'*D_f;

N_overlapRatio=0.5;
N_overlap=ceil(N*N_overlapRatio);
Q=20;
N_total=Q*(N-N_overlap)+N_overlap;

f_0=200;
t_col_total=(0:N_total-1).'*D_t;
x=sin(2*pi*f_0*t_col_total);
y=x*2;
x_noisy=x+rand(N_total,1);
y_noisy=y+rand(N_total,1);

%Welch (averaged) cpsd estimation
window=hann(N);
R_XX_noisy_our=our_cpsd(x_noisy,x_noisy,window,N_overlap);
R_YY_noisy_our=our_cpsd(y_noisy,y_noisy,window,N_overlap);
R_XY_noisy_our=our_cpsd(x_noisy,y_noisy,window,N_overlap);

R_XX_noisy_Matlab=cpsd(x_noisy,x_noisy,window,N_overlap,N,f_s,'twosided')*f_s;
R_YY_noisy_Matlab=cpsd(y_noisy,y_noisy,window,N_overlap,N,f_s,'twosided')*f_s;
R_XY_noisy_Matlab=cpsd(y_noisy,x_noisy,window,N_overlap,N,f_s,'twosided')*f_s;

Error=max(abs(R_XX_noisy_our-R_XX_noisy_Matlab))
Error=max(abs(R_YY_noisy_our-R_YY_noisy_Matlab))
Error=max(abs(R_XY_noisy_our-R_XY_noisy_Matlab))

subplot(3,1,1)
semilogy(f_col(1:floor(N/2)),R_XX_noisy_Matlab(1:floor(N/2)))
grid
title(['$R_{XX}^{\textrm{Welch}}$ for $N=',int2str(N),'$ and $Q=',int2str(Q),'$'], 'interpreter', 'latex')
subplot(3,1,2)
semilogy(f_col(1:floor(N/2)),R_YY_noisy_Matlab(1:floor(N/2)))
grid
title(['$R_{YY}^{\textrm{Welch}}$ for $N=',int2str(N),'$ and $Q=',int2str(Q),'$'], 'interpreter', 'latex')
subplot(3,1,3)
semilogy(f_col(1:floor(N/2)),abs(R_XY_noisy_Matlab(1:floor(N/2))))
grid
title(['$\left|R_{XY}^{\textrm{Welch}}\right|$ for $N=',int2str(N),'$ and $Q=',int2str(Q),'$'], 'interpreter', 'latex')
xlabel('$f$ (Hz)', 'interpreter', 'latex')

export_figure(gcf,'||',{'WelchSpectralDensity'})