function Example_10_3()
clc
close all

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultLineLineWidth',1);

delta=1;
f_s=20;
T=1.1;  %T=1.1, 2, 5, 50;
[D_t,K,D_f]=samplingParameters_T_fs(T,f_s);
Nt=1000*K;

rng(0);  
x_vec=randn(1,Nt+delta/D_t); 
y_vec=x_vec(1:length(x_vec)-delta/D_t);
x_vec=x_vec(delta/D_t+1:end);

[R_XY, f]=cpsd(y_vec,x_vec, rectwin(K),0, 1000, f_s, 'twosided');
R_XY=R_XY*f_s;

figure
[ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(f,[R_XY,ones(size(R_XY))]);
%plot(f,abs(R_XY), f, ones(size(f)))
%xlabel('$f$ (Hz)','interpreter','latex')
legend(ax_mag_h,{'$R_{\hat{X}\hat{Y}}(f)$','$R_{XY}(f)$'},'interpreter','latex')
%ylabel('Estimate of |\itS_x_y\rm(\itf\rm)| (linear scale)')
ylim(ax_mag_h,[0 1.1])
xlim(ax_mag_h,[0 f_s/2])
xlim(ax_phase_h,[0 f_s/2])

figure
plot(f,unwrap(angle(R_XY)), [0 f_s/2], [0 -2*pi*10*delta])
xlabel('Frequency (Hz)')
ylabel('Estimate of arg\itS_x_y\rm(\itf\rm) (rad)')
axis([0 f_s/2 -65 0])

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove');