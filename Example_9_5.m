function Example_9_5()
clc
close all

set(groot,'DefaultAxesColorOrder',[0,0,1;0,0,0;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultLineLineWidth',1);

load beam_experiment

f_s=256;
T=4;
[D_t,K,D_f]=samplingParameters_T_fs(T,f_s);
N=K;
f_col=(0:N-1).'*D_f;
%f_col_half=(0:floor(N/2)).'*D_f;

R_XX=cpsd(x,x,hanning(K),K/2,K,f_s,'twosided');
R_YY=cpsd(y,y,hanning(K),K/2,K,f_s,'twosided');
R_XY=cpsd(y,x,hanning(K),K/2,K,f_s,'twosided');
R_YX=conj(R_XY);

figure;
subplot(8,1,1:2);
semilogy(f_col,R_XX)
xlim([0,f_s/2])
grid
set(gca,'XTickLabel',[]);
ylabel('$R_{XX}(f)$','interpreter','latex')

gamma_YX_2=abs(R_XY).^2./R_XX./R_YY;
subplot(8,1,3:4);
plot(f_col,gamma_YX_2)
axis([0,f_s/2,0,1.1])
grid
set(gca,'XTickLabel',[]);
ylabel('$\gamma_{YX}^{2}(f)$','interpreter','latex')

H_1=R_XY./R_XX; 
H_2=R_YY./R_YX;
H_T=(R_YY-R_XX + sqrt((R_XX-R_YY).^2 + 4*abs(R_XY).^2))./(2*R_YX);

ax_mag_h=subplot(8,1,5:7);
ax_phase_h=subplot(8,1,8);
plot_FRF_mag_phase(f_col,[H_1,H_2,H_T],false,ax_mag_h,ax_phase_h);
legend(ax_mag_h,["$H_{1}$","$H_{2}$","$H_{\mathrm{T}}$"],'interpreter','latex')
xlim(ax_mag_h,[0,f_s/2])
xlim(ax_phase_h,[0,f_s/2])

%Figure zoom
figure
subplot(1,2,1)
semilogy(f_col,abs([H_1,H_2,H_T]))
grid
title('zoom on 1st resonance')
xlabel('$f$ (Hz)','interpreter','latex')
ylabel('$\left|H\right|$','interpreter','latex')
legend(["$H_{1}$","$H_{2}$","$H_{\mathrm{T}}$"],'interpreter','latex')
xlim([9.5,12])

subplot(1,2,2)
semilogy(f_col,abs([H_1,H_2,H_T]))
grid
title('zoom on 2nd resonance')
xlabel('$f$ (Hz)','interpreter','latex')
ylabel('$\left|H\right|$','interpreter','latex')
legend(["$H_{1}$","$H_{2}$","$H_{\mathrm{T}}$"],'interpreter','latex')
xlim([61,64])

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove');

export_figure(1,'||',"H_estimators_experimental")
export_figure(2,'',"H_estimators_experimental_zoom")