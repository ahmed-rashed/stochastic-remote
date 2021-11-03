clc
close all
clearvars

rng(0); %Generate the same random noise each time this script is executed

A=1; f_0=1;
T_0=1/f_0;

f_s=200*f_0;

SNR=0.5;  %SNR_db=-3
f_c_by_f_Nyquist=0.5;

T_vec=[2,10,100]*T_0;
N_T=length(T_vec);

T_max=max(T_vec);
[D_t,K_max]=samplingParameters_T_fs(T_max,f_s);
t_row=(0:K_max-1)*D_t;
K_min=floor(min(T_vec)/D_t);
tau_1_row=(-(K_min-1):K_min-1)*D_t;

x_row=A*sin(2*pi*f_0*t_row);
x_hat_row=addNoise(x_row,SNR,f_c_by_f_Nyquist);

figure
ax1=subplot(N_T+1,1,1);
plot(t_row(1:K_min),[x_hat_row(1:K_min);x_row(1:K_min)])
xlabel('$t$ (sec.)','interpreter','latex');
legend({'$\hat{x}(t)=x(t)+n(t)$','$x(t)$'},'interpreter','latex')
title("$\mathrm{SNR}="+SNR+'$','interpreter','latex')

%Additional optimization of the axes for correct comparison with the correlation curves
pos=get(ax1,'Position');
pos(1)=pos(1)+pos(3)/2;
pos(3)=pos(3)/2;
set(ax1,'Position',pos,'XGrid','on')

n=1;
for T=T_vec
    K=floor(T/D_t);
    r_xx_hat_part=xcorr(x_hat_row(1:K),K_min-1,'unbiased');   %This is part of r_xx_hat. r_xx_hat extends over tau_row=(-(K-1):K-1)*D_t
    
    ax=subplot(N_T+1,1,n+1);
    plot(tau_1_row,r_xx_hat_part)
    xlabel('$\tau$ (sec.)','interpreter','latex');
    ylabel("$r_{\hat{x}\hat{x}}(\tau); \quad T="+T+'$ sec','interpreter','latex');
    set(ax,'XGrid','on')

    n=n+1;
end

export_figure(1,'||',{'s1'})
