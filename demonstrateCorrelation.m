clc
close all
clearvars

%Noise Parameters
rng(0); %Generate the same random noise each time this script is executed
SNR=0.5;
f_c_by_f_Nyq=0.8;

%Sampling parameters
f_s=1000;
K_min=500;
K_vec=[1,5,10]*K_min;
N_K=length(K_vec);

K_max=max(K_vec);
[D_f_min,T_max,D_t]=samplingParameters_fs_N(f_s,K_max);

t_tot_row=(0:K_max-1)*D_t;
tau_min_row=(-(K_min-1):K_min-1)*D_t;

%Signal Parameters
A=1;
f_0=(floor(.1*K_min)+.5)*D_f_min; %f_0 yields worst-case leakage

%Generate long signal
x_tot_row=A*square(2*pi*f_0*t_tot_row,50-D_t/10);
x_tot_hat_row=addNoise(x_tot_row,SNR,f_c_by_f_Nyq);

figure
ax1=subplot(N_K+1,1,1);
plot(t_tot_row(1:K_min),[x_tot_hat_row(1:K_min);x_tot_row(1:K_min)])
xlabel('$t$ (sec.)','interpreter','latex');
legend(["$\hat{x}(t)=x(t)+n(t)$","$x(t)$"],'interpreter','latex')
title("$T_0="+1/f_0+"$ \& SNR="+SNR,'interpreter','latex')

%Additional optimization of the axes for correct comparison with the correlation curves
pos=get(ax1,'Position');
pos(1)=pos(1)+pos(3)/2;
pos(3)=pos(3)/2;
set(ax1,'Position',pos,'XGrid','on')

n=1;
for K=K_vec
    r_xx_hat_part=xcorr(x_tot_hat_row(1:K),K_min-1)/K;  % This is part of r_xx_hat. r_xx_hat extends over tau_row=(-(K-1):K-1)*D_t
                                                        % For preiodic and random signals, correlation is scalled with K
    ax=subplot(N_K+1,1,n+1);
    plot(tau_min_row,r_xx_hat_part)
    ylabel("$T="+K*D_t+'$','interpreter','latex');
    set(ax,'XGrid','on')

    if n==1
        title("$r_{\hat{x}\hat{x}}(\tau)$",'interpreter','latex');
    end

    if n~=N_K
        set(ax,'XTickLabel',[]);
    end

    if n==N_K
        xlabel('$\tau$ (sec.)','interpreter','latex');
    end

    n=n+1;
end

% export_figure(1,'||',"s1")