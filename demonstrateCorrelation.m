function demonstrateCorrelation(signal_fn,SNR,f_max_by_f_Nyq,bestCaseLeakage,linearCorrelation)
rng(0); %Generate the same random noise each time this script is executed

%Signal Parameters
A=1;    f_0=1; T_0=1/f_0;

%Sampling parameters
f_s=100*f_0;  %To be able to generate both the best and worst leakage cases, f_s/f_0 must be even integer
if rem(f_s/f_0,2)~=0,error('To be able to generate both the best and worst leakage cases, f_s/f_0 must be even integer.'),end
if bestCaseLeakage
    T_min=5*T_0;  %Best case leakage
else
    T_min=5.5*T_0;  %Worst case leakage
end

[D_t,K_min,D_f_max]=samplingParameters_T_fs(T_min,f_s);

K_vec=[1,5,16]*K_min;
N_K=length(K_vec);
K_max=max(K_vec);

t_tot_row=(0:K_max-1)*D_t;

%Generate long signal
x_tot_row=A*square(2*pi*f_0*t_tot_row,50-D_t/10);
% x_tot_row=A*signal_fn(2*pi*f_0*t_tot_row);
x_tot_hat_row=addNoise(x_tot_row,SNR,f_max_by_f_Nyq);

figure
ax1=subplot(N_K+1,2,1);
plot(t_tot_row(1:K_min)/T_min,[x_tot_hat_row(1:K_min);x_tot_row(1:K_min)])
xlabel('$t/T_{\min}$','interpreter','latex');
legend(["$\hat{x}(t)=x(t)+n(t)$","$x(t)$"],'interpreter','latex')
title("SNR="+SNR,'interpreter','latex')
grid

%Additional optimization of the axes for correct comparison with the correlation curves
pos=get(ax1,'Position');
pos(1)=pos(1)+pos(3)/2;
pos(3)=pos(3)/2;
set(ax1,'Position',pos)

n=1;
for K=K_vec
    T=K*D_t;

    %Auto correlation
    r_xx_hat_biased=xcorr(x_tot_hat_row(1:K),'biased');  % For preiodic and random signals, correlation is scalled with K
    r_xx_hat_unbiased=xcorr(x_tot_hat_row(1:K),'unbiased');  % For preiodic and random signals, correlation is scalled with K
    tau_row=(-(K-1):K-1)*D_t;
    ax=subplot(N_K+1,2,2*(n+1)-1);
    plot(tau_row/T,[r_xx_hat_unbiased;r_xx_hat_biased])
    xlim([-1,1])
    grid

    % Fourier transform
    [D_f,~,~,N_no_fold]=samplingParameters_fs_N(f_s,K);
    f_row=(0:N_no_fold-1)*D_f;
    ax1=subplot(N_K+1,2,2*(n+1));
    DFT_r_xx_hat_biased=real(fft(ifftshift(r_xx_hat_biased)));
    DFT_r_xx_hat_unbiased=real(fft(ifftshift(r_xx_hat_unbiased)));
    semilogy(f_row/f_s,[DFT_r_xx_hat_biased(1:N_no_fold);DFT_r_xx_hat_unbiased(1:N_no_fold)]/K)
    set(ax1,'YAxisLocation','right');
    ylim tight
    grid

    if n==1
        ylabel(ax,"$T_{\min}="+T+'$ (sec.)','interpreter','latex');
%         ylabel(ax1,"$\Delta f=\Delta f_{\max}=1/T_{\min}$",'interpreter','latex');
        legend(ax,["unbiased","biased"])

        title(ax,"$\hat{r}_{xx}(\tau)$",'interpreter','latex');
        title(ax1,"$\mathrm{DFT}\left( \hat{r}_{xx}(\tau) \right)/N \neq \hat{R}_{XX}(f)/N$",'interpreter','latex');
    else
        ylabel(ax,"$T="+K/K_min+'T_{\min}'+'$','interpreter','latex');
%         ylabel(ax1,"$\Delta f=\Delta f_{\max}/"+D_f_max/D_f+'$','interpreter','latex');
    end

    if n~=N_K
        set([ax,ax1],'XTickLabel',[]);
    end

    if n==N_K
        xlabel(ax,'$\tau/T$','interpreter','latex');
        xlabel(ax1,'$f/f_{\mathrm{s}}$','interpreter','latex');
    end

    n=n+1;
end

% export_figure(1,'||',"s1")