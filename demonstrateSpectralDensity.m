function demonstrateSpectralDensity(signal_fn,SNR,i_win,linearCorrelation,f_max_by_f_Nyq,bestCaseLeakage,unbiasedCorr)
rng(0); %Generate the same random noise each time this script is executed

%Signal Parameters
A=1;    f_0=1; T_0=1/f_0;

%Sampling parameters
f_s=100*f_0;  %To be able to generate both the best and worst leakage cases, f_s/f_0 must be even integer
if rem(f_s/f_0,2)~=0,error('To be able to generate both the best and worst leakage cases, f_s/f_0 must be even integer.'),end
if bestCaseLeakage
    T=5*T_0;  %Best case leakage
else
    T=5.5*T_0;  %Worst case leakage
end

[D_t,K]=samplingParameters_T_fs(T,f_s);

%Welch parapeters
P_min=1;   %No of Welch averages
P_vec=[1,10,1000]*P_min;
N_P=length(P_vec);
P_max=max(P_vec);
K_o=0;  %Overlap points
K_tot=P_max*(K-K_o)+K_o;
t_tot_row=(0:K_tot-1)*D_t;

% Windows
fn_temp=@(w)(w(1:end-1));
window_fn_cvec={@(K) (window(@rectwin,K)),@(K) (window(@hann,K,'periodic')),@(K) (window(@hamming,K,'periodic')),@(K) (fn_temp(window(@kaiser,K+1,9.44683)))};
windowName_vec=["rectangular","Hann","Hamming","Kaiser-Bessel"];

%Generate long signal
x_tot_row=A*square(2*pi*f_0*t_tot_row,50-D_t/10);
% x_tot_row=A*signal_fn(2*pi*f_0*t_tot_row);
x_tot_hat_row=addNoise(x_tot_row,SNR,f_max_by_f_Nyq);

%Plot one record of the signal
win_col=window_fn_cvec{i_win}(K);
figure;
ax=subplot(N_P+1,2,1);
P_plot=1;
ind_vec=(P_plot-1)*K+(1:K);
plot(t_tot_row(ind_vec)/T,[x_tot_hat_row(ind_vec);x_tot_row(ind_vec)].*win_col.')
xlabel('$t/T$','interpreter','latex');
legend(["$\hat{x}(t)w(t);\;\hat{x}(t)=x(t)+n(t)$","$x(t)w(t)$"],'interpreter','latex')
title("Record "+P_plot+" of $x(t)$ with "+windowName_vec(i_win)+" window and SNR="+SNR,'interpreter','latex');
grid

%Additional optimization of the axes for correct comparison with the correlation curves
pos=get(ax,'Position');
pos(1)=pos(1)+pos(3)/2;
pos(3)=pos(3)/2;
set(ax,'Position',pos)

if linearCorrelation    %Linear correlation
    K_z=2*K-1; %Perform linear (instead of circular) CSD by zero padding
    tau_row=(-(K-1):K-1)*D_t;
else                    %Circular correlation
    K_z=K;
    tau_row=(-ceil((K-1)/2):floor((K-1)/2))*D_t;
end
N_z=K_z;
f_s_z=f_s;
[D_f_z,~,~,N_z_no_fold]=samplingParameters_fs_N(f_s_z,N_z);
f_z_row=(0:N_z_no_fold-1)*D_f_z;

n=1;
for P=P_vec
    K_tot_i=P*(K-K_o)+K_o;
    R_XX_Welch=cpsd(x_tot_hat_row(1:K_tot_i),x_tot_hat_row(1:K_tot_i),win_col,K_o,K_z,'twosided')*2*pi;   % For preiodic and random signals, correlation is scalled with K

    r_xx_Welch=fftshift(ifft(R_XX_Welch));
    ax=subplot(N_P+1,2,2*(n+1)-1);
    plot(tau_row/T,r_xx_Welch)
%     xlim([-1,1]*0.5)
    grid

    ax1=subplot(N_P+1,2,2*(n+1));
    semilogy(f_z_row/f_s,R_XX_Welch(1:N_z_no_fold))
    set(ax1,'YAxisLocation','right');
    ylabel("$P="+P+'$','interpreter','latex');
    ylim tight
    grid
%     set(gca,'PlotBoxAspectRatio',[2,1,1])
%     legend_str=["$\hat{R}_{XX}^{\mathrm{raw}}$","$\hat{R}_{\bar{X}\bar{X}}^{\mathrm{raw}}$ using Hann window","$\hat{R}_{XX}^{\mathrm{Welch}}$ using Hann window \& $P="+P_min+'$'];
%     legend(legend_str,'interpreter','latex','Location','southeast')

    if n==1
        title(ax,"$\mathrm{IDFT}\left( \hat{R}_{XX}^{\mathrm{Welch}} \right) \neq \hat{r}_{xx}(\tau)$",'interpreter','latex');
        title(ax1,"$\hat{R}_{XX}^{\mathrm{Welch}}$ using $K_{\mathrm{o}}="+K_o+'$','interpreter','latex');
    end

    if n~=N_P
        set([ax,ax1],'XTickLabel',[]);
    end

    if n==N_P
        xlabel(ax,'$\tau/T$','interpreter','latex');
        xlabel(ax1,'$f/f_{\mathrm{s}}$','interpreter','latex');
    end

    n=n+1;
end


% %Periodogram Estimator
% R_XX_periodogram=periodogram(x_tot_hat_row(1:K),[],N,f_s,'twosided')*f_s;
% 
% %Modified Periodogram Estimator
% %win_col=ones(K,1);
% R_XX_m_periodogram=periodogram(x_tot_hat_row(1:K),win_col,N,f_s,'twosided')*f_s;
% 
% %Welch (averaged) cpsd estimation
% R_XX_Welch=cpsd(x_tot_hat_row,x_tot_hat_row,win_col,K_o,N,f_s,'twosided')*f_s;
% 
% figure
% subplot(3,1,1)
% semilogy(f_row/f_s,[periodogram(x_tot_row(1:K),[],N,f_s,'twosided')*f_s,periodogram(x_tot_row(1:K),win_col,N,f_s,'twosided')*f_s,cpsd(x_tot_row,x_tot_row,win_col,K_o,N,f_s,'twosided')*f_s])
% xlim([0,.5])
% 
% subplot(3,1,2)
% semilogy(f_row/f_s,[periodogram(noisee(1:K),[],N,f_s,'twosided')*f_s,periodogram(noisee(1:K),win_col,N,f_s,'twosided')*f_s,cpsd(noisee,noisee,win_col,K_o,N,f_s,'twosided')*f_s])
% xlim([0,.5])
% 
% subplot(3,1,3)
% semilogy(f_row/f_s,[R_XX_periodogram,R_XX_m_periodogram,R_XX_Welch])
% xlim([0,.5])
% grid on
% % set(gca,'PlotBoxAspectRatio',[2,1,1])
% xlabel('$f/f_{\mathrm{s}}$','interpreter','latex');
% legend_str=["$\hat{R}_{XX}^{\mathrm{raw}}$","$\hat{R}_{\bar{X}\bar{X}}^{\mathrm{raw}}$ using Hann window","$\hat{R}_{XX}^{\mathrm{Welch}}$ using Hann window \& $P="+P_min+'$'];
% legend(legend_str,'interpreter','latex','Location','southeast')


% export_figure(gcf,'',"NonParametricSpectralEstimators");