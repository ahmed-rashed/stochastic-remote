clc
close all
clearvars

set(groot,'DefaultAxesColorOrder',[1,0,0;0,0,1;0,0.5,0;1,0,1;0,0,0])
set(groot,'DefaultLineLineWidth',1);

%DFT parameters
T=4;
K=2^10; %K=1024
[Delta_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_vec=(0:K-1)*Delta_t;
tau_vec=(-(K-1):K-1)*Delta_t;

%s(t) parameters
T_burst=T/8;
f_0=(f_s/2)/10;

%x(t) and y(t) parameters
a=1;b=.8; c=.75;
T_1=T/4;
T_2=3/4*T;

s_fn=@(t) (sin(2*pi*f_0*t));

x_vec=a*signalPulse(t_vec,T_burst,s_fn);
y_vec=b*signalPulse(t_vec-T_1,T_burst,s_fn)+c*signalPulse(t_vec-T_2,T_burst,s_fn);

SNR_vec=[1/2,inf];
N_SNR=length(SNR_vec);
legend_str=cell(N_SNR,1);

x_ax=subplot(5,1,1);hold on
set(x_ax,'XGrid','on','YLimSpec','tight','XTickLabel',[])
ylabel(x_ax,'$\hat{x}(t)$','interpreter','latex')
%Additional optimization of the axes for correct comparison with the correlation curves
pos=get(x_ax,'Position');
pos(1)=pos(1)+pos(3)/2;
pos(3)=pos(3)/2;
set(x_ax,'Position',pos)

y_ax=subplot(5,1,2);hold on
set(y_ax,'XGrid','on','YLimSpec','tight','XTickLabel',[])
ylabel(y_ax,'$\hat{y}(t)$','interpreter','latex')
xlabel(y_ax,'$t$ (s)','interpreter','latex');
%Additional optimization of the axes for correct comparison with the correlation curves
pos=get(y_ax,'Position');
pos(1)=pos(1)+pos(3)/2;
pos(3)=pos(3)/2;
set(y_ax,'Position',pos)

r_xx_hat_ax=subplot(5,1,3);hold on;set(r_xx_hat_ax,'XGrid','on')
set(r_xx_hat_ax,'XTickLabel',[],'YLimSpec','tight');
ylabel(r_xx_hat_ax,'$r_{\hat{x}\hat{x}}(\tau)$','interpreter','latex')

r_yy_hat_ax=subplot(5,1,4);hold on;set(r_yy_hat_ax,'XGrid','on')
set(r_yy_hat_ax,'YLimSpec','tight','XTickLabel',[])
ylabel(r_yy_hat_ax,'$r_{\hat{y}\hat{y}}(\tau)$','interpreter','latex')

r_xy_hat_ax=subplot(5,1,5);hold on;set(r_xy_hat_ax,'XGrid','on')
set(r_xy_hat_ax,'YLimSpec','tight')
ylabel(r_xy_hat_ax,'$r_{\hat{x}\hat{y}}(\tau)$','interpreter','latex')
xlabel(r_xy_hat_ax,'$\tau$ (s)','interpreter','latex')

rng(0);
for ii=1:N_SNR
    if SNR_vec(ii)==inf
        legend_str{ii}='$\mathrm{SNR}=\infty$';
        line_width=get(groot,'DefaultLineLineWidth')*1.5;
    else
        legend_str{ii}=['$\mathrm{SNR}=',num2str(SNR_vec(ii)),'$'];
        line_width=get(groot,'DefaultLineLineWidth');
    end

    x_hat_vec=addNoise(x_vec,SNR_vec(ii));
    y_hat_vec=addNoise(y_vec,SNR_vec(ii));
    r_xx_hat=xcorr(x_hat_vec,x_hat_vec); %unscaled linear auto-correlation for transient signals
    r_yy_hat=xcorr(y_hat_vec,y_hat_vec); %unscaled linear auto-correlation for transient signals
    r_xy_hat=xcorr(y_hat_vec,x_hat_vec); %unscaled linear cross-correlation for transient signals

    plot(x_ax,t_vec,x_hat_vec,'LineWidth',line_width)
    plot(y_ax,t_vec,y_hat_vec,'LineWidth',line_width)
    plot(r_xx_hat_ax,tau_vec,r_xx_hat,'LineWidth',line_width)
    plot(r_yy_hat_ax,tau_vec,r_yy_hat,'LineWidth',line_width)
    plot(r_xy_hat_ax,tau_vec,r_xy_hat,'LineWidth',line_width)
end
legend(x_ax,legend_str,'Location','northeast','interpreter','latex')

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove');

export_figure(gcf,'==',{'Echo_Correlation'})