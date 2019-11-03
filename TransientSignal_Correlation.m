clc
close all
clearvars

set(groot,'DefaultAxesColorOrder',[1,0,0;0,0,1;0,0.5,0;1,0,1;0,0,0])
set(groot,'DefaultLineLineWidth',1);

%DFT parameters
T=4;
K=2^10; %K=1024
[D_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_vec=(0:K-1)*D_t;
tau=(-(K-1):K-1)*D_t;

%s(t) parameters
T_burst=T/8;
f_0=(f_s/2)/10;

%x(t) and y(t) parameters
a=1;b=.8; c=.75;
D_1=T/4;D_2=3/4*T;

x_vec=a*burst_tone(t_vec,f_0,T_burst);
y_vec=b*burst_tone(t_vec-D_1,f_0,T_burst)+c*burst_tone(t_vec-D_2,f_0,T_burst);

SNR_vec=[1/sqrt(2),inf];
N_SNR=length(SNR_vec);
legend_str=cell(N_SNR,1);

x_ax=subplot(5,1,1);hold on
set(x_ax,'XGrid','on','YLimSpec','tight','XTickLabel',[])
ylabel(x_ax,'$\hat{x}(t)$', 'interpreter', 'latex')
%Additional optimization of the axes for correct comparison with the correlation curves
pos=get(x_ax,'Position');
pos(1)=pos(1)+pos(3)/2;
pos(3)=pos(3)/2;
set(x_ax,'Position',pos)

y_ax=subplot(5,1,2);hold on
set(y_ax,'XGrid','on','YLimSpec','tight','XTickLabel',[])
ylabel(y_ax,'$\hat{y}(t)$', 'interpreter', 'latex')
xlabel(y_ax,'$t$ (sec.)', 'interpreter', 'latex');
%Additional optimization of the axes for correct comparison with the correlation curves
pos=get(y_ax,'Position');
pos(1)=pos(1)+pos(3)/2;
pos(3)=pos(3)/2;
set(y_ax,'Position',pos)

r_xx_hat_ax=subplot(5,1,3);hold on;set(r_xx_hat_ax,'XGrid','on')
set(r_xx_hat_ax,'XTickLabel',[],'YLimSpec','tight');
ylabel(r_xx_hat_ax,'$\hat{r}_{xx}(\tau)$', 'interpreter', 'latex')

r_yy_hat_ax=subplot(5,1,4);hold on;set(r_yy_hat_ax,'XGrid','on')
set(r_yy_hat_ax,'YLimSpec','tight','XTickLabel',[])
ylabel(r_yy_hat_ax,'$\hat{r}_{yy}(\tau)$', 'interpreter', 'latex')

r_xy_hat_ax=subplot(5,1,5);hold on;set(r_xy_hat_ax,'XGrid','on')
set(r_xy_hat_ax,'YLimSpec','tight')
ylabel(r_xy_hat_ax,'$\hat{r}_{xy}(\tau)$', 'interpreter', 'latex')
xlabel(r_xy_hat_ax,'$\tau$ (sec.)', 'interpreter', 'latex')

rng(0);
for ii=1:N_SNR
    if SNR_vec(ii)==inf
        legend_str{ii}='$\mathrm{SNR}=\infty$';
        line_width=get(groot,'DefaultLineLineWidth')*1.5;
    else
        legend_str{ii}=['$\mathrm{SNR}=',num2str(SNR_vec(ii)),'$'];
        line_width=get(groot,'DefaultLineLineWidth');
    end

    x_vec_hat=x_vec+std(x_vec)*randn(1,K)/SNR_vec(ii);
    y_vec_hat=y_vec+std(y_vec)*randn(1,K)/SNR_vec(ii);
    r_xx_hat=xcorr(x_vec_hat,x_vec_hat,'none'); %unscaled inear auto-correlation for transient signals
    r_yy_hat=xcorr(y_vec_hat,y_vec_hat,'none'); %unscaled linear auto-correlation for transient signals
    r_xy_hat=xcorr(y_vec_hat,x_vec_hat,'none'); %unscaled linear cross-correlation for transient signals
    %tau=(-(K-1):K-1)*D_t;

    plot(x_ax,t_vec,x_vec_hat,'LineWidth',line_width)
    plot(y_ax,t_vec,y_vec_hat,'LineWidth',line_width)
    plot(r_xx_hat_ax,tau,r_xx_hat,'LineWidth',line_width)
    plot(r_yy_hat_ax,tau,r_yy_hat,'LineWidth',line_width)
    plot(r_xy_hat_ax,tau,r_xy_hat,'LineWidth',line_width)
end
legend(r_xy_hat_ax,legend_str,'Location','southwest', 'interpreter', 'latex')

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove');

export_figure(gcf,'==',{'Echo_Correlation'})

function s_vec=burst_tone(t_vec,f_0,T_burst)
s_vec=zeros(size(t_vec));
mask=t_vec<=T_burst & t_vec>=0;
s_vec(mask)=sin(2*pi*f_0*t_vec(mask));

end

function s_vec=baseBandRandom(t_vec,t_before,t_after,f_c_by_f_s_2)

end
