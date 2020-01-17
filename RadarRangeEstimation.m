clc
close all
clearvars

set(groot,'DefaultAxesColorOrder',[1,0,0;0,0,1;0,0.5,0;1,0,1;0,0,0])

%DFT parameters
T=4;
K=2^10; %K=1024
[Delta_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_vec=(0:K-1)*Delta_t;
tau_vec=(-(K-1):K-1)*Delta_t;

%x(t) parameters
T_burst=T/8;

%y(t) parameters
alpha=.2;
T_echo=2*T_burst;

%Tone and Chirp signal parameters
f_0=f_s/24;
T_1=T_burst;
f_1=f_s/4;

%Random signal parameters
K_burst=T_burst/Delta_t;
t_rnd_vec=(0:K_burst-1)*Delta_t;
if rem(K_burst,1)~=0,error('T_burst should be selected as integer multiples of Delta_t'),end
rng(11);
x_rnd=randn(1,K_burst);
f_c=100;    %cut off frequency
[b,a]=butter(9,f_c/(f_s/2));
x_rnd=filtfilt(b,a,x_rnd); 
x_rnd=x_rnd-mean(x_rnd);
x_rnd=x_rnd/std(x_rnd); % Makes mean(s)=0 & std(s)=1;

s_fn_cvec={@(t) (cos(2*pi*f_0*t)),@(t) (cos(2*pi*f_0*t-pi).*(1+cos(2*pi*t/T_burst-pi))),@(t) (chirp(t,f_0,T_1,f_1)),@(t) (interp1(t_rnd_vec,x_rnd,t,'linear',0))};

SNR_vec=[1/2,inf];

N_SNR=length(SNR_vec);
legend_str=cell(N_SNR,1);

for s_fn=s_fn_cvec
    figure
    x_ax=subplot(3,1,1);
    hold on
    set(x_ax,'XGrid','on','YLimSpec','tight','XTickLabel',[])
    ylabel(x_ax,'$\hat{x}(t)$','interpreter','latex')
    %Additional optimization of the axes for correct comparison with the correlation curves
    pos=get(x_ax,'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(x_ax,'Position',pos)

    y_ax=subplot(3,1,2);
    hold on
    set(y_ax,'XGrid','on','YLimSpec','tight','XTickLabel',[])
    ylabel(y_ax,'$\hat{y}(t)$','interpreter','latex')
    xlabel(y_ax,'$t$ (s)','interpreter','latex');
    %Additional optimization of the axes for correct comparison with the correlation curves
    pos=get(y_ax,'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(y_ax,'Position',pos)

    r_xy_hat_ax=subplot(3,1,3);
    hold on;
    set(r_xy_hat_ax,'XGrid','on')
    set(r_xy_hat_ax,'YLimSpec','tight')
    ylabel(r_xy_hat_ax,'$r_{x\hat{y}}(\tau)$','interpreter','latex')
    xlabel(r_xy_hat_ax,'$\tau$ (s)','interpreter','latex')

    for ii=1:N_SNR
        if SNR_vec(ii)==inf
            legend_str{ii}='$\mathrm{SNR}=\infty$';
            line_width=1.5;
        else
            legend_str{ii}=['$\mathrm{SNR}=',num2str(SNR_vec(ii)),'$'];
            line_width=1;
        end
        
        x_vec=signalPulse(t_vec,T_burst,s_fn{1});
        y_vec=x_vec+alpha*signalPulse(t_vec-T_echo,T_burst,s_fn{1});

        rng(1);
        x_hat_vec=addNoise(x_vec,SNR_vec(ii));
        y_hat_vec=addNoise(y_vec,SNR_vec(ii));
        r_x_y_hat=xcorr(y_hat_vec,x_vec); %unscaled linear cross-correlation

        plot(x_ax,t_vec,x_hat_vec,'LineWidth',line_width)
        plot(y_ax,t_vec,y_hat_vec,'LineWidth',line_width)
        plot(r_xy_hat_ax,tau_vec,r_x_y_hat,'LineWidth',line_width)
    end
    legend(x_ax,legend_str,'Location','southeast','interpreter','latex')
end

set(groot,'DefaultAxesColorOrder','remove')

% export_figure(gcf,'==',{'Echo_Correlation'})