clc
close all
clearvars

set(groot,'DefaultAxesColorOrder',[1,0,0;0,0,1;0,0.5,0;1,0,1;0,0,0])

%DFT parameters
T=4;
K=2^10; %K=1024
[Delta_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_row=(0:K-1)*Delta_t;
tau_row=(-(K-1):K-1)*Delta_t;

%x(t) parameters
T_burst=T/8;

%y(t) parameters
alpha_col=[.2;.1;.2;.1];
T_echo_col=[2;4;6;6.5]*T_burst;

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
s_title_cvec={'Tone pulse','Hann-weighted-tone pulse','LFM pulse','Random pulse'};

SNR_vec=[1/2,inf];

N_SNR=length(SNR_vec);
legend_str=cell(N_SNR,1);

iii=1;
for s_fn=s_fn_cvec
    x_row=signalPulse(t_row,T_burst,s_fn{1});
    
    figure
    x_ax=subplot(3,1,1);
    plot(t_row/T,x_row,'LineWidth',1.5)
    set(x_ax,'XGrid','on','XTick',0:.2:1,'XTickLabel',[],'XLim',[0,1],'YLimSpec','tight')
    title(s_title_cvec{iii})
    ylabel(x_ax,'$x(t)$','interpreter','latex')
    %Additional optimization of the axes for correct comparison with the correlation curves
    pos=get(x_ax,'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(x_ax,'Position',pos)

    y_ax=subplot(3,1,2);
    hold on
    set(y_ax,'XGrid','on','XTick',0:.2:1,'XTickLabel',[],'XLim',[0,1],'YLimSpec','tight')
    ylabel(y_ax,'$\hat{y}(t)$','interpreter','latex')
    xlabel(y_ax,'$t/T$','interpreter','latex');
    %Additional optimization of the axes for correct comparison with the correlation curves
    pos=get(y_ax,'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(y_ax,'Position',pos)

    r_xy_hat_ax=subplot(3,1,3);
    hold on;
    set(r_xy_hat_ax,'XGrid','on')
    set(r_xy_hat_ax,'XTick',-1:.2:1,'XLim',[-1,1],'YLimSpec','tight')
    ylabel(r_xy_hat_ax,'$r_{x\hat{y}}(\tau)$','interpreter','latex')
    xlabel(r_xy_hat_ax,'$\tau/T$','interpreter','latex')
    
    y_row=x_row+sum(alpha_col.*signalPulse(t_row-T_echo_col,T_burst,s_fn{1}),1);

    for ii=1:N_SNR
        if SNR_vec(ii)==inf
            legend_str{ii}='$\mathrm{SNR}=\infty$';
            line_width=1.5;
        else
            legend_str{ii}=['$\mathrm{SNR}=',num2str(SNR_vec(ii)),'$'];
            line_width=1;
        end

        rng(1);
        y_hat_row=addNoise(y_row,SNR_vec(ii));
        plot(y_ax,t_row/T,y_hat_row,'LineWidth',line_width)
        
        r_x_y_hat=xcorr(y_hat_row,x_row); %unscaled linear cross-correlation
        plot(r_xy_hat_ax,tau_row/T,r_x_y_hat,'LineWidth',line_width)
    end
    legend(y_ax,legend_str,'Location','southeast','interpreter','latex')
    
    iii=iii+1;
end

set(groot,'DefaultAxesColorOrder','remove')

export_figure(1:4,'==',{'HarmonicPulseRadar','WeightedHarmonicPulseRadar','LFM_PulseCompressionRadar','RandomPulseRadar'})