clearvars
% close all
figure
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;0,0,0;1,0,0;0,0.5,0;1,0,1])
%set(groot,'DefaultLineLineWidth',1);

P=6;
K=2^7; %K=1024
T=1;
[D_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_col=(0:K-1).'*D_t;

a_win=[1,0,0,0,0;1,-1,0,0,0;1,-1.24,.244,-.00305,0;1,-1.933,1.286,-.388,.032];
win_name={'Rectangular','Hann','Kaiser-Bessel','Flat top'};
i_win=2;  %Hann window
win_col=(a_win(i_win,:)*cos((0:2:8).'*pi*t_col.'/T)).';
% win=hann(K+1,'symmetric');

alpha_vec=[0,0.5,2/3,0.75];
ii=1;
for alpha=alpha_vec
    K_o=round(alpha*K+1);
    if K_o~=alpha*K+1
        warning('check here'),
    end

    K_tot=P*(K-K_o+1)+K_o-1;
    win_2_eff_col=zeros(K_tot,1);
    t_tot_col=(0:K_tot-1).'*D_t;

    subplot(length(alpha_vec),2,2*ii-1)
    hold on
    grid
    for p=1:P
        k_start=(p-1)*(K-K_o+1)+1;
        p_sig_ind=(0:K-1)+k_start;
        t_col=t_tot_col(p_sig_ind);
        plot(t_col/T,win_col)

        win_2_eff_col(p_sig_ind)=win_2_eff_col(p_sig_ind)+win_col.^2;
    end
    win_2_eff_col=win_2_eff_col/P;
    if ii<length(alpha_vec)
        set(gca,'XTickLabel',[]);
    else
        xlabel('$t/T$', 'interpreter', 'latex')
    end
    
    if ii==1
        title('Overlapped windows, $w(t)$', 'interpreter', 'latex')
    end
    ylabel(['$\alpha=',num2str(alpha*100),'\%$'], 'interpreter', 'latex')
    
    subplot(length(alpha_vec),2,2*ii)
    plot(t_tot_col/T,win_2_eff_col)
    grid
    if ii<length(alpha_vec)
        set(gca,'XTickLabel',[]);
    else
        xlabel('$t/T$', 'interpreter', 'latex')
    end
    
    if ii==1
        title(['$w_{\mathrm{eff}}^{2}(t)$ for ',win_name{i_win},' window'], 'interpreter', 'latex')
    end
%     ylabel(['$\alpha=',num2str(alpha*100),'\%$'], 'interpreter', 'latex')

    ii=ii+1;
end

%Additional formatting of the x-axis
subplot(length(alpha_vec),2,1)
v=xlim;
for ii=2:length(alpha_vec)
    subplot(length(alpha_vec),2,2*ii-1)
    v=max(v,xlim);
end

for ii=1:length(alpha_vec)*2
    subplot(length(alpha_vec),2,ii)
    xlim(v);
end


set(groot,'DefaultAxesColorOrder','remove')
%set(groot,'DefaultLineLineWidth','remove');

export_figure(gcf,'==',{'OverlapEffectivPowerWeighting'})