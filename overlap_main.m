clearvars
close all
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1;0,0,0])
% set(groot,'DefaultAxesClipping','off')    %This is needed if you commented "win_2_eff_col=win_2_eff_col./N_vec;"


K=1200;
T=1;
[D_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_col=(0:K-1).'*D_t;
T_max=8*T;

win_col=window(@hann,K,'periodic')*2;
win_name='Hann';

alpha_vec=[0,0.5,2/3,0.75];
N_alpha_vec=length(alpha_vec);
legend_str_vec=strings(N_alpha_vec,1);
subplot(N_alpha_vec+1,1,N_alpha_vec+1);hold on;grid
for ii=1:N_alpha_vec
    if alpha_vec(ii)==0
        legend_str_vec(ii)="$\alpha=0$";
    else
        [numerator,denomerator]=rat(alpha_vec(ii));   %Display the floating alpha_vec(ii) as a fraction
        legend_str_vec(ii)="$\alpha=\frac{"+numerator+'}{'+denomerator+'}$';
    end
    
    subplot(N_alpha_vec+1,1,ii)
    if ii==1
        title(['Overlapped ',win_name,' windows'],'interpreter','latex')
    end
    if ii<=N_alpha_vec
        set(gca,'XTickLabel',[]);
    end

    yyaxis left
    set(gca,'LineStyleOrder','-');hold on;grid
    
    K_o=alpha_vec(ii)*K+1;
    if mod(K_o,1)~=0
        error('Overlap fraction cannot be realized. Consider changing K.'),
    end
    P=floor((T_max-alpha_vec(ii)*T)/(1-alpha_vec(ii))/T);
    K_tot=P*(K-K_o+1)+K_o-1;
    win_2_eff_col=zeros(K_tot,1);
    N_vec=zeros(size(win_2_eff_col));
    t_tot_col=(0:K_tot-1).'*D_t;
    for p=1:P
        k_start=(p-1)*(K-K_o+1)+1;
        p_sig_ind=(0:K-1)+k_start;
        plot(t_tot_col(p_sig_ind)/T,win_col)

        win_2_eff_col(p_sig_ind)=win_2_eff_col(p_sig_ind)+win_col.^2;
        N_vec(p_sig_ind)=N_vec(p_sig_ind)+ones(K,1);
    end
    ylabel('$w(t)$','interpreter','latex')
    
    win_2_eff_col=win_2_eff_col./N_vec;
    yyaxis right
    plot(t_tot_col/T,sqrt(win_2_eff_col))
    ylabel('$w_{\mathrm{eff}}(t)$','interpreter','latex')

    subplot(N_alpha_vec+1,1,N_alpha_vec+1)
    plot(t_tot_col/T,N_vec)
end

subplot(N_alpha_vec+1,1,N_alpha_vec+1)
xlabel('$t/T$','interpreter','latex')
ylabel('$N(t)$','interpreter','latex')
h_leg=legend(legend_str_vec,'interpreter','latex','Orientation','horizontal','Units','normalized');

%Improve figure display
for ii=1:N_alpha_vec
    subplot(N_alpha_vec+1,1,ii)
    xlim([0,T_max]);
    
    yyaxis left
    ylims=ylim;
    
    yyaxis right
    ylim(ylims)
    
    text(.985,.9,legend_str_vec(ii),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','EdgeColor',[0,0,0],'BackgroundColor',[1,1,1],'interpreter','latex')
end
subplot(N_alpha_vec+1,1,N_alpha_vec+1)
xlim([0,T_max]);
ylim([0,inf]);   %fix the lower ylim to 0

pos=get(h_leg,'Position');
pos=[1-pos(3),pos(4),pos(3:4)];
set(h_leg,'Position',pos)

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesClipping','remove')

export_figure(gcf,'==',"OverlapEffectivPowerWeighting")
