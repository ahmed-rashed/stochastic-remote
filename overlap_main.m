clearvars
close all
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1;0,0,0])
set(groot,'DefaultAxesClipping','off')

P=8;
K=1200;
T=1;
[D_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_col=(0:K-1).'*D_t;

a_win=[1,0,0,0,0;1,-1,0,0,0;1,-1.24,.244,-.00305,0;1,-1.933,1.286,-.388,.032];
win_name={'Rectangular','Hann','Kaiser-Bessel','Flat top'};
i_win=2;  %Hann window
win_col=(a_win(i_win,:)*cos((0:2:8).'*pi*t_col.'/T)).';
% win=hann(K+1,'symmetric');

alpha_vec=[0,0.5,2/3,0.75];
N_alpha_vec=length(alpha_vec);
legend_str_vec=cell(N_alpha_vec,1);
subplot(N_alpha_vec+1,1,N_alpha_vec+1);hold on;grid
for ii=1:N_alpha_vec
    legend_str_vec{ii}=['$\alpha=',num2str(alpha_vec(ii)*100),'\%$'];
    
    subplot(N_alpha_vec+1,1,ii)
    if ii==1
        title(['Overlapped ',win_name{i_win},' windows'], 'interpreter', 'latex')
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
    ylabel('$w(t)$', 'interpreter', 'latex')
    
    win_2_eff_col=win_2_eff_col./N_vec;
    yyaxis right
    plot(t_tot_col/T,win_2_eff_col)
    ylabel('$w_{\mathrm{eff}}^{2}(t)$', 'interpreter', 'latex')

    subplot(N_alpha_vec+1,1,N_alpha_vec+1)
    plot(t_tot_col/T,N_vec)
end

subplot(N_alpha_vec+1,1,N_alpha_vec+1)
xlabel('$t/T$', 'interpreter', 'latex')
ylabel('$N(t)$', 'interpreter', 'latex')
legend(legend_str_vec, 'interpreter', 'latex','Orientation','horizontal')

%Improve curves display
for ii=1:N_alpha_vec
    subplot(N_alpha_vec+1,1,ii)
    xlim([0,P*T]);
    yyaxis left
    vy=ylim;
    ylim(vy.^2)
    
    yyaxis right
    ylim(vy.^2)
    
    text(.985,.9,legend_str_vec{ii},'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', 'EdgeColor',[0,0,0], 'BackgroundColor',[1,1,1], 'interpreter', 'latex')
end
subplot(N_alpha_vec+1,1,N_alpha_vec+1)
xlim([0,P*T]);


set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesClipping','remove')

export_figure(gcf,'==',{'OverlapEffectivPowerWeighting'})