function overlap_main()
close all
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;0,0,0;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultAxesLineStyleOrder','-|--|-.')
set(groot,'DefaultLineMarkerSize',5);
set(groot,'DefaultLineLineWidth',1);

Q=4;
K=1024;
T=1;
[D_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_vec=(0:K).'*D_t;

aa=[1,0,0,0,0;1,1,0,0,0;1,1.24,.244,.00305,0;1,1.933,1.286,.388,.032];
a=aa(2,:);
win=a(1)-a(2)*cos(2*pi*t_vec/T)+a(3)*cos(4*pi*t_vec/T)-a(4)*cos(6*pi*t_vec/T)+a(5)*cos(8*pi*t_vec/T);
% win=hann(K+1,'symmetric');

figure
alpha_vec=[0,0.5,2/3,0.75];
ii=1;
for alpha=alpha_vec
    K_overlap=round(alpha*K+1);
    if K_overlap~=alpha*K+1,warning('check here'),end

    K_total=Q*(K-K_overlap+1)+K_overlap;
    win_2_effective=zeros(K_total,1);
    t_total_vec=(0:K_total-1)*D_t;

    subplot(length(alpha_vec),2,2*ii-1)
    hold on
    for q=1:Q
        k_start=(q-1)*(K-K_overlap+1)+1;
        t_vec=((k_start:k_start+K)-1)*D_t;

        subplot(length(alpha_vec),2,2*ii-1)
        plot(t_vec/T,win)

        win_2_effective(k_start:k_start+K)=win_2_effective(k_start:k_start+K)+win.^2;
    end
    win_2_effective=win_2_effective/Q;
%     if ii<length(OL_vec)
%         set(gca,'XTickLabel',[]);
%     else
        xlabel('$t/T$', 'interpreter', 'latex')
%     end
    
    if ii==1
        title('Overlapped weighting functions, $w(t)$', 'interpreter', 'latex')
    end
    ylabel([num2str(alpha*100),'% overlap'])
    
    subplot(length(alpha_vec),2,2*ii)
    plot(t_total_vec/T,win_2_effective)
    grid
%     if ii<length(OL_vec)
%         set(gca,'XTickLabel',[]);
%     else
        xlabel('$t/T$', 'interpreter', 'latex')
%     end
    
    if ii==1
        title('Effective power weighting, $w_{\textrm{eff}}^{2}(t)$', 'interpreter', 'latex')
    end
    ylabel([num2str(alpha*100),'% overlap'])

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
set(groot,'DefaultAxesLineStyleOrder','remove')
set(groot,'DefaultLineMarkerSize','remove');
set(groot,'DefaultLineLineWidth','remove');

export_figure(gcf,'==',{'OverlapEffectivPowerWeighting'})