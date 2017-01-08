function overlap_main()
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;0,0,0;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultAxesLineStyleOrder','-|--|-.')
set(groot,'DefaultLineMarkerSize',5);
set(groot,'DefaultLineLineWidth',1);

Q=4;
K=1024;
T=1;
[D_t,f_s,D_f]=samplingParameters_T_N(T,K);

win=hann(K+1,'symmetric');

figure
OL_vec=[0,0.5,2/3,0.75];
ii=1;
for OL=OL_vec
    K_overlap=ceil(OL*(K+1))

    K_total=Q*(K-K_overlap+1)+K_overlap;
    win_2_effective=zeros(K_total,1);
    t_total_vec=(0:K_total-1)*D_t;

    subplot(length(OL_vec),2,2*ii-1)
    hold on
    for q=1:Q
        k_start=(q-1)*(K-K_overlap+1)+1;
        t_vec=((k_start:k_start+K)-1)*D_t;

        subplot(length(OL_vec),2,2*ii-1)
        plot(t_vec/T,win)

        win_2_effective(k_start:k_start+K)=win_2_effective(k_start:k_start+K)+win.^2;
    end
%     if ii<length(OL_vec)
%         set(gca,'XTickLabel',[]);
%     else
        xlabel('$t/T$', 'interpreter', 'latex')
%     end
    
    if ii==1
        title('Overlapped weighting functions')
    end
    ylabel([num2str(OL*100),'% overlap'])
    
    subplot(length(OL_vec),2,2*ii)
    plot(t_total_vec/T,win_2_effective)
%     if ii<length(OL_vec)
%         set(gca,'XTickLabel',[]);
%     else
        xlabel('$t/T$', 'interpreter', 'latex')
%     end
    
    if ii==1
        title('Effective power weighting')
    end
    ylabel([num2str(OL*100),'% overlap'])

    ii=ii+1;
end

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesLineStyleOrder','remove')
set(groot,'DefaultLineMarkerSize','remove');
set(groot,'DefaultLineLineWidth','remove');