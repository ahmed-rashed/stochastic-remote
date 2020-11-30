function Idependence_vs_LinearUncorrelation()
close all
clc

N=100;

rng(0);
x_row=2*rand(1,N)-1;
mu_x=mean(x_row);
StdDev_x=std(x_row);    %Unbiased covariance

noise_row=randn(1,N); % Uncorrelated noise

a_col=[2:-1:-2].';
N_a=length(a_col);

figure
for ii=1:3
    for i_a=1:N_a
        if ii==1
            y_row=a_col(i_a)*x_row;
        elseif ii==2
            y_row=a_col(i_a)*x_row+noise_row;
        else
            y_row=x_row.^(i_a+1)/10;
        end
        
        mu_y=mean(y_row);
        StdDev_y=std(y_row);    %Unbiased covariance

        c_xy=sum((x_row-mu_x).*(y_row-mu_y))/(N-1);  %Unbiased covariance

        CorrCoeff_xy=c_xy/(StdDev_x*StdDev_y);

        subplot(3,N_a,(ii-1)*N_a+i_a)
        plot(x_row,y_row,'.')
        set(gca,'XTick',[],'YTick',[],'Box','Off')
        axis equal
        axis tight
        xlabel('$x$','interpreter','latex','FontSize',8);
        ylabel('$y$','interpreter','latex','FontSize',8);
        title(['$\rho_{xy}=$',num2str(CorrCoeff_xy)],'interpreter','latex','FontSize',8)
    end
end

export_figure(gcf,'==',{'CorrelationCoeff'})