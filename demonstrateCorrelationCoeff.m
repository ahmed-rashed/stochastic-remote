function demonstrateCorrelationCoeff()
close all
clc

N=100;

rng(0);
x_row=randn(1,N);
mu_x=mean(x_row);
StdDev_x=std(x_row);    %Unbiased covariance

noise_row=randn(1,N); % Uncorrelated noise

a_vec=[2:-1:-2].';
b_vec=[2;3;4;5;6];
if any(size(a_vec)~=size(b_vec)),error('a_col & b_col must have the same size'),end
N_a=length(a_vec);

figure
for i_row=1:3
    for i_a=1:N_a
        if i_row==1
            y_row=a_vec(i_a)*x_row;
            if a_vec(i_a)==1
                ylabel_text='$y=x$';
            elseif a_vec(i_a)==-1
                ylabel_text='$y=-x$';
            elseif a_vec(i_a)==0
                ylabel_text='$y=0$';
            else
                ylabel_text=['$y=',num2str(a_vec(i_a)),'x$'];
            end
        elseif i_row==2
            y_row=a_vec(i_a)*x_row+noise_row;
            if a_vec(i_a)==1
                ylabel_text='$y=x+\mathrm{noise}$';
            elseif a_vec(i_a)==-1
                ylabel_text='$y=-x+\mathrm{noise}$';
            elseif a_vec(i_a)==0
                ylabel_text='$y=\mathrm{noise}$';
            else
                ylabel_text=['$y=',num2str(a_vec(i_a)),'x+\mathrm{noise}$'];
            end
        else
            y_row=x_row.^(b_vec(i_a))/10;
            if b_vec(i_a)~=1
                ylabel_text=['$y=0.1x^{',num2str(b_vec(i_a)),'}$'];
            else
                ylabel_text='$y=0.1x$';
            end
        end
        
        mu_y=mean(y_row);
        StdDev_y=std(y_row);    %Unbiased covariance

        c_xy=sum((x_row-mu_x).*(y_row-mu_y))/(N-1);  %Unbiased covariance

        CorrCoeff_xy=c_xy/(StdDev_x*StdDev_y);

        subplot(3,N_a,(i_row-1)*N_a+i_a)
        plot(x_row,y_row,'.')
        set(gca,'XTick',[],'YTick',[],'Box','Off')
        axis equal
        axis tight
        xlabel('$x$', 'interpreter', 'latex', 'FontSize',10);
        ylabel(ylabel_text, 'interpreter', 'latex', 'FontSize',10);
        title(['$\rho_{xy}=$',num2str(CorrCoeff_xy)], 'interpreter', 'latex', 'FontSize',10)
    end
end

export_figure(gcf,'==',{'CorrelationCoeff'})