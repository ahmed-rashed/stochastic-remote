clc
close all
clearvars

N=100;

rng(1); % seeds the random number generator to obtain the same sequence of random numbers.
x_row=randn(1,N);
mu_x=mean(x_row);
sigma_x=std(x_row);    %Unbiased covariance

noise_row=randn(1,N); % Uncorrelated noise

a_vec=(2:-1:-2).';
b_vec=[2;3;4;5;6];
if any(size(a_vec)~=size(b_vec)),error('a_vec & b_vec must have the same size'),end
N_a=length(a_vec);

figure
for i_row=1:3
    for i_a=1:N_a
        if i_row==1
            y_row=a_vec(i_a)*x_row;
            if a_vec(i_a)==1
                ylabel_text="$y=x$";
            elseif a_vec(i_a)==-1
                ylabel_text="$y=-x$";
            elseif a_vec(i_a)==0
                ylabel_text="$y=0$";
            else
                ylabel_text="$y="+a_vec(i_a)+'x$';
            end
        elseif i_row==2
            y_row=a_vec(i_a)*x_row+noise_row;
            if a_vec(i_a)==1
                ylabel_text="$y=x+n$";
            elseif a_vec(i_a)==-1
                ylabel_text="$y=-x+n$";
            elseif a_vec(i_a)==0
                ylabel_text="$y=n$";
            else
                ylabel_text="$y="+a_vec(i_a)+'x+n$';
            end
        else
            y_row=.1*x_row.^(b_vec(i_a));
            if b_vec(i_a)~=1
                ylabel_text="$y=0.1x^{"+b_vec(i_a)+'}$';
            else
                ylabel_text="$y=0.1x$";
            end
        end
        
        mu_y=mean(y_row);
        sigma_y=std(y_row);    %Unbiased covariance
        c_xy=sum((x_row-mu_x).*(y_row-mu_y))/(N-1);  %Unbiased covariance
        CorrCoeff_xy=c_xy/(sigma_x*sigma_y);
        CorrCoeff_xy_mat=corrcoef(x_row,y_row);
%         CorrCoeff_xy_1=CorrCoeff_xy_mat(1,2);
%         error=CorrCoeff_xy-CorrCoeff_xy_1

        ax=subplot(3,N_a,(i_row-1)*N_a+i_a);
        plot(x_row,y_row,'.')
        set(ax,'Box','Off')
        axis equal
        axis tight
        xlabel('$x$','interpreter','latex');
        ylabel(ylabel_text,'interpreter','latex');
        title("$\rho_{xy}="+CorrCoeff_xy+'$','interpreter','latex')
    end
end

export_figure(gcf,'==',"CorrelationCoeff")