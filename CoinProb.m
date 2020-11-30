clc
clearvars
close all

N_max=3e2;

rng(0);
%x_row=round(rand(1,N_max)); % 1: head_row,0: tail_row
x_row=randi(2,1,N_max)-1; % 1: head_row,0: tail_row

N_vec=1:N_max;
head_row= x_row==1;
tail_row= x_row==0;

n_head_vec=cumsum(head_row);
n_tail_vec=cumsum(tail_row);

prob_H_hat=n_head_vec./N_vec;
prob_T_hat=n_tail_vec./N_vec;

subplot(3,1,1)
plot(x_row,'.-')
ylabel('Coin Face','interpreter','latex')
ylim([-.1 1.1])
set(gca,'XTickLabel',[],'YTick',[0,1],'YTickLabel',{'H','T'});
grid

subplot(3,1,2)
plot(prob_H_hat)
ylabel('$\widehat{\mathrm{Prob}}\left(\mathrm{H}\right)$','interpreter','latex')
ylim([0 1])
set(gca,'YTick',[0:.1:1])
set(gca,'XTickLabel',[]);
grid

subplot(3,1,3)
plot(prob_T_hat)
xlabel('$N$; Number of trials','interpreter','latex')
ylabel('$\widehat{\mathrm{Prob}}\left(\mathrm{T}\right)$','interpreter','latex')
ylim([0 1])
set(gca,'YTick',[0:.1:1])
grid

export_figure(1,'||',{'RelativeFrequency'})