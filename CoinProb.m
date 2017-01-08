function CoinProb()

rng(0);
K_max=3e2;

%x=round(rand(1,K_max)); % 1: head, 0: tail
x=randi(2,1,K_max)-1; % 1: head, 0: tail

K_vec=1:K_max;
head= x==1;
tail= x==0;

n_head_vec=cumsum(head);
n_tail_vec=cumsum(tail);

prob_H_hat=n_head_vec./K_vec;
prob_T_hat=n_tail_vec./K_vec;

subplot(3,1,1)
plot(x,'.-')
ylabel('Coin Face', 'interpreter', 'latex')
ylim([-.1 1.1])
set(gca,'XTickLabel',[],'YTick',[0,1],'YTickLabel',{'H','T'});
grid

subplot(3,1,2)
plot(prob_H_hat)
ylabel('$\widehat{\mathrm{Prob}}\left(\textrm{H}\right)$', 'interpreter', 'latex')
ylim([0 1])
set(gca,'YTick',[0:.1:1])
set(gca,'XTickLabel',[]);
grid

subplot(3,1,3)
plot(prob_T_hat)
xlabel('$K$; Number of trials', 'interpreter', 'latex')
ylabel('$\widehat{\mathrm{Prob}}\left(\textrm{T}\right)$', 'interpreter', 'latex')
ylim([0 1])
set(gca,'YTick',[0:.1:1])
grid

export_figure(1,'||',{'RelativeFrequency'})