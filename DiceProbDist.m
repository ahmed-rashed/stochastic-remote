function DiceProbDist()
clc

x=1:13;
p=[[0:6],[5:-1:0]]/36;
P=cumsum(p);
stairs(x,P)

Y_Ax_TickLabels=strings(7,1);
Y_Ax_TickLabels(1)="0";
for ii=1:5
    Y_Ax_TickLabels(ii+1)=(ii*6)+"/36";
end
Y_Ax_TickLabels(7)="1";

axis([1,13,-.05,1.05])
set(gca,'XTick',[2:12],'YTick',[0:6:36]/36,'YTickLabel',Y_Ax_TickLabels)
grid on
xlabel('$x$','interpreter','latex')
title('$P(x)$','interpreter','latex')
export_figure(gcf,'',"Dice_ProbDist")