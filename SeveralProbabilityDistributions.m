function SeveralProbabilityDistributions()
close all
clc

N=2000;

rng(0);
samples_rows=[randn([1,N])      %Normaly distributed noise (White noise)
                rand([1,N])];   %Uniformly distributed noise
Titles={'Normally distributed noise','Uniformly distributed noise'};
x_lims=[-3,3];
N_bins=37;

equiSpacedintervals=linspace(x_lims(1),x_lims(2),N_bins);
midIntervals=(equiSpacedintervals(2:end)+equiSpacedintervals(1:end-1))/2;

N_signals=size(samples_rows,1);
for ii=1:N_signals
    subplot(3,N_signals,ii);
    plot(samples_rows(ii,:),1:N);
    title(Titles{ii})
    ylabel('samples')
    xlim(x_lims)
    set(gca,'XTick',equiSpacedintervals,'XGrid','On');
    set(gca,'XTickLabel',[]);

    [Prob,ProbDensity,ProbDist]=CalculateProbabilityFunctions(samples_rows(ii,:),equiSpacedintervals);    
    subplot(3,N_signals,2*N_signals-N_signals+ii);
    bar(midIntervals,ProbDensity,'hist');
    ylabel('$p(x)$', 'interpreter', 'latex')
    xlim(x_lims)
    set(gca,'XTick',equiSpacedintervals);
    set(gca,'XTickLabel',[]);

    subplot(3,N_signals,3*N_signals-N_signals+ii);
    stairs([x_lims(1),midIntervals],[0,ProbDist]);
    ylabel('$P(x)$', 'interpreter', 'latex')
    xlabel('$x$', 'interpreter', 'latex')
    xlim(x_lims)
end

export_figure(1,'==',{'WhiteNoise'})
