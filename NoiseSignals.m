clearvars
close all
clc

N=2000;

rng(0);
signal_rows=[randn([1,N])      %Normaly distributed noise (White noise)
                rand([1,N])];   %Uniformly distributed noise
Titles={'Normally distributed noise','Uniformly distributed noise'};
y_lims=[-3,3];
N_bins=37;

equiSpacedintervals=linspace(y_lims(1),y_lims(2),N_bins);
midIntervals=(equiSpacedintervals(2:end)+equiSpacedintervals(1:end-1))/2;

N_signals=size(signal_rows,1);
for n=1:N_signals
    ax=subplot(N_signals,4,(n-1)*4+(1:3));
    plot(signal_rows(n,:));
    title(Titles{n})
    xlabel('samples')
    ylim(y_lims);
    ax.YAxis.MinorTickValues=equiSpacedintervals;
    ax.YAxis.MinorTick='on';
    ax.YMinorGrid='on';

    [Prob,ProbDensity]=CalculateProbabilityFunctions(signal_rows(n,:),equiSpacedintervals);    
    
    ax=subplot(N_signals,4,(n-1)*4+4);
    barh(midIntervals,ProbDensity,'hist');
    title('$p(x)$', 'interpreter', 'latex')
    xlim([0,1.1]);
    ylim(y_lims);
    ax.YAxis.MinorTickValues=equiSpacedintervals;
    ax.YAxis.MinorTick='on';
end

export_figure(1,'==',{'WhiteNoise'})
