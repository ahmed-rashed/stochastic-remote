function NoisySignalAveraging()
clc
close all

N=2^9
n_f_max=floor(N/2);
T=5
[dt,f_s,df]=samplingParameters_T_N(T,N);
f0=10*df;

f_c_by_f_s=.25;    %Low pass filter cut off frequency ratio
SNR_db=-20;

N_avg_vec=[1, 1000, 10000];
%%%%%%%%%%%%%%%%%%%%%%%
t=(0:N-1).'*dt;
f=(0:N-1).'*df;

N_pad=2^nextpow2(2*N-1)
T_pad=N_pad/N*T;
[dt_pad,f_s_pad,df_pad]=samplingParameters_T_N(T_pad,N_pad);
t_pad=(0:N_pad-1).'*dt_pad;
f_pad=(0:N_pad-1).'*df_pad;

f_c=f_c_by_f_s*f_s;   %Low pass filter cut off frequency
[b,a]=butter(9, f_c/f_s);    %designs a ninth-order low-pass digital Butterworth filter (IIR), where b is a vector containing coefficients of a moving average part and a is a vector containing coefficients of an auto-regressive part of the transfer function [Shin; Equation (6.12)].

N_N_avg=length(N_avg_vec);
signal_expectation_current=zeros(N,1);
signal_expectation_cols=zeros(N,N_N_avg);
r_ss_padded_current=zeros(N_pad,1);
r_ss_padded_cols=zeros(N_pad,N_N_avg);
nn=1;
rng(10);
for N_avg=1:N_avg_vec(end)
    singnalPure_i=sin(2*pi*f0*t+randn*T);   %Pure sine but with random start phase
    noise_i=filtfilt(b,a,randn(size(t)));   %performs zero-phase digital filtering. The resulting sequence N_avg is the band-limited (zero to f_c) white noise.
    noise_i=std(singnalPure_i)/std(noise_i)/sqrt(db2lin(SNR_db))*noise_i;  %Scales the noise so that the SNR is as given
    signal_i=singnalPure_i+noise_i;
    
    r_nn_i_padded=xcorr_xspec_linear(signal_i,signal_i,'unbiased');

    signal_expectation_current=(1-1/N_avg)*signal_expectation_current+signal_i/N_avg;  %equivalent to sum(signal_i)/N_avg, for i =1 to i=N_avg
    
    r_ss_padded_current=(1-1/N_avg)*r_ss_padded_current+r_nn_i_padded/N_avg;        %equivalent to sum(r_nn_i_padded)/N_avg, for i =1 to i=N_avg

    if N_avg==N_avg_vec(nn)
        signal_expectation_cols(:,nn)=signal_expectation_current;
        r_ss_padded_cols(:,nn)=r_ss_padded_current;

        nn=nn+1;
    end
end











%The following is for plotting and annotation
%You need not understand this

f_spectrum=figure;maximizeFigure(f_spectrum);
x_xlabel_Latex='$t/T$';
X_xlabel_Latex='$f/f_{\textrm{s}}$';
x_Title_Latex=cell(N_N_avg,1);
X_Title_Latex=cell(N_N_avg,1);
for ii=1:N_N_avg
    x_Title_Latex{ii}=['$E[x(t)]$,: $N = ',int2str(N),'$ \& $N_{\textrm{avg}} = ',int2str(N_avg_vec(ii)), '$'];
    X_Title_Latex{ii}=['$|E[X(f)]|$,: $N = ',int2str(N),'$ \& $N_{\textrm{avg}} = ',int2str(N_avg_vec(ii)), '$'];
end
X_VeticalLines=(f_c/2)/f_s_pad;
VerticalLinesTextLatex={'${\displaystyle f_{\textrm{c}} / f_{\textrm{s}}}$'};
plotFourierTransformPair(t/T,signal_expectation_cols,f(1:n_f_max)/f_s,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,[],[],[],[],[],[],X_VeticalLines,VerticalLinesTextLatex,true);

f_autoCorr1=figure;maximizeFigure(f_autoCorr1);
x_xlabel_Latex='$\tau/T^{\textrm{pad}}$ ,:$T^{\textrm{pad}}=N^{\textrm{pad}} / N \times T$';
X_xlabel_Latex='$f/f_{\textrm{s}}^{\textrm{pad}}$ ,:$f_{\textrm{s}}^{\textrm{pad}}=f_{\textrm{s}}$';
for ii=1:N_N_avg
    x_Title_Latex{ii}=['$r_{xx}^{\textrm{pad}}(\tau)$ (@ $\Delta t^{\textrm{pad}} = \Delta t$), $N^{\textrm{pad}} = ',int2str(N_pad),'$ \& $N_{\textrm{avg}} = ',int2str(N_avg_vec(ii)), '$'];
    X_Title_Latex{ii}=['$R_{XX}^{\textrm{pad}}(f)$ (@$\Delta f^{\textrm{pad}}=N/N^{\textrm{pad}} \times \Delta f$), $N^{\textrm{pad}} = ',int2str(N_pad),'$ \& $N_{\textrm{avg}} = ',int2str(N_avg_vec(ii)), '$'];
end
%signal_variance_row=r_ss_padded_cols(1,:);
y_HorizontalLines=[];
horizontalLinesTextLatex=[];
x_VeticalLines=T/T_pad;
verticalLinesTextLatex={'$T/T^{\textrm{pad}}$'};
Y_HorizontalLines=[];
HorizontalLinesTextLatex=[];
X_VeticalLines=[(f_c/2)/f_s_pad;(f_s/2)/f_s_pad];
VerticalLinesTextLatex={'${\displaystyle f_{\textrm{c}} / f_{\textrm{s}}^{\textrm{pad}}}$';'${\displaystyle (f_{\textrm{s}}^{\textrm{pad}} / 2) /f_{\textrm{s}}^{\textrm{pad}}= (f_{\textrm{s}} / 2) / f_{\textrm{s}}}$'};
plotFourierTransformPair(t_pad/T_pad,r_ss_padded_cols,f_pad/f_s_pad,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,y_HorizontalLines,horizontalLinesTextLatex,x_VeticalLines,verticalLinesTextLatex,Y_HorizontalLines,HorizontalLinesTextLatex,X_VeticalLines,VerticalLinesTextLatex,true);

f_autoCorr2=figure;maximizeFigure(f_autoCorr2);
x_xlabel_Latex='$\tau/T$  ,:$T^{\textrm{pad}}=N^{\textrm{pad}} / N \times T$';
X_xlabel_Latex='$f/f_{\textrm{s}}$  ,:$f_{\textrm{s}}^{\textrm{pad}}=f_{\textrm{s}}$';
for ii=1:N_N_avg
    x_Title_Latex{ii}=['$r_{xx}^{\textrm{pad}}(\tau)$ (@ $\Delta t^{\textrm{pad}} = \Delta t$), $N = ',int2str(N),'$ \& $N_{\textrm{avg}} = ',int2str(N_avg_vec(ii)), '$'];
    X_Title_Latex{ii}=['$R_{XX}^{\textrm{pad}}(f)$ (@$\Delta f^{\textrm{pad}}=N/N^{\textrm{pad}} \times \Delta f$), $N = ',int2str(N),'$ \& $N_{\textrm{avg}} = ',int2str(N_avg_vec(ii)), '$'];
end
x_VeticalLines=[];
verticalLinesTextLatex={};
X_VeticalLines=X_VeticalLines(1);
VerticalLinesTextLatex={'${\displaystyle f_{\textrm{c}} / f_{\textrm{s}}}$'};
plotFourierTransformPair(t/T,r_ss_padded_cols,f_pad(1:N)/f_s,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,y_HorizontalLines,horizontalLinesTextLatex,x_VeticalLines,verticalLinesTextLatex,Y_HorizontalLines,HorizontalLinesTextLatex,X_VeticalLines,VerticalLinesTextLatex,true);

%export_figure([f_spectrum,f_autoCorr1,f_autoCorr2],'==',{'SpectrumAvg','CorrAvg-details','CorrAvg'})