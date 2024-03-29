clc
close all
clearvars

N=2^9
n_f_max=floor(N/2);
T=5
f_c_by_f_s=.25;    %Low pass filter cut off frequency ratio
N_avg_vec=[1,10,1000];
%%%%%%%%%%%%%%%%%%%%%%%
[dt,f_s,df]=samplingParameters_T_N(T,N);
t=(0:N-1).'*dt;
f=(0:N-1).'*df;

N_pad=2^nextpow2(2*N-1)
T_pad=N_pad/N*T;
[dt_pad,f_s_pad,df_pad]=samplingParameters_T_N(T_pad,N_pad);
t_pad=(0:N_pad-1).'*dt_pad;
f_pad=(0:N_pad-1).'*df_pad;

f_c=f_c_by_f_s*f_s;   %Low pass filter cut off frequency
[b,a]=butter(9,f_c/f_s);    %designs a ninth-order low-pass digital Butterworth filter (IIR), where b is a vector containing coefficients of a moving average part and a is a vector containing coefficients of an auto-regressive part of the transfer function [Shin; Equation (6.12)].

N_N_avg=length(N_avg_vec);
noise_expectation_current=zeros(N,1);
noise_expectation_cols=zeros(N,N_N_avg);
r_nn_padded_current=zeros(N_pad,1);
r_nn_padded_cols=zeros(N_pad,N_N_avg);
nn=1;
for N_avg=1:N_avg_vec(end)
    noise_i=randn(size(t));
    noise_i=filtfilt(b,a,noise_i);    %performs zero-phase digital filtering. The resulting sequence N_avg is the band-limited (zero to f_c) white noise.
    
    r_nn_i_padded=xcorr_xspec_linear(noise_i,noise_i,'unbiased');

    noise_expectation_current=(1-1/N_avg)*noise_expectation_current+noise_i/N_avg;  %equivalent to sum(noise_i)/N_avg, for i =1 to i=N_avg
    
    r_nn_padded_current=(1-1/N_avg)*r_nn_padded_current+r_nn_i_padded/N_avg;        %equivalent to sum(r_nn_i_padded)/N_avg, for i =1 to i=N_avg

    if N_avg==N_avg_vec(nn)
        noise_expectation_cols(:,nn)=noise_expectation_current;
        r_nn_padded_cols(:,nn)=r_nn_padded_current;

        nn=nn+1;
    end
end











%The following is for plotting and annotation
%You need not understand this

f_spectrum=figure;maximizeFigure(f_spectrum);
x_xlabel_Latex='$t/T$';
X_xlabel_Latex='$f/f_{\mathrm{s}}$';
x_Title_Latex=strings(N_N_avg,1);
X_Title_Latex=strings(N_N_avg,1);
for ii=1:N_N_avg
    x_Title_Latex(ii)="$E[x(t)]$,: $N="+N+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
    X_Title_Latex(ii)="$|E[X(f)]|$,: $N="+N+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
end
X_VeticalLines=(f_c/2)/f_s_pad;
VerticalLinesTextLatex="${\displaystyle f_{\mathrm{c}} / f_{\mathrm{s}}}$";
plotFourierTransformPair(t/T,noise_expectation_cols,f(1:n_f_max)/f_s,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,[],[],[],[],[],[],X_VeticalLines,VerticalLinesTextLatex,true);

f_autoCorr1=figure;maximizeFigure(f_autoCorr1);
x_xlabel_Latex='$\tau/T^{\mathrm{pad}}$ ,:$T^{\mathrm{pad}}=N^{\mathrm{pad}} / N \times T$';
X_xlabel_Latex='$f/f_{\mathrm{s}}^{\mathrm{pad}}$ ,:$f_{\mathrm{s}}^{\mathrm{pad}}=f_{\mathrm{s}}$';
for ii=1:N_N_avg
    x_Title_Latex(ii)="$r_{xx}^{\mathrm{pad}}(\tau)$ (@ $\Delta t^{\mathrm{pad}}=\Delta t$), $N^{\mathrm{pad}}="+N_pad+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
    X_Title_Latex(ii)="$R_{XX}^{\mathrm{pad}}(f)$ (@$\Delta f^{\mathrm{pad}}=N/N^{\mathrm{pad}} \times \Delta f$), $N^{\mathrm{pad}}="+N_pad+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
end
signal_variance_row=r_nn_padded_cols(1,:);
y_HorizontalLines=signal_variance_row;
horizontalLinesTextLatex="$\sigma_{x}^{2}$";
x_VeticalLines=T/T_pad;
verticalLinesTextLatex="$T/T^{\mathrm{pad}}$";
Y_HorizontalLines=signal_variance_row*f_s_pad/(f_c/2)/2;
HorizontalLinesTextLatex="${\displaystyle \frac{\sigma_{x}^{2}}{2 f_{\mathrm{c}} /f_{\mathrm{s}}^{\mathrm{pad}}}}$";
X_VeticalLines=[(f_c/2)/f_s_pad;(f_s/2)/f_s_pad];
VerticalLinesTextLatex=["${\displaystyle f_{\mathrm{c}} / f_{\mathrm{s}}^{\mathrm{pad}}}$";"${\displaystyle (f_{\mathrm{s}}^{\mathrm{pad}} / 2) /f_{\mathrm{s}}^{\mathrm{pad}}= (f_{\mathrm{s}} / 2) / f_{\mathrm{s}}}$"];
plotFourierTransformPair(t_pad/T_pad,r_nn_padded_cols,f_pad/f_s_pad,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,y_HorizontalLines,horizontalLinesTextLatex,x_VeticalLines,verticalLinesTextLatex,Y_HorizontalLines,HorizontalLinesTextLatex,X_VeticalLines,VerticalLinesTextLatex,true);

f_autoCorr2=figure;maximizeFigure(f_autoCorr2);
x_xlabel_Latex='$\tau/T$  ,:$T^{\mathrm{pad}}=N^{\mathrm{pad}} / N \times T$';
X_xlabel_Latex='$f/f_{\mathrm{s}}$  ,:$f_{\mathrm{s}}^{\mathrm{pad}}=f_{\mathrm{s}}$';
for ii=1:N_N_avg
    x_Title_Latex(ii)="$r_{xx}^{\mathrm{pad}}(\tau)$ (@ $\Delta t^{\mathrm{pad}}=\Delta t$), $N="+N+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
    X_Title_Latex(ii)="$R_{XX}^{\mathrm{pad}}(f)$ (@$\Delta f^{\mathrm{pad}}=N/N^{\mathrm{pad}} \times \Delta f$),$N="+N+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
end
x_VeticalLines=[];
verticalLinesTextLatex=[];
X_VeticalLines=X_VeticalLines(1);
VerticalLinesTextLatex="${\displaystyle f_{\mathrm{c}} / f_{\mathrm{s}}}$";
plotFourierTransformPair(t/T,r_nn_padded_cols,f_pad(1:N)/f_s,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,y_HorizontalLines,horizontalLinesTextLatex,x_VeticalLines,verticalLinesTextLatex,Y_HorizontalLines,HorizontalLinesTextLatex,X_VeticalLines,VerticalLinesTextLatex,true);

%export_figure([f_spectrum,f_autoCorr1,f_autoCorr2],'==',["WhiteNoiseVsN_Avg","WhiteNoiseVsAvg-details","WhiteNoiseVsAvg"])
