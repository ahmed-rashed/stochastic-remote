clc
close all
clearvars

N=2^9
N_Nyq=floor(N/2);
T=5
[dt,f_s,df]=samplingParameters_T_N(T,N);
f0=10*df;

SNR=0.5;  %SNR_db=-3
f_c_by_f_Nyq=.8;

N_avg_vec=[1,1000,10000];
%%%%%%%%%%%%%%%%%%%%%%%
t_col=(0:N-1).'*dt;
f_row=(0:N-1)*df;

N_pad=2^nextpow2(2*N-1)
T_pad=N_pad/N*T;
[dt_pad,f_s_pad,df_pad]=samplingParameters_T_N(T_pad,N_pad);
t_pad=(0:N_pad-1).'*dt_pad;
f_pad=(0:N_pad-1).'*df_pad;

N_N_avg=length(N_avg_vec);
signal_expectation_current=zeros(N,1);
signal_expectation_cols=zeros(N,N_N_avg);
r_ss_padded_current=zeros(N_pad,1);
r_ss_padded_cols=zeros(N_pad,N_N_avg);
nn=1;
rng(10);
for N_avg=1:N_avg_vec(end)
    x_col_i=sin(2*pi*f0*t_col+randn*T);   %Pure sine but with random start phase
    [x_hat_col_i,n_col_i]=addNoise(x_col_i,SNR,f_c_by_f_Nyq);
    
    r_nn_i_padded=xcorr_xspec_linear(x_hat_col_i,x_hat_col_i,'unbiased');

    signal_expectation_current=(1-1/N_avg)*signal_expectation_current+x_hat_col_i/N_avg;  %equivalent to sum(x_hat_col_i)/N_avg,for i =1 to i=N_avg
    
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
X_xlabel_Latex='$f/f_{\mathrm{s}}$';
x_Title_Latex=strings(N_N_avg,1);
X_Title_Latex=strings(N_N_avg,1);
for ii=1:N_N_avg
    x_Title_Latex(ii)="$E[x(t)]$,: $N="+N+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
    X_Title_Latex(ii)="$|E[X(f)]|$,: $N="+N+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
end
f_c=f_c_by_f_Nyq*(f_s/2);
X_VeticalLines=(f_c/2)/f_s_pad;
VerticalLinesTextLatex="${\displaystyle f_{\mathrm{c}} / f_{\mathrm{s}}}$";
plotFourierTransformPair(t_col/T,signal_expectation_cols,f_row(1:N_Nyq)/f_s,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,[],[],[],[],[],[],X_VeticalLines,VerticalLinesTextLatex,true);

f_autoCorr1=figure;maximizeFigure(f_autoCorr1);
x_xlabel_Latex='$\tau/T^{\mathrm{pad}}$ ,:$T^{\mathrm{pad}}=N^{\mathrm{pad}} / N \times T$';
X_xlabel_Latex='$f/f_{\mathrm{s}}^{\mathrm{pad}}$ ,:$f_{\mathrm{s}}^{\mathrm{pad}}=f_{\mathrm{s}}$';
for ii=1:N_N_avg
    x_Title_Latex(ii)="$r_{xx}^{\mathrm{pad}}(\tau)$ (@ $\Delta t^{\mathrm{pad}}=\Delta t$), $N^{\mathrm{pad}}="+N_pad+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
    X_Title_Latex(ii)="$R_{XX}^{\mathrm{pad}}(f)$ (@$\Delta f^{\mathrm{pad}}=N/N^{\mathrm{pad}} \times \Delta f$), $N^{\mathrm{pad}}="+N_pad+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
end
%signal_variance_row=r_ss_padded_cols(1,:);
y_HorizontalLines=[];
horizontalLinesTextLatex=[];
x_VeticalLines=T/T_pad;
verticalLinesTextLatex="$T/T^{\mathrm{pad}}$";
Y_HorizontalLines=[];
HorizontalLinesTextLatex=[];
X_VeticalLines=[(f_c/2)/f_s_pad;(f_s/2)/f_s_pad];
VerticalLinesTextLatex=["${\displaystyle f_{\mathrm{c}} / f_{\mathrm{s}}^{\mathrm{pad}}}$";"${\displaystyle (f_{\mathrm{s}}^{\mathrm{pad}} / 2) /f_{\mathrm{s}}^{\mathrm{pad}}= (f_{\mathrm{s}} / 2) / f_{\mathrm{s}}}$"];
plotFourierTransformPair(t_pad/T_pad,r_ss_padded_cols,f_pad/f_s_pad,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,y_HorizontalLines,horizontalLinesTextLatex,x_VeticalLines,verticalLinesTextLatex,Y_HorizontalLines,HorizontalLinesTextLatex,X_VeticalLines,VerticalLinesTextLatex,true);

f_autoCorr2=figure;maximizeFigure(f_autoCorr2);
x_xlabel_Latex='$\tau/T$  ,:$T^{\mathrm{pad}}=N^{\mathrm{pad}} / N \times T$';
X_xlabel_Latex='$f/f_{\mathrm{s}}$  ,:$f_{\mathrm{s}}^{\mathrm{pad}}=f_{\mathrm{s}}$';
for ii=1:N_N_avg
    x_Title_Latex(ii)="$r_{xx}^{\mathrm{pad}}(\tau)$ (@ $\Delta t^{\mathrm{pad}}=\Delta t$), $N="+N+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
    X_Title_Latex(ii)="$R_{XX}^{\mathrm{pad}}(f)$ (@$\Delta f^{\mathrm{pad}}=N/N^{\mathrm{pad}} \times \Delta f$), $N="+N+'$ \& $N_{\mathrm{avg}}='+N_avg_vec(ii)+'$';
end
x_VeticalLines=[];
verticalLinesTextLatex=[];
X_VeticalLines=X_VeticalLines(1);
VerticalLinesTextLatex="${\displaystyle f_{\mathrm{c}} / f_{\mathrm{s}}}$";
plotFourierTransformPair(t_col/T,r_ss_padded_cols,f_pad(1:N)/f_s,x_xlabel_Latex,x_Title_Latex,X_xlabel_Latex,X_Title_Latex,y_HorizontalLines,horizontalLinesTextLatex,x_VeticalLines,verticalLinesTextLatex,Y_HorizontalLines,HorizontalLinesTextLatex,X_VeticalLines,VerticalLinesTextLatex,true);

%export_figure([f_spectrum,f_autoCorr1,f_autoCorr2],'==',["SpectrumAvg","CorrAvg-details","CorrAvg"])
