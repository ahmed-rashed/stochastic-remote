clc
close all
clearvars

%Noise Parameters
SNR_vec=[
%             inf     %Deterministic signal
%             10      %Deterministic signal buried in noise
            1
%             0.5
%             0
            ].';   %White noise

f_max_by_f_Nyq=0.8;   %BAndwidth of the band-limited white noise

bestCaseLeakage=false;

i_win=1;    %windowName_vec=["rectangular","Hann","Hamming","Kaiser-Bessel"];

linearCorrelation=true;

signal_fn=@(t) (sawtooth(t));

for SNR=SNR_vec
    demonstrateCorrelation(signal_fn,SNR,f_max_by_f_Nyq,bestCaseLeakage,linearCorrelation)
    demonstrateSpectralDensity(signal_fn,SNR,i_win,linearCorrelation,f_max_by_f_Nyq,bestCaseLeakage)
end