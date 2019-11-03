function M_coh=MultipleCoherehence_Ahmed(fs,x,y,n,nfft)
% Calculation of multiple coherence.
%
%This is based on equation (5.126) of [Random Data: Analysis and measurements procedures, Bendat and Piersol]

% x :inputs (#time samples x #inputs)
% y :outputs (#time samples x #outputs)
% index : selected output 
% nfft : number of fft
% fs : sampling rate
% window : you can change window function if your data are not periodic
%
%
% Sample usage:
%   mCoh=MultipleCoherehence(2048,inputs,Datas,5,2048);
%
%
%
%  For other Modal Analysis Routines please contact :
%  mmbicak@mtu.edu
%  http://www.me.mtu.edu/~mmbicak
noverlap=[];
M=length(x(1,:));% # inputs
N=length(y(1,:));% # outputs
if (n<=0) || (n>N) 
    error('index error, check index value (index)') 
end
window =[];% assuming pure harmonic signal,otherwise use a window
%window=hann(nfft);
NN=floor(nfft/2)+1;
R_XX=nan(M,M,NN);
for m1=1:M 
    for m2=1:M
      R_XX(m1,m2,:)=cpsd(x(:,m1),x(:,m2),window, noverlap,nfft,fs);
    end
end

%R_XY=cell(M,1);
R_YX=nan(M,1,NN);
for m1=1:M 
%     R_XY{m1,1}=cpsd(y(:,n),x(:,m1),window, noverlap,nfft,fs);
     R_YX(m1,1,:)=cpsd(x(:,m1),y(:,n),window, noverlap,nfft,fs);   %By Le Fristoker
end

R_YY=cpsd(y(:,n),y(:,n),window, noverlap,nfft,fs);
%GY=[R_YY;R_YX]; %Le Fristoker claims this is a bug
%GY=[R_YY;R_YX]; %Le Fristoker claims this is a solution to the bug
%R_YXX=[[R_YY;R_YX] [R_YX.';R_XX]];

M_coh=nan(1,NN);
for ii= 1:NN
    R_YXX_f=[[R_YY(ii);R_YX(:,:,ii)] [R_YX(:,:,ii).';R_XX(:,:,ii)]];
    M_coh(ii)=1-det(R_YXX_f)/R_YY(ii)/det(R_XX(:,:,ii));
end
