function estimate_H()
clc

N_avg=1000;

N=2^10;
j_inputs=[1,2];
i_outputs=[1,2,3];

T=15;
[dt,f_s,df]=samplingParameters_T_N(T,N);
t_row=(0:N-1)*dt;
f_row=(0:N-1)*df;
n_f_max=floor(N/2);
f_max=n_f_max*df;

n_inputs=length(j_inputs);
n_outputs=length(i_outputs);


H=zeros(2,3);
R_XX=zeros(n_inputs,n_inputs,n_f_max);
R_XY=zeros(n_outputs,n_inputs,n_f_max);
j_inputs=[1,2];
i_outputs=[1,2,3];
for ii=1:N_avg
    [x_rows,y_rows]=MDOF_System(i_outputs,j_inputs,N,T);
    X_rows=fft(x_rows,[],2);
    Y_rows=fft(y_rows,[],2);
    for nf=1:n_f_max
        R_XX_temp=X_rows(:,nf)*X_rows(:,nf)'/N;
        R_XY_temp=Y_rows(:,nf)*X_rows(:,nf)'/N;
        
        R_XX(:,:,nf)=(1-1/N_avg)*R_XX(:,:,nf)+R_XX_temp/N_avg;        %equivalent to sum(R_XX(:,:,nf))/N_avg,for i =1 to i=N_avg
        R_XY(:,:,nf)=(1-1/N_avg)*R_XY(:,:,nf)+R_XY_temp/N_avg;        %equivalent to sum(R_XY(:,:,nf))/N_avg,for i =1 to i=N_avg
    end
end

H=ones(n_outputs,n_inputs,n_f_max);
for nf=1:n_f_max
    H(:,:,nf)=R_XY(:,:,nf)/R_XX(:,:,nf);
end

figure(7);clf
subplot(6,1,1:2)
semilogy(f_row(1:n_f_max),reshape(abs(H(:,1,:)),[n_outputs,n_f_max]))
subplot(6,1,3)
plot(f_row(1:n_f_max),reshape(unwrap(angle(H(:,1,:))),[n_outputs,n_f_max]))

subplot(6,1,4:5)
semilogy(f_row(1:n_f_max),reshape(abs(H(:,2,:)),[n_outputs,n_f_max]))
subplot(6,1,6)
plot(f_row(1:n_f_max),reshape(unwrap(angle(H(:,2,:))),[n_outputs,n_f_max]))
