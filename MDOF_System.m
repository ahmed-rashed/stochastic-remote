function [x_rows,y_rows]=MDOF_System(i_outputs,j_inputs,N,T)

[dt,f_s,df]=samplingParameters_T_N(T,N);
t_column=(0:N-1).'*dt;
f_column=(0:N-1).'*df;
n_f_max=floor(N/2);
f_max=n_f_max*df;

m=[1,1.2,1.1];
k=80*[100,150,100,120];
c=2*ones(1,N+1);
[M,C,K]=N_DOF_sys(m,c,k);

[EigVectors_Normalized, EigValues_vec]=MDOF_Eig_Visc(M, C, K);

%IRF
n_inputs=length(j_inputs);
n_outputs=length(i_outputs);
x_rows=randn(n_inputs,N);
y_rows=zeros(n_outputs,N);
ii_row=i_outputs;
for n_in=1:n_inputs
    jj_row=j_inputs(n_in)*ones(1,n_outputs);
    h_mat=MDOF_IRF_Visc(EigValues_vec, EigVectors_Normalized, t_column, ii_row, jj_row);

    for n_out=1:n_outputs
        %y_rows(n_out,:)=y_rows(n_out,:)+ifft(fft(h_mat(:,n_out)).'.*fft(x_rows(n_in,:)));
        y_rows(n_out,:)=y_rows(n_out,:)+filter(h_mat(:,n_out).',1,x_rows(n_in,:),[],2);
    end
end

x_rows=x_rows+0.1*randn(size(x_rows));
y_rows=y_rows+0.1*randn(size(y_rows));