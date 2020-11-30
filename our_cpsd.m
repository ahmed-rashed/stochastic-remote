function R_XY_vec=our_cpsd(x_tot_vec,y_tot_vec ,win_vec ,K_o)

if any(size(x_tot_vec)~=size(y_tot_vec)),error('x_tot_vec and y_tot_vec must have the same size'),end
if mod(K_o,1)~=0,error('N_o must be integer'),end

K_tot=length(x_tot_vec);
K=length(win_vec);
if K>K_tot,error('win_vec length must be smaller than or equals to x_tot_vec length'),end

if K_o>(K+1),error('N_o cannot exceed win_vec length + 1'),end

ind_1=find(size(x_tot_vec)==1,1);
if isempty(ind_1),error('x_tot_vec must be a row or column vector');end
if size(win_vec,ind_1)~=1,error('Both win_vec and x_tot_vec must be rows or columns'),end

P=floor((K_tot-K_o+1)/(K-K_o+1));
R_XY_vec=zeros(size(win_vec));
for p=1:P
    k_start=(p-1)*(K-K_o)+1;
    p_sig_ind=(0:K-1)+k_start;
    R_XY_vec=R_XY_vec+ourPeriodogram(x_tot_vec(p_sig_ind),y_tot_vec(p_sig_ind),win_vec);
end
R_XY_vec=R_XY_vec/P;