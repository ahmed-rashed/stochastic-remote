function s=signalPulse(t_mat,T_burst,s_fn)

if ~isa(s_fn,'function_handle')
    error('s_fn must be a function handle');
end

s=zeros(size(t_mat));
i_vec=find(t_mat>=0 & t_mat<=T_burst);
s(i_vec)=s_fn(t_mat(i_vec));