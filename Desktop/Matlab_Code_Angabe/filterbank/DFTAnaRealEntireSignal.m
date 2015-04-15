function X_out=DFTAnaRealEntireSignal(x_in,K,N,p)
Lp = length(p);
Lx = length(x_in);
dim = size(x_in);
if min(dim) ~= 1
    error('input must be a vector');
end
n_blocks = ceil(Lx/N);

x_in = x_in(:).';
x_buffer = buffer([zeros(1,N-1) x_in],Lp,Lp-N);
x_buffer = x_buffer(Lp:-1:1,1:n_blocks);

U = reshape(bsxfun(@times, x_buffer, p.'),K,ceil(Lp/K),n_blocks);

V = squeeze(fft(sum(U,2),[],1));
X_out = V(1:K/2+1,:);