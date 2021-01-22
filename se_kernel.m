% Standard Squared Exponential (SE) kernel
% Input
%     X: training points
%     hyp: hyperparameters
%     XX: testing points (optional)
% Output
%     K: dense kernel matrix
%     dKhyp: kernel matrix derivatives w.r.t. hyperparameters

function [K, dKhyp] = se_kernel(X, hyp, XX)
if nargin == 2
    XX = X;
end

ell = exp(hyp.ell); 
siqma= exp(hyp.siqma);
D = pdist2(XX, X);
K = siqma^2 * exp(-D.^2/(2*ell^2));
% 
% sprintf('size of D: %d * %d \n', size(D))
% sprintf('size of K: %d * %d \n', size(K))
if nargout == 2
    dKhyp = [ 2*K; 1/ell^2 * (D.^2 .* K)];
end
end