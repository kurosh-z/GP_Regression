function [K] = per_kernel(X, hyp, XX)
if nargin == 2
    XX = X;
end

l = exp(hyp.l);
p = exp(hyp.p);
siqma= exp(hyp.siqma);
D = pdist2(XX, X);
K = siqma^2 * exp(-sin((pi/p)*D).^2/(l^2));

disp('finished')
end