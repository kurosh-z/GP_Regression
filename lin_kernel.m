
%% linear kernal of form:
%  siqma_b^2 + siqma_a^2(x-c)(x'-c)
%
function [K] = lin_kernel(X, hyp, XX)
if nargin == 2
    XX = X;
end

c = exp(hyp.c);
siqma_b = exp(hyp.siqma_b);
siqma_a= exp(hyp.siqma_a);

n = size(X,1);
m = size(XX,1);

K = zeros(m,n);
d = zeros(m,n);
for i=1:m
    for j=1:n
        K(i,j) = siqma_b^2 + siqma_a^2* dot((X(j,:)-c), (XX(i,:)-c));
%         d(i,j) = -norm((XX(i,:)-c))-norm((XX(j,:)-c));
    end
    
end

% if nargout == 2
%     dKhyp = [2*siqma_b^2*ones(n); 2*(K-siqma_b^2);siqma_a^2*c*d ];
% end


end


