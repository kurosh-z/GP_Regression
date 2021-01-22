% Gaussian process Regression

function [mpost, vpost] = GPRegression(data,kernel,  hyper)
xTrain= data.xTrain';
yTrain= data.yTrain;
xTest = data.xTest';
yTest = data.yTest;


K = kernel(xTrain, hyper);
Kx= kernel(xTrain, hyper, xTest);
Kxx = kernel(xTest, hyper);

N = size(yTrain, 1);
Nt = size(yTest, 1);
mu = mean(yTrain)*ones(N,1);
mux = mean(yTest)*ones(Nt,1);

% calculating the inverse of Kx^T(K + siqma^2I) is computationally expensive,
% so what we do instead is to use cholansky decomposition
siqma = .2;
A = K + siqma^2*eye(N);
% G = Kx* inv(A)
G = chol_solve(A', Kx');

% posterior mean and variance:
mpost =mux+ G * (yTrain -  mu);
vpost = Kxx - G * Kx';

end


function [x] = chol_solve(A, b)
[cho_fac, flag] = chol(A);

x = cho_fac\(cho_fac'\b);
x=x';


end

