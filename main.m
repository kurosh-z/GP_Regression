close all, clear all,
clc

% gaussian kernal hyperparameters
hyper.siqma =log(3);
hyper.ell = log(1.3);



% periodic kernal hyperparameters
% hyper.l = log(.5);
% hyper.siqma = log(1.2);
% hyper.p = log(.5);

% linear kernel hyperparameters
hyper.siqma_a = log(.2);
hyper.siqma_b = log(.1);
hyper.c = log(2.5);
data = load('tempData.mat')

disp("-----------------------------------------------------")
disp("Calculating the GP with linear kernel and hyperparametrs:")
sprintf('sigma_a: %f, sigma_b: %f and c: %f \n', exp(hyper.siqma_a), exp(hyper.siqma_b), exp(hyper.c))
[mpost_lin, vpost_lin]=GPRegression(data, @lin_kernel, hyper);
disp("Linear Kernel GP finished the posterior mean and variance are saved in variables mpost_lin and vpost_lin ")
fig1= plot_gp_2d(data.xTest, data.yTest ,'Test Data', true);
fig2=plot_gp_2d(data.xTest, mpost_lin ,'Predictions GP with Linear Kernel', true);

disp("-----------------------------------------------------")
disp("Calculating the GP with gaussian kernel and hyperparametrs:")
sprintf('sigma: %f and l: %f \n', exp(hyper.siqma), exp(hyper.ell))
[mpost_se, vpost_se]=GPRegression(data, @se_kernel, hyper);
plot_gp_2d(data.xTest, mpost_se ,'Predictions GP with se-kernel', true);


%% optimizing the the hyperparameters for se_kernel
%% Notice: to make it faster to run, I just used 3500 measurements to
% optimized the parameters , if you want to test it with full measurements 
%simply set num_meas to N !
rng(1);
N = size(data.xTrain, 2);
num_meas = 3500;
randIndices = randperm(N, num_meas);

data2.xTrain = data.xTrain(:, randIndices);
data2.yTrain = data.yTrain(randIndices);
data2.xTest = data.xTest(:, 1:5);
data2.yTest = data.yTest(1:5);

P.objective = @(hyper) loglike(hyper, data2.xTrain', data2.yTrain);
P.options = optimoptions('fmincon','SpecifyObjectiveGradient',true, 'Display','iter', 'Algorithm','sqp', 'TolFun',1e-5,'TolX',1e-3);
P.x0 = [hyper.siqma, hyper.ell];
P.A = [];
P.b = [];
P.Aeq = [];
P.beq = [];
P.lb = [];
P.ub = [];
P.nonlcon = [];
P.solver= 'fmincon';
tic;
[optTheta, objVal,exitFlag,output] = fmincon(P);

sprintf('optimization finished with parametrs theta=[siqma, ell] = [%f, %f]:', exp(optTheta(1)), exp(optTheta(2)))
optHyper.siqma = optTheta(1);
optHyper.ell = optTheta(2);

disp("-----------------------------------------------------")
disp("Calculating the GP with gaussian kernel and optimized hyperparameter:")
sprintf('sigma: %f and l: %f \n', exp(optHyper.siqma), exp(optHyper.ell))
[mpost_seOpt, vpost_seOpt]=GPRegression(data, @se_kernel, optHyper);
plot_gp_2d(data.xTest, mpost_seOpt ,'Predictions with optimized hyperparameter', true);

disp("-------------------------------------------------------------")
dOpt=mpost_seOpt-data.yTest;
errOpt = sum(dOpt.^2);
d=mpost_se-data.yTest;
err = sum(d.^2);
sprintf('squared error of predictions for optimized parameter is %f \n in comparision with none-optimized parameters: %f', errOpt, err)






%% Sub-Routines
% computes the log likelihood of the generative model on the training data,
%         as a function of the hyperparameters, with derivative.
%         Input:
%         hypers — log hyperparameters, as defined for the kernel

function [loglik, dloglik]=loglike(hyper, X, y)

se_hyper.siqma = hyper(1);
se_hyper.ell = hyper(2);

% lin_hyper.siqma_a = hyper(1);
% lin_hyper.siqma_b = hyper(2);
% lin_hyper.c = hyper(3);
[K, dK] = se_kernel(X, se_hyper);
% [K, dK] = lin_kernel(X, lin_hyper);

% adding noise to the kernel
siqma = .2;
N = size(X,1);
G = K + siqma^2*eye(N);
ld = logdet(G);
%a = G\y;
a = linsolve(G,y);%a simple trick to solve G\y faster!
loglik = dot(a, y) +ld; % (Y / G) * Y + log |G|


num_hyper = length(hyper);
dloglik = zeros(num_hyper,1);
for i= 1: num_hyper
    b= dK(N*(i-1)+1:(i)*N,:); % slicing dk to get dk for i-th parameter
    %dloglik(i)=-dot(a, b*a) + trace(G\b);
    dloglik(i)=dot(-a, b*a) + trace(linsolve(G, b)); %solving linear eq instead of G\b makes it faster to run!
end

end

%log determinant of a positive definite matrix:
function lgdet = logdet(A)
% L = chol(M);
% lgdet = 2*sum(log(diag(L)));

[L, U, P] = lu(A);
du = diag(U);
c = det(P) * prod(sign(du));
lgdet = log(c) + sum(log(abs(du)));
end