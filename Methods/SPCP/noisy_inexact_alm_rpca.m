function [A_hat E_hat iter value_F rankA] = noisy_inexact_alm_rpca(D, sigma, lambda, tol, maxIter)

% Apr 2010
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust PCA with noise.
%
% min \|A\|_* + lambda * \|E\|_1   subj to  D = A + E + G, \|G\|_F <= sigma  
%
%   
% This matlab code is based on the following implementation of the inexact
% augmented Lagrange multiplier method for Robust PCA.
% 
%   "The Augmented Lagrange Multiplier Method for Exact Recovery of Corrupted
%   Low-Rank Matrices", Z. Lin, M. Chen, L. Wu, and Y. Ma
% 
% D - m x n matrix of observations/data (required input)
%
% sigma - Frobenius norm of noise (required input) 
%
% lambda - weight on sparse error term in the cost function
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
%
% Zihan Zhou (zzhou7@illinois.edu)

addpath PROPACK;

[m n] = size(D);

if nargin < 3
    lambda = 1 / sqrt(m);
end

if nargin < 4
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end
% 

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;
% 
% A_hat = zeros( m, n);
% E_hat = zeros( m, n);
% G_hat = zeros( m, n);
A_hat = randn( m, n);
E_hat = randn( m, n);
G_hat = randn( m, n);
mu = 1.25/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.5;          % this one can be tuned
d_norm = norm(D, 'fro');

%mu = 0.1;

value_F = 0;
iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = round(min(m,n)/10);
while ~converged       
    iter = iter + 1;
%     
%     old_A = A_hat;
%     old_E = E_hat;
    
    temp_T = D - A_hat - E_hat + (1/mu)*Y;
    scale = min(1,sigma / norm(temp_T,'fro'));
    G_hat = temp_T * scale;
    
    temp_T = D - A_hat - G_hat + (1/mu)*Y;
    E_hat = max(temp_T - lambda/mu, 0);
    E_hat = E_hat+min(temp_T + lambda/mu, 0);

    if choosvd(n, sv) == 1
        [U S V] = lansvd(D - E_hat - G_hat + (1/mu)*Y, sv, 'L');
    else
        [U S V] = svd(D - E_hat - G_hat + (1/mu)*Y, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    svp_up = diagS(1:svp) - 1/mu;
    A_hat = U(:, 1:svp) * diag(svp_up) * V(:, 1:svp)';    

    total_svd = total_svd + 1;
    
    Z = D - A_hat - E_hat - G_hat;
    
    Y = Y + mu*Z;
    
    mu = min(mu*rho, mu_bar);
%     rankA(iter) = rank(A_hat);
    %+ 1/(2*norm(A_hat + E_hat -D, 'fro')^2)
%     value_F(iter) = mu_obj*(sum(svp_up)+lambda*sum(abs(E_hat(:)))+(1/(2*mu_obj))*norm(A_hat + E_hat -D, 'fro')^2);    
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
    
%     A_obj = A_hat + epsV; %low rank
%     E_obj = abs(E_hat) + epsV; %sparse
    %norm(x,'fro') is sqrt(x_i ^2) Value_F(iter)
   
    
%     if mod( total_svd, 1) == 0
%         disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
%             ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
%             ' stopCriterion ' num2str(stopCriterion)]);
%         disp(['A_real: ' num2str(norm(A_hat-A,'fro') / norm(A,'fro')) ...
%             ' E_real: ' num2str(norm(E_hat-E,'fro') / norm(E,'fro'))]);
%         disp(['A_diff: ' num2str(norm(A_hat-old_A,'fro') / norm(A_hat,'fro')) ...
%             ' E_diff: ' num2str(norm(E_hat-old_E,'fro') / norm(E_hat,'fro'))]);
%         disp(['Objective: ' num2str(sum(diagS(1:svp) - 1/mu) + lambda * norm(E_hat(:),1))]);
%     end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end
%  r = min(m,n);
%  if size(svp_up) > r
%        svp_up = svp_up(1:r);
%  end
 