function [A E iter obj rankA] = pqpca(D, lambda_1, lambda_2, p, q, A_init, E_init)

% function: lpq-norm stable pca
% Created by Jing Wang (jw998@rutgers.edu) since 2013/9/06
% Wang, Jing and Wang, Meng and Hu, Xuegang and Yan, Shuicheng, Visual data denoising with a unified Schatten-p norm and Lq norm regularized principal component pursuit},
% Pattern Recognition, 48(10), 3135--3144, 2015.
% Note that there are many parameters to be tuned for better performance.
%         The stop condition of the iteration can be controlled in other ways
% Input
% D : m x n matrix of observations/data (required input)%
% lambda_1 : weight on low rank  error term in the cost function
% lambda_2 : weight on sparse error term in the cost function
%                   the parameter lambda_1 and lambda_2 should both be small
% p : Lp
% q : Lq
% A_init : initialized A
% E_init: initialized E
% Output
% A, E, iter (number of iterations), obj (value of objective function), rankA (rank of A)

[m n] = size(D);
if nargin < 2
    lambda_1 = 1/5;
end
if nargin < 3
    lambda_2 = 1/100;
end
if nargin < 4
    p = 0.2;
end
if nargin < 5
    q = 0.3;
end
maxIter = 400;%1000;

% initialize
% A = randn(m, n);
% E = randn(m, n);
% A = zeros( m, n);
% E = zeros( m, n);
A = A_init;
E = E_init;
% r = min(m,n);
mu = 2.1 ; %1.1  mu > 2, could be tuned
iter = 0;
%the value of the objective function
converged = false;
epsion = 0 ; %0.1 could be tuned.0.1,1e-1
rho = 1.1; %1.1
tol1 = 1e-3; % synthetic data 1e-6 1e-1
tol2 = 1e-3; % 1e-1 for video problem
w = ones(m, 1); % initialize the weight of low rank part, vector
M = ones(m, n); % initialize the weight of sparse part, matrix
par_A = lambda_1/mu; % the parameter of the low rank part
par_E = lambda_2/mu; % the parameter of the sparse part

while ~converged
    iter = iter + 1;
    if mod(iter,100) == 1
        fprintf(' %d ', iter);
    end
        
    % Calculate the sparse part
    E_pre = E ;
    M = q./((abs(E) + epsion).^(1-q));  %update the weight of sparse part
    %  M = 1./((abs(E) + epsion).^1); 
    %  the M_ij is the weight of E_ij, thus ./ is performed on each element of the matrix.    
    Temp = E - (1/mu)*(A + E - D);
    E = max(Temp - par_E*M, 0) + min(Temp + par_E*M, 0);    
     
    % Calculate the low rank part   
    A_pre = A;
    Temp = A - (1/mu)*(A + E_pre - D) ;
    [U S V] = svd(Temp,'econ');    
    diagS = diag(S);
    w = p./((diagS + epsion).^(1-p));
    % w=1./((diagS + epsion).^1);
    % svp is the index of the singular values bigger than lambda * W
    svp = find(diagS > par_A * w);    
    svp_up = diagS(svp) - par_A*w(svp);
    A = U(:, svp) * diag(svp_up) * V(:, svp)';    
   
   % caclualte the objective function
   % obj(iter) = lambda_1 * sum(A_obj.^p) + lambda_2 * sum(E_obj(:).^q) + 1/2*norm(A + E - D, 'fro');
   
    obj(iter) = lambda_1 * sum((svp_up + epsion).^p) + lambda_2 * sum((abs(E(:)) + epsion).^q) + (1/2)*(norm(A + E - D, 'fro')^2);
    
    % stop Criterion
    rankA(iter) = rank(A) ;
%      disp(['#r(A) ' num2str(rank(A))...
%             ' |E|_0 ' num2str(length(find(abs(E)>0)))
%              ]);
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;
    end
    
    % when E is all zero, then the first item is inf and > epsion
    if norm(E_pre-E,'fro')/norm(E,'fro') < tol1 && norm(A_pre - A,'fro')/norm(A, 'fro') < tol2
% if iter>=2
%     if abs(obj(iter)-obj(iter-1)) < tol1
        disp('converged');
        converged = 1;
    end     
% end
    epsion = epsion / rho ;  
end
%  figure
% subplot(1,2,1)
% plot(rankA)
% subplot(1,2,2)
% plot(obj);

