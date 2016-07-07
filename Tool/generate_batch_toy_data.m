function generate_batch_toy_data( sigma, c, r, num_iter)
% GENERATE_BATCH_TOY_DATA
for i = 1:num_iter
    U = c * randn(m,r) ; V = c * randn(n,r) ;
    A = U*V' ;
    save_fold = '../Data/';
    if ~exist(save_fold, 'dir')
        mkdir(save_fold);
    end
    temp = randperm(m*n) ;
    numCorruptedEntries = round(ps*m*n) ;
    corruptedPositions = temp(1:numCorruptedEntries) ;
    E = zeros(m,n) ;
    E(corruptedPositions) = 10*(rand(numCorruptedEntries,1)-0.5) ;
    
    filename = ['data_s' num2str(sigma) '_d_' num2str(m) '_ps_' num2str(ps) '_r_' num2str(r) '_itr_' num2str(i) ];
    save([save_fold filename '.mat'],'A','E','G');
end

