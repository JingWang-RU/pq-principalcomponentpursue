

sigma = 0.1; %0.01, 0.2;

m = 200;
ps = 0.2;

for i = 1:10

%load data
filename = ['data_s' num2str(sigma) '_d_' num2str(m) '_ps_' num2str(ps) '_r_' num2str(r) '_itr_' num2str(i) ];
load (['.\Data\' filename '.mat'],'A','E','G');

[X,S] = nsa(D,stdev,tol,0,optimal_X,optimal_S);

save filename    ;
        
end