function [ I, G, DG, r ] = gaussian_image( inputImage, sigma,p )
%GAUSSIAN_IMAGE Summary of this function goes here
% Input:
% inputImage: image, three channel
% sigma: noise
% p: add noise to p % of the image maxtrix

% Output:
% I : the clean data, three channel;
% G : the gaussian noise;
% DG : the occluded image.
% r : the rank of the image.

I = double(inputImage)/256;
r = rank(I(:,:,1));

[m,n] = size(I(:,:,1));
noise = sigma*randn(1,floor(p*m*n));

temp = randperm(m*n) ;
numCorruptedEntries = floor(p*m*n) ;%round(0.5*m*n) ;
corruptedPositions = temp(1:numCorruptedEntries) ;
G = zeros(m,n);
G(corruptedPositions) = noise;
for i = 1:3
    Tp = I(:,:,i);
    
    Tp = Tp + G;
    D{i} = Tp;
    DG(:,:,i) = D{i};
end

end

