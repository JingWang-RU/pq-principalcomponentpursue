function generate_batch_gaussian_image( file_fold, sigma, p, num_img )
% function: add noise to a batch of image pictures and save into matlab format
% Input:
% file_fold: images fold
% sigma: noise parameter
% p: the percent of the occluded pixels
% num_img: number of images in the file_fold

% example:
% sigma = 0.2;
% p = 0.5;
% num_img = 50;
% file_fold = ['.\Data\picture\'];
% generate_batch_gaussian_image( file_fold, sigma, p, num_img );
if ~exist(file_fold, 'dir')
    error('no such fold \n');
end
org_fold = [file_fold 'orgImage\'];
noise_fold = [file_fold  'gaussImage\'];
if ~exist(org_fold, 'dir')
    mkdir(org_fold);
    mkdir(noise_fold);
end
for i = 1: num_img
    inputFile = [org_fold  num2str(i) '.jpg'];
    inputImage = imread(inputFile);
    [ I, G, DG, r ] = gaussian_image( inputImage, sigma, p);
    save([noise_fold 'pic_s_' num2str(sigma) '_p_' num2str(i) '.mat'],'I','G','DG','r');
end
