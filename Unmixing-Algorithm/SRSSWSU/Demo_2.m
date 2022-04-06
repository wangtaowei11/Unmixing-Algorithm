%%
% Simulation datasets generation part code 
% cite by : A Fast Multiscale Spatial Regularization
% for Sparse Hyperspectral Unmixing

close all
clear all
clc
p = 9;  
SNR = 30; %10;20;30;40
bandwidth = 10000; % 10000 == iid noise
load spatial2.mat

nl = size(Xim,1);
nc = size(Xim,2);
np = nl*nc;     % number of pixels
% buid the dictionary 
load USGS_1995_Library.mat
%  order bands by increasing wavelength
[dummy index] = sort(datalib(:,1));
A =  datalib(index,4:end);
names = names(4:end,:);
min_angle = 4.44;       
[A, index] = prune_library2(A,min_angle); % 240  signature 
names = names(index',:);

% order  the columns of A by decreasing angles 
[A, index, angles] = sort_library_by_angle(A);
names = names(index',:);
namesStr = char(names);
supp = [2 3 4 5 6 7 8 9 10]; % dont take 2 Jarosites

M = A(:,supp);
[L,p] = size(M);  % L = number of bands; p = number of material

%---------------------------------
% generate  the observed  data X
%---------------------------------
% set noise standard deviation
sigma = sqrt(sum(sum((M*X).^2))/np/L/10^(SNR/10));
% generate Gaussian iid noise
noise = sigma*randn(L,np);
% make noise correlated by low pass filtering
% low pass filter (Gaussian)
filter_coef = exp(-(0:L-1).^2/2/bandwidth.^2)';
scale = sqrt(L/sum(filter_coef.^2));
filter_coef = scale*filter_coef;
noise = idct(dct(noise).*repmat(filter_coef,1,np));

%  observed spectral vector
Y = M*X + noise;

n = size(A,2);   % 
N = nl*nc;       % 
XT = zeros(n,N);
XT(supp,:) = X;  % 

%%
%parameter settings
if SNR == 20
    % 40db
    lambda1 = 0.001;
    lambda2 = 0.005;
    beta       = 2e-1;
    slic_size  = 5;
    slic_reg   = 0.0001;
    
elseif SNR == 30
    % 30db
    lambda1 = 0.003;
    lambda2 = 0.03;
    beta       = 3e-3;
    slic_size  = 6;
    slic_reg   = 0.00125;
end
%% Superpixels segmentation
Y2 = reshape(Y', nl, nc, L);   
Y2a = Y2;  
% reorder and rescale data into 2-D array
[numRows,numCols,numSpectra] = size(Y2);
norm_y = mean(reshape(sqrt(sum(Y2.^2,3)), numRows*numCols, 1));
Y2 = Y2./norm_y;

imgVec = reshape(Y2, [numRows*numCols numSpectra]);
spSegs = vl_slic(single(Y2), slic_size, slic_reg);
numSuperpixels = double(max(spSegs(:)))+1; 
%%
Y3 = zeros(size(Y2));
avg_superpx = zeros(1, numSuperpixels+1, L);

for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);%find()function
    
    for j=1:length(rowi)
        % Averages all pixels inside each superpixel
        if j == 1
            avg_superpx(1,i+1,:) = (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        else
            avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        end
    end
    
    % This is optional (for visualization)
    for j=1:length(rowi)
        
        Y3(rowi(j),coli(j),:) = avg_superpx(1,i+1,:);
    end
end
%%
%Spectral angular distance 
Y4 = reshape(Y3, [numRows*numCols numSpectra])';
Y5 = Y4./norm_y;
Y6 = Y./norm_y;
x_1=diag(Y6'*Y5)';
W_11=abs(acos(x_1));
beta = 10*W_11.*W_11;
beta = repmat(beta,240,1);

%% unmixing

[X_C,res_d,res_p]= SRSSWSU_1(squeeze(avg_superpx)',A,lambda1,0.5);
temp = zeros(size(Y2,1), size(Y2,2), n);
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = X_C(:,i+1);
    end
end
X_hat_C = reshape(temp, [size(Y2,1)*size(Y2,2) n])';

%%
[X_hat_SRSSWSU]= SRSSWSU_2(Y,A,lambda1,beta,0.5,X_hat_C);
%SREvalue
SRE_SRSSWSU = 20*log10(norm(XT,'fro')/norm(X_hat_SRSSWSU-XT,'fro'));

%%