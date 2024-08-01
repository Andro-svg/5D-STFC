function Y = myWhite(X,DIM)
%MYWHITE - The whitening function.
%   To calculate the white matex of input matrix X and 
%   the result after X being whitened. 
%   
%   Res = myWhite(X,DIM)
% 
%   Input - 
%   X: a N*M matrix containing M datas with N dimensions;
%   DIM: specifies a dimension DIM to arrange X.
%       DIM = 1: X(N*M)
%       DIM = 2: X(M*N)
%       DIM = otherwisw: error
%   Output - 
%   Y  : result matrix of X after being whitened;
%       Y.PCAW: PCA whiten result;
%       Y.ZCAW: ZCA whiten result.
% 
%   Copyright (c) 2018 CHEN Tianyang
%   more info contact: tychen@whu.edu.cn


% ZCA白化（Zero-phase Component Analysis Whitening）是一种常用的数据预处理方法，用于降低数据特征之间的相关性，并将数据转换为具有单位方差和零均值的形式。
% 
% ZCA白化的主要步骤如下：
% 
% 1. 计算协方差矩阵：对数据进行中心化，即减去数据的均值，然后计算数据的协方差矩阵。
% 
% 2. 特征值分解：对协方差矩阵进行特征值分解，得到特征值和特征向量。
% 
% 3. 白化矩阵计算：将特征值进行平方根操作，然后取其倒数，得到白化矩阵。
% 
% 4. 白化操作：将数据乘以白化矩阵，即可完成ZCA白化操作。
% 
% ZCA白化的作用是降低数据特征之间的相关性，使得数据的特征更加独立，同时保持数据的原始分布形态。这有助于提高机器学习算法的性能，尤其在特征之间存在强相关性的情况下。
% 
% ZCA白化可以应用于图像处理、语音识别、信号处理等领域。在图像处理中，ZCA白化可以帮助去除图像中的冗余信息，减少噪声干扰，提高图像识别和分类的准确性。
% 
% 需要注意的是，ZCA白化可能会导致数据的尺度变化，因此在应用ZCA白化之前，通常需要对数据进行标准化或归一化处理，以保持数据的尺度一致性。

% 白化可以分为PCA白化和ZCA白化两种。PCA白化是通过PCA获得新的特征空间中的坐标，
% 并对这些坐标进行标准差归一化处理。ZCA白化则是在PCA的基础上进一步处理，将数据变换回原来的坐标系下的坐标。

%% parameter test
if nargin < 2
    DIM = 1;
end
if DIM == 2
    X = X';
elseif DIM~=1 && DIM~=2
    error('Error! Parameter DIM should be either 1 or 2.');
end
[~,M] = size(X);

%% whitening
% step1 PCA pre-processing
% X = X - repmat(mean(X,2),1,M);        % de-mean  
C = 1/M*X*(X');                       % calculate cov(X), or: C = cov((X)')
[eigrnvector,eigenvalue] = eig(C);    % calculate eigenvalue, eigrnvector
% TEST NOW: eigrnvector*(eigrnvector)' should be identity matrix.
% step2 PCA whitening
if all(diag(eigenvalue))    % no zero eigenvalue
    Xpcaw = eigenvalue^(-1/2) * (eigrnvector)' * X;
else
    vari = 1./sqrt(diag(eigenvalue)+1e-6);
    Xpcaw = diag(vari) * (eigrnvector)' * X;
end
% Xpczw = (eigenvalue+diag(ones(size(X,1),1)*(1e-5)))^(-1/2)*(eigrnvector)'*X;    % 数据正则化
% step3 ZCA whitening
Xzcaw = eigrnvector*Xpcaw;
% TEST NOW: cov((Xpczw)') and cov((Xzcaw)') should be identity matrix.

%% result output
Y.PCAW = Xpcaw;
Y.ZCAW = Xzcaw; %(mat2gray(Xzcaw)).^1.5;%.*;

end