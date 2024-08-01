function [patchTen, patchNumber, patchPosition] = construct_patch_ten(img, patchSize, slideStep)
if ~exist('patchSize', 'var')
    patchSize = 50;
end

if ~exist('slideStep', 'var')
    slideStep = 10;
end

% img = reshape(1:9, [3 3])
% img = reshape(1:12, [3 4])
% patchSize = 2;
% slideStep = 1;
[imgHei, imgWid] = size(img);

rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;
colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;
rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];

%% arrayfun version, identical to the following for-loop version
% [meshCols, meshRows] = meshgrid(colPosArr, rowPosArr);
% idx_fun = @(row,col) img(row : row + patchSize - 1, col : col + patchSize - 1);
% patchCell = arrayfun(idx_fun, meshRows, meshCols, 'UniformOutput', false);
% patchTen = cat(3, patchCell{:});

%% for-loop version
% 常规的排队法：
%  1  2 
%  3  4  5
%  6  7  8  9  10
patchTen = zeros(patchSize/3, patchSize/3, 9, rowPatchNum * colPatchNum);
patchPosition = zeros(1,2,9,rowPatchNum*colPatchNum);
k = 0;
for col = colPosArr
    for row = rowPosArr
        k = k + 1;
        tmp_patch = img(row : row + patchSize - 1, col : col + patchSize - 1); % tmp_patch 来构造一个small_patch * small_patch * 9 的小张量
        l=0;
        for subcol = 1:patchSize/3:patchSize
            for subrow = 1:patchSize/3:patchSize
                l = l + 1;
                small_tmp_patch = tmp_patch( subrow  : subrow + patchSize/3 - 1, subcol : subcol + patchSize/3 - 1); % 构造small_patch * small_patch * 9 的小张量
                smallpatchTen(:, :, l) = small_tmp_patch;
                subpatchPosition(:,:,l) = [row+subrow-1 , col+subcol-1];
            end
        end
        patchTen(:, :, :, k) = smallpatchTen;
        patchPosition(:,:,:,k) = subpatchPosition;
        smallpatchTen = zeros(size(smallpatchTen));
    end
end
patchNumber=k;

