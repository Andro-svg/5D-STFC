function target_detection(readPath, savePath, tuneopts)

if isfield(tuneopts, 'temporal_step');   temporal_step = tuneopts.temporal_step;    end
if isfield(tuneopts, 'patchSize');               patchSize = 3*tuneopts.patchSize;    end
if isfield(tuneopts, 'lambdaL');                 lambdaL = tuneopts.lambdaL;        end
if isfield(tuneopts, 'mu');                      mu = tuneopts.mu;                  end
if isfield(tuneopts, 'x');                      x = tuneopts.x;                  end

slideStep = patchSize;


%% Get all image file names, please make sure that the image file order is correct by this reading way.
filesdir = dir([readPath '/*.jpg']);
if isempty( filesdir )
    filesdir = dir( [readPath '/*.bmp'] );
end
if isempty( filesdir )
    filesdir = dir([readPath '/*.png']);
end
if isempty( filesdir )
    fprintf('\n There is no any image in the folder of %s', readPath);
    return;
end
% get all image file names into a cell array;
files = { filesdir.name };
files = sort_nat(files);

iteration = 0;

%% begin to process images using mog based detection method.
sliding_step = temporal_step;
if mod(length(files),sliding_step)==0
    t_list =  1 : sliding_step : length(files)-temporal_step + 1;
else
    t_list =  [1 : sliding_step : length(files)-sliding_step + 1, length(files) - temporal_step + 1];
end

for t = t_list
    iteration = iteration + 1;
    disp('========================================');
    fprintf('%s %d%s\n','Starting', iteration, '-th loop');
    spat_temp_ten = []; % [n1,n2,n3,tensorNumber]
    X=[];

    %% read images and construct the patch image
    for tt = 1 : temporal_step
        img = imread([readPath '/' files{tt+t-1}]);
        if size(img, 3) > 1
            img = rgb2gray( img );
        end
        [imgHei, imgWid] = size(img);

        img = double(img);

        imwrite(mat2gray(img), [savePath '/' files{tt+t-1}]);
        X(:,:,tt) = img;
    end
    
    sz=size(X); 
    X = reshape(X, [sz(3)  prod(sz(1:2))]);

    %% whitening
    Y = myWhite(X,1);
    Z = Y.ZCAW;
    Z = reshape(Z, [sz(1) sz(2)  sz(3)]);
    for ttt = 1 : temporal_step
        img = Z(:,:,ttt);
        %% construct patch tensor
        [tenF, patchNumber, patchPosition] = construct_patch_ten(img, patchSize, slideStep);
        spat_temp_ten(:,:,:,:,ttt) = tenF;
        Position(:,:,:,:,ttt) = patchPosition;

    end
    Nway = size(spat_temp_ten);

    %% Independent component space projection
    for i = 1:Nway(5)
        spat_temp_ten(:,:,:,:,i) = 255*mat2gray(spat_temp_ten(:,:,:,:,i));
    end

    %% The proposed model
    opts=[];
    opts.tol =0.0001;
    opts.max_iter = 400;
    opts.max_mu = 1e5;
    opts.gamma = 1.2;
    opts.lambda1 = lambdaL/((max([Nway(1),Nway(2),Nway(3),Nway(4)])) * Nway(5));
    opts.lambda2 =  100*opts.lambda1; 
    opts.mu = mu;
    opts.R1 =    3 * [0,  1,  1, 1;
                      0,  0,  1, 1;
                      0,  0,  0, 1;
                      0,  0,  0, 0];
    opts.R2 =    3 * [0,  1,  1, 1;
                      0,  0,  1, 1;
                      0,  0,  0, 1;
                      0,  0,  0, 0];
    opts.R3 =    3 * [0,  1,  1, 1, 1;
                      0,  0,  1, 1, 1;
                      0,  0,  0, 1, 1;
                      0,  0,  0, 0, 1;
                      0,  0,  0, 0, 0];
    opts.epsilon = 1e-4;
    opts.max_R1 = opts.R1 /3 * 5;
    opts.max_R2 = opts.R2 /3 * 5;
    opts.max_R3 = opts.R3 /3 * 5;
    opts.x =  x;

    tenT=[];
    tenB=[];
    tenN=[];
    [tenB, tenT, tenN] = LRSD(spat_temp_ten,opts);

    %% reconstrcut images
    for kk = 1:temporal_step 
        tarImg = zeros([imgHei, imgWid]);
        BKGImg = zeros([imgHei, imgWid]);
        noiseImg = zeros([imgHei, imgWid]);
        tarposition = Position(:,:,:,:,kk);
        tar_tmp = tenT(:,:,:,:,kk);
        BKG_tmp = tenB(:,:,:,:,kk);
        noise_tmp = tenN(:,:,:,:,kk);
        for ss = 1:patchNumber
            tar_patch = tar_tmp(:,:,:,ss);
            noise_patch = noise_tmp(:,:,:,ss);
            BKG_patch = BKG_tmp(:,:,:,ss);
            tmp_position = tarposition(:,:,:,ss);
            for jj=1:9
                position = tmp_position(:,:,jj);
                row = position(1);
                col = position(2);
                tarImg(row: row + patchSize/3-1, col:col + patchSize/3-1) = tar_patch(:,:,jj);
                noiseImg(row: row + patchSize/3-1, col:col + patchSize/3-1) = noise_patch(:,:,jj);
                BKGImg(row: row + patchSize/3-1, col:col + patchSize/3-1) = BKG_patch(:,:,jj);
            end
        end
        imwrite(mat2gray(tarImg), [savePath '/' strtok([files{kk+t-1}],'.') '_tar.jpg']);
    end
end
end