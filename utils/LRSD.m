function [B,T,Noise] = LRSD(D, opts)

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');             tol       = opts.tol;             end
if isfield(opts, 'max_iter');        max_iter  = opts.max_iter;        end
if isfield(opts, 'max_mu');          max_mu    = opts.max_mu;          end
if isfield(opts, 'gamma');           gamma     = opts.gamma;           end
if isfield(opts, 'lambda1');         lambda1   = opts.lambda1;         end
if isfield(opts, 'lambda2');         lambda2   = opts.lambda2;         end
if isfield(opts, 'R1');              R1        = opts.R1;              end
if isfield(opts, 'R2');              R2        = opts.R2;              end
if isfield(opts, 'R3');              R3        = opts.R3;              end
if isfield(opts, 'epsilon');         epsilon   = opts.epsilon;         end
if isfield(opts, 'max_R1');          max_R1    = opts.max_R1;          end
if isfield(opts, 'max_R2');          max_R2    = opts.max_R2;          end
if isfield(opts, 'max_R3');          max_R3    = opts.max_R3;          end
if isfield(opts, 'mu');              mu        = opts.mu;              end
if isfield(opts, 'x');               x         = opts.x;               end

%% Initialization
Nway = size(D);
TT = Nway(5);
SS = Nway(4);
N = ndims(D);

B = zeros(Nway);
T = zeros(Nway);
Noise = zeros(Nway);

M1 = zeros(Nway);
Gt = cell(TT,4);
Qt = cell(TT,4);
M2 = cell(TT,4);
Gtway = [Nway(1),Nway(2),Nway(3),Nway(4)];
Gtdim = diag(Gtway)+R1+R1';
max_Gtdim = diag(Gtway)+max_R1+max_R1';
for j = 1:Nway(5)
    for i = 1:4
        Gt{j,i} = ones(Gtdim(i,:));
        Qt{j,i} = ones(Gtdim(i,:));
        M2{j,i} = zeros(Gtdim(i,:));
    end
end

Gs = cell(SS,4);
Qs = cell(SS,4);
M3 = cell(SS,4);
Gsway = [Nway(1),Nway(2),Nway(3),Nway(5)];
Gsdim = diag(Gsway)+R2+R2';
max_Gsdim = diag(Gsway)+max_R2+max_R2';
for j = 1:Nway(4)
    for i = 1:4
        Gs{j,i} = ones(Gsdim(i,:));
        Qs{j,i} = ones(Gsdim(i,:));
        M3{j,i} = zeros(Gsdim(i,:));
    end
end

G = cell(1,N);
Q = cell(1,N);
M4 = cell(1,N);
Gway = [Nway(1),Nway(2),Nway(3),Nway(4),Nway(5)];
Gdim = diag(Gway)+R3+R3';
max_Gdim = diag(Gway)+max_R3+max_R3';
for i = 1:N
    G{i} = ones(Gdim(i,:));
    Q{i} = ones(Gdim(i,:));
    M4{i} = zeros(Gdim(i,:));
end

M5 = cell(1,Nway(4));
for i = 1:Nway(4)
    M5{i} = zeros([Nway(1) Nway(2) Nway(3) Nway(5)]);
end

M6 = cell(1,Nway(5));
for i = 1:Nway(5)
    M6{i} = zeros([Nway(1) Nway(2) Nway(3) Nway(4)]);
end

M7 = zeros(Nway);
W = ones(size(T));
preNumT = numel(T);

for iter = 1 : max_iter
    %% update Gs
    for i = 1:Nway(4)
        for j = 1:4
            V = squeeze(B(:,:,:,i,:)) - M5{i};
            Vi = my_Unfold(V,size(V),j);
            Qtii = my_Unfold(Qs{i,j},size(Qs{i,j}),j);
            M3ii = my_Unfold(M3{i,j},size(M3{i,j}),j);
            Gi = my_Unfold(Gs{i,j},Gsdim(j,:),j);
            Girest = tnreshape(tnprod_rest(Gs(i,:),j),4,j);
            tempC = (Qtii-M3ii)+Vi*Girest';
            tempA = (Girest*Girest')+eye(size(Gi,2));
            Gs{i,j}  = real(my_Fold(tempC*pinv(tempA),Gsdim(j,:),j));
        end
    end

    %% update Gt
    for i = 1:Nway(5)
        for j = 1:4
            V = squeeze(B(:,:,:,:,i))-M6{i};
            Vi = my_Unfold(V,size(V),j);
            Qtii = my_Unfold(Qt{i,j},size(Qt{i,j}),j);
            M2ii = my_Unfold(M2{i,j},size(M2{i,j}),j);
            Gi = my_Unfold(Gt{i,j},Gtdim(j,:),j);
            Girest = tnreshape(tnprod_rest(Gt(i,:),j),4,j);
            tempC = (Qtii-M2ii)+Vi*Girest';
            tempA = (Girest*Girest')+eye(size(Gi,2));
            Gt{i,j}  = real(my_Fold(tempC*pinv(tempA),Gtdim(j,:),j));
        end
    end

    %% update G
    for j = 1:N
        V = B-M7;
        Vi = my_Unfold(V,size(V),j);
        Qii = my_Unfold(Q{j},size(Q{j}),j);
        M4ii = my_Unfold(M4{j},size(M4{j}),j);
        Gi = my_Unfold(G{j},Gdim(j,:),j);
        Girest = tnreshape(tnprod_rest(G,j),N,j);
        tempC = (Qii-M4ii)+Vi*Girest';
        tempA = (Girest*Girest')+eye(size(Gi,2));
        G{j}  = real(my_Fold(tempC*pinv(tempA),Gdim(j,:),j));
    end

    %% update Qs
    alpha2 = 0.25;
    beta2 = 0.25;
    tau2 = alpha2 * beta2 / mu;
    for i = 1:Nway(4)
        for j=1:4
            Qs{i,j}  = prox_tnn_my( Qs{i,j}, Gs{i,j} + M3{i,j}, tau2,  epsilon);
        end
    end

    %% update Qt
    alpha1 = 0.25;
    beta1 = 0.25;
    tau1 = alpha1 * beta1 / mu;
    for i = 1:Nway(5)
        for j=1:4
            Qt{i,j}  = prox_tnn_my( Qt{i,j}, Gt{i,j} + M2{i,j}, tau1,  epsilon);
        end
    end

    %% update Q
    alpha3 = 0.5;
    beta3 = 0.2;
    tau3 = alpha3 * beta3 / mu;
    for i = 1:N 
        Q{i}  = prox_tnn_my( Q{i}, G{i} + M4{i}, tau3,  epsilon);
    end

    %% update T
    T = (D-B-Noise-M1) .* max(W .* T - 2*sqrt(x),0) + (D-B-Noise-M1)./(lambda1 .*  W .*  W ./(mu*x) + ones(size(T))) .* max(2*sqrt(x)- W .* T,0);
    
    %% update B_4D_s
    for i = 1:Nway(4)    
        tmp = tnprod(Gs(i,:));
        tempGs = reshape(tmp, size(tmp, 1), size(tmp, 2), size(tmp, 3), 1, size(tmp, 4));
        M5i = reshape(M5{i}, size(M5{i}, 1), size(M5{i}, 2), size(M5{i}, 3), 1, size(M5{i}, 4));
        B(:,:,:,i,:) = double((D(:,:,:,i,:)-T(:,:,:,i,:)-Noise(:,:,:,i,:)+ tempGs - M1(:,:,:,i,:)+M5i)/2);
    end

    %% update B_4D_t
    for i = 1:Nway(5)
        B(:,:,:,:,i) = double((D(:,:,:,:,i)-T(:,:,:,:,i)-Noise(:,:,:,:,i)+tnprod(Gt(i,:))-M1(:,:,:,:,i)+M6{i})/2);
    end

    %% update B_5D
    B = double((D-T-Noise+tnprod(G)-M1-M7)/2);

    %% update W
    W = 1 ./ ((abs(T))+ 0.01);

    %% update Noise
    Noise = mu*(D-B-T-M1)/(2*lambda2+mu);
        
    %% check the convergence
    currNumT = sum(T(:) > 0); 
    chg =norm(D(:)-T(:)-B(:)-Noise(:))/norm(D(:));
    fprintf('iter = %d   res=%.10f  \n', iter, chg);
    if (chg < tol) || (currNumT == preNumT)
        break;
    end
    preNumT = currNumT;

    %% update Lagrange multipliers M
    M1 = M1 - (D-B-T-Noise);
    for j = 1:Nway(5)
        for i = 1:4
            M2{j,i} = M2{j,i} - (Qt{j,i} - Gt{j,i});
        end
    end
    for j = 1:Nway(4)
        for i = 1:4
            M3{j,i} = M3{j,i} - (Qs{j,i} - Gs{j,i});
        end
    end
    for i = 1:N
        M4{i} = M4{i}-(Q{i}-G{i});
    end
    for i = 1:Nway(4)
        M5{i} = M5{i} - (squeeze(B(:,:,:,i,:))-tnprod(Gs(i,:)));
    end
    for i = 1:Nway(5)
        M6{i} = M6{i} - (B(:,:,:,:,i)-tnprod(Gt(i,:)));
    end
    M7 = M7 - (B-tnprod(G));

    %% update penalty parameter mu
    mu = min(gamma*mu,max_mu); 

    %% update the estimated rank
    Gt_rank_inc=double(Gtdim<max_Gtdim);
    if sum(Gt_rank_inc(:))~=0
        for j = 1:Nway(5)
            Gt(j,:) = rank_inc_adaptive(Gt(j,:),Gt_rank_inc,4);
            Qt(j,:) = rank_inc_adaptive(Qt(j,:),Gt_rank_inc,4);
            M2(j,:) = rank_inc_adaptive(M2(j,:),Gt_rank_inc,4);
        end
        Gtdim = Gtdim+Gt_rank_inc;
    end
    
    Gs_rank_inc=double(Gsdim<max_Gsdim);
    if sum(Gs_rank_inc(:))~=0
        for j = 1:Nway(4)
            Gs(j,:) = rank_inc_adaptive(Gs(j,:),Gs_rank_inc,4);
            Qs(j,:) = rank_inc_adaptive(Qs(j,:),Gs_rank_inc,4);
            M3(j,:) = rank_inc_adaptive(M3(j,:),Gs_rank_inc,4);
        end
        Gsdim = Gsdim+Gs_rank_inc;
    end
    
    G_rank_inc=double(Gdim<max_Gdim);
    if sum(G_rank_inc(:))~=0
        G = rank_inc_adaptive(G,G_rank_inc,N);
        Q = rank_inc_adaptive(Q,G_rank_inc,N);
        M4 = rank_inc_adaptive(M4,G_rank_inc,N);
        Gdim = Gdim+G_rank_inc;
    end
   
end
end

function [G]=rank_inc_adaptive(G,rank_inc,N)
    for j = 1:N
        G{j} = padarray(G{j},rank_inc(j,:),1,'post');
    end
end