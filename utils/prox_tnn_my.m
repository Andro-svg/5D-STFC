function [X] = prox_tnn_my(Z, YY, rho, epsilon)

dim = ndims(YY);
if dim == 4
    [n1, n2, n3, n4] = size(YY);
    for k = 1:n4
        Y = YY(:,:,:,k);
        n12 = min(n1, n2);
        Yf = fft(Y, [], dim);
        Uf = zeros(n1, n12, n3);
        Vf = zeros(n2, n12, n3);
        Sf = zeros(n12,n12, n3);
        
        % Yf(isnan(Yf)) = 0;
        % Yf(isinf(Yf)) = 0;
        
        Zf = fft(Z, [], dim);
        U = zeros(n1, n12, n3);
        V = zeros(n2, n12, n3);
        S = zeros(n12,n12, n3);
        
        % Z(isnan(Z)) = 0;
        % Z(isinf(Z)) = 0;
        
        trank = 0;
        endValue = n3/2 + 1;
        for i = 1 : endValue
        %     sf=[];
        %     s=[];
            [Uf(:,:,i), Sf(:,:,i), Vf(:,:,i)] = svd(Yf(:,:,i), 'econ');
            [U(:,:,i), S(:,:,i), V(:,:,i)] = svd(Zf(:,:,i), 'econ');
            sf = diag(Sf(:, :, i));
            s = diag(S(:,:,i));
        %     Sf_(:, :, i) = max(sf-rho/epsilon*exp(-1/epsilon*s),0);
            Sf(:, :, i) = diag(max(sf-rho/epsilon*exp(-s/epsilon),0));
            temp = length(find(Sf(:, :, i)>0));
            trank = max(temp, trank);
        
        end
        
        for j =n3:-1:endValue+1
            Uf(:,:,j) = conj(Uf(:,:,n3-j+2));
            Vf(:,:,j) = conj(Vf(:,:,n3-j+2));
            Sf(:,:,j) = Sf(:,:,n3-j+2);
        end
        
        Uf = Uf(:, 1:trank, :);
        Vf = Vf(:, 1:trank, :);
        Sf = Sf(1:trank, 1:trank, :);
        
        U = ifft(Uf, [], dim);
        S = ifft(Sf, [], dim);
        V = ifft(Vf, [], dim);
        
        X(:,:,:,k) = tprod( tprod(U,S), tran(V) );
    end
else
    [n1, n2, n3, n4, n5] = size(YY);
    for m = 1:n4
        for n= 1:n5
            Y = YY(:,:,:,m,n);
            n12 = min(n1, n2);
            Yf = fft(Y, [], dim);
            Uf = zeros(n1, n12, n3);
            Vf = zeros(n2, n12, n3);
            Sf = zeros(n12,n12, n3);
            
            % Yf(isnan(Yf)) = 0;
            % Yf(isinf(Yf)) = 0;
            
            Zf = fft(Z, [], dim);
            U = zeros(n1, n12, n3);
            V = zeros(n2, n12, n3);
            S = zeros(n12,n12, n3);
            
            % Z(isnan(Z)) = 0;
            % Z(isinf(Z)) = 0;
            
            trank = 0;
            endValue = n3/2 + 1;
            for i = 1 : endValue
            %     sf=[];
            %     s=[];
                [Uf(:,:,i), Sf(:,:,i), Vf(:,:,i)] = svd(Yf(:,:,i), 'econ');
                [U(:,:,i), S(:,:,i), V(:,:,i)] = svd(Zf(:,:,i), 'econ');
                sf = diag(Sf(:, :, i));
                s = diag(S(:,:,i));
            %     Sf_(:, :, i) = max(sf-rho/epsilon*exp(-1/epsilon*s),0);
                Sf(:, :, i) = diag(max(sf-rho/epsilon*exp(-s/epsilon),0));
                temp = length(find(Sf(:, :, i)>0));
                trank = max(temp, trank);            
            end
            
            for j =n3:-1:endValue+1
                Uf(:,:,j) = conj(Uf(:,:,n3-j+2));
                Vf(:,:,j) = conj(Vf(:,:,n3-j+2));
                Sf(:,:,j) = Sf(:,:,n3-j+2);
            end
            
            Uf = Uf(:, 1:trank, :);
            Vf = Vf(:, 1:trank, :);
            Sf = Sf(1:trank, 1:trank, :);
            
            U = ifft(Uf, [], dim);
            S = ifft(Sf, [], dim);
            V = ifft(Vf, [], dim);
            
            X(:,:,:,m,n) = tprod( tprod(U,S), tran(V) );
        end
    end
end
   
end
