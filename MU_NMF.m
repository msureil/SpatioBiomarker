function [U, V] = MU_NMF(A, K, nruns, niters)
    [N, D] = size(A); % A is N x D
    Ubest = zeros(N, K);
    Vbest = zeros(D, K);
    resBest = Inf;

    for run = 1:nruns
        
        U = rand(N, K);
        V = rand(D, K);

        for iter = 1:niters
            
            U = U .* (A * V) ./ (U * (V' * V) + 1e-9);
            V = V .* (A' * U) ./ (V * (U' * U) + 1e-9);

            
            U(U < 0) = 0;
            V(V < 0) = 0;
        end

        
        res = sum((A - U * V').^2, 'all') / (N * D);
        if res < resBest
            Ubest = U;
            Vbest = V;
            resBest = res;
        end
    end

    U = Ubest;
    V = Vbest;
end
