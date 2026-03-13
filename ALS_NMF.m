function [U, V] = ALS_NMF(A, K, nruns, niters)
    [N, D] = size(A);
    Ubest = zeros(N, K);
    Vbest = zeros(D, K);
    resBest = Inf;

    for run = 1:nruns
        U = rand(N, K);
        V = rand(D, K);

        for iter = 1:niters
            for i = 1:N
                U(i, :) = lsqnonneg(V, A(i, :)');
            end

            for j = 1:D
                V(j, :) = lsqnonneg(U, A(:, j));
            end

 
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
