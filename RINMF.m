% Alternating Direction Method of Multipliers (ADMM) based NMF which
% enforce integer value costrains. 
function [U,V] = RINMF(A,K,nruns,niters,rho,flag)
N = size(A,1);
D = size(A,2);

resBest = Inf;
Ubest = zeros(N,K);
Vbest = zeros(D,K);
rhoI = rho*eye(K);
% rng(s); % initialize rng if you want the same numbers again
for run = 1:nruns    
    % random initialization 
    V = .01*normrnd(0,1,[D,K]); 
  
    V(V<0) = 0; % use of non-negative initial is not required
    

    Y1 = zeros(N,K);
    Y2 = zeros(D,K);
    Z1 = Y1;
    Z2 = Y2;

    for iter = 1:niters
        U = (A*V+rho*(Z1-Y1))/(V'*V+rhoI);
        V = (A'*U+rho*(Z2-Y2))/(U'*U+rhoI);
        
        Z1 = U + Y1;
        Z1(Z1 < 0) = 0;
        if flag(1)
            Z1 = round(Z1);
        end
        
        Z2 = V + Y2;
        Z2(Z2 < 0) = 0;
        if flag(2)
            Z2 = round(Z2);
        end
        Y1 = Y1+U-Z1;
        Y2 = Y2+V-Z2;
    end
    clear iter
    
    % for ADMM, non-negative is only guaranteed in the limit of convergence
    % therefore negative values in output is a sign of non-convergence.
    % However, tiny negative values will always occur because convergence
    % is asymptotic. If we are confident of convergence, we can zero out
    % the small non-negative values
    U(U<0) = 0;
    V(V<0) = 0;
    % likewise, we can polish the integer contraints if we want to
    if flag(1)
        U = round(U);
    end
    if flag(2)
        V = round(V);
    end
   
    % compute residual
    res = sum((A -U*V').^2,'all')/(N*D);
    if res < resBest
        Ubest = U;
        Vbest = V;
        resBest = res;
    end
end
clear run
U = Ubest;
V = Vbest;
res = resBest;


end

