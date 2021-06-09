function L = makeLaplacian(x,N)
% This function takes the decision variables and formats them into a 
% graph laplacian.
    upperL = tril(ones(N),-1);
    upperL(upperL==1) =x(1:N*(N-1)/2);
    L = upperL + upperL';

    diagcomps = L*ones(N,1);
    for i =1:N
        L(i,i) = -diagcomps(i);
    end
end

