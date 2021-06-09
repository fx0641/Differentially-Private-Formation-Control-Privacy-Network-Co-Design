function [c,ceq] = constraints(x,N,eps_min, eps_max, max_error, gamma, delta, adjparam, edge_budget, gammai, taui, rhoi, givenL,d)
    %First lets make a laplacian
    L = makeLaplacian(x,N);
    
    eigvals = sort(real(eig(L)));

    lambda2 = eigvals(2);
    
    %c(1) = 2*abs(min(x(1:N*(N-1)/2)))*(1-cos(pi/N)) - lambda2; %make sure L is connected
    c(1) = .3-lambda2;
    c(2:1+N*(N-1)/2) = x(1:N*(N-1)/2);%make sure off diagonals are negative
    c(2+N*(N-1)/2) = .5*trace(L) - edge_budget;
    
    %Epsilon stuff, this is for the homogenous case
    epsilon = x(N*(N-1)/2+1:end);

    c(3+N*(N-1)/2:2+N*(N-1)/2+N) = epsilon-eps_max;
    
    c(3+N*(N-1)/2+N:2+N*(N-1)/2+2*N) = gammai.*diag(L) + taui.*epsilon - rhoi;
    
    
    %The gnarly constraint   
    eSS = calculateSSerror(min(epsilon), delta, gamma, N, L, adjparam,d);

    c(3+N*(N-1)/2+2*N) = eSS - max_error;
    
    c(4+N*(N-1)/2+2*N:3+N*(N-1)/2+3*N) = eps_min-epsilon;
    
    % Zero out the terms that have no edge in the given L
    idx = logical(tril(ones(size(givenL)), -1));
    x_from_givenL = givenL(idx);
    location_of_zeros = find(~x_from_givenL);
    ceq(1:length(location_of_zeros)) = x(location_of_zeros)+.00001;
end