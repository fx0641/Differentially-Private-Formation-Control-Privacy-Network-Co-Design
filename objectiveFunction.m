function obj = objectiveFunction(x,N)
    %First lets make a laplacian
    L = makeLaplacian(x,N);
    
    epsilon = x(N*(N-1)/2+1:end);
    %obj = x(1:N*(N-1)/2)'*ones(N*(N-1)/2,1) + sum(1./epsilon.^2) + .5*norm(x(1:N*(N-1)/2),1);
    obj = trace(L) + sum(1./epsilon.^2);
end