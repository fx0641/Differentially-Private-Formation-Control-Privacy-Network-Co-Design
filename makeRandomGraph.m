function [L, eps_min, eps_max, gammai, taui, rhoi] = makeRandomGraph(N)
    probability = .30;
    noOfEdges = round(N*(N-1)/2*probability); 
    
    x = [ones(1, noOfEdges) zeros(1,N*(N-1)/2-noOfEdges)];
    x = -x(randperm(N*(N-1)/2));
    
    L = makeLaplacian(x,N);
    
    eps_min = abs(.1*randn(N,1) +.1);
    eps_max = abs(.2*randn(N,1) +.6);
    
    gammai = abs(.4*randn(N,1));
    taui = abs(.2*randn(N,1) +.25);
    rhoi = abs(.2*randn(N,1) +15);
end

