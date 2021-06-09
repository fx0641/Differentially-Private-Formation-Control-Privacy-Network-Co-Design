function eSS = calculateSSerror(epsilon, delta, gamma, N, L, adjparam,d)
    eigvec = sort(real(eig(L)));
    lambda2 = eigvec(2);

    Kdelta = qfuncinv(delta);
    c = gamma*(N-1)^2*adjparam^2*d;
    kappa_squared = ((Kdelta +sqrt(Kdelta^2 +2*epsilon))/(2*epsilon))^2;
    g = lambda2*(2-gamma*lambda2);
    eSS = c*kappa_squared/g;
end

