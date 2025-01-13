function nodes = cheb_nodes(n, smin, smax)
    % CHEB_NODES Compute Chebyshev nodes
    N = n:-1:1;
    nodes = (smin+smax)/2 + (smax-smin)/2 * cos((2*N-1)*pi/(2*n));
    
end