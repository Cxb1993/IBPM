function [W, Wsq] = build_W_new( parms )

nb = parms.nb;

%# of x-vel (flux) points
nu = get_velx_ind( parms.m-1, parms.n, parms.mg, parms );
%# of y-vel (flux) points
nv = get_vely_ind( parms.m, parms.n-1, parms.mg, parms );
%Total # of vel (flux) points
nq = nu + nv;

h = parms.len / parms.m;




% nq = 4; nb = 4;
% h = 1;




n = nq + 4 * nb;


W = sparse( n, n );

W = W + sparse( 1 : nq, 1 : nq, ones(size(1 : nq)), n, n );

%--Structural part

%x-displacements

    %diag
    ind1 = 1 : nb; ind2 = ind1;
    W = W - sparse( nq + ind1, nq  + ind2, -2 * ...
        ones(size( ind1 ) )/ (h^2), n, n ) * parms.Kb;

    %superdiag
    ind1 = 1 : nb-1; ind2 = 2 : nb;
    W = W - sparse( nq + ind1, nq + ind2, 1 * ...
        ones(size( ind1 ) ) / (h^2), n, n ) * parms.Kb;

    %subdiag
    ind1 = 2 : nb; ind2 = 1 : nb-1;
    W = W - sparse( nq + ind1, nq + ind2, 1 * ...
        ones(size( ind1 ) ) / (h^2), n, n ) * parms.Kb;
    
%y-displacements

    %diag
    ind1 = 1 : nb; ind2 = ind1;
    W = W - sparse( nq + nb + ind1, nq + nb + ind2, -2 * ...
        ones(size( ind1 ) )/ (h^2), n, n ) * parms.Kb;
    
    %superdiag
    ind1 = 1 : nb-1; ind2 = 2 : nb;
    W = W - sparse( nq + nb + ind1, nq + nb + ind2, 1 * ...
        ones(size( ind1 ) ) / (h^2), n, n ) * parms.Kb;
    
    %subdiag
    ind1 = 2 : nb; ind2 = 1 : nb-1;
    W = W - sparse( nq + nb + ind1, nq + nb + ind2, 1 * ...
        ones(size( ind1 ) ) / (h^2), n, n ) * parms.Kb;
    
% velocities
    ind = 1 : 2*nb;
    W = W + parms.M_rho * sparse( nq + 2*nb + ind, nq + 2*nb + ind, ...
        ones(size( ind ) ), n, n ); 

%--

Wsq = W * W;

% max(max( abs( W - W' ) ) )
% 
% max(max(abs( Wsq - Wsq' ) ) )
% 
% full( W ) 
% % max(max( abs( Wsqrt*Wsqrt - W ) ) )
% % max(max( abs( Wsqrt - Wsqrt' ) ) )
% [V,D] = eig( full(W) );
% diag(D)




