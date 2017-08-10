function sfn_tr = trunc_fort_sfn( sfn, parms )

%Takes streamfunction from the fortran code (which lives on multiple
%grids)and truncates the parts the overlap with finer grids to return a
%streamfunction that is consistent with the matlab base flow solver.

mg = parms.mg; m = parms.m; n = parms.n;

%size of streamfunction on truncated grid used by Matlab base flow solver
sfn_tr_sz = get_vort_ind( m-1, n-1, mg, parms );
sfn_tr = zeros( sfn_tr_sz, 1 );
sfn_tr( 1 : (m-1)*(n-1) ) = sfn( :, 1);
%1st grid level is unchanged


for lev = 2 : mg

    %--bottom part of streamfunction (no overlap with finer grid)

        %indices for truncated grid:
        ind_tr = get_vort_ind( 1, 1, lev, parms) : ...
            get_vort_ind( m-1, n/4, lev, parms);
        
        sfn_tr( ind_tr ) = sfn( 1 : (n/4)*(m-1), lev );

    %--
    
    %--top part (no overlap with finer grid)
    
        %indices for truncated grid
        ind_tr = get_vort_ind( 1, 3*n/4, lev, parms ) : ...
            get_vort_ind( m-1, n-1, lev, parms );
        sfn_tr( ind_tr ) = sfn( (3*n/4-1)*(m-1) + 1 : (n-1)*(m-1), lev );
    
    %--
    
    %--middle part (overlap here we have to truncate)
    
        %indices for truncated grid
        ind_tr = get_vort_ind( 1, n/4+1, lev,parms ) : ...
            get_vort_ind( m-1, 3*n/4-1, lev, parms) ;
        
        %indices for big grid
            %first row
            ind1 = n/4*(m-1) + [1 : m/4, 3*m/4 : m-1];
            Ind1 = repmat( ind1, [1, n/2-1] );
            
            %indices to add for remaining rows
            indadd = repelem( 0 : (m-1) : (n/2-2) ...
                * (m-1), length(ind1) );
            
            %compile this into a set of indices for sfn on big grid
            indbg = Ind1 + indadd;
            
        %assign values to sfn
        sfn_tr( ind_tr ) = sfn( indbg, lev );
            
            
    %--

    
end
