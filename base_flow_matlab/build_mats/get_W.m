function W = get_W( grid_parms )

%Build averaging operator W = [Wx; Wy]. The x-velocity block (Wx) takes 
%vorticity and averages it onto the x-velocity edges, and the y-velocity
%block (Wy) averages it onto y-velocity edges.

%Note that W is not the full W used in the code. The W in the code is M *
%the W computed here.

%Inputs: grid_parms -- data structure containing m (number of points in x
%dirn), n (number of points in y dirn), and mg (number of grid levels)

m = grid_parms.m; n = grid_parms.n; mg = grid_parms.mg; 

%Get size of W
nrows = get_velx_ind( m-1, n, mg, grid_parms ) + ...
    get_vely_ind( m, n-1, mg, grid_parms );
ncols = get_vort_ind( m-1, n-1, mg, grid_parms );

grid_parms.nrows = nrows;
grid_parms.ncols = ncols;

%Initialize
W = sparse( nrows, ncols );

%Loop through gridlevels
for glev = 1 : mg
    
    %Get main blocks (not accounting for BC's)
    W = get_W_main( W, glev, grid_parms );
    
    %BC contributions from coarser grid on current grid
    if mg > 1 & glev < mg
        
        W = get_W_fineedge_BCs( W, glev, grid_parms );
        
    end
    
end

grid_parms = rmfield( grid_parms, 'nrows'); 
grid_parms = rmfield( grid_parms, 'ncols');
