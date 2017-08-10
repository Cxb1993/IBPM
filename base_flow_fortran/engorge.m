function q = engorge( qsm, params )

%Take velocity matrix of size nq x mg as given by fortran code and truncate
%to eliminate redundant information on overlapping grids. Also, convert the
%velocity as given by the code into physical units.

nq = params.nq; mg = params.mg ;  len = params.len;
m = params.m; n = params.n; inx = params.inx; iny = params.iny;
bndryx = params.bndryx; bndryy = params.bndryy; 
indbgx = params.indbgx; indbgy = params.indbgy;


q = zeros( mg * nq, 1);


%Grid 1:

indx = indbgx;
indy = indbgy;
indx(bndryx) = 0;
indy(bndryy) = 0;

qtr = qsm( 1 : params.nmg1 );

%Scale to get fluxes from vels:
fac = 1;
dlt = len ./ m *fac;
qtr = qtr .* dlt;

qxtr = qtr( 1 : sum(indx~=0));
qytr = qtr( sum(indx~=0) + 1 : end );

qx = zeros( (m+1)*n, 1 );
qy = zeros( (n+1)*m, 1);

qx( indx~= 0 ) = qxtr;
qy( indy~= 0 ) = qytr;

q( 1 : nq ) = [qx; qy];


if mg > 1
    
    for j = 2 : mg 
        
        qtr = qsm( params.nmg1 + 1 + (j-2)*params.nmid : ...
            params.nmg1 + (j-1)*params.nmid );
        
        %Scale to get velocities from fluxes:
        fac = 2^(j-1);
        dlt = len ./ m *fac;
        qtr = qtr .* dlt;

   
        indx = indbgx;
        indx(bndryx) = 0;
        indx(inx) = 0;
        indy = indbgy;
        indy(bndryy) = 0;
        indy(iny) = 0;
        
        qxtr = qtr( 1 : sum(indx~=0));
        qytr = qtr( sum(indx~=0) + 1 : end );
        
        qx = zeros( (m+1)*n, 1 );
        qy = zeros( (n+1)*m, 1);

        qx( indx~= 0 ) = qxtr;
        qy( indy~= 0 ) = qytr;

        q( (j-1)*nq + (1 : nq) ) = [qx; qy];
        
        
    end
    
        
%     qtr = qsm( params.nmg1 + 1 + (mg-2)*params.nmid : ...
%             params.nmg1 + (mg-2)*params.nmid + params.nmglst );
%         
%         %Scale to get velocities from fluxes:
%         fac = 2^(mg-1);
%         dlt = len ./ m *fac;
%         qtr = qtr .* dlt;
% 
%    
%         indx = indbgx;
%         indx(inx) = 0;
%         indy = indbgy;
%         indy(iny) = 0;
%         
%         qxtr = qtr( 1 : sum(indx~=0));
%         qytr = qtr( sum(indx~=0) + 1 : end );
%         
%         qx = zeros( (m+1)*n, 1 );
%         qy = zeros( (n+1)*m, 1);
%         
%         qx( indx~= 0 ) = qxtr;
%         qy( indy~= 0 ) = qytr;
% 
%         q( (mg-1)*nq + (1 : nq) ) = [qx; qy];
%     
    
    
end






