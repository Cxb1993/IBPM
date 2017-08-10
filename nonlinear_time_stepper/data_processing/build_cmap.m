function cmap = build_cmap( levs )


cmap = zeros( levs , 3);
indh = round( levs/2 );

cmap(1:indh-2, 1) = linspace(0,1,indh-2);
cmap(1:indh-2,2) = linspace(0,1,indh-2);
cmap(1:indh-2,3) = 1;
cmap(indh-1 : indh,:) = 1;
cmap(indh +1:end,1) = 1;
cmap(indh +1:end,2) = linspace(1,0,levs - indh);
cmap(indh +1:end,3) = linspace(1,0,levs - indh);
