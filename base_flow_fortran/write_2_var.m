function [] = write_2_var( V, params  )

nq = params.nq;
mg = params.mg; 
nb = params.nb;


u0sm = V(1 : params.nvel_tot); 

u0bg = engorge( u0sm, params );



u0 = reshape( u0bg, nq, mg );


[~,~] = system('rm output/*.chd');

[~,~] = system('rm output/q_init.var');

fileID = fopen('output/q_init.var','w','b');
for j = 1 : nq
    fwrite(fileID,u0(j,:),'real*8','b');
end
fclose(fileID);


uib0 = V( params.nvel_tot + (1 : 3*nb ) );


[~,~] = system('rm output/uib_init.var');

fileID = fopen('output/uib_init.var','w','b');
for j = 1 : 3 * nb
    fwrite(fileID,uib0(j),'real*8','b');
end
fclose(fileID);


udib0 = V( params.nvel_tot + 3*nb + (1 : 3*nb ) );


[~,~] = system('rm output/udib_init.var');

fileID = fopen('output/udib_init.var','w','b');
for j = 1 : 3 * nb
    fwrite(fileID,udib0(j),'real*8','b');
end
fclose(fileID);


fb0 = V( params.nvel_tot + 6*nb + ( 1 : 2*nb ) );

fb0 = fb0 * params.dt;


[~,~] = system('rm output/fb_init.var');

fileID = fopen('output/fb_init.var','w','b');
for j = 1 : 2 * nb
    fwrite(fileID,fb0(j),'real*8','b');
end
fclose(fileID);


