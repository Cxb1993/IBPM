function [] = build_inp(stp,isv,params)

%Build input file:

fileID = fopen('input/ib.inp','w');
fprintf(fileID, '%-20s \n','&READ_PARAMETERS');
fprintf(fileID, '%-11s \n','ISTART = 0,');
fprintf(fileID, '%-8s %-1d%-1s \n','ISTOP = ', stp, ',');
fprintf(fileID, '%-8s %-1d%-1s \n','ISAVE = ',isv, ',');
fprintf(fileID, '%-8s %-1d%-1s \n','M = ',params.m,',');
fprintf(fileID, '%-8s %-1d%-1s \n','N = ',params.n,',');
fprintf(fileID, '%-8s %-1f%-1s \n','DT = ',params.dt,',');
fprintf(fileID, '%-8s %-1f%-1s \n','Re = ', params.Re,',');
fprintf(fileID, '%-20s \n','CGTOL = 1.E-8, ');
fprintf(fileID, '%-20s \n','CG_MAX_ITER = 3000,');
fprintf(fileID, '%-8s %-1f%-1s \n','LEN = ', params.len,',');
fprintf(fileID, '%-8s %-1f%-1s \n','OFFSETX = ',params.offx,',');
fprintf(fileID, '%-8s %-1f%-1s \n','OFFSETY = ',params.offy,',');
fprintf(fileID, '%-8s %-1d%-1s \n','MGRIDLEV = ', params.mg , ',');
fprintf(fileID, '%-20s \n','COMPUTE_PRESSURE = F,');
fprintf(fileID, '%-8s %-1s%-1s \n', 'STANDARD = ', params.standard,',');
fprintf(fileID, '%-8s %-1s%-1s \n', 'PINNED = ', params.pinned,',');
fprintf(fileID, '%-8s %-1.15f%-1s \n', 'R_rho = ', params.R_rho,',');
fprintf(fileID, '%-8s %-1.15f%-1s \n', 'R_E = ', params.R_E,',');
fprintf(fileID, '%-8s %-1.15f%-1s \n', 'R_sh = ', params.R_sh,',');
fprintf(fileID, '%-8s %-1.15f%-1s \n', 'R_th = ', params.R_th,',');
fprintf(fileID, '%-20s \n','/');
fclose(fileID);


