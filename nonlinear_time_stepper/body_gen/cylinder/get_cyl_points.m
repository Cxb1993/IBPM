function n_pts = get_cyl_points(circum,h)

%Get the number of points required for a cylinder of circumference circum
%to have the same grid spacing as h

n_pts = floor(circum/h) / 2;