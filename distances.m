function [ xDistance, yDistance, magDist, angDist ] = distances( x, y )
%DISTANCES calculates pairwise distances between a set of points
%   Calculates distance along x-axis, distance along y-axis, total
%   distance, and angular distance

dim = length(x);

xDistance = zeros(dim, dim);
yDistance = zeros(dim, dim);
<<<<<<< HEAD
magDist = zeros(dim, dim);
=======
>>>>>>> 75910a67982bd62f51b4532644c8c042fd7c58c6
angDist = zeros(dim, dim);

for i=1:dim
    for j=1:dim
        
        xDistance(i, j) = x(i) - x(j);
        yDistance(i, j) = y(i) - y(j);
        
        angDist(i, j) = atan2(y(j), x(j)) - atan2(y(i), x(i));
        
    end
end

magDist = sqrt(xDistance.^2 + yDistance.^2);
angDist = mod(angDist, 2*pi);

end
<<<<<<< HEAD

=======
>>>>>>> 75910a67982bd62f51b4532644c8c042fd7c58c6
