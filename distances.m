function [ xDistance, yDistance, magDist, angDist ] = distances( x, y )
%DISTANCES calculates pairwise distances between a set of points
%   Calculates distance along x-axis, distance along y-axis, total
%   distance, and angular distance

dim = length(x);

xDistance = zeros(dim, dim);
yDistance = zeros(dim, dim);
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
