function [ xDistance, yDistance, magDist, angDist ] = distances( x, y )
%DISTANCES calculates pairwise distances between a set of points
%   Calculates distance along x-axis, distance along y-axis, total
%   distance, and angular distance

dim = length(x);

xDistance = zeros(dim, dim);
yDistance = zeros(dim, dim);
magDist = zeros(dim, dim);
angDist = zeros(dim, dim);

for i=1:dim
    for j=1:dim
        
        xDistance(i, j) = x(i) - x(j);
        yDistance(i, j) = y(i) - y(j);
        angle1 = normalizeAngle(atan(y(i)/x(i)));
        angle2 = normalizeAngle(atan(y(j)/x(j)));
        angDist(i, j) = normalizeAngle(angle2 - angle1, 0);
        
    end
end

magDist = sqrt(xDistance.^2 + yDistance.^2);

end

