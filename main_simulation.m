%% DEFINE PARAMETERS
% 

% number of timesteps
numSteps = 1000000;

% number of yield and growth strategists for seeding the simulation
numYS = 10; % cell type "0"
numGS = 10; % cell type "1"
numCells = numYS + numGS;
initCellTypes = [zeros(1, numYS) ones(1, numGS)];

% maximum dimensions of simulation grid and cell size
xMax = 25;
yMax = 25;
cellDiameter = 1;

% define grid size and initialize nutrient grid
nutrientGridSize = 1;
nutrientGrid = zeros(xMax/nutrientGridSize, yMax/nutrientGridSize);

% tracks newest ID of cells to be created (first "numCells" cells already initialized at beginning of simulation)
nextCellID = numCells + 1;

%% INITIALIZE GRID/OBJECTS
%

% randomly generate initial positions of cells on simulation grid, then initialize cell objects in same loop
for i=1:numCells
    
    % sequentially place/initialize each cell--try until no collision
    placed = false;
    
    while (placed == false)
        
        % randomly initialize positions of all particles (no wall collisions allowed)
        xPos(i) = 0.5*cellDiameter + rand()*(xMax - cellDiameter);
        yPos(i) = 0.5*cellDiameter + rand()*(yMax - cellDiameter);

        % check for collisions with other particles/wall--if none, continue to next particle
        if (sum(pdist([xPos(1:i); yPos(1:i)]') < cellDiameter) == 0)
            placed = true;
        end
        
    end
    
    % initialize cell objects based on randomly generated initial positions
    activeCells(i) = cell_obj(i, initCellTypes(i), xPos(i), yPos(i));
    
end

% draw grid with all active cells (YS strategists blue, GS strategists red)
h = zeros(1, length(activeCells));

for i=1:length(activeCells)
    
    if (activeCells(i).cellType == 0)
        h(i) = rectangle('Position', [(activeCells(i).xCoor - cellDiameter/2) (activeCells(i).yCoor - cellDiameter/2) cellDiameter cellDiameter], 'Curvature', [1, 1], 'edgecolor', 'b');
    elseif (activeCells(i).cellType == 1)
        h(i) = rectangle('Position', [(activeCells(i).xCoor - cellDiameter/2) (activeCells(i).yCoor - cellDiameter/2) cellDiameter cellDiameter], 'Curvature', [1, 1], 'edgecolor', 'r');
    end
    
end

axis([0 xMax 0 yMax]);

%% PROCESS SIMULATION
%

% start simulation
for i=1:numSteps
    
    i
    
    % retrieve rectangular coordinates and velocities of all active cells at start of current time step
    xPos = zeros(1, length(activeCells));
    yPos = zeros(1, length(activeCells));
    xVel = zeros(1, length(activeCells));
    yVel = zeros(1, length(activeCells));
    
    for j=1:length(activeCells)
        
        xPos(j) = activeCells(j).xCoor;
        yPos(j) = activeCells(j).yCoor;
        xVel(j) = cos(activeCells(j).velocityAng)*activeCells(j).velocityMag;
        yVel(j) = sin(activeCells(j).velocityAng)*activeCells(j).velocityMag;
        
    end
    
    % process cell-wall collisions
    for j=1:length(activeCells)
        
        if (xPos(j) - cellDiameter/2 < 0)
            xVel(j) = abs(xVel(j));
        elseif (xPos(j) + cellDiameter/2 > xMax)
            xVel(j) = -abs(xVel(j));
        end
        
        if (yPos(j) - cellDiameter/2 < 0)
            yVel(j) = abs(yVel(j));
        elseif (yPos(j) + cellDiameter/2 > yMax)
            yVel(j) = -abs(yVel(j));
        end
        
    end
    
    % process cell-cell collisions
    [xDist, yDist, magDist, angDist] = distances(xPos, yPos);
    collisions = (magDist < cellDiameter).*(magDist > 0);
    [a b] = find(triu(collisions));
    numCollisions = length(a);
    
    for j=1:length(activeCells)
        
        % use rotation matrix to calculate elastic collisions, taking into account mass of clusters
        cell1 = a(j);
        cell2 = b(j);
        m1 = 1;
        m2 = 1;
        
        % calculate collision angle and rotate frame of reference
        collAngle = angDist(cell1, cell2);
        v_rot = [cos(-collAngle) -sin(-collAngle); sin(-collAngle) cos(-collAngle)]*[xVel(cell1) xVel(cell2); yVel(cell1) yVel(cell2)];
        v1 = v_rot(1, 1);
        v2 = v_rot(1, 2);
        
        %adjust velocities according to elastic energy transfer
        v_rot(1, 1) = (v1*(m1 - m2) + 2*m2*v2)/(m1 + m2);
        v_rot(1, 2) = (v2*(m2 - m1) + 2*m1*v1)/(m1 + m2);
        
        %rotate back to original frame of reference
        v_new = [cos(collAngle) -sin(collAngle); sin(collAngle) cos(collAngle)]*v_rot;
        xVel(cell1) = v_new(1, 1);
        xVel(cell2) = v_new(1, 2);
        yVel(cell1) = v_new(2, 1);
        yVel(cell2) = v_new(2, 2);
        
    end
    
    % update positions of all cells, update nutrients, then process chemotaxis to update velocities for next time step
    for j=1:length(activeCells)
        
        activeCells(j).xCoor = xPos(j) + xVel(j);
        activeCells(j).yCoor = yPos(j) + yVel(j);

        %activeCells(j).update_nutrients(1); % need to implement code for detecting nutrient concentration from underlying grid
        %activeCells(j).update_velocityAng();
        
    end
    
    % process nutrient consumption/cell divisions/death
    
    % update positions of active cells based on events of most recent time step
    for j=1:length(activeCells)
        set(h(j), 'Position', [(activeCells(j).xCoor - cellDiameter/2) (activeCells(j).yCoor - cellDiameter/2) cellDiameter cellDiameter]);
    end
    drawnow
    
end