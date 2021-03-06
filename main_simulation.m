%% DEFINE PARAMETERS
% 

clear all

% number of timesteps
numSteps = 1000;

% figures on?
fig_on = 0;

% number of yield and growth strategists for seeding the simulation
numYS = 50; % cell type "0"
numGS = 50; % cell type "1"
numCells = numYS + numGS;
initCellTypes = [zeros(1, numYS) ones(1, numGS)];

% maximum dimensions of simulation grid and cell size
xMax = 25;
yMax = 25;
cellDiameter = 1;

% some variables to make the simulation more robust
velocityScalingFactor = 0.25; % controls ratio of cell velocity to cell diameter
collisionBuffer = 1.1; % makes collisions more robust by artifically inflating cell diameter by a small percentage

% define grid size and initialize nutrient grid
totalStartingNutrients = 625*0.01; % this is the total amount distributed uniformly over the entire grid
nutrientReplenishmentRate = 0.025*625; % this is the total amount that is added uniformly over the entire grid
nutrientGridSize = 0.5;
nutrientGrid = zeros(xMax/nutrientGridSize, yMax/nutrientGridSize);
numGridCells = size(nutrientGrid, 1)*size(nutrientGrid, 2);
nutrientGrid(:, :) = totalStartingNutrients/numGridCells/nutrientGridSize;
%nutrientGrid(:, :) = totalStartingNutrients/numGridCells;
nutrientConsumptionRate = 0.1; % should be the same as defined in object

% to implement nutrient gradient at start of simulation
%for j=1:size(nutrientGrid, 1)
%    nutrientGrid(j,:) = j;
%end

% code to report populations of two species over time
popYS = [];
popGS = [];

% tracks newest ID of cells to be created (first "numCells" cells already initialized at beginning of simulation)
nextCellID = numCells + 1;

% seed random number generator
%seed(311);

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

if fig_on == 1
    % draw grid with all active cells (YS strategists blue, GS strategists red)
    fig1 = figure();
    h = zeros(1, length(activeCells));

    for i=1:length(activeCells)

        if (activeCells(i).cellType == 0)
            h(i) = rectangle('Position', [(activeCells(i).xCoor - cellDiameter/2) (activeCells(i).yCoor - cellDiameter/2) cellDiameter cellDiameter], 'Curvature', [1, 1], 'edgecolor', 'b');
        elseif (activeCells(i).cellType == 1)
            h(i) = rectangle('Position', [(activeCells(i).xCoor - cellDiameter/2) (activeCells(i).yCoor - cellDiameter/2) cellDiameter cellDiameter], 'Curvature', [1, 1], 'edgecolor', 'r');
        end

    end

    axis([0 xMax 0 yMax]);
    title('simulation grid')

    % plot nutrient grid as well
    fig2 = figure();
    surf = surface(rot90(nutrientGrid)); % we flip the nutrient grid since the origin is at the top left versus the bottom left for the simulation grid
    colorbar
    axis([0 26/nutrientGridSize 0 26/nutrientGridSize])
    title('nutrient grid')
    hold on
end
%% PROCESS SIMULATION
%

% start simulation
for i=1:numSteps

    % retrieve rectangular coordinates and velocities of all active cells at start of current time step
    xPos = zeros(1, length(activeCells));
    yPos = zeros(1, length(activeCells));
    xVel = zeros(1, length(activeCells));
    yVel = zeros(1, length(activeCells));
    
    for j=1:length(activeCells)
        
        xPos(j) = activeCells(j).xCoor;
        yPos(j) = activeCells(j).yCoor;
        xVel(j) = cos(activeCells(j).velocityAng)*activeCells(j).velocityMag*velocityScalingFactor;
        yVel(j) = sin(activeCells(j).velocityAng)*activeCells(j).velocityMag*velocityScalingFactor;
        
    end
    
    % process cell-wall collisions
    xWallBounce = (xPos <= cellDiameter/2) + (xPos >= (xMax - cellDiameter/2));
    yWallBounce = (yPos <= cellDiameter/2) + (yPos >= (yMax - cellDiameter/2));
    
    if any(xWallBounce)
        xVel = xVel - 2*xVel.*xWallBounce;
    end
    
    if any(yWallBounce)
        yVel = yVel - 2*yVel.*yWallBounce;
    end
    
    % ensure that cells don't get stuck to walls
    xWallStuck = (xPos <= (cellDiameter/2 - abs(xVel))) + (xPos >= (xMax - cellDiameter/2 + abs(xVel)));
    yWallStuck = (yPos <= (cellDiameter/2 - abs(yVel))) + (yPos >= (yMax - cellDiameter/2 + abs(yVel)));
    xPos = xPos.*(1 - xWallStuck) + ceil(xPos).*(xWallStuck);
    yPos = yPos.*(1 - yWallStuck) + ceil(yPos).*(yWallStuck);
    
    % process cell-cell collisions
    [xDist, yDist, magDist, angDist] = distances(xPos, yPos);
    collisions = (magDist <= cellDiameter*collisionBuffer).*(magDist > 0);
    [a b] = find(triu(collisions));
    numCollisions = length(a);
    %numCollisions = 0;
    
    for j=1:numCollisions
        
        % use rotation matrix to calculate elastic collisions, taking into account mass of clusters
        cell1 = a(j);
        cell2 = b(j);
        
        % calculate collision angle and rotate frame of reference
        collAngle = angDist(cell1, cell2);
        v_rot = [cos(-collAngle) -sin(-collAngle); sin(-collAngle) cos(-collAngle)]*[xVel(cell1) xVel(cell2); yVel(cell1) yVel(cell2)];
        v1 = v_rot(1, 1);
        v2 = v_rot(1, 2);
        
        %adjust velocities according to elastic energy transfer
        v_rot(1, 1) = v2;
        v_rot(1, 2) = v1;
        
        %rotate back to original frame of reference
        v_new = [cos(collAngle) -sin(collAngle); sin(collAngle) cos(collAngle)]*v_rot;
        xVel(cell1) = v_new(1, 1);
        xVel(cell2) = v_new(1, 2);
        yVel(cell1) = v_new(2, 1);
        yVel(cell2) = v_new(2, 2);
        
    end

    % A cell array to keep track of which cells are closest to which grid positions in the nutrient array.
    closestCells = cell(xMax/nutrientGridSize, yMax/nutrientGridSize);
    
    % update positions of all cells, update nutrients, then process chemotaxis to update velocities for next time step
    for j=1:length(activeCells)
        
        % update positions based on collision adjusted velocities
        activeCells(j).xCoor = xPos(j) + xVel(j);
        activeCells(j).yCoor = yPos(j) + yVel(j);
        
        % code that updates velocities (circumvents random walk)
        %activeCells(j).velocityMag = sqrt(xVel(j)^2 + yVel(j)^2)/velocityScalingFactor;
        %activeCells(j).velocityAng = atan2(yVel(j), xVel(j));

        % find nearest nutrient grid space (by center) to the current cell
        nearestNutrientPosX = ceil(activeCells(j).xCoor/nutrientGridSize);
        nearestNutrientPosY = ceil(activeCells(j).yCoor/nutrientGridSize);
        if (nearestNutrientPosX > size(nutrientGrid, 1))
            nearestNutrientPosX = 25/nutrientGridSize;
        elseif (nearestNutrientPosX < 1)
            nearestNutrientPosX = 1;
        end
        if (nearestNutrientPosY > size(nutrientGrid, 2))
            nearestNutrientPosY = 25/nutrientGridSize;
        elseif (nearestNutrientPosY < 1)
            nearestNutrientPosY = 1;
        end
        closestCells{nearestNutrientPosX, nearestNutrientPosY} = [closestCells{nearestNutrientPosX, nearestNutrientPosY} j];
        
    end
    
    % now that we've identified which cells correspond to which nutrient grid spaces, process nutrient consumption accordingly (update nutrient grid as well as cell objects)
    for j=1:size(closestCells, 1)
        for k=1:size(closestCells, 2)
            
            numCorrespondingCells = length(closestCells{j, k});
            
            % if at least one cell corresponds to the current nutrient grid space, then process nutrient consumption for that grid space
            if (numCorrespondingCells > 0)
                
                neededNutrients = nutrientConsumptionRate*numCorrespondingCells;
                
                % if there aren't enough nutrients to be fully consumed by all the cells, then consume whatever's available and update accordingly
                if (neededNutrients < nutrientGrid(j, k))
                    
                    for index=1:numCorrespondingCells
                        
                        activeCells(closestCells{j, k}(index)) = activeCells(closestCells{j, k}(index)).update_nutrients(nutrientGrid(j, k), nutrientConsumptionRate);
                        activeCells(closestCells{j, k}(index)) = activeCells(closestCells{j, k}(index)).update_velocityAng();
                        %[nutrientGrid(j,k) nutrientConsumptionRate]
                    end
                    
                    nutrientGrid(j, k) = nutrientGrid(j, k) - neededNutrients;
                    
                else
                    
                    for index=1:numCorrespondingCells
                        activeCells(closestCells{j, k}(index)) = activeCells(closestCells{j, k}(index)).update_nutrients(nutrientGrid(j, k), nutrientConsumptionRate*nutrientGrid(j, k)/neededNutrients);
                        activeCells(closestCells{j, k}(index)) = activeCells(closestCells{j, k}(index)).update_velocityAng();
                    end
                    
                    nutrientGrid(j, k) = 0;
                    
                end
                
            end
            
        end
    end
    
    % replenish nutrients uniformly across entire grid
    nutrientGrid(:, :) = nutrientGrid(:, :) + nutrientReplenishmentRate/numGridCells;
    
    % process cell division/death
    cellsToRemove = [];
    newID = length(activeCells)+1;
    
    for j=1:length(activeCells)        
        
        % ***nutrient "concentration" must inputed depending on the closest
        % nutrient source. We might need to change the update_nutrient
        % function a little bit depending on how the nutrient detection is
        % implemented.***
        %activeCells(j) = update_nutrients(activeCells(j), concentration);
        % MS Note: nutrient update routine was implemented in above loop.
        % I slightly changed the update_nutrients() function to not only
        % accepted the currently sensed concentration, but also how much
        % the cell object consumed in that given time step.  Previously it
        % seemed to be consuming a proportional amount of the total
        % nutrient concentration at that grid space, so for example if the
        % nutrient consumption rate was 1/5 then if more than 5 cells were
        % to try to consume from that space, there would be a shortage.
        % Not sure if that was the intended behavior but I changed it so
        % that only a fixed amount of nutrients is consumed at each step.
        
        [activeCells(j), boolDivision] = activeCells(j).check_division();
        [activeCells(j), boolDeath] = activeCells(j).check_death();
        
        % diagnostic code
        %[activeCells(j).prevSensedConc activeCells(j).currSensedConc]
        %[activeCells(j).nutrientsConsumed activeCells(j).timeStepsNutComplete activeCells(j).age boolDivision boolDeath]
        
        if boolDeath == 1
            
            % remove cell from activeCells list
            cellsToRemove = [cellsToRemove j];
            
        elseif boolDivision == 1
            
            % randomly place cell close to the original cell
            activeCells = [activeCells cell_obj(nextCellID, activeCells(j).cellType, activeCells(j).xCoor, activeCells(j).yCoor)];
            nextCellID = nextCellID + 1;
            
            if fig_on == 1
                % update plot handles based on cell type
                if (activeCells(j).cellType == 0)
                    new_h = rectangle('Position', [(activeCells(end).xCoor - cellDiameter/2) (activeCells(end).yCoor - cellDiameter/2) cellDiameter cellDiameter], 'Curvature', [1, 1], 'edgecolor', 'b');
                    h(end+1) = new_h;
                elseif (activeCells(j).cellType == 1)
                    new_h = rectangle('Position', [(activeCells(end).xCoor - cellDiameter/2) (activeCells(end).yCoor - cellDiameter/2) cellDiameter cellDiameter], 'Curvature', [1, 1], 'edgecolor', 'r');
                    h(end+1) = new_h;
                end
            end
            
        end
    end
    
    % delete any cells marked for deletion
    activeCells(cellsToRemove) = [];
    if fig_on == 1
        % update positions of active cells based on events of most recent time step
        figure(fig1)
        clf
        axis([0 xMax 0 yMax]);
        title(i)
    end
    
    %for j=1:length(activeCells)
    %    set(h(j), 'Position', [(activeCells(j).xCoor - cellDiameter/2) (activeCells(j).yCoor - cellDiameter/2) cellDiameter cellDiameter]);
    %end
    
    % need code to count populations too
    countYS = 0;
    countGS = 0;
    
    if fig_on == 1
        h = zeros(1, length(activeCells));
    end
    
    for j=1:length(activeCells)

        if (activeCells(j).cellType == 0)
            if fig_on == 1
                h(j) = rectangle('Position', [(activeCells(j).xCoor - cellDiameter/2) (activeCells(j).yCoor - cellDiameter/2) cellDiameter cellDiameter], 'Curvature', [1, 1], 'edgecolor', 'b');
            end
            countYS = countYS + 1;
        elseif (activeCells(j).cellType == 1)
            if fig_on == 1
                h(j) = rectangle('Position', [(activeCells(j).xCoor - cellDiameter/2) (activeCells(j).yCoor - cellDiameter/2) cellDiameter cellDiameter], 'Curvature', [1, 1], 'edgecolor', 'r');
            end
            countGS = countGS + 1;
        end

    end
    
    popYS = [popYS countYS];
    popGS = [popGS countGS];
    
    if fig_on == 1
        if mod(i,10) == 0
            % update cell plot
            figure(fig1)
            drawnow
            title(i)
            print(['cell' num2str(i)],'-dpng')

            % update plot of nutrient grid
            figure(fig2)
            delete(surf);
            surf = surface(flip(rot90(nutrientGrid, 1)));
            print(['grid' num2str(i)],'-dpng')
        end
    end
    
    if mod(i,10) == 0
        i
    end
end
%%
figure()
t = 1:numSteps;
plot(t, popYS, 'color', [51 51 255]./255, 'linewidth', 3)
hold all
plot(t, popGS, 'color', [255 51 51]./255, 'linewidth', 3)
plot(t, popYS + popGS, 'color', [102 204 0]./255, 'linewidth', 3)

xlabel('Time Steps')
ylabel('Cell Number')

% title('')
legend('YS (blue)', 'GS (red)', 'YS + GS', 'location', 'best')

set(gca, 'fontsize', 16)
set(gcf, 'color', 'w')

grid on

