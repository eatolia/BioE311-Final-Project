classdef cell_obj
    properties
        ID
        cellType
        
        velocityMag
        velocityAng
        growthRate
        growthYield
        deathRate
        nutrientConsumRate
        
        xCoor 
        yCoor
        nutrientsConsumed
        prevSensedConc
        currSensedConc
        
        timeStepsNutComplete % timesteps since the amount of nutrients needed to divde has been reached
    end
    
    methods
        function obj = cell_obj(ID, cellType, xCoor, yCoor)
            obj.ID = ID;
            obj.cellType = cellType;
            
            if cellType == 0
                obj.growthRate = 1/5; % divides every 5 time steps 
                obj.growthYield = 1/10; % can divide when 10 nutrient particles have been consumed
                obj.deathRate = 1/100; % every 100 time steps without suficient nutrients
            elseif cellType == 1
                obj.growthRate = 1/2; % divides every 5 time steps 
                obj.growthYield = 1/5; % can divide when 10 nutrient particles have been consumed
                obj.deathRate = 1/100; % every 100 time steps without suficient nutrients            
            end   
            
            obj.velocityMag = 1;
            obj.velocityAng = 2*pi*rand(1, 1);
            
%             vx = cos(obj.veocityAng).*obj.velocityMag;
%             vy = sin(obj.veocityAng).*obj.velocityMag;
            
            obj.nutrientConsumRate = 1/5; % per timestep
            
            obj.xCoor = xCoor;
            obj.yCoor = yCoor;
            obj.nutrientsConsumed = 0;
            obj.prevSensedConc = 0;
            obj.currSensedConc = 0;
            
            obj.timeStepsNutComplete = 0;
        end
        
        % updating the current concetration of nutrients and also keeping
        % track of time till next cell division
        function obj = update_nutrients(obj, concentration)
            obj.nutrientsConsumed = obj.nutrientsConsumed + concentration*obj.nutrientConsumRate;
            obj.prevSensedConc = obj.currSensedConc;
            obj.currSensedConc = concentration;
            
            if obj.nutrientsConsumed >= 1/obj.growthYield && obj.timeStepsNutComplete >= 1/obj.growthRate
                obj = obj.divide_cell();
            elseif obj.nutrientsConsumed >= 1/obj.growthYield && obj.timeStepsNutComplete < 1/obj.growthRate
                obj.timeStepsNutComplete = obj.timeStepsNutComplete + 1;
            end
                
        end
        
        % implementation of biased random walk
        function obj = update_velocityAng(obj)
            if obj.prevSensedConc >= obj.currSensedConc
                if rand(1,1) < .2
                    obj.velocityAng = 2*pi*rand(1, 1);
                end
            elseif obj.prevSensedConc < obj.currSensedConc
                obj.velocityAng = 2*pi*rand(1, 1);
            end
        end
        
        % implementation of cell division
        function obj = divide_cell(obj)
            obj.timeStepsNutComplete = 0;
            obj.nutrientsConsumed = 0;
        end
        
        
    end
        
        
end

