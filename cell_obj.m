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
        
        age
        timeStepsNutComplete % timesteps since the amount of nutrients needed to divde has been reached

    end
    
    methods
        function obj = cell_obj(ID, cellType, xCoor, yCoor, varargin)
            obj.ID = ID;
            obj.cellType = cellType;
            obj.deathRate = 1/200; % every 100 time steps without suficient nutrients            
            
            if cellType == 0 % yield strategists
                obj.growthRate = 0.1; % can divide every 10 time steps 
                obj.growthYield = 0.16; % can divide when 5 nutrient particles have been consumed
            elseif cellType == 1 % growth strategists
                obj.growthRate = 1; % can divide every 5 time steps 
                obj.growthYield = 0.1; % can divide when 10 nutrient particles have been consumed
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
            
            obj.age = 0;
        end
        
        % updating the current concetration of nutrients and also keeping
        % track of time till next cell division
        function obj = update_nutrients(obj, concentration, consumed)
            %obj.nutrientsConsumed = obj.nutrientsConsumed + concentration*obj.nutrientConsumRate;
            obj.nutrientsConsumed = obj.nutrientsConsumed + consumed;
            obj.prevSensedConc = obj.currSensedConc;
            obj.currSensedConc = concentration;
            
        end
        
        % implementation of biased random walk
        function obj = update_velocityAng(obj)
            if obj.prevSensedConc < obj.currSensedConc
                if rand(1,1) < .2
                    obj.velocityAng = 2*pi*rand(1, 1);
                end
            elseif obj.prevSensedConc >= obj.currSensedConc
                obj.velocityAng = 2*pi*rand(1, 1);
            end
        end
        
        % implementation of cell division
        function obj = divide_cell(obj)
            obj.timeStepsNutComplete = 0;
            obj.nutrientsConsumed = 0;
            obj.age = 0;
        end
        
        function [obj, boolDivision] = check_division(obj)
            if obj.nutrientsConsumed >= 1/obj.growthYield && obj.age >= 1/obj.growthRate
                boolDivision = 1;
                obj = obj.divide_cell();
            elseif obj.nutrientsConsumed >= 1/obj.growthYield && obj.age < 1/obj.growthRate
                boolDivision = 0;
                obj.timeStepsNutComplete = obj.timeStepsNutComplete + 1;
                obj.age = obj.age + 1;
            else
                boolDivision = 0;
                obj.age = obj.age + 1;
            end
        end
        
        function [obj, boolDeath] = check_death(obj)
            if obj.age > 1/obj.deathRate 
                boolDeath = 1;
            else
                boolDeath = 0;
            end
        end  
    end
end

