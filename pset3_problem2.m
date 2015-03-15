% code for question 2, see also function connected_from_bound.m
close all;

% set parameters
numStep=1000000;
numBlue=10;
numRed=10;
numCells=numBlue+numRed;
simSize=12;
vMax=0.01;
bind_dist=0.1;

% % 3 different k_on and k_off values now
% k_on_bb=0.0000001;
% k_on_br=0.001;
% k_on_rr=0.3;
% k_off_bb=0.3;
% k_off_br=0.000001;
% k_off_rr=0.000001;

% % k_on and off now are matrices instead of single elements. Top left corner
% % correspond to blue-blue interactions, bottom right corresponds to red-red,
% % off-axis corresponds to blue-red.
% k_on=[k_on_bb*ones(numBlue) k_on_br*ones(numBlue, numRed); k_on_br*ones(numRed, numBlue) k_on_rr*ones(numRed)];
% k_off=[k_off_bb*ones(numBlue) k_off_br*ones(numBlue, numRed); k_off_br*ones(numRed, numBlue)
% k_off_rr*ones(numRed)];

D=1; %diameter of ball
sim=figure(1);
axis([0 simSize 0 simSize]);
axis square;

% pick random starting coordinates, make sure they don't intersect
conflict=1;

while conflict
    x=(simSize-D)*rand(numCells, 1);
    y=(simSize-D)*rand(numCells, 1);
    [xDist, yDist, magDist, angDist]=distances(x, y);
    conflict=sum(sum((magDist<D).*(magDist>0)));
end

%set initial velocities
vMag=vMax*rand(numCells, 1);
vAng=2*pi*rand(numCells, 1);
vx=cos(vAng).*vMag;
vy=sin(vAng).*vMag;

%create circles, first set are blue, the rest are red
for index=1:numCells
    if index<=numBlue
    h(index)=rectangle('Position',[x(index)-.5*D y(index)-.5*DDD],'Curvature',[1 1],'edgecolor','b');
    else
    h(index)=rectangle('Position',[x(index)-.5*D y(index)-.5*DDD],'Curvature',[1 1],'edgecolor','r');
    end
end

bound=zeros(numCells);
connected=zeros(numCells);

%loop through simulation, add velocity to position
for timestep=1:numStep
    %calculate binding, give close unbound cells a chance to bind, and
    %bound cells a chance to unbind based on k_on, k_off matrices calculated above
    %Note this implements inelastic collisions
    to_bind=((triu(magDist<D+bind_dist).*(magDist>D)-bound)).*(rand(numCells)<k_on);
    un_bind=triu(bound).*(rand(numCells)<k_off);
    bound=triu(bound)+to_bind-un_bind;
    bound=bound+bound';
    
    % calculate connected clusters based on new bound matrix, use to
    % calculate mass
    [connected cluster_list numClusters] = connected_from_bound(bound);
    mass=1+sum(connected)';
    
    % ball wall collisions
    wallbounce_x=(x<=0)+(x>=simSize-D);
    wallbounce_y=(y<=0)+(y>=simSize-D);
    if any(wallbounce_x)
        wallbounce_x=wallbounce_x+sum(connected(:, find(wallbounce_x)), 2);
        vx=vx-2*vx.*wallbounce_x;
    end
    
    if any(wallbounce_y)
        wallbounce_y=wallbounce_y+sum(connected(:, find(wallbounce_y)), 2);
        vy=vy-2*vy.*wallbounce_y;
    end
    
%     % the following code ensures balls don't stick to the wall
%     wallstuck_x=(x<=-abs(vx)) + (x>=simSize-D+abs(vx));
%     wallstuck_y=(y<=-abs(vy)) + (y>=simSize-D+abs(vy));
%     x=x.*(1-wallstuck_x)+round(x).*(wallstuck_x);
%     y=y.*(1-wallstuck_y)+round(y).*(wallstuck_y);
    
    % ball-ball collisions, elastic
    [xDist, yDist, magDist, angDist]=distances(x, y);
    collisions=(magDist<D).*(magDist>0); %collisions if distances are less than diameter
    [a b]=find(triu(collisions)); %extract individual collisions from pairwise collision matrix
    numCollisions=length(a);
    
    for index=1:numCollisions
        % use rotation matrix to calculate elastic collisions, taking into account mass of clusters
        ball1=a(index);
        ball2=b(index);
        m1=mass(ball1);
        m2=mass(ball2);
        
        %calculate collision angle and rotate frame of reference
        collAngle=angDist(ball1, ball2);
        v_rot=[cos(-collAngle) -sin(-collAngle);sin(-collAngle) cos(-collAngle)]*[vx(ball1) vx(ball2); vy(ball1) vy(ball2)];
        v1=v_rot(1, 1);
        v2=v_rot(1, 2);
        
        %adjust velocities according to elastic energy transfer
        v_rot(1, 1)=(v1*(m1-m2)+2*m2*v2)/(m1+m2);
        v_rot(1, 2)=(v2*(m2-m1)+2*m1*v1)/(m1+m2);
        
        %rotate back to original frame of reference
        v_new=[cos(collAngle) -sin(collAngle);sin(collAngle) cos(collAngle)]*v_rot;
        vx(ball1)=v_new(1, 1);
        vx(ball2)=v_new(1, 2);
        vy(ball1)=v_new(2, 1);
        vy(ball2)=v_new(2, 2);
        
        %all balls belonging to cluster of either of these cells need to have velocity adjusted
        connected1=connected(:, ball1);
        connected2=connected(:, ball2);
        vx(find(connected1))=v_new(1, 1);
        vx(find(connected2))=v_new(1, 2);
        vy(find(connected1))=v_new(2, 1);
        vy(find(connected2))=v_new(2, 2);
    end
    
    %for each cluster, set all velocities equal to average velocity
    for index=1:numClusters
        cluster=cluster_list(:, index);
        vx_av=sum(vx.*cluster)/sum(cluster);
        vy_av=sum(vy.*cluster)/sum(cluster);
        vx(find(cluster))=vx_av;
        vy(find(cluster))=vy_av;
    end
    
    %update positions and plot
    x=x+vx;
    y=y+vy;
    
    for index=1:numCells
        set(h(index),'Position',[x(index) y(index) D D]);
    end
    drawnow
end

