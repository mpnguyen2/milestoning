%% Function to calculate the coordinate of the hitting point starting
% from an initial pt until it reach the top or bottom boundary or the
% adjacent milestones. The pt is assume to follow the Brownian dynamics:
% dXt = V(Xt)dt + sigma dBt
% This function is different from first version since it uses random
% walk instead of (sum of) normal random vars

% Args:
% init: y-coordinates of the initial point
% x, y: x and y-coordinates of the lowest point of the middle milestone,
% the milestone we are considering
% V: drift fct
% sigma: constant variance describing the Brownian dynamics
% ms_dist: equidistance btw milestones
% N: the number of pts of milestone
% vert_dist: the distance btw pts on a milestone
% num_traj: number of trajectory initiated from [x, init]
% max_step: maximum number of steps allowed to run to avoid too long 
% trajectory

% Do normal Gaussian approach first and then come back to random walk approach
% later
% (n: the integer so that h = n * step_size (step size of the discretized
% random walk)

% Return:
% l_hit: an array where ith element indicate percentage of trajectories 
% hitting ith pt on the left milestone
% r_hit: an array where ith element indicate percentage of trajectories
% hitting ith pt on the right milestone

function [l_hit, r_hit] =shootTraj2(init, x, y, ...
        V, sigma, vert_dist, ms_dist, N, n, num_traj, max_step)

if num_traj == 0
    l_hit = zeros(1,N); r_hit = zeros(1,N);
    return;
end
pts = zeros(2, num_traj);
pts(1,:) = x * ones(1, num_traj);
pts(2,:) = init * ones(1,num_traj);
left = x - ms_dist; right = x + ms_dist; 
top = y + (N-1) * vert_dist; bottom = y;

i = 1;
Dt = 0.005;

%{
min((vert_dist/sigma)^2, ...
    vert_dist/(max(abs(v(1)), abs(v(2))) + 1));
%}
ind = 1:num_traj;
h = vert_dist/n;
moves = [[h 0],[-h 0], [0, h], [0,-h]];
while ~isempty(ind)
    % for other index, we can start processing left or right
    % for l_hit or r_hit
    ind2 = ind(find(pts(1, ind) < left)); 
    ind3 = ind(find(pts(1, ind) > right)); 
    ind4 = ind(find(pts(2, ind) > top | pts(2, ind) < bottom)); 
    % Now update the indices of trajectories haven't hit milestones yet
    ind = ind(find(pts(1, ind) >= left & pts(1, ind) <= right ...
        & pts(2, ind) >= bottom & pts(2, ind) <= top));
    % Assign new value to done trajectories
    pts(1, ind2) = 1*ones(1, length(ind2)); %left
    pts(1, ind3) = -1*ones(1, length(ind3)); % right
    pts(1, ind4) = 0*ones(1, length(ind4)); % top or bottom
   
    % Continue other trajectory based on the SDE given above
    V_vec = V(pts(:,ind));
    pts(:,ind) = pts(:,ind) + V_vec(:,:) * Dt + ...
        sigma * moves(ceil(4*rand(1, length(ind))));
    %disp(pts(:,2));
    % Update number of steps done
    i = i+1;
end

% Now use trajectories hitting info to find l_hit and r_hit
l_hit = zeros(1, N); r_hit = zeros(1,N);
ind = find(pts(1,:) == 1 | pts(1,:) == -1);
pts(2,ind) = max(min(round((pts(2,ind) - y)/vert_dist), N-1),0) + 1;

for t = 1:length(ind)
    if pts(1,ind(t)) == 1
        l_hit(pts(2,ind(t))) = l_hit(pts(2,ind(t)))+1;
    else
        r_hit(pts(2,ind(t))) = r_hit(pts(2,ind(t)))+1;
    end
end

% r_hit = r_hit/num_traj;
% l_hit = l_hit/num_traj
end