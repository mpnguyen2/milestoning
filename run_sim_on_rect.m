%% Function for running simulation to compare PDE vs milestone approach
% on a 2D rectangle
% Args
% num_m: Number of milestones
% ms_dist: distance btw milestones
% lower: the coordinate of the bottom left corner of the domain
% N: Number of point per milestones 
% vert_dist: distance btw pt on milestone
% sigma: standard deviation or variance of the process

% ADDITIONAL ARGUMENTS:

% PDE approach:
% n: step size multiple used in finite difference for PDE approach
% V: drift function for diffusion (on a single variable)
% div_V: divergence of the drift function

% Milestone approach:
% V_arr: drift functions but applies on an entire array and 
% also return an array of results
% Num_traj: number of traj used 
% max_step: maximum number of steps a single traj can take to avoid 
% too long trajectory
% big_num: is used to rescale the result since in this case, 
% the result is too small

% Initial data: our initial data will be uniformly concentrated on
% (start_ms) milestone to (end_ms) milestone

% Animation:
% file_name: Name of the video file to save the simulation
% rescale: Since the density is small and is hard to see in the image, the
% density data must be made bigger by rescale time to be displayed properly
% in the video
% factor: The number of image frames may be too small with rapid changes
% The factor is used to increase this number to provide a smooth
% interpolation and transition

function run_sim_on_rect(num_m, ms_dist, lower, N, vert_dist, sigma, ...
    n, V, V_arr, div_V, Num_traj, max_step, ...
    big_num, start_ms, end_ms, file_name, rescale, factor)

%% DATA SET UP
% initial data uniformly concentrated on milestone start_ms to end_ms
init_data = zeros(num_m, N);
temp = end_ms - start_ms + 1;
init_data(start_ms:end_ms,:) = 1/(N*temp) * ones(temp,N);
% Print info
fprintf("SETUP:\n");
fprintf("Our domain is a rectangle whose lower left corner ");
fprintf("has coordinate (%.2f, %.2f)\n", lower(1), lower(2));
fprintf("There are %d milestones, each has %d points\n", num_m, N);
fprintf("The distance between milestone is %1.2f, ", ms_dist)
fprintf("and the distance between pts is %1.2f \n", vert_dist);
fprintf("We use the constant variance %2.2f \n", sigma);
fprintf("Initial density is uniformly concentrated on")
fprintf(" milestones %d to %d \n", start_ms, end_ms);
fprintf("The number of trajectories used for milestone method is %d \n", Num_traj);
fprintf("Each trajectory can only be allowed to run at most %d steps \n", max_step);
fprintf("In the PDE approach, each PDE is solve on the domain with step");
fprintf(" size equals 1/%d the distance btw pts on milestone\n\n", n);

%% Milestone approach
[data_ms, b1_data, b2_data] = milestone(Num_traj, max_step, num_m, ...
    ms_dist, N, vert_dist, lower, init_data, sigma, V_arr, big_num);
fprintf("\nMilestone methods results:\n");
fprintf("Density (after scale %d times bigger) on the left boundary is: \n", big_num);
disp(b1_data);
fprintf("Density (after scale %d times bigger) on the right boundary is: \n", big_num);
disp(b2_data);

%% PDE approach
[data_de, b1_de, b2_de] = pdeMilestoneV2(n, num_m, ...
    ms_dist, N, vert_dist, lower, init_data, sigma, V, V_arr, div_V, big_num);
fprintf("\nPDE methods results:\n");
fprintf("Density (after scale 10000 times bigger) on the left boundary is: \n");
disp(b1_de);
fprintf("Density (after scale 10000 times bigger) on the right boundary is: \n");
disp(b2_de);

%% Drawing the animation
fprintf("\nFor the animation we rescale the density by %d times \n", rescale);
fprintf("\nThe interpolation factor in the animation is %d \n", factor);
draw_anim(file_name, data_de, data_ms, N, num_m, rescale, factor);

end