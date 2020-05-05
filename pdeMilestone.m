%% Function use pde approach to output data over time on milestones
% Args:
% There're num_m milestones which are parallel edges of a rectangle
% the distance btw milestones is ms_dist
% the lower left corner have coordinate 'lower'
% N is the number of point per milestones and distance btw pt on milestone
% sigma: variance (standard deviation) for the process
% init_data: initial data on milestone
% big_num: A big number to help rescale milestone

% Return:
% outputs for simulation and result to display 
% (zoom out by a factor of big_num)

function [total_dat, b1_data, b2_data] = pdeMilestone(n, num_m, ms_dist, N, ...
    vert_dist, lower, init_data, sigma, V, V_arr, div_V, big_num)

total_dat(1,:) = init_data(:);
data = init_data;

% Specify the parallel pool
p = parpool(num_m);
        
tic
%% Iterative algorithm
std_err = 1e-10; err = 1;
cpy_data = zeros(num_m,N);
b1_data = zeros(1,N); b2_data = zeros(1,N);
while err > std_err
    cpy_data(:) = data(:);
    cpy_data(1,:) = zeros(1,N);
    cpy_data(num_m,:) = zeros(1,N);
    parfor i = 1:num_m
        if i ~= 1 && i ~= num_m
            data(i,:) = solvePDE(lower + [(i-1)* ms_dist, 0], ...
                vert_dist, ms_dist, N, n, V, div_V, sigma,...
                cpy_data(i-1,:), cpy_data(i+1,:));
        elseif i == 1
           added = zeros(1,N);
           num_traj = 10000;
           for j = 1:N
               [l, ~] = shootTraj(lower(2)+(j-1)*vert_dist, ...
                   lower(1)+(i-1)*ms_dist, lower(2), V_arr, sigma, ...
                    vert_dist, ms_dist, N, num_traj, 400);
                added = added + cpy_data(i+1,j)*l/num_traj;
           end
           b1_data = b1_data + big_num * added * err;
        elseif i == num_m
           added = zeros(1,N);
           num_traj = 10000;
           for j = 1:N
               [~, r] = shootTraj(lower(2)+(j-1)*vert_dist, ...
                   lower(1)+(i-1)*ms_dist, lower(2), V_arr, sigma, ...
                   vert_dist, ms_dist, N, num_traj, 200);
                added = added + (cpy_data(i-1,j)/num_traj)*r;
           end
           b2_data = b2_data + big_num * added * err;
        end
    end
    temp = [b1_data/big_num; data(2:num_m-1,:) * err; ...
        b2_data/big_num];
    total_dat(end+1,:) = temp(:);
    my_sum = sum(data(:));
    data = data/my_sum;
    err = err * my_sum;  
    %fprintf("Error is: %3.6f \n", err);
end

fprintf("The time this PDE approach takes is: %3.6f seconds\n", toc);

delete(p);
end