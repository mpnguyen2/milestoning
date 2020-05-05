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

function [total_dat, b1_data, b2_data] = pdeMilestoneV2(n, num_m, ms_dist, N, ...
    vert_dist, lower, init_data, sigma, V, V_arr, div_V, big_num)

total_dat(1,:) = init_data(:);
data = init_data;

% Specify the parallel pool
p = parpool(num_m-2);
        
tic
%% Iterative algorithm
std_err = 1e-10; err = 1;
% r_data, l_data is induced density to the right/left
l_data = zeros(num_m-2, N); r_data = zeros(num_m-2, N);
b1_data = zeros(1,N); b2_data = zeros(1,N);
while err > std_err
    parfor i = 2:num_m-1
        % Solve Poisson eq on each of num_m - 2 inner rectangles
        [l_data(i-1,:), r_data(i-1,:)] = solvePoisson(lower + [(i-1)* ms_dist, 0], ...
            vert_dist, ms_dist, N, n, V, sigma,...
            data(i,:));
    end
    %disp(size(l_data(1,:) ~= 0));
    data = [l_data; zeros(2,N)] + [zeros(2,N); r_data];
    b1_data = b1_data + big_num * err * l_data(1,:);
    b2_data = b2_data + big_num * err * r_data(end,:);
    temp = [b1_data/big_num; data(2:num_m-1,:) * err; ...
        b2_data/big_num];
    total_dat(end+1,:) = temp(:);
    data(1,:) = zeros(1,N); data(end,:) = zeros(1,N);
    my_sum = sum(data(:));
    data = data/my_sum;
    err = err * my_sum;  
    %fprintf("Error is: %3.6f \n", err);
end

fprintf("The time this PDE approach takes is: %3.6f seconds\n", toc);

delete(p);
end