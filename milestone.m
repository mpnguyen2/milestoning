%% Function use milestone approach to output data over time on milestones
% Args:
% Num_traj: Use a total of Num_traj trajectories for simulation at each
% iteration
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

function [total_dat, b1_data, b2_data] = milestone(Num_traj, max_step, num_m, ...
    ms_dist, N, vert_dist, lower, init_data, sigma, V_arr, big_num)

% Matrix A to calculate w in the algorithm
% A = zeros(num_m, 2);
data = init_data;
total_dat(1,:) = init_data(:);

% Specify the parallel pool
p = parpool(num_m);
     
tic
%% Iterative algorithm
std_err = 1e-10; err = 1;
b1_data = zeros(1,N); b2_data = zeros(1,N);
while err > std_err
    % Simulate to get transition prob matrix between milestones
    % transition matrix where each row contain two density distribution on 
    % 2 neighboring milestones: 1 to N left ms and N+1 to 2N right ms
    tran_mat = zeros(num_m, 2*N);
    s = zeros(1, num_m);
    parfor i = 2:num_m-1      
        l = zeros(1,N); r = zeros(1,N);
        s(i) = sum(data(i,:));
        % Shoot trajectory from 1 milestone to get update transition matrix
        for j = 1:N
            if s(i) ~= 0
                temp_num_traj = round(Num_traj * data(i,j));
            else
                temp_num_traj = 0;
            end
            [l_j, r_j] = shootTraj(lower(2)+(j-1)*vert_dist, lower(1) + ...
                ms_dist*(i-1), lower(2), V_arr, sigma, vert_dist, ms_dist, N, ...
                temp_num_traj, max_step);
            l = l+l_j; r = r+ r_j;
        end
        l = l/Num_traj; r = r/Num_traj;
        
        tran_mat(i,:) = [l r]; % 1 by 2N vector 
    end
    % Update new distribution (only on milestones)
    new_dat = zeros(num_m, N);
    parfor i = 1:num_m
        if i == 1
            b1_data = b1_data + big_num * tran_mat(i+1,1:N) * err;
        elseif i == num_m
            b2_data = b2_data + big_num * tran_mat(i-1, (N+1):2*N) * err;
        else
            new_dat(i,:) = tran_mat(i-1, (N+1):2*N) + ...
                tran_mat(i+1, 1:N);
        end
    end
    % Normalize new distribution and find next error
    my_sum = sum(s);
    temp = [b1_data/big_num; err*data(2:num_m-1,:); b2_data/big_num];
    total_dat(end+1,:) = temp(:);
    err = err * my_sum;
    if my_sum ~= 0
        data = new_dat/my_sum;
    end
    % DEBUG
    %disp(data * 10000);
    %fprintf("The error is %3.6f \n", err);
end
fprintf("The time this milestone approach takes is: %3.6f seconds\n", toc);

delete(p);
end