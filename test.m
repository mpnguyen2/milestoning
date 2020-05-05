function test
clear;clc;

%% ARGUMENTS SET UP
% Domain setting
num_m = 10; ms_dist = 2; lower = [0 0]; N = 30; vert_dist = 0.2; 
sigma = 1; % constant variance of diffusion process
n = 20; % PDE approach step size
Num_traj = 200000; max_step = 1000; % Milestone setting
% big_num is used to rescale the result (since it is small)
big_num = 10000;
% Animation
file_name = 'sim4.avi';
rescale = 100;
factor = 3;
start_ms = 5; end_ms = 8;

run_sim_on_rect(num_m, ms_dist, lower, N, vert_dist, sigma, ...
    n, @V, @V_arr, @div_V, Num_traj, max_step, ...
    big_num, start_ms, end_ms, file_name, rescale, factor);

end
    
%% Argument functions
% Drift function with a single argument
function v = V(x)
v = zeros(2, 1);
%{
x1 = x(1); x2 = x(2);
% drift fct for Muller-Brown potential
v(1) = -200*exp(- (x1 - 1)^2 - 10*x2^2)*(2*x1 - 2) - 200*x1*exp(- 10*(x2 - 1/2)^2 - x1^2) - ...
     15*exp((7*(x1 + 1)^2)/10 + (7*(x2 - 1)^2)/10 + ((3*x1)/5 + 3/5)*(x2 - 1))*((7*x1)/5 + ...
     (3*x2)/5 + 4/5) - 170*exp((11*x1 + 11/2)*(x2 - 3/2) - ...
     (13*(x2 - 3/2)^2)/2 - (13*(x1 + 1/2)^2)/2)*(13*x1 - 11*x2 + 23);
v(2) = 170*exp((11*x1 + 11/2)*(x2 - 3/2) - (13*(x2 - 3/2)^2)/2 - ...
    (13*(x1 + 1/2)^2)/2)*(11*x1 - 13*x2 + 25) - ...
    4000*x2*exp(- (x1 - 1)^2 - 10*x2^2) - 15*exp((7*(x1 + 1)^2)/10 + ...
    (7*(x2 - 1)^2)/10 + ((3*x1)/5 + 3/5)*(x2 - 1))*((3*x1)/5 + (7*x2)/5 - 4/5) ...
    - 100*exp(- 10*(x2 - 1/2)^2 - x1^2)*(20*x2 - 10);
%}
v(1) = -x(1); v(2) = -x(2);
end
% Drift function with an array/matrix argument
function v_arr = V_arr(X)
v_arr = zeros(2, size(X,2));
for i = 1:size(X,2)
    v_arr(:, i) = V(X(:,i));
end
end

% Divergence of drift function (used later)
function div = div_V(x)
%x1 = x(1); x2 = x(2);

div = 0;
%{
div = 170*exp((11*x1 + 11/2)*(x2 - 3/2) - (13*(x2 - 3/2)^2)/2 - ...
    (13*(x1 + 1/2)^2)/2)*(13*x1 - 11*x2 + 23)^2 - 4420*exp((11*x1 + 11/2)*(x2 - 3/2)...
    - (13*(x2 - 3/2)^2)/2 - (13*(x1 + 1/2)^2)/2) - 4400*exp(- (x1 - 1)^2 - 10*x2^2) - ...
    2200*exp(- 10*(x2 - 1/2)^2 - x1^2) - 15*exp((7*(x1 + 1)^2)/10 + (7*(x2 - 1)^2)/10 + ...
    ((3*x1)/5 + 3/5)*(x2 - 1))*((3*x1)/5 + (7*x2)/5 - 4/5)^2 - 15*exp((7*(x1 + 1)^2)/10 + ...
    (7*(x2 - 1)^2)/10 + ((3*x1)/5 + 3/5)*(x2 - 1))*((7*x1)/5 + (3*x2)/5 + 4/5)^2 - ...
    42*exp((7*(x1 + 1)^2)/10 + (7*(x2 - 1)^2)/10 + ((3*x1)/5 + 3/5)*(x2 - 1)) + ...
    170*exp((11*x1 + 11/2)*(x2 - 3/2) - (13*(x2 - 3/2)^2)/2 - ...
    (13*(x1 + 1/2)^2)/2)*(11*x1 - 13*x2 + 25)^2 + 200*exp(- (x1 - 1)^2 - 10*x2^2)*(2*x1 - 2)^2 + ...
    100*exp(- 10*(x2 - 1/2)^2 - x1^2)*(20*x2 - 10)^2 + 80000*x2^2*exp(- (x1 - 1)^2 - 10*x2^2) + ...
    400*x1^2*exp(- 10*(x2 - 1/2)^2 - x1^2);
%}
end