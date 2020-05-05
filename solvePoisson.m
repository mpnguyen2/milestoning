%% Function to solve Poisson equation to find induced density 
% by finite difference method on a small rectangle in 2D
%
% We will solve the Poisson equation with differential operator being the of the 
% infinitesimal generator of the actual process:
% L^T(p) = -div(V(x)p(x)) + sigma^2/2 * \Delta(p(x))
% Solve PDE L^T u = f with f is the initial density to find potential u
% Then take partial derivative wrt normal of u to find left/right boundary
% induced density

% Args:
% t_xy: (2d) coordinates of the lowest point on the target milestone
% ms_dist: the equal distance between 2 (1D) milestones
% vert_dist: distance btw pts on one milestone
% N: number of pts on a milestone
% n: an integer so that vert_dist = n * grid_size (the grid size used to 
% discretize and solve for the PDE)

% Note that top and bottom boundary values are zero

% Return:
% l_b: the left boundary density
% r_b: the right boundary density induced by the electric potential
% that is the solution of the possion equation
function [lb, rb] = solvePoisson(t_xy, vert_dist, ...
    ms_dist, N, n, V, sigma, middle)

% Build coefficient matrix for the discretization

% grid_size
h = vert_dist/n; N_1 = n*(N-1)+1;
% M = number of pts on horizontal direction in our small domain
M = 2*round(ms_dist/h) + 1;
% Some more extra helper vars
mat_is = M*N_1*5 + 4 - 8*(M+N_1);
i = zeros(1, mat_is); j = zeros(1, mat_is); val = zeros(1, mat_is);
cur_ind = 1;

% Find the coefficient matrix corresponding to forward operator
for t = 1:(M*N_1)
    if mod(t-1,M) == 0 || mod(t-1,M) == M-1 || ...
            floor((t-1)/M) == 0 || floor((t-1)/M) == N_1-1
        j(cur_ind) = t; i(cur_ind) = t; val(cur_ind) = 1;
        cur_ind = cur_ind + 1;
        continue;
    end
    x_t = t_xy + [(mod(t-1, M)-(M-1)/2)*h, (t-1)/M * h];
    v = V(x_t);
    % Grid size is h = dx and we need to take the transpose of the
    % infinitesimal generator
    %{
    if v(1) + v(2) ~= 0
        temp = 2 * sigma^2/((v(1)+v(2))^2) + h/(v(1)+v(2));
        if temp > 0
            temp = sqrt(temp);
        else
            temp = 0;
        end
        dt = (temp - sqrt(2)*sigma/(v(1)+v(2)))^2;
        
    else
        p1 = 1/4; p2 = 1/4; p3 = 1/4; p4 = 1/4;
    end
    %}
    dt = h^2/(2 * sigma^2);
    p1 = v(1)*dt/(2*h) + 1/4;
    p2 = -v(1)*dt/(2*h) + 1/4;
    p3 = v(2)*dt/(2*h) + 1/4;
    p4 = 1 - p1 - p2 - p3;
    ind = cur_ind:cur_ind+4;
    % the coef matrix needs to be transpose of prob transition matrix
    j(ind) = [t, t, t, t, t];
    i(ind) = [t, t-1, t+1, t-M, t+M];
    val(ind) = [1, -p2, -p1, -p4, -p3];
    cur_ind = cur_ind + 5;
end
Coef = sparse(i,j,val,M*N_1,M*N_1); 
% Build data b based on the middle milestone.
% Lu = f and u = 0 on the boundary
b = zeros(M*N_1,1);
ind = (0:(N_1-1)) * M + (M+1)/2;
b(ind) = interp1(1:N, middle, 1:1/n:N);
% Solve linear equation to get sol
sol = Coef\b;
% disp(sol(1:5));
% Calculate the left and right value by taking the directional derivative
% wrt normal direction of sol
% rb and lb has size (1, N)
lb = [mean(reshape(sol(1 + M * (0:(N_1-2))), n, N-1)), ...
   mean(sol(1+M*((N_1-n):(N_1-1))))];
rb = [mean(reshape(sol(M + M * (0:(N_1-2))), n, N-1)), ...
   mean(sol(M*((N_1-n+1):N_1)))];
end