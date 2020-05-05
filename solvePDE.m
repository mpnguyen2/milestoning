%% Function to solve PDE equation by finite difference method 
% on a small rectangle in 2D
%
% We will solve the PDE with differential operator being the of the 
% infinitesimal generator of the actual process:
% L^T(p) = -div(V(x)p(x)) + sigma^2/2 * \Delta(p(x))

% Args:
% t_xy: (2d) coordinates of the lowest point on the target milestone
% ms_dist: the equal distance between 2 (1D) milestones
% vert_dist: distance btw pts on one milestone
% N: number of pts on a milestone
% n: an integer so that vert_dist = n * grid_size (the grid size used to 
% discretize and solve for the PDE)
% l_b: the left boundary condition
% r_b: the right boundary condition
% Note that top and bottom boundary values are zero

% Return:
% left, middle, right: the solution of the PDE at 
% left/right boundary milestones and at middle milestone
function middle = solvePDE(t_xy, vert_dist, ...
    ms_dist, N, n, V, div_V, sigma, l_b, r_b)

% Build coefficient matrix for the discretization

% grid_size
h = vert_dist/n; N_1 = n*(N-1)+1;
% M = number of pts on horizontal direction in our small domain
M = 2*round(ms_dist/h) + 1;
% Some more extra helper vars
mat_is = M*N_1*5 + 4 - 8*(M+N_1);
i = zeros(1, mat_is); j = zeros(1, mat_is); v = zeros(1, mat_is);
cur_ind = 1;

% Find the coefficient matrix
for t = 1:(M*N_1)
    if mod(t-1,M) == 0 || mod(t-1,M) == M-1 || ...
            floor((t-1)/M) == 0 || floor((t-1)/M) == N_1-1
        i(cur_ind) = t; j(cur_ind) = t; v(cur_ind) = 1;
        cur_ind = cur_ind + 1;
        continue;
    end
    x_t = t_xy + [(mod(t-1, M)-(M-1)/2)*h, (t-1)/M * h];
    v_t = V(x_t)/(2*h);
    temp = sigma^2/(2*h^2);
    ind = cur_ind:cur_ind+4;
    i(ind) = [t, t, t, t, t];
    j(ind) = [t, t-1, t+1, t-M, t+M];
    v(ind) = [-div_V(x_t)- 4*temp, v_t(1) + temp, ...
        -v_t(1) + temp, v_t(2) + temp, -v_t(2) + temp];
    cur_ind = cur_ind + 5;
end
Coef = sparse(i,j,v,M*N_1,M*N_1); 
% Build data b based on boundary condition
b = zeros(M*N_1,1);
ind = (0:(N_1-1)) * M;
b(1+ ind) = interp1(1:N, l_b, 1:1/n:N);
b(M+ ind) = interp1(1:N, r_b, 1:1/n:N);
% Solve linear equation to get sol
sol = Coef\b;

% Extract target value (value on middle milestone) from sol
%left = zeros(N, 1); middle = zeros(N,1); right = zeros(N,1);
%left = sol(1+ M*n*(0:(N-1)));
middle = (sol((M+1)/2 +M*n*(0:(N-1))))';
%right = sol(M + M*n*(0:(N-1)));
%left = left/n; middle = middle/n; right = right/n;

end