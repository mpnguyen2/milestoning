Files description:

solvePDE.m: fct to solve Fokker-Planck (Kolmogorov forward) using given boundary condition
shootTraj.m: fct to calculate the hitting density of a diffusion process from a specific point
milestone.m: fct calculate the density using milestone approach
pdeMilestone.m: fct calculate the density using pde approach
draw_anim.m: fct to draw the animation given the data from two approaches
run_sim_on_rect.m: run the simulation comparing two approach (on rectangle) using the functions draw_anim, milestone, and pdeMilestone
test.m: testing function on various data using run_sim_on_rect

3 Tests implemented:

First test: num_m = 5, ms_dist = 2, N = 20, vert_dist = 0.2
Second test: num_m = 10, ms_dist = 2, N = 20, vert_dist = 0.2
Third test: num_m = 10, ms_dist = 2, N = 30, vert_dist = 0.2

All have the same lower corner (0, 0)