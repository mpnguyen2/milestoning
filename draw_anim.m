%% Function to draw and save the animation for simulation for comparing
% two approaches PDE vs Milestone
% Args:
% file_name: name of the file to be saved
% data_de: Data on all iteration step of the PDE approach; each row 
% contains num_m * N density values on each iteration
% data_ms: Data on all iteration step of the milestone approach; 
% each row contains num_m * N density values on each iteration
% num_m: number of milestone
% N: number of (discrete) pts on each milestone
% rescale: rescale factor to make density values more visible in the
% grayscale image
% factor: interpolation factor to make the image transition smoother
function draw_anim(file_name, data_de, data_ms, N, num_m, rescale, factor)

num_pics = size(data_de,1) * size(data_ms,1) * factor;
step_de = num_pics/size(data_de,1);
step_ms = num_pics/size(data_ms, 1);

axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';

%num_row = 300; num_col = floor(num_row * ...
%    ms_dist *(num_m-1)/((vert_dist)*(N-1)));
vid = VideoWriter(file_name);
open(vid);

for t = 0:num_pics-1
    % additional helper vars
    ind_de = floor(t/step_de); ind_ms = floor(t/step_ms);
    fact_de = (t-ind_de*step_de)/step_de;
    fact_ms = (t-ind_ms*step_ms)/step_ms;
    
    % Interpolate and reshape data for pde approach
    if ind_de == size(data_de,1)-1
        data_odd = data_de(ind_de + 1,:);
    else
        data_odd = (1-fact_de) * data_de(ind_de + 1,:) + ...
            fact_de * data_de(ind_de+2,:);
    end
    data_odd = rescale * reshape(data_odd, num_m, N);
    data_odd = imcomplement(data_odd);
    
    % Interpolate and reshape data for milestone approach
    if ind_ms == size(data_ms,1)-1
        data_even = data_ms(ind_ms + 1,:);
    else
        data_even = (1-fact_ms) * data_ms(ind_ms + 1,:) + ...
            fact_ms * data_ms(ind_ms+2,:);
    end
    data_even = rescale * reshape(data_even, num_m, N);
    data_even = imcomplement(data_even);
    
    % Draw 2 figures iteratively (to create an animation)
    subplot(1,2,1);
    imshow((data_odd)', 'InitialMagnification', 'fit');
    title('PDE method');
    subplot(1,2,2);
    imshow((data_even)', 'InitialMagnification', 'fit');
    title('Milestone method');
    drawnow
    frame = getframe(gcf);
    writeVideo(vid, frame);
end
close(vid);

end