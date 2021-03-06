function [direction_indices] = isolate_direction(interp_behav,direction)
%ISOLATE_DIRECTION Summary of this function goes here
%   Detailed explanation goes here

z_delta_x = diff(interp_behav);
z_delta_x(end+1) = 0;
z_delta_x(isnan(z_delta_x)) = 0;
z_delta_x = zscore(z_delta_x);

switch direction
    case 'right'
    direction_indices = z_delta_x > 0.2;
    case 'left'
    direction_indices = z_delta_x < -0.2;
    otherwise
    error('Please indicate a valid direction. The input variable direction should be either "right" or "left"');
end

end