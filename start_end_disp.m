%% Start/End Dispatch
% Current Version 6/30/2020
% Used for any non-negative (currently, hourly) time-series data.
% Finds the start and end time for each non-negative dispatch sequence 
% in the data.
% Non-negative dispatch sequence = start from 0, positive, return to 0

% Input ts_data is an Nx1 vector.
function [ts_pos_ind] = start_end_disp(ts_pos)

    % Starting point of each discharge:
    ts_pos_start = find(ts_pos(1:end-1)==0 & ts_pos(2:end) > 0);

    % Ending point of each discharge:
    ts_pos_end = find(ts_pos(1:end-1) > 0 & ts_pos(2:end) == 0) + 1;
    ts_pos_ind = [ts_pos_start ts_pos_end];

end