%% Rest Lengths
% Function that returns vector of lengths when time-series data supplied is
% 0.
% Example: [0 0 0 0 0 1 2 3 2 0 0 3 5 1 0 32 1000 02 0 0]
% Returns: [5 2 1 2]
function duration = rest_lengths(timeseries_data)
    zero_idx = find(timeseries_data == 0);
    a = diff(zero_idx);
    b = find([a; inf] > 1);
    duration = diff([0; b]);
end