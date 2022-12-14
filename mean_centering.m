%% Mean-centering
% Current Version 6/23/2020
% Used for any (currently, hourly) time-series data.
% General idea is to separate discharge (positive) and charge (negative)
% components, and create new vectors where each individual discharge/charge
% instance is inserted along with a flipped (reversed sign) instance.

% Input ts_data is an 8760x1 vector.
% TODO: Add support for other time-series data (e.g. 15-min data)
function [ts_pos_nz, ts_neg_nz] = mean_centering(ts_data)

    % Threshold for considering zero dispatch:
    zero_thr = 1e-5;
    
    % Need to add leading and trailing zeros to ensure the dispatch of the
    % beginning and ending of ts_data is properly captured as a charge or
    % discharge.

    % Add a 0 to the beginning if the beginning value is not 0
    if ts_data(1) ~=0
        ts_data = [0; ts_data];
    end
    
    % Add a 0 to the end if the ending value is not 0
    if ts_data(end) ~=0
        ts_data = [ts_data; 0];
    end
    
    % ---- Discharge ---- %
    % Set negative values to 0:
    ts_pos = max(ts_data,0);

    ts_pos_ind = start_end_disp(ts_pos);

    % Insert negative + reversed copy of discharge into vector
    ts_pos2 = [];
    for i = 1:size(ts_pos_ind, 1)
        st = ts_pos_ind(i,1);
        en = ts_pos_ind(i,2);
        ts_pos2 = [ts_pos2; ts_pos(st:en); -flip(ts_pos(st:en))];
    end

    % Remove zero entries:
    ts_pos_nz = ts_pos2;
    ts_pos_nz(abs(ts_pos2) < zero_thr) = [];
    
    % ---- Charge ---- %
    % Set initially positive values to 0 and then make all values positive:
    ts_neg = -min(ts_data,0);

    ts_neg_ind = start_end_disp(ts_neg);

    % Insert negative + reversed copy of discharge into vector
    ts_neg2 = [];
    for i = 1:size(ts_neg_ind, 1)
        st = ts_neg_ind(i,1);
        en = ts_neg_ind(i,2);
        ts_neg2 = [ts_neg2; ts_neg(st:en); -flip(ts_neg(st:en))];
    end

    % Remove zero entries:
    ts_neg_nz = ts_neg2;
    ts_neg_nz(abs(ts_neg2) < zero_thr) = [];
end