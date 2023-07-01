%% Read data from VED
% Kevin Moy, 12/15/21
clearvars
close all
clc
set(0,'defaultTextInterpreter','latex');

city_1_dat = readtable("City_1_CellPowerProfile.csv");
city_2_dat = readtable("City_2_CellPowerProfile.csv");
Q_nom = 5; % Ah, assumed
freq = 10; % Hz

% Assume current from front pack
C1_curr = city_1_dat.FrontCurrent_A_./Q_nom;
C2_curr = city_2_dat.FrontCurrent_A_./Q_nom;

% Get velocity
C1_vel = city_1_dat.Velocity_m_s_;
C2_vel = city_2_dat.Velocity_m_s_;

% Get time
C1_time = city_1_dat.Time_s_;
C2_time = city_2_dat.Time_s_;

%% Conversions
m_s_kph = 3.6; % m/s per km/hr


%% Read data

csvfiles = dir("VED_DynamicData_Part1/*.csv");

tmp_9 = table;
tmp_455 = table;
tmp_541 = table;
for file = csvfiles'
    fprintf(1, 'Read %s.\n', file.name)
    tmp = readtable(strcat("VED_DynamicData_Part1/", file.name));
    idx_9 = tmp.VehId == 9;
    idx_455 = tmp.VehId == 455;
    idx_541 = tmp.VehId == 541;

    tmp_9 = [tmp_9; tmp(idx_9, :)];
    tmp_455 = [tmp_455; tmp(idx_455, :)];
    tmp_541 = [tmp_541; tmp(idx_541, :)];
end

% csvfiles = dir("VED_DynamicData_Part2/*.csv");
% 
% tmp_9 = table;
% tmp_455 = table;
% tmp_541 = table;
% for file = csvfiles'
%     fprintf(1, 'Read %s.\n', file.name)
%     tmp = readtable(strcat("VED_DynamicData_Part2/", file.name));
%     idx_9 = tmp.VehId == 9;
%     idx_455 = tmp.VehId == 455;
%     idx_541 = tmp.VehId == 541;
% 
%     tmp_9 = [tmp_9; tmp(idx_9, :)];
%     tmp_455 = [tmp_455; tmp(idx_455, :)];
%     tmp_541 = [tmp_541; tmp(idx_541, :)];
% end

% 
%% Get individual trip data
tripIDs = unique(tmp_455.Trip);
SOC_diff = [];
Q_packs = [];
hFig = figure(100);
set(hFig, 'Position', [100 100 800 800])
hold on
for trip = 1:length(tripIDs)
    idx_tmp = tmp_455.Trip == tripIDs(trip);
    time_tmp = table2array(tmp_455(idx_tmp, 'Timestamp_ms_'))./ 1000 ./60;
    socs_tmp = table2array(tmp_455(idx_tmp, 'HVBatterySOC___'));
    SOC_diff(trip) = socs_tmp(1) - socs_tmp(end);
    plot(time_tmp, socs_tmp, 'LineWidth', 2)
end
ylim([0 100])
xlabel("Time [min]")
ylabel("SOC [\%]")
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
%% Estimate Ah capacity of pack
% Find index of largest SOC difference
[maxdiff, idxdiff] = max(SOC_diff);
idx_tmp = tmp_455.Trip == tripIDs(idxdiff);
time_tmp = table2array(tmp_455(idx_tmp, 'Timestamp_ms_'))./ 1000;
curr_tmp = table2array(tmp_455(idx_tmp, 'HVBatteryCurrent_A_'));

% Get Ah-throughput of current during this trip
max_Ah = trapz(time_tmp, -curr_tmp)./3600;

% Divide by delta SOC to get Ah-capacity
Q_pack = max_Ah ./ (maxdiff/100); % value: 42.7797 Ah

plot(time_tmp, curr_tmp./Q_pack)
%%
metrics_455 = [];
hFig = figure(1);
set(hFig, 'Position', [100 100 800 1200])
hold on
for trip = 1:length(tripIDs)
    idx_tmp = tmp_455.Trip == tripIDs(trip);
    time_tmp = table2array(tmp_455(idx_tmp, 'Timestamp_ms_'))./ 1000;
    velo_tmp = table2array(tmp_455(idx_tmp, 'VehicleSpeed_km_h_')) ./ m_s_kph;
    curr_tmp = table2array(tmp_455(idx_tmp, 'HVBatteryCurrent_A_'))./ Q_pack;
    socs_tmp = table2array(tmp_455(idx_tmp, 'HVBatterySOC___'));

    % Resample to 1 second
    time_new = [0:1:floor(time_tmp(end))]';
    curr_new = interp1(time_tmp, -curr_tmp, time_new);
    velo_new = interp1(time_tmp, velo_tmp, time_new);
    socs_new = interp1(time_tmp, socs_tmp, time_new);

    % Compute metrics
    metrics_tmp = metrics(curr_new, velo_new, 1);

    % Save to matrix
    metrics_455 = [metrics_455, metrics_tmp];
    
    subplot(3,1,1)
    hold on
    plot(time_new./60, curr_new, 'LineWidth', 2)
    subplot(3,1,2)
    hold on
    plot(time_new./60, velo_new, 'LineWidth', 2)
    subplot(3,1,3)
    hold on
    plot(time_new./60, socs_new, 'LineWidth', 2)
    
end

subplot(3,1,1)
title("Trip Current")
% xlabel("Time [s]")
ylabel("C-rate [-]")
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
subplot(3,1,2)
title("Trip Velocity")
% xlabel("Time [s]")
ylabel("Velocity [m/s]")
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
subplot(3,1,3)
title("Trip SOC")
ylim([0 100])
xlabel("Time [min]")
ylabel("SOC [\%]")
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)

%% Compute metrics for drive cycles
metrics_C1 = metricsdc(C1_curr, C1_vel, freq);
metrics_C2 = metricsdc(C2_curr, C2_vel, freq);



%% Save to table
T1 = array2table(metrics_455, ...
    'RowNames',{'Peak Positive C-rate','Peak Negative C-rate', ...
    'Average Positive C-rate', 'Average Negative C-rate', 'Average C-rate', 'Peak Frequency, Positive C-rate [Hz]', ...
    'Peak Frequency, Negative C-rate [Hz]', 'Variance of Positive C-rate', 'Variance of Negative C-rate', ...
    'Distance Travelled [m]', 'Average Velocity [m/s]', 'Peak Velocity [m/s]', 'Variance of velocity [m^2/s^2]', ...
    'Number of start/stops'});
writetable(T1, 'Summary_455_drive_cycles.csv','WriteRowNames',true)   

% Save to table
T2 = array2table(metrics_C1', ...
    'RowNames',{'Peak Positive C-rate','Peak Negative C-rate', ...
    'Average Positive C-rate', 'Average Negative C-rate', 'Average C-rate', 'Peak Frequency, Positive C-rate [Hz]', ...
    'Peak Frequency, Negative C-rate [Hz]', 'Variance of Positive C-rate', 'Variance of Negative C-rate', ...
    'Distance Travelled [m]', 'Average Velocity [m/s]', 'Peak Velocity [m/s]', 'Variance of velocity [m^2/s^2]', ...
    'Number of start/stops'});
writetable(T2, 'Summary_C1_drive_cycles.csv','WriteRowNames',true)   

% Save to table
T3 = array2table(metrics_C2', ...
    'RowNames',{'Peak Positive C-rate','Peak Negative C-rate', ...
    'Average Positive C-rate', 'Average Negative C-rate', 'Average C-rate', 'Peak Frequency, Positive C-rate [Hz]', ...
    'Peak Frequency, Negative C-rate [Hz]', 'Variance of Positive C-rate', 'Variance of Negative C-rate', ...
    'Distance Travelled [m]', 'Average Velocity [m/s]', 'Peak Velocity [m/s]', 'Variance of velocity [m^2/s^2]', ...
    'Number of start/stops'});
writetable(T3, 'Summary_C2_drive_cycles.csv','WriteRowNames',true)   

%% scratch
idx_tmp = tmp_455.Trip == 1184;
time_tmp = table2array(tmp_455(idx_tmp, 'Timestamp_ms_'))./ 1000;
velo_tmp = table2array(tmp_455(idx_tmp, 'VehicleSpeed_km_h_')) .* m_s_kph;
curr_tmp = table2array(tmp_455(idx_tmp, 'HVBatteryCurrent_A_'));
socs_tmp = table2array(tmp_455(idx_tmp, 'HVBatterySOC___'));

figure; hold on
yyaxis left
plot(time_tmp, velo_tmp)
yyaxis right
plot(time_tmp, curr_tmp)
plot(time_tmp, socs_tmp)
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 26)

%% Functions

% Function to get metrics for a given duty cycle

function metric_matrix = metrics(cycl_d, cyclv_d, freq)
    % - number of occurrences
    % - hours spent (duration of) discharging
    % - peak discharge power
    % - average discharge power
    % - hour starting discharge
    % - hour ending discharge
    % - SOC starting discharge
    % - SOC ending discharge
    % - Average SOC during discharge
    % - peak frequency of discharge? (Maybe largest 3 peaks?)

    % Number of discharging periods
    disp_pos = max(cycl_d,0);
    dischg_start = find(disp_pos(1:end-1)==0 & disp_pos(2:end) > 0); % Starting point of each discharge:
    num_dischg = length(dischg_start);

    % Number of charging periods
    disp_neg = min(cycl_d,0);
    chg_start = find(disp_neg(1:end-1)==0 & disp_neg(2:end) < 0); % Starting point of each charge:
    num_chg = length(chg_start);

    % Find when battery is actively discharging:
    idx_pos = find(cycl_d > 0);

    % Find when battery is actively charging:
    idx_neg = find(cycl_d < 0);

    % Time spent discharging:
    duration_pos = length(idx_pos);

    % Time spent charging:
    duration_neg = length(idx_neg);

    % Distance travelled over the trip
    dist = sum(cyclv_d)/freq;

    % Average velocity over the trip 
    avg_vel = mean(cyclv_d);

    %peak velocity over the trip % -- remove due to strong correlation with avg vel
    pk_vel = max(cyclv_d);

    % Variance of velocity over the trip -- remove due to strong correlation with avg vel
    var_vel = var(cyclv_d);

    % # of times start/stop (start only) within trip
    % Remove noise by thresholding
    v_d = cyclv_d;
    v_d(v_d < 0.1) = 0;
    diff_v = diff(v_d);
    sts = find(diff_v(1:end-1)==0 & diff_v(2:end) > 0); % point of each start
    num_sts = length(sts);

    % Peak discharging power
    pk_pow_pos = max(disp_pos);

    % Peak charging power (magnitude)
    pk_pow_neg = -min(disp_neg);

    % Average discharging power
    avg_pow_pos = mean(disp_pos(disp_pos~=0));

    % Average charging power (magnitude)
    avg_pow_neg = -mean(disp_neg(disp_neg~=0));

    % Average power
    avg_pow = mean(cycl_d);

    % Find peak frequency of discharge and charge

    P_day_rep = repmat(cycl_d,100,1);

    [disp_pos, disp_neg] = mean_centering(P_day_rep);

    N = length(disp_pos);
    y = fft(disp_pos);
    ft_dischg_freq_drive_cycle = (0:N-1)'*(freq/N);       % frequency range
    ft_dischg_power_drive_cycle = abs(y).^2/N/freq;      % periodogram of the DFT
    N = length(disp_neg);
    y = fft(disp_neg);
    ft_chg_freq_drive_cycle = (0:N-1)'*(freq/N);       % frequency range
    ft_chg_power_drive_cycle = abs(y).^2/N/freq;      % periodogram of the DFT

    if ~isempty(ft_dischg_freq_drive_cycle)
        [~,locs_dischg,~,~] = findpeaks(ft_dischg_power_drive_cycle(1:(length(ft_dischg_freq_drive_cycle)/2 - 1)), 'SortStr','descend','NPeaks',1);
        freq_peak_pos = ft_dischg_freq_drive_cycle(locs_dischg);
    end

    if ~isempty(ft_chg_freq_drive_cycle)
        [~,locs_chg,~,~] = findpeaks(ft_chg_power_drive_cycle(1:(length(ft_chg_freq_drive_cycle)/2 - 1)), 'SortStr','descend','NPeaks', 1);
        freq_peak_neg = ft_chg_freq_drive_cycle(locs_chg);
    end
    
    var_disp_pos = var(cycl_d(cycl_d > 0));
    var_disp_neg = var(cycl_d(cycl_d < 0));

    % Set metrics to 0 if there is no charge or discharge over the day
    if isempty(idx_neg)
        pk_pow_neg = 0;
        avg_pow_neg = 0;
        freq_peak_neg = 0;
        var_disp_neg = 0;

    end
    if isempty(idx_pos)
        pk_pow_pos = 0;
        avg_pow_pos = 0;
        freq_peak_pos = 0;
        var_disp_pos = 0;
    end
    if isempty(freq_peak_pos)
        freq_peak_pos = 0;
    end

    if isempty(freq_peak_neg)
        freq_peak_neg = 0;
    end
    % Construct metrics summary for the duty cycle

    metric_matrix = [pk_pow_pos, pk_pow_neg, avg_pow_pos, avg_pow_neg, avg_pow, ...
                    freq_peak_pos, freq_peak_neg,...
                    var_disp_pos, var_disp_neg, ...
                    dist, avg_vel, pk_vel, var_vel, num_sts]';
end

function metrics_dc = metricsdc(drive_cycle_crate, drive_cycle_vel, freq)

    diff_thresh = 1e-5; % Threshold to consider zero velocity
    
    twomin = 30*freq; % Length of 30 seconds per stop (approx.)
    
    % Complex diff threshold method
    dtest = diff(drive_cycle_vel);
    nf = find(abs(dtest) > 0);
    dnf = diff(nf);   
    flats = find(dnf > twomin);    
    flats_nf = nf(flats);
    cyc_st = [nf(1)];
    cyc_end = [];
    for i = 1:length(flats)
        cyc_end = [cyc_end; nf(flats(i))+1];
        cyc_st = [cyc_st; nf(flats(i)+1)];
    end
    cyc_end = [cyc_end; nf(end)];

    % --- DAILY METRICS ---
    % Compute PCA Metrics matrix
    % Each row is one nonzero dispatch day in the year
    % Each column is a different metric
    % Split into two metric matrices: urban and highway driving
    
    % Set threshold for highway driving
    hwy_thresh = 20; % m/s
    
    close all
    metrics_dc = [];
    metrics_hwy = [];
    cyc = {};
    cyc_v = {};
    cyc_hwy = {};
    cyc_v_hwy = {};
    for d = 1:length(cyc_st) 
        % - number of occurrences
        % - hours spent (duration of) discharging
        % - peak discharge power
        % - average discharge power
        % - hour starting discharge
        % - hour ending discharge
        % - SOC starting discharge
        % - SOC ending discharge
        % - Average SOC during discharge
        % - peak frequency of discharge? (Maybe largest 3 peaks?)
        
        cycl_d = drive_cycle_crate(cyc_st(d):cyc_end(d));
        cyclv_d = drive_cycle_vel(cyc_st(d):cyc_end(d));
    
        % Number of discharging periods
        disp_pos = max(cycl_d,0);
        dischg_start = find(disp_pos(1:end-1)==0 & disp_pos(2:end) > 0); % Starting point of each discharge:
        num_dischg = length(dischg_start);
    
        % Number of charging periods
        disp_neg = min(cycl_d,0);
        chg_start = find(disp_neg(1:end-1)==0 & disp_neg(2:end) < 0); % Starting point of each charge:
        num_chg = length(chg_start);
    
        % Find when battery is actively discharging:
        idx_pos = find(cycl_d > 0);
    
        % Find when battery is actively charging:
        idx_neg = find(cycl_d < 0);
    
        % Time spent discharging:
        duration_pos = length(idx_pos);
    
        % Time spent charging:
        duration_neg = length(idx_neg);
    
        % Distance travelled over the trip
        dist = sum(cyclv_d)/freq;
    
        % Average velocity over the trip 
        avg_vel = mean(cyclv_d);
    
        %peak velocity over the trip % -- remove due to strong correlation with avg vel
        pk_vel = max(cyclv_d);
    
        % Variance of velocity over the trip -- remove due to strong correlation with avg vel
        var_vel = var(cyclv_d);
    
        % # of times start/stop (start only) within trip
        % Remove noise by thresholding
        v_d = cyclv_d;
        v_d(v_d < 0.1) = 0;
        diff_v = diff(v_d);
        sts = find(diff_v(1:end-1)==0 & diff_v(2:end) > 0); % point of each start
        num_sts = length(sts);
    
        % Peak discharging power
        pk_pow_pos = max(disp_pos);
    
        % Peak charging power (magnitude)
        pk_pow_neg = -min(disp_neg);
    
        % Average discharging power
        avg_pow_pos = mean(disp_pos(disp_pos~=0));
    
        % Average charging power (magnitude)
        avg_pow_neg = -mean(disp_neg(disp_neg~=0));
    
        % Average power
        avg_pow = mean(cycl_d);
    
        % Find peak frequency of discharge and charge
    
        P_day_rep = repmat(cycl_d,100,1);
    
        [disp_pos, disp_neg] = mean_centering(P_day_rep);
    
        N = length(disp_pos);
        y = fft(disp_pos);
        ft_dischg_freq_drive_cycle = (0:N-1)'*(freq/N);       % frequency range
        ft_dischg_power_drive_cycle = abs(y).^2/N/freq;      % periodogram of the DFT
        N = length(disp_neg);
        y = fft(disp_neg);
        ft_chg_freq_drive_cycle = (0:N-1)'*(freq/N);       % frequency range
        ft_chg_power_drive_cycle = abs(y).^2/N/freq;      % periodogram of the DFT
    
        if ~isempty(ft_dischg_freq_drive_cycle)
            [~,locs_dischg,~,~] = findpeaks(ft_dischg_power_drive_cycle(1:(length(ft_dischg_freq_drive_cycle)/2 - 1)), 'SortStr','descend','NPeaks',1);
            freq_peak_pos = ft_dischg_freq_drive_cycle(locs_dischg);
        end
    
        if ~isempty(ft_chg_freq_drive_cycle)
            [~,locs_chg,~,~] = findpeaks(ft_chg_power_drive_cycle(1:(length(ft_chg_freq_drive_cycle)/2 - 1)), 'SortStr','descend','NPeaks', 1);
            freq_peak_neg = ft_chg_freq_drive_cycle(locs_chg);
        end
        
        var_disp_pos = var(cycl_d(cycl_d > 0));
        var_disp_neg = var(cycl_d(cycl_d < 0));
    
        % Set metrics to 0 if there is no charge or discharge over the day
        if isempty(idx_neg)
            pk_pow_neg = 0;
            avg_pow_neg = 0;
            freq_peak_neg = 0;
            var_disp_neg = 0;
    
        end
        if isempty(idx_pos)
            pk_pow_pos = 0;
            avg_pow_pos = 0;
            freq_peak_pos = 0;
            var_disp_pos = 0;
        end
        if isempty(freq_peak_pos)
            freq_peak_pos = 0;
        end
    
        if isempty(freq_peak_neg)
            freq_peak_neg = 0;
        end
        % Construct metrics summary for the day
                    
    %     metrics(d,:) = [ num_dischg, num_chg, duration_pos, duration_neg, ...
    %                     pk_pow_pos, pk_pow_neg, avg_pow_pos, avg_pow_neg, ...
    %                     freq_peak_pos, freq_peak_neg,...
    %                     var_disp_pos, var_disp_neg, ...
    %                     dist, avg_vel, pk_vel, var_vel, num_sts];    
    
        metrics_d = [pk_pow_pos, pk_pow_neg, avg_pow_pos, avg_pow_neg, avg_pow, ...
                        freq_peak_pos, freq_peak_neg,...
                        var_disp_pos, var_disp_neg, ...
                        dist, avg_vel, pk_vel, var_vel, num_sts];
    
    %         metrics_d = [pk_pow_pos, pk_pow_neg, avg_pow_pos, avg_pow_neg, ...
    %                     var_disp_pos, var_disp_neg, ...
    %                     dist, avg_vel, num_sts];    
            
            % Sort into urban and highway matrices
            if avg_vel < hwy_thresh
                metrics_dc = [metrics_dc; metrics_d];
                cyc = [cyc; cycl_d];
                cyc_v = [cyc_v; cyclv_d]; 
            else
                metrics_hwy = [metrics_hwy; metrics_d];
                cyc_hwy = [cyc_hwy; cycl_d];
                cyc_v_hwy = [cyc_v_hwy; cyclv_d]; 
            end
    
    end
    
    % Remove microcycle 46 (highway driving)? Get metrics separately for 46
    % metrics(46,:) = [];
    % metrics(251,:) = [];
    % metrics(731,:) = [];
    
    urban_vel = vertcat(cyc_v{:});
    urban_crt = vertcat(cyc{:});
    highw_vel = vertcat(cyc_v_hwy{:});
    highw_crt = vertcat(cyc_hwy{:});
    
    disp("Computation of metrics finished!")
end
