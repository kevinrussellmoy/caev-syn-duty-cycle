%% Synthetic duty cycle generation using PCA and k-means clustering
% Publications:
% https://doi.org/10.1016/j.adapen.2021.100065
% https://doi.org/10.1115/1.4050192
% Kevin Moy
% Latest revisions: Dec 12 2022

clearvars
close all
clc
set(0,'defaultTextInterpreter','latex'); 

% ----- Load in Drive Cycle data -----
% Will be removed once public -- will not publish confidential data
city_1_dat = readtable("City_1_CellPowerProfile.csv");
city_2_dat = readtable("City_2_CellPowerProfile.csv");
Q_nom = 5; % Ah, assumed
freq = 1; % Hz

% Obtain current from front pack
city_1_curr = downsample(city_1_dat.FrontCurrent_A_, 10);
city_2_curr = downsample(city_2_dat.FrontCurrent_A_, 10);

% Vehicle velocity
city_1_vel = downsample(city_1_dat.Velocity_m_s_, 10);
city_2_vel = downsample(city_2_dat.Velocity_m_s_, 10);

% Combine both into cell arrays for data selection (below)

city_curr = {city_1_curr, city_2_curr};
city_velo = {city_1_vel, city_2_vel};

%% Select drive cycle velocity, crate data

% SELECT:
% 1 = City 1
% 2 = City 2
% 3 = City 1 and City 2 combined
city_select = 1;

drive_cycle_velo = [];
drive_cycle_crat = [];
drive_cycle_name = [];

% Get drive cycle velocity
if city_select == 1 || city_select == 2
    drive_cycle_velo = city_velo{city_select};
    drive_cycle_crat = city_curr{city_select}/Q_nom;
    drive_cycle_name = ['City ' num2str(city_select)];
    drive_cycle_fn   = ['City_' num2str(city_select)];
elseif city_select == 3
    drive_cycle_velo = [city_velo{1}; city_velo{2}];
    drive_cycle_crat = [city_curr{1}; city_curr{2}]/Q_nom;
    drive_cycle_name = 'City 1+2 combined';
    drive_cycle_fn   = 'City_1_2_combined';
end

% % drive_cycle_velo = city_1_vel;
% drive_cycle_velo = city_2_vel;
% % drive_cycle_velo = [city_1_vel; city_2_vel];
% 
% % Convert current to C-rate
% % drive_cycle_crat = city_1_curr/Q_nom;
% drive_cycle_crat = city_2_curr/Q_nom;
% % drive_cycle_crat = [city_1_curr; city_2_curr]/Q_nom;

%% Get individual "microtrips"
% Assume each microtrip is 30 seconds of rest (zero velocity) apart

diff_thresh = 1e-5; % Threshold to consider zero velocity
stop_len = 30*freq; % Length of 30 seconds per stop (assumed)

% Diff threshold method for stops
dtest = diff(drive_cycle_velo);
nf = find(abs(dtest) > 0);
dnf = diff(nf);
flats = find(dnf > stop_len);
flats_nf = nf(flats);

cyc_st = [nf(1)];
cyc_end = [];
for i = 1:length(flats)
    cyc_end = [cyc_end; nf(flats(i))+1];
    cyc_st = [cyc_st; nf(flats(i)+1)];
end
cyc_end = [cyc_end; nf(end)];

hFig = figure(2);
set(hFig, 'Position', [100 100 800 800])
hold on
plot((1:length(drive_cycle_velo))/freq/3600, drive_cycle_velo, 'LineWidth', 1)
scatter(cyc_st/freq/3600, drive_cycle_velo(cyc_st), 100, 'LineWidth', 2)
scatter(cyc_end/freq/3600, drive_cycle_velo(cyc_end), 100, 'LineWidth', 2)
ylabel('Velocity [m/s]', 'interpreter', 'latex')
xlabel('Time [hr]', 'interpreter', 'latex')
legend('Velocity', 'Cycle Start', 'Cycle End', 'Location', 'best', 'interpreter', 'latex')
title({[drive_cycle_name ' Drive Cycle Velocity Data']})
set(gca, 'FontSize', 24)
set(gca, 'TickLabelInterpreter','latex')
hold off
box on

hFig = figure(3);
set(hFig, 'Position', [100 100 800 800])
hold on
for i = 1:length(cyc_st)
    dc_vel_i = drive_cycle_velo(cyc_st(i):cyc_end(i));
    plot((1:length(dc_vel_i))/freq/60, dc_vel_i)
end
ylabel('Velocity [m/s]', 'interpreter', 'latex')
xlabel('Time [min]', 'interpreter', 'latex')
xlim([0 40])
ylim([0 30])
title({[drive_cycle_name ' Drive Cycle Duty Cycles']})
set(gca, 'FontSize', 24)
set(gca, 'TickLabelInterpreter','latex')
hold off
box on

hFig = figure(4);
set(hFig, 'Position', [100 100 800 800])
hold on
for i = 1:length(cyc_st)
    dc_curr_i = drive_cycle_crat(cyc_st(i):cyc_end(i));
    plot((1:length(dc_curr_i))/freq/60, dc_curr_i)
end
xlim([0 40])
ylim([-2 2])
ylabel('C-rate [-]', 'interpreter', 'latex')
xlabel('Time [min]', 'interpreter', 'latex')
title({[drive_cycle_name ' Drive Cycle Duty Cycles']})
set(gca, 'FontSize', 24)
set(gca, 'TickLabelInterpreter','latex')
hold off
box on

%% --- DAILY METRICS ---
% Compute PCA Metrics matrix
% Each row is one nonzero dispatch day in the year
% Each column is a different metric
% Split into two metric matrices: urban and highway driving

% Set threshold for highway driving
hwy_thresh = 20; % m/s

close all
metrics = [];
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
    
    cycl_d = drive_cycle_crat(cyc_st(d):cyc_end(d));
    cyclv_d = drive_cycle_velo(cyc_st(d):cyc_end(d));

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

    % Construct metrics summary for the microtrip
    metrics_d = [pk_pow_pos, pk_pow_neg, avg_pow_pos, avg_pow_neg, avg_pow, ...
                    freq_peak_pos, freq_peak_neg,...
                    var_disp_pos, var_disp_neg, ...
                    dist, avg_vel, pk_vel, var_vel, num_sts]; 
        
    % Sort into urban and highway matrices
    if avg_vel < hwy_thresh
        metrics = [metrics; metrics_d];
        cyc = [cyc; cycl_d];
        cyc_v = [cyc_v; cyclv_d]; 
    else
        metrics_hwy = [metrics_hwy; metrics_d];
        cyc_hwy = [cyc_hwy; cycl_d];
        cyc_v_hwy = [cyc_v_hwy; cyclv_d]; 
    end

end

urban_drive_cycle_velocity = vertcat(cyc_v{:});
urban_drive_cycle_crate = vertcat(cyc{:});
highway_drive_cycle_velocity = vertcat(cyc_v_hwy{:});
highway_drive_cycle_crate = vertcat(cyc_hwy{:});

disp("Computation of metrics finished!")


%% Apply PCA/k-means algorithm
% Note position 3 of cyc in pca_k_means func call, since we don't have any
% state data
% it's not actually being returned anyway, but if we had relevant state
% data we would return it in position 2 of the function return

[char_crate, ~, char_vel, char_interval] = pca_k_means(metrics, cyc, cyc, cyc_v, true);

%% Plotting + Average current
% close all
% linestyle = ["-.", "-", ":"];

dctable = [];
dcrate = urban_drive_cycle_crate;
num_clust = length(char_crate);

% TODO: Add max/mean discharge, charge current separately; charge
% throughput for discharge and charge separately
dctable(1,1) = max(dcrate(dcrate > 0));
dctable(1,2) = min(dcrate(dcrate < 0));
dctable(1,3) = mean(dcrate(dcrate > 0));
dctable(1,4) = mean(dcrate(dcrate < 0));
dctable(1,5) = length(dcrate(dcrate > dctable(1,3)))/length(dcrate);
dctable(1,6) = length(dcrate(dcrate < dctable(1,3)))/length(dcrate);

fprintf('\n')
hFig = figure(12);
set(hFig, 'Position', [100 100 1000 500])
set(0,'DefaultAxesLineStyleOrder',{'-',':'});
hold on
for i = 1:num_clust
    syn_duty_cycle = char_crate{i};
    subplot(1,num_clust,i)
    plot((1:length(syn_duty_cycle))./freq/60, syn_duty_cycle, 'LineWidth', 1.5)
    set(gca, 'TickLabelInterpreter','latex')
    set(gca,'FontSize', 20)
    xlabel('Time [min]')
    ylabel('C-rate [-]')
%     xlim([0, 60])
    ylim([-1, 1])
    title(sprintf('Characteristic Duty Cycle %d', i))
    dctable(i+1,1) = max(syn_duty_cycle(syn_duty_cycle > 0));
    dctable(i+1,2) = min(syn_duty_cycle(syn_duty_cycle < 0));
    dctable(i+1,3) = mean(syn_duty_cycle(syn_duty_cycle > 0));
    dctable(i+1,4) = mean(syn_duty_cycle(syn_duty_cycle < 0));
    dctable(i+1,5) = length(syn_duty_cycle(syn_duty_cycle > dctable(1,3)))/length(syn_duty_cycle);
    dctable(i+1,6) = length(syn_duty_cycle(syn_duty_cycle < dctable(1,3)))/length(syn_duty_cycle);
    fprintf('Characteristic Duty Cycle %d, average discharge current %4.4f \n', i, dctable(i+1,3))
    fprintf('Characteristic Duty Cycle %d, average charge current %4.4f \n', i, dctable(i+1,4))
    fprintf('Characteristic Duty Cycle %d, average current %4.4f \n', i, mean(syn_duty_cycle))
end
box on
fprintf('Average Discharge Current for entire drive cycle: %4.4f \n', dctable(1,3))
fprintf('Average Charge Current for entire drive cycle: %4.4f \n', dctable(1,4))
fprintf('Average Current for entire drive cycle: %4.4f \n', mean(dcrate))

hFig = figure(13);
set(hFig, 'Position', [100 100 1000 500])
set(0,'DefaultAxesLineStyleOrder',{'-',':'});
hold on
for i = 1:num_clust
    syn_duty_cycle_v = char_vel{i};
    subplot(1,num_clust,i)
    plot((1:length(syn_duty_cycle_v))./freq/60, syn_duty_cycle_v, 'LineWidth', 1.5)
    set(gca, 'TickLabelInterpreter','latex')
    set(gca,'FontSize', 20)
    xlabel('Time [min]')
    ylabel('Velocity [m/s]')
    title(sprintf('Characteristic Duty Cycle %d', i))
end
box on

fn = [drive_cycle_fn '.mat'];

save(fn, 'char_crate', 'char_vel')


%%
T = array2table(dctable, ...
    'VariableNames', {'Maximum Discharge Current', ...
    'Maximum Charge Current' 'Average Discharge Current', ...
    'Average Charge Current', 'Percentage of Time above Average Discharge Current', ...
    'Percentage of Time Below Average Discharge Current'}, ...
    'RowNames', {'Drive Cycle', 'Characteristic Duty Cycle 1', ...
    'Characteristic Duty Cycle 2', 'Characteristic Duty Cycle 3'});

%% Write to table
% writetable(T, 'City_1_Duty_Cycle_Comparison.xls', 'WriteRowNames', true)
% writetable(T, 'City_2_Duty_Cycle_Comparison.xls', 'WriteRowNames', true)

%% HIGHWAY DRIVE CYCLE

hFig = figure(1); 
set(hFig, 'Position', [100 100 800 1000])
plot((1:length(highway_drive_cycle_crate))/freq/60, highway_drive_cycle_crate)
ylabel('C-rate [-]')
xlabel('Time [min]')
title('Highway Drive Cycle')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)

avg_crate = 1/10;

% Determine length of rest s.t. avg. dischg, chg current = that of
% the drive cycle
hwy_rest_len = floor(sum(highway_drive_cycle_crate)/avg_crate - ...
    length(highway_drive_cycle_crate));

% Construct synthetic duty cycle 1a
syn_duty_cycle_1a = [ highway_drive_cycle_crate; zeros(hwy_rest_len, 1 ) ];
syn_duty_cycle_1a_vel = [ highway_drive_cycle_veloocity; zeros(hwy_rest_len, 1 ) ];

% Construct synthetic duty cycle 1b
% Parameter for number of replicates in a row
reps = 4;
syn_duty_cycle_1b = [repmat(highway_drive_cycle_crate, reps, 1); repmat(zeros(hwy_rest_len, 1), reps, 1)];
syn_duty_cycle_1b_vel = [repmat(highway_drive_cycle_veloocity, reps, 1); repmat(zeros(hwy_rest_len, 1), reps, 1)];

hFig = figure(12);
set(hFig, 'Position', [100 100 1000 600])
tiledlayout(1,2);
nexttile
plot((1:length(syn_duty_cycle_1a))./freq/60, syn_duty_cycle_1a, 'LineWidth', 1)
xlim([0, length(syn_duty_cycle_1b)/freq/60])
xlabel('Time [min]')
ylabel('C-rate [-]')
title('Syn. Duty Cycle 1a')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
nexttile
plot((1:length(syn_duty_cycle_1b))./freq/60, syn_duty_cycle_1b, 'LineWidth', 1)
xlim([0, length(syn_duty_cycle_1b)/freq/60])
xlabel('Time [min]')
ylabel('C-rate [-]')
title('Syn. Duty Cycle 1b')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)

t = sgtitle('Highway Duty Cycles');
t.Interpreter = 'latex';
t.FontSize = 24;

hFig = figure(13);
set(hFig, 'Position', [100 100 1000 600])
tiledlayout(1,2);
nexttile
plot((1:length(syn_duty_cycle_1a_vel))./freq/60, syn_duty_cycle_1a_vel, 'LineWidth', 1)
xlim([0, length(syn_duty_cycle_1a_vel)/freq/60])
xlabel('Time [min]')
ylabel('Velocity [m/s]')
title('Syn. Duty Cycle 1a')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
nexttile
plot((1:length(syn_duty_cycle_1b_vel))./freq/60, syn_duty_cycle_1b_vel, 'LineWidth', 1)
xlim([0, length(syn_duty_cycle_1b_vel)/freq/60])
xlabel('Time [min]')
ylabel('Velocity [m/s]')
title('Syn. Duty Cycle 1b')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)

t = sgtitle('Highway Duty Cycles');
t.Interpreter = 'latex';
t.FontSize = 24;

syn_duty_1ab = {syn_duty_cycle_1a, syn_duty_cycle_1b};
syn_duty_vel_1ab = {syn_duty_cycle_1a_vel, syn_duty_cycle_1b_vel}; 

fn = [drive_cycle_fn '_hwy.mat'];

save(fn, 'syn_duty_1ab', 'syn_duty_vel_1ab')
