%% Produce synthetic duty cycle summary tables and package them for use 
% in cyclers

% Kevin Moy, 12/15/21
clearvars
close all
clc
set(0,'defaultTextInterpreter','latex');

sf_dat = readtable("Confidential - Not subject to FOIA - 202201_Urban_CellPowerProfile.csv");
lv_dat = readtable("Confidential - Not subject to FOIA - 202201_Las Vegas Urban+_CellPowerProfile.csv");
Q_nom = 5; % Ah... from Zoox slides??
freq = 10; % Hz

% Assume current from front pack
sf_curr = sf_dat.FrontCurrent_A_/Q_nom;
lv_curr = lv_dat.FrontCurrent_A_/Q_nom;

% Get velocity
sf_vel = sf_dat.Velocity_m_s_;
lv_vel = lv_dat.Velocity_m_s_;
%% Load all characteristic duty cycles

% SF Synthetic Duty Cycles
load clust_char_SF.mat
load clust_char_SF_vel.mat

% LV Synthetic Duty Cycles
load clust_char_LV.mat
load clust_char_LV_vel.mat

% SF + LV Synthetic Duty Cycles
load clust_char_SF+LV.mat
load clust_char_SF+LV_vel.mat

% Highway Synthetic Duty Cycles
load hwy_syn_duty_cycles.mat
load hwy_syn_duty_cycles_vel.mat

% Highway Duty Cycles
load hwy_drive_cycle.mat
load hwy_drive_cycle_vel.mat

%% COMBINE ALL DUTY CYCLES INTO SYNTHETIC DUTY CYCLES

syn_duty_cycle_2a = [clust_char_SF{1}; clust_char_SF{2}];
syn_duty_cycle_2a_vel = [clust_vel_SF{1}; clust_vel_SF{2}];

syn_duty_cycle_2b = [clust_char_LV{1}; clust_char_LV{2}];
syn_duty_cycle_2b_vel = [clust_vel_LV{1}; clust_vel_LV{2}];

syn_duty_cycle_2c = [clust_char_SF_LV{1}; clust_char_SF_LV{2}];
syn_duty_cycle_2c_vel = [clust_vel_SF_LV{1}; clust_vel_SF_LV{2}];

syn_duty_cycle_3 = [clust_char_SF{1}; clust_char_SF{2}; ...
    clust_char_LV{1}; clust_char_LV{2};... 
    syn_duty_1ab{1}];

syn_duty_cycle_3_vel = [clust_vel_SF{1}; clust_vel_SF{2}; ...
    clust_vel_LV{1}; clust_vel_LV{2};...
    syn_duty_vel_1ab{1}];

duty_cycle_metrics_1a = metrics(syn_duty_1ab{1}, syn_duty_vel_1ab{1}, freq);
duty_cycle_metrics_1b = metrics(syn_duty_1ab{2}, syn_duty_vel_1ab{2}, freq);


duty_cycle_metrics_2a = metrics(syn_duty_cycle_2a, syn_duty_cycle_2a_vel, freq);
duty_cycle_metrics_2b = metrics(syn_duty_cycle_2b, syn_duty_cycle_2b_vel, freq);
duty_cycle_metrics_2c = metrics(syn_duty_cycle_2c, syn_duty_cycle_2c_vel, freq);
duty_cycle_metrics_3 = metrics(syn_duty_cycle_3, syn_duty_cycle_3_vel, freq);


%%
T2 = array2table([duty_cycle_metrics_1a, duty_cycle_metrics_1b, ...
    duty_cycle_metrics_2a, duty_cycle_metrics_2b, ...
    duty_cycle_metrics_2c, duty_cycle_metrics_3], ...
    'VariableNames',{'Synthetic Duty Cycle 1a', 'Synthetic Duty Cycle 1b', ...
    'Synthetic Duty Cycle 2a', 'Synthetic Duty Cycle 2b', ...
    'Synthetic Duty Cycle 2c', 'Synthetic Duty Cycle 3'}, ...
    'RowNames',{'Peak Positive C-rate','Peak Negative C-rate', ...
    'Average Positive C-rate', 'Average Negative C-rate', 'Average C-rate', 'Peak Frequency, Positive C-rate [Hz]', ...
    'Peak Frequency, Negative C-rate [Hz]', 'Variance of Positive C-rate', 'Variance of Negative C-rate', ...
    'Distance Travelled [m]', 'Average Velocity [m/s]', 'Peak Velocity [m/s]', 'Variance of velocity [m^2/s^2]', ...
    'Number of start/stops'});

% T2 = array2table([SF_drive_cycle_metrics, LV_drive_cycle_metrics, SF_LV_drive_cycle_metrics, SF_syn_cycle_metrics, LV_syn_cycle_metrics, SF_LV_syn_cycle_metrics, ...
%     hwy_duty_cycle_metrics, SF_all_syn_cycle_metrics], ...
%     'VariableNames',{'SF Drive Cycle', 'LV Drive Cycle', ...
%     'SF+LV Drive Cycle', 'SF Synthetic Duty Cycle', 'LV synthetic Duty Cycle', ...
%     'SF+LV Synthetic Duty Cycle', 'Highway Synthetic Duty Cycle 1a', ...
%     'Highway Synthetic Duty Cycle 1b', 'Combined Synthetic Duty Cycle'}, ...
%     'RowNames',{'Peak Positive C-rate','Peak Negative C-rate', ...
%     'Average Positive C-rate', 'Average Negative C-rate', 'Average C-rate', 'Peak Frequency, Positive C-rate [Hz]', ...
%     'Peak Frequency, Negative C-rate [Hz]', 'Variance of Positive C-rate', 'Variance of Negative C-rate', ...
%     'Distance Travelled [m]', 'Average Velocity [m/s]', 'Peak Velocity [m/s]', 'Variance of velocity [m^2/s^2]', ...
%     'Number of start/stops'});
%%
writetable(T2, 'Summary_Syn_Duty_Cycles_v3.csv','WriteRowNames',true)   

%% TODO: PLOT ALL SYNTHETIC DUTY CYCLES!!
% Highway
hFig = figure(1);
set(hFig, 'Position', [100 100 1000 600])
tiledlayout(1,2);
nexttile
plot((1:length(syn_duty_1ab{1}))./freq/60, syn_duty_1ab{1}, 'LineWidth', 1)
xlim([0, length(syn_duty_1ab{2})/freq/60])
xlabel('Time [min]')
ylabel('C-rate [-]')
title('Syn. Duty Cycle 1a')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
nexttile
plot((1:length(syn_duty_1ab{2}))./freq/60, syn_duty_1ab{2}, 'LineWidth', 1)
xlim([0, length(syn_duty_1ab{2})/freq/60])
xlabel('Time [min]')
ylabel('C-rate [-]')
title('Syn. Duty Cycle 1b')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
t = sgtitle('Highway Duty Cycles');
t.Interpreter = 'latex';
t.FontSize = 24;

% Non-highway
hFig = figure(2);
set(hFig, 'Position', [100 100 1000 600])
tiledlayout(1,2);
nexttile
plot((1:length(syn_duty_cycle_2a))./freq/60, syn_duty_cycle_2a, 'LineWidth', 1)
xlim([0, length(syn_duty_cycle_2a)/freq/60])
xlabel('Time [min]')
ylabel('C-rate [-]')
title('Syn. Duty Cycle 2a')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
nexttile
plot((1:length(syn_duty_cycle_2b))./freq/60, syn_duty_cycle_2b, 'LineWidth', 1)
xlim([0, length(syn_duty_cycle_2b)/freq/60])
xlabel('Time [min]')
ylabel('C-rate [-]')
title('Syn. Duty Cycle 2b')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
t = sgtitle('Non-Highway Duty Cycles');
t.Interpreter = 'latex';
t.FontSize = 24;

hFig = figure();
set(hFig, 'Position', [100 100 600 500])
plot((1:length(syn_duty_cycle_2a))./freq/60, syn_duty_cycle_2a*2, 'color', rgb('pine green'), 'LineWidth', 2)
title({'City 1 Synthetic Duty Cycle,' 'C/5 Average C-Rate'})
xlabel('Time [min]')
ylabel('Dispatch C-rate [-]')
set(gca,'FontSize', 26)

hFig = figure();
set(hFig, 'Position', [100 100 600 500])
plot((1:length(syn_duty_cycle_2b))./freq/60, syn_duty_cycle_2b*5, 'color', rgb('pine green'), 'LineWidth', 2)
title({'City 2 Synthetic Duty Cycle,' 'C/2 Average C-Rate'})
xlabel('Time [min]')
ylabel('Dispatch C-rate [-]')
set(gca,'FontSize', 26)

% hFig = figure(2);
% set(hFig, 'Position', [100 100 1500 600])
% tiledlayout(1,3);
% nexttile
% plot((1:length(syn_duty_cycle_SF))./freq/60, syn_duty_cycle_SF, 'LineWidth', 1)
% xlim([0, length(syn_duty_cycle_SF)/freq/60])
% xlabel('Time [min]')
% ylabel('C-rate [-]')
% title('Syn. Duty Cycle 2a')
% set(gca, 'TickLabelInterpreter','latex')
% set(gca,'FontSize', 24)
% nexttile
% plot((1:length(syn_duty_cycle_LV))./freq/60, syn_duty_cycle_LV, 'LineWidth', 1)
% xlim([0, length(syn_duty_cycle_LV)/freq/60])
% xlabel('Time [min]')
% ylabel('C-rate [-]')
% title('Syn. Duty Cycle 2b')
% set(gca, 'TickLabelInterpreter','latex')
% set(gca,'FontSize', 24)
% nexttile
% plot((1:length(syn_duty_cycle_SF_LV))./freq/60, syn_duty_cycle_SF_LV, 'LineWidth', 1)
% xlim([0, length(syn_duty_cycle_SF_LV)/freq/60])
% xlabel('Time [min]')
% ylabel('C-rate [-]')
% title('Syn. Duty Cycle 2c')
% set(gca, 'TickLabelInterpreter','latex')
% set(gca,'FontSize', 24)
% t = sgtitle('Non-Highway Duty Cycles');
% t.Interpreter = 'latex';
% t.FontSize = 24;

hFig = figure(3);
set(hFig, 'Position', [100 100 800 600])
tiledlayout(1,1);
nexttile
plot((1:length(syn_duty_cycle_3))./freq/60, syn_duty_cycle_3, 'LineWidth', 1)
xlim([0, length(syn_duty_cycle_3)/freq/60])
xlabel('Time [min]')
ylabel('C-rate [-]')
title('Syn. Duty Cycle 3')
set(gca, 'TickLabelInterpreter','latex')
set(gca,'FontSize', 24)
t = sgtitle('Combined Highway + Non-Highway Duty Cycle');
t.Interpreter = 'latex';
t.FontSize = 24;

%% Boxplot of summary table
summ_arr = table2array(T2)';
hFig = figure(4); 
set(hFig, 'Position', [100 100 1000 600])
h = boxplot(summ_arr./max(summ_arr, [], 1), 'Labels',{'Peak Pos. C-rate','Peak Neg. C-rate', ...
    'Avg. Pos. C-rate', 'Avg. Neg. C-rate', 'Avg. C-rate', 'Peak Freq., Pos. C-rate [Hz]', ...
    'Peak Freq., Neg. C-rate [Hz]', 'Variance of Pos. C-rate', 'Variance of Neg. C-rate', ...
    'Distance Travelled [m]', 'Avg. Velocity [m/s]', 'Peak Velocity [m/s]', 'Variance of velocity [m$^2$/s$^2$]', ...
    'Number of start/stops'});
set(gca,'FontSize', 18);
set(gca,'TickLabelInterpreter','latex');
ylim([0 1.05])
xlabel('Metric', 'FontSize', 24)
ylabel('Normalized Quantity', 'FontSize', 24)
set(h,{'linew'},{1})

%% Barplot of summary table
summ_arr2 = table2array(T2);
hFig = figure(5); 
set(hFig, 'Position', [100 100 1500 700])
X = categorical({'Peak Pos. C-rate','Peak Neg. C-rate', ...
    'Avg. Pos. C-rate', 'Avg. Neg. C-rate', 'Avg. C-rate', 'Peak Freq., Pos. C-rate [Hz]', ...
    'Peak Freq., Neg. C-rate [Hz]', 'Variance of Pos. C-rate', 'Variance of Neg. C-rate', ...
    'Distance Travelled [m]', 'Avg. Velocity [m/s]', 'Peak Velocity [m/s]', 'Variance of velocity [m$^2$/s$^2$]', ...
    'Number of start/stops'});
X = reordercats(X,{'Peak Pos. C-rate','Peak Neg. C-rate', ...
    'Avg. Pos. C-rate', 'Avg. Neg. C-rate', 'Avg. C-rate', 'Peak Freq., Pos. C-rate [Hz]', ...
    'Peak Freq., Neg. C-rate [Hz]', 'Variance of Pos. C-rate', 'Variance of Neg. C-rate', ...
    'Distance Travelled [m]', 'Avg. Velocity [m/s]', 'Peak Velocity [m/s]', 'Variance of velocity [m$^2$/s$^2$]', ...
    'Number of start/stops'});
h = bar(X, summ_arr2./max(summ_arr2, [], 2));
set(gca,'FontSize', 18);
set(gca,'TickLabelInterpreter','latex');
ylim([0 1.05])
xlabel('Metric', 'FontSize', 24)
ylabel('Normalized Quantity', 'FontSize', 24)
legend('SF Drive Cycle', 'LV Drive Cycle', 'SF+LV Drive Cycle', ...
    'Syn. Duty Cycle 1a', 'Syn. Duty Cycle 1b', 'Syn. Duty Cycle 2a', ...
    'Syn. Duty Cycle 2b', 'Syn. Duty Cycle 2c', 'Syn. Duty Cycle 3', 'Location', 'eastoutside', 'Interpreter', 'latex')
% set(h,{'linew'},{1})



%% Plots

reps_lv1 = floor(length(drive_cycle)/length(SF_lv1));

hFig = figure(1); 
set(hFig, 'Position', [100 100 1000 1000])
set(0,'DefaultAxesLineStyleOrder',{'-',':'});
tiledlayout(2,2);

nexttile
hold on
yyaxis left
plot(repmat(SF_lv1, reps_lv1, 1), 'LineWidth', 1)
ylabel('Current [A]')
ylim([-3 4.5])
yyaxis right
plot((Q_nom + cumtrapz(-repmat(SF_lv1, reps_lv1, 1))*0.1/3600)/Q_nom, 'LineWidth', 2)
ylabel('SOC [-]')
hold off
xlabel('Time [s]')
xlim([0 length(repmat(SF_lv1, reps_lv1, 1))])
ylim([0 1])
title('Level 1')
set(gca,'FontSize', 20)

nexttile
hold on
yyaxis left
plot(repmat(SF_lv2, reps, 1), 'LineWidth', 1)
ylabel('Current [A]')
ylim([-3 4.5])
yyaxis right
plot((Q_nom + cumtrapz(-repmat(SF_lv2, reps, 1))*0.1/3600)/Q_nom, 'LineWidth', 2)
ylabel('SOC [-]')
hold off
xlabel('Time [s]')
xlim([0 length(repmat(SF_lv1, reps_lv1, 1))])
ylim([0 1])
title('Level 2')
set(gca,'FontSize', 20)

nexttile
hold on
yyaxis left
plot(SF_lv2_rests, 'LineWidth', 1)
ylabel('Current [A]')
ylim([-3 4.5])
yyaxis right
plot((Q_nom + cumtrapz(-SF_lv2_rests)*0.1/3600)/Q_nom, 'LineWidth', 2)
ylabel('SOC [-]')
hold off
xlabel('Time [s]')
xlim([0 length(repmat(SF_lv1, reps_lv1, 1))])
ylim([0 1])
title('Level 2 with Rests')
set(gca,'FontSize', 20)

nexttile
hold on
yyaxis left
plot(drive_cycle, 'LineWidth', 1)
ylabel('Current [A]')
ylim([-3 4.5])
yyaxis right
plot((Q_nom + cumtrapz(-drive_cycle)*0.1/3600)/Q_nom, 'LineWidth', 2)
ylabel('SOC [-]')
hold off
xlabel('Time [s]')
xlim([0 length(repmat(SF_lv1, reps_lv1, 1))])
ylim([0 1])
title('Full Drive Cycle')
set(gca,'FontSize', 20)

(Q_nom + trapz(-repmat(SF_lv1, reps_lv1, 1))*0.1/3600)/Q_nom
(Q_nom + trapz(-repmat(SF_lv2, reps, 1))*0.1/3600)/Q_nom
(Q_nom + trapz(-SF_lv2_rests)*0.1/3600)/Q_nom
(Q_nom + trapz(-drive_cycle)*0.1/3600)/Q_nom

hFig = figure(2); 
set(hFig, 'Position', [100 100 1000 1000])
set(0,'DefaultAxesLineStyleOrder',{'-',':'});
tiledlayout(2,2);

nexttile
hold on
plot(SF_lv1, 'LineWidth', 1)
ylabel('Current [A]')
ylim([-3 4.5])
hold off
xlabel('Time [s]')
xlim([0 length(SF_lv1)])
title('Level 1')
set(gca,'FontSize', 20)

nexttile
hold on
plot(SF_lv2, 'LineWidth', 1)
ylabel('Current [A]')
ylim([-3 4.5])
hold off
xlabel('Time [s]')
xlim([0 length(SF_lv2)])
title('Level 2')
set(gca,'FontSize', 20)

nexttile
hold on
plot(SF_lv2_rests, 'LineWidth', 1)
ylabel('Current [A]')
ylim([-3 4.5])
hold off
xlabel('Time [s]')
xlim([0 length(SF_lv2_rests)])
title('Level 2 with Rests')
set(gca,'FontSize', 20)


%% Collect all synthetic duty cycles at each C-rate

C_16_avg_current = Q_nom/16;
C_8_avg_current = Q_nom/8;
C_2_avg_current = Q_nom/2;

SF_lv1_C_16 =  C_16_avg_current/mean(SF_lv1) * SF_lv1;
SF_lv1_C_8 =  C_8_avg_current/mean(SF_lv1) * SF_lv1;
SF_lv1_C_2 =  C_2_avg_current/mean(SF_lv1) * SF_lv1;

SF_lv2_C_16 =  C_16_avg_current/mean(SF_lv2) * SF_lv2;
SF_lv2_C_8 =  C_8_avg_current/mean(SF_lv2) * SF_lv2;
SF_lv2_C_2 =  C_2_avg_current/mean(SF_lv2) * SF_lv2;

SF_lv2_rests_C_16 =  C_16_avg_current/mean(SF_lv2_rests) * SF_lv2_rests;
SF_lv2_rests_C_8 =  C_8_avg_current/mean(SF_lv2_rests) * SF_lv2_rests;
SF_lv2_rests_C_2 =  C_2_avg_current/mean(SF_lv2_rests) * SF_lv2_rests;

LV_lv1_C_16 =  C_16_avg_current/mean(LV_lv1) * LV_lv1;
LV_lv1_C_8 =  C_8_avg_current/mean(LV_lv1) * LV_lv1;
LV_lv1_C_2 =  C_2_avg_current/mean(LV_lv1) * LV_lv1;

hFig = figure(2); 
set(hFig, 'Position', [100 100 800 800])
set(0,'DefaultAxesLineStyleOrder',{'-',':'});
tiledlayout(3,1);

nexttile
hold on
plot(SF_lv1_C_2, 'LineWidth', 1)
plot(SF_lv1_C_8, 'LineWidth', 1)
plot(SF_lv1_C_16, 'LineWidth', 1)
ylabel('Current [A]')
hold off
xlabel('Time [s]')
xlim([0 length(SF_lv1)])
title('Level 1')
set(gca,'FontSize', 20)

nexttile
hold on
plot(SF_lv2_C_2, 'LineWidth', 1)
plot(SF_lv2_C_8, 'LineWidth', 1)
plot(SF_lv2_C_16, 'LineWidth', 1)
ylabel('Current [A]')
hold off
xlabel('Time [s]')
xlim([0 length(SF_lv2)])
title('Level 2')
set(gca,'FontSize', 20)

nexttile
hold on
plot(SF_lv2_rests_C_2, 'LineWidth', 1)
plot(SF_lv2_rests_C_8, 'LineWidth', 1)
plot(SF_lv2_rests_C_16, 'LineWidth', 1)
ylabel('Current [A]')
hold off
xlabel('Time [s]')
xlim([0 length(SF_lv2_rests)])
title('Level 2 with Rests')
set(gca,'FontSize', 20)

t = sgtitle('SF Duty Cycles');
t.Interpreter = 'latex';
t.FontSize = 24;

lg  = legend('C/2', 'C/8', 'C/16','Orientation','Vertical'); 
lg.Layout.Tile = 'East'; % <-- Legend placement with tiled layout

hFig = figure(3); 
set(hFig, 'Position', [100 100 600 500])
set(0,'DefaultAxesLineStyleOrder',{'-',':'});
hold on
plot(LV_lv1_C_2, 'LineWidth', 1)
plot(LV_lv1_C_8, 'LineWidth', 1)
plot(LV_lv1_C_16, 'LineWidth', 1)
ylabel('Current [A]')
hold off
xlabel('Time [s]')
xlim([0 length(SF_lv1)])
title('LV Duty Cycle')
set(gca,'FontSize', 20)
lg  = legend('C/2', 'C/8', 'C/16','Orientation','Vertical', 'Location', 'eastoutside'); 

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
