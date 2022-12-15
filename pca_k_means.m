%% Function for PCA + k-means Characteristic Duty Cycle Algorithm
% Kevin Moy
% See publications for details:
% https://doi.org/10.1016/j.adapen.2021.100065
% https://doi.org/10.1115/1.4050192

function [char_disp, char_state, char_aux, char_interval] = pca_k_means(metrics, dispatch, state, auxdata, showplots)
    % Inputs: 
    % metrics -- Interval dispatch metrics matrix, as 2d matrix (number of intervals x number of metrics)
    % dispatch -- intervals of power or C-rate or current of dispatch, as cell array
    % state -- intervals of SOC or SOE of dispatch, as cell array
    % auxdata -- intervals e.g. temperature, velocity, etc., as cell array
    % showplots -- boolean of true/false whether to produce plots as well

    % Outputs:
    % char_disp -- characteristic duty cycles from dispatch, as cell array
    % char_state -- corresponding state for char. duty cycles, as cell array
    % char_aux -- corresponding auxiliary data for char. duty cycles, as cell array
    % char_intervals -- intervals that correspond to the characteristic
    % duty cycles

    %% Principal Component Analysis
    % Normalize metrics (columns)
    % metrics_normalized = metrics./(std(metrics, 1, 1, 'omitnan'));
    metrics_normalized = metrics./(max(abs(metrics)));
    
    % PCA on normalized metrics matrix
    [coeff,~,~,~,explained,~] = pca(metrics_normalized, 'Centered', 'on');
    
    % Determine number of principal components to use
    % (minimum number needed to explain 90% of variance)
    num_comp = find(cumsum(explained) > 90, 1);
    
    % Map to reduced-order PCA space
    reduce = coeff(:,1:num_comp); % reduced-order PCA space
    reconst_daily = metrics_normalized*reduce;
    
    %% k-means Clustering
    % Determine cluster number
    % Maximum k to try
    maxk = 30;
    
    % Use 'Replicates' tag, which runs kmeans n times and keeps instance with
    % lowest sumd (= Within-cluster sums of point-to-centroid distances)
    avg_sumd = [];
    sumd_sums = [];
    clust_dist = [];
    pdist_clust = [1, 0];
    for d = 1:maxk
        [idx, clust, sumd] = kmeans(reconst_daily, d, 'Replicates', 100);
        num_in_clust = [];
        idxj = [];
        for j = 1:d % for each cluster
            idxj(j) = length(find(idx == j));
        end
        avg_sumd = [avg_sumd; repmat(d,d,1), sumd./idxj']; % Average distance from cluster member to centroid
        sumd_sums(d) = max(sumd./idxj');
        if d > 1 % Pairwise distances need at least 2 clusters
            pdists = pdist(clust);
            clust_dist(d) = min(pdists);
            pdist_clust = [pdist_clust; repmat(d, length(pdists),1) , pdists'];
        end
    end
    
    clust_dist_min = [];
    for d = 1:maxk
        idx = pdist_clust(:,1) == d;
        clust_dist_min(d) = min(pdist_clust(idx,2));
    end
    
    sumd_max = [];
    for d = 1:maxk
        idx = avg_sumd(:,1) == d;
        sumd_max(d) = max(avg_sumd(idx,2));
    end
    
    % Goal: Keep intracluster distances small, 
    % while intercluster distances are large (close clusters, far apart)
    % k_metric = (sumd_max(2:end)./clust_dist_min(2:end));
    k_metric = clust_dist_min(2:end)-sumd_max(2:end);
    [~, k_opt] = max(k_metric); 
    k_opt = k_opt + 1;
    num_clust = k_opt;
    
    % Given optimal number of clusters num_clust, obtain clusters and
    % centroids
    [idx_cl, clusters] = kmeans(reconst_daily, num_clust, 'Replicates', 100);
    
    % Allocate arrays
    char_disp = cell(num_clust,1); % initialize cell array for characteristic duty cycles
    char_state = cell(num_clust,1); % initialize cell array for characteristic duty cycle states
    char_aux = cell(num_clust,1); % initialize cell array for characteristic duty cycle aux data
    char_interval = zeros(num_clust,1); % days closest to cluster centroids in PCA subspace
    
    % Collect all characteristic duty cycle data
    for i = 1:num_clust
        [~, char_interval(i)] = min(pdist2(reconst_daily, clusters(i,:)));
        char_disp{i} = dispatch{char_interval(i)};
        char_state{i} = state{char_interval(i)};
        char_aux{i} = auxdata{char_interval(i)};
    end


%% Plot figures (if showplots is true)
    if showplots
    
        hFig = figure();
        set(hFig, 'Position', [100 100 600 800])
        h = heatmap(metrics_normalized, 'ColorMap', parula, ...
        'XLabel', 'Metric', 'YLabel', 'Interval');
        set(gca,'FontSize', 20)
    
        hFig = figure();
        set(hFig, 'Position', [100 100 600 800])
        h = heatmap(metrics, 'ColorMap', parula, ...
        'XLabel', 'Metric', 'YLabel', 'Interval');
        set(gca,'FontSize', 20)

        figure()
        heatmap(coeff)
        set(gca,'FontSize', 20)
        xlabel('Principal Components')
        ylabel('Metrics')
        
        hFig = figure();
        set(hFig, 'Position', [100 100 600 500])
        hold on
        plot(cumsum(explained)/100, 'r.', 'MarkerSize', 30)
        plot(0:length(explained), 0.9*ones(1+length(explained),1), 'k-.', 'LineWidth', 2)
        set(gca,'FontSize', 20)
        xlabel('Number of Principal Components')
        ylabel('Fraction of Retained Variation')
        ylim([min(cumsum(explained)/100) 1])
        xlim([1, length(explained)])
        box on
       
        hFig = figure();
        set(hFig, 'Position', [100 100 600 600])
        h = heatmap(reconst_daily, 'ColorMap', parula, ...
        'XLabel', 'Principal Component', 'YLabel', 'Interval');
        set(gca,'FontSize', 26)
        
        C = cov(metrics_normalized, 'omitrows');
        figure()
        heatmap(C)
        
        hFig = figure();
        set(hFig, 'Position', [100 100 600 600])
        h = heatmap(reduce, 'ColorMap', parula, ...
        'XLabel', 'Principal Component', 'YLabel', 'Metric', 'CellLabelColor','none');
        set(gca,'FontSize', 26)

        hFig = figure();
        set(hFig, 'Position', [100 100 600 200])
        h = heatmap(string(1:size(reconst_daily,2)), string(1:num_clust), clusters, 'ColorMap', parula, ...
            'XLabel', 'Principal Component', 'YLabel', 'Cluster', 'CellLabelColor','none');
        set(gca,'FontSize', 26)

        hFig = figure();
        set(hFig, 'Position', [100 100 600 600])
        hold on
        plot(2:maxk,k_metric, 'color', rgb('magenta'), 'LineWidth', 2)
        set(gca,'FontSize', 26)
        xlabel('Cluster Number, k')
        ylabel('$dist_{between}(k) - dist_{within}(k)$, [-]')
        ylim([0 max(k_metric)])
        box on
        hold off

        % Plot idx for non-zero dispatch days
        hFig = figure(21);
        set(hFig, 'Position', [100 100 600 500])
        plot(idx_cl, 'k.', 'MarkerSize',20)
        set(gca,'FontSize', 26)
        title(['Cluster Assignment, $N_c$ = ' ,num2str(num_clust)])
        xlabel('Interval')
        ylabel('Cluster')
        ylim([0,num_clust+1])
        yticks(1:num_clust)

    end
end