function [point_idx, loc] = assignPointIndexSimple(loc, tim, time_threshold)

% labels = dbscan(loc, 5, 5);
% valid_indices = labels ~= -1;
% loc1=loc
% loc = loc(valid_indices, :);
% tim = tim(valid_indices);
% valid_labels=labels(valid_indices);
% % 计算聚类的数量
% num_clusters = max(valid_labels);
% colors = rand(num_clusters, 3);
% scatter_colors = colors(valid_labels, :);
% figure
% scatter(loc1(:,1)*3, loc1(:,2)*3, 20, 'k', 'filled');
% hold on
% scatter(loc(:,1)*3, loc (:,2)*3, 15, scatter_colors, 'filled');

    % 参数检查
    % if nargin < 3 || isempty(time_threshold)
    %     time_threshold = 0.004;
    % end
    
    n = length(tim);
    
    % 验证时间向量是否已排序
    if ~issorted(tim)
        warning('时间向量未排序！建议先对数据进行时间排序。');
        % 可以选择在此处添加排序逻辑，或者继续处理
    end
    
    % 计算连续点的时间差
    time_diffs = diff(tim);
    
    % 找到时间差大于阈值的位置（新组的开始）
    group_starts = [true; time_diffs > time_threshold];
    
    % 使用cumsum分配索引（一行代码解决）
    point_idx = cumsum(group_starts);
    
    % 计算分组信息
    group_info = calculateGroupInfoSimple(point_idx, loc, tim, time_threshold);
    
   
end

function group_info = calculateGroupInfoSimple(point_idx, loc, tim, time_threshold)
% 计算分组统计信息（简化版）
    unique_ids = unique(point_idx);
    num_groups = length(unique_ids);
    
    group_info.num_groups = num_groups;
    group_info.group_ids = unique_ids;
    group_info.group_sizes = zeros(num_groups, 1);
    group_info.group_time_ranges = zeros(num_groups, 2);
    group_info.group_durations = zeros(num_groups, 1);
    
    for i = 1:num_groups
        group_id = unique_ids(i);
        mask = (point_idx == group_id);
        
        group_info.group_sizes(i) = sum(mask);
        
        group_times = tim(mask);
        group_info.group_time_ranges(i, :) = [min(group_times), max(group_times)];
        group_info.group_durations(i) = max(group_times) - min(group_times);
    end
    
    group_info.time_threshold = time_threshold;
    group_info.total_points = length(tim);

% 计算聚类的数量
num_clusters = max(point_idx);
colors = rand(num_clusters, 3);

% 创建颜色映射：点数≤3的cluster用灰色，其他的用随机颜色
scatter_colors = zeros(length(point_idx), 3);
for i = 1:num_clusters
    cluster_mask = (point_idx == i);
    cluster_size = sum(cluster_mask);
    if cluster_size <= 3
        % 点数≤3的cluster用灰色
        scatter_colors(cluster_mask, :) = repmat([0.9, 0.9, 0.9], cluster_size, 1);
    else
        % 点数>3的cluster用随机颜色
        scatter_colors(cluster_mask, :) = repmat(colors(i,:), cluster_size, 1);
    end
end

% figure
scatter(loc(:,1)*3, loc(:,2)*3, 15, scatter_colors, 'filled');
set(gcf,'color','white'); %窗口背景白色
colordef white; %2D/3D图背景黑色
grid on;
axis equal tight; 
xlabel('X (nm)'); ylabel('Y (nm)');
% 保存为.fig文件
savefig('frames_all.fig');

%% 更新后的绘图
% 确保保存目录存在
if ~exist('frames', 'dir')
    mkdir('frames');
end
if ~exist('density', 'dir')
    mkdir('density');
end

% 第一组图：每个cluster单独绘制，但只保存点数>3的cluster
for i = 1:num_clusters
    % 检查cluster大小，点数≤3的跳过不保存
    cluster_mask = (point_idx == i);
    cluster_size = sum(cluster_mask);
    
    if cluster_size <= 3
        continue; % 跳过点数≤3的cluster
    end
    
    % 获取当前cluster对应的时间索引
    current_time_indices = find(cluster_mask);
    
    % 计算时间范围（小数点后三位）
    if ~isempty(current_time_indices)
        min_time = min(tim(current_time_indices));
        max_time = max(tim(current_time_indices));
        time_range_str = sprintf('Time range: %.3f - %.3f', min_time, max_time);
    else
        time_range_str = 'Time range: No data';
    end
    
    % % 绘制所有点（其他点用浅灰色）
    scatter(loc(~cluster_mask,1)*3, loc(~cluster_mask,2)*3, 15, [0.8 0.8 0.8], 'filled');
    hold on;

    % 高亮当前cluster的点，使用相同的颜色
    scatter(loc(cluster_mask,1)*3, loc(cluster_mask,2)*3, 15, colors(i,:), 'filled');

    set(gcf,'color','white');
    colordef white;
    grid on;
    xlabel('X (nm)'); ylabel('Y (nm)');
    % 设置标题，只有一个换行符
    title(sprintf('Frame %d\n%s (s)', i, time_range_str), 'Interpreter', 'none');

    axis equal tight; 

    % 保存单个cluster图
    savefig(sprintf('frames/frame%d.fig', i));
    close;
end

% 第二组图：每个cluster单独计算密度并绘制，只处理点数>3的cluster
points_nm = loc * 3; % 将原始loc数据乘以3
R = 15;

% 为每个cluster单独计算和绘制密度图
all_avg_densities=[];
for i = 1:num_clusters
    % 检查cluster大小，点数≤3的跳过不保存
    cluster_mask = (point_idx == i);
    cluster_size = sum(cluster_mask);
    
    if cluster_size <= 3
        continue; % 跳过点数≤3的cluster
    end
    
    % 获取当前cluster的点
    current_points = points_nm(cluster_mask, :);
    num_current_points = size(current_points, 1);
    
    % 为当前cluster单独计算密度
    current_densities = zeros(num_current_points, 1);
    
    for j = 1:num_current_points
        distances = sqrt(sum((current_points - current_points(j,:)).^2, 2));
        Ninside = sum(distances <= R);
        current_densities(j) = Ninside / (pi * R^2);
    end
    
    % 计算平均密度
    avg_density = mean(current_densities);
    all_avg_densities = [all_avg_densities; avg_density]; 

    
    scatter(current_points(:,1), current_points(:,2), 15, current_densities, 'filled');
    colorbar;
    set(gcf,'color','white');
    colordef white;
    colormap(jet); % 设置colormap为jet
    clim([0 0.010]); % 设置显色范围为0到0.005
    title(sprintf('Local Density - Frame %d\nAverage Density: %.3f', i, avg_density));
    axis equal tight; 
    grid on;
    xlabel('X (nm)'); ylabel('Y (nm)');

    savefig(sprintf('density/frame%d_density.fig', i));
    close;
end

save('density_data.mat','all_avg_densities')
end
