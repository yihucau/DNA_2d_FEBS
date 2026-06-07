function createCombinedFramePlot12_5_timeOrder()
    %% ==================== 参数设置 ====================
    NN_THRESHOLD = 5;          % 最近邻距离阈值（nm）
    SPLINE_SAMPLES = 800;       % 样条采样密度
    MIN_POINTS = 5;             % 最少点数要求（与旧代码一致）5
    
    % 点大小设置
    INDIVIDUAL_POINT_SIZE = 40;  % 单帧图的点大小
    COMBINED_POINT_SIZE = 30;    % 总图的点大小

    indi_cylinder_size=0.4;
   

    %% ==================== 主函数 ====================
    % Load .mat file containing multiple frames
    [filename, pathname] = uigetfile('*.mat', 'Select .mat file containing frame data');
    if isequal(filename, 0)
        fprintf('User cancelled file selection\n');
        return;
    end
    
    filepath = fullfile(pathname, filename);
    data = load(filepath);
    
    % Get all variable names starting with 'frame'
    var_names = fieldnames(data);
    frame_vars = var_names(startsWith(var_names, 'frame'));
    
    if isempty(frame_vars)
        fprintf('No variables starting with "frame" found in the file\n');
        return;
    end
    
    fprintf('Found %d frame datasets:\n', length(frame_vars));
    
    % Create directories
    model_dir = 'model';
    if ~exist(model_dir,'dir'), mkdir(model_dir); end
    output_dir = 'combined_analysis_3d';
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    
    % Create main 3D figure
    fig_3d = figure('Position',[100,100,1200,900],'Name','3D Chromatin Conformation','Visible','off');
    
    % Initialize summary
    summary_data = []; % [frame_number, total_points, red_points, blue_points, compaction]
    
    % Initialize min_distance collection
    all_min_distances = [];

    %% ==================== 逐帧处理 ====================
    for frame_idx = 1:length(frame_vars)
        current_frame = data.(frame_vars{frame_idx});
        frame_number = str2double(regexp(frame_vars{frame_idx}, '\d+', 'match'));
        if isempty(frame_number), frame_number = frame_idx; end
        
        if size(current_frame,2) >= 4
            % x, y, z, time/order
            scale_xy = 3; scale_z = 1;
            a_2d = current_frame(:,1:2) * scale_xy;
            a_3d = [a_2d, current_frame(:,3)*scale_z, current_frame(:,4)];
            num_points = size(a_2d,1);
            
            if num_points < MIN_POINTS
                fprintf('  Skipping frame %d: insufficient points\n', frame_number);
                continue;
            end
            
            % Classify points by nearest neighbor and collect min_distances
            [point_colors, frame_min_distances] = classifyPointsByNN(a_2d, NN_THRESHOLD, frame_number);
            
            % Add frame_min_distances to all_min_distances
            if ~isempty(frame_min_distances)
                frame_info = [ones(size(frame_min_distances)) * frame_number, frame_min_distances];
                all_min_distances = [all_min_distances; frame_info];
            end
            
            % Count red and blue points
            red_points = sum(point_colors == 1);
            blue_points = sum(point_colors == 0);
            compaction = red_points / num_points;
            
            % Individual figure - 使用大点
            fig_individual = createIndividualFrameFigure3D(a_3d, frame_number, compaction);
            drawTimeOrderedPath_individual(a_3d, point_colors, INDIVIDUAL_POINT_SIZE,indi_cylinder_size);
            finalizeIndividualFigure3D(fig_individual, frame_number, compaction, num_points);
            saveIndividualFigure3D(fig_individual, frame_number, model_dir);
            close(fig_individual);
            
            % Plot on main figure - 使用小点
            figure(fig_3d);
            drawTimeOrderedPath_combined(a_3d, point_colors, COMBINED_POINT_SIZE);
            % 设置背景淡灰色
             set(gcf,'color','black'); %窗口背景白色
             colordef black; %2D/3D图背景黑色
            set(gca, ...
                'Box', 'on', ...
                'TickDir', 'in', ...
                'LineWidth', 0.5, ...
                'FontSize', 12, ...
                'GridLineStyle', '--');
            % 开启网格
           grid on;
           set(gca, 'GridColor', [0.2 0.2 0.2]);  % 白色网格 

            % 调整网格为虚线
            set(gca, 'GridLineStyle', '--');
            
            % Summary
            summary_data = [summary_data; frame_number, num_points, red_points, blue_points, compaction];
            fprintf('  Frame %d: Points=%d, Red Points=%d (%.2f%%), Blue Points=%d (%.2f%%), Compaction=%.3f\n', ...
                frame_number, num_points, red_points, 100*red_points/num_points, blue_points, 100*blue_points/num_points, compaction);
        else
            fprintf('  Skipping frame %d: insufficient columns\n', frame_number);
        end
    end
    
    %% ==================== 保存统计数据 ====================
    saveLinearitySummary3D_timeOrder(summary_data, output_dir);
    
    %% ==================== 保存最小距离数据 ====================
    saveMinDistanceData(all_min_distances, output_dir);
    
    % Finalize main figure
    finalize3DFigure(fig_3d);
    save3DFigure(fig_3d);
    close(fig_3d);
    
    fprintf('All individual frame figures saved to %s directory\n', model_dir);
    fprintf('3D analysis summary saved to %s/\n', output_dir);
end

%% ==================== 2D 最近邻点分类 ====================
function [point_colors, min_distances] = classifyPointsByNN(points_nm, threshold, frame_number)
    num_points = size(points_nm, 1);
    point_colors = zeros(num_points, 1); % 0 - blue, 1 - red
    min_distances = zeros(num_points, 1); % 存储每个点的最小距离
    
    for i = 1:num_points
        % 计算2D距离
        distances = sqrt(sum((points_nm - points_nm(i,:)).^2, 2)); 
        
        % 排除自身的距离（将其设为一个大值，避免误判为最近邻）
        distances(i) = Inf;
        
        % 找到最小距离
        min_distance = min(distances);
        min_distances(i) = min_distance; % 记录最小距离
        
        % 判断分类
        if min_distance < threshold
            point_colors(i) = 1; % 红色
        else
            point_colors(i) = 0; % 蓝色
        end
    end
    
    % 打印统计信息
    fprintf('  Frame %d: Min distance stats - Mean=%.3f, Std=%.3f, Min=%.3f, Max=%.3f\n', ...
        frame_number, mean(min_distances), std(min_distances), min(min_distances), max(min_distances));
end

%% ==================== 保存最小距离数据 ====================
function saveMinDistanceData(all_min_distances, output_dir)
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % 保存为MAT文件
    save(fullfile(output_dir, 'all_min_distances.mat'), 'all_min_distances');
    fprintf('Minimum distances saved to: %s\n', fullfile(output_dir, 'all_min_distances.mat'));
    
    % 保存为文本文件
    txt_filename = fullfile(output_dir, 'min_distances_statistics.txt');
    fid = fopen(txt_filename, 'w');
    
    if ~isempty(all_min_distances)
        % 写入标题
        fprintf(fid, 'Frame\tMin_Distance(nm)\n');
        
        % 写入数据
        for i = 1:size(all_min_distances, 1)
            fprintf(fid, '%d\t%.4f\n', all_min_distances(i, 1), all_min_distances(i, 2));
        end
        
        % 写入统计信息
        fprintf(fid, '\n====================================\n');
        fprintf(fid, 'STATISTICAL SUMMARY:\n');
        fprintf(fid, '====================================\n');
        fprintf(fid, 'Total points analyzed: %d\n', size(all_min_distances, 1));
        fprintf(fid, 'Mean minimum distance: %.4f nm\n', mean(all_min_distances(:, 2)));
        fprintf(fid, 'Median minimum distance: %.4f nm\n', median(all_min_distances(:, 2)));
        fprintf(fid, 'Standard deviation: %.4f nm\n', std(all_min_distances(:, 2)));
        fprintf(fid, 'Minimum distance: %.4f nm\n', min(all_min_distances(:, 2)));
        fprintf(fid, 'Maximum distance: %.4f nm\n', max(all_min_distances(:, 2)));
        
        % 按帧统计
        fprintf(fid, '\n====================================\n');
        fprintf(fid, 'PER-FRAME STATISTICS:\n');
        fprintf(fid, '====================================\n');
        
        unique_frames = unique(all_min_distances(:, 1));
        for frame = unique_frames'
            frame_indices = all_min_distances(:, 1) == frame;
            frame_distances = all_min_distances(frame_indices, 2);
            
            fprintf(fid, 'Frame %d:\n', frame);
            fprintf(fid, '  Number of points: %d\n', sum(frame_indices));
            fprintf(fid, '  Mean: %.4f nm\n', mean(frame_distances));
            fprintf(fid, '  Median: %.4f nm\n', median(frame_distances));
            fprintf(fid, '  Std: %.4f nm\n', std(frame_distances));
            fprintf(fid, '  Min: %.4f nm\n', min(frame_distances));
            fprintf(fid, '  Max: %.4f nm\n', max(frame_distances));
            fprintf(fid, '  ---\n');
        end
    else
        fprintf(fid, 'No distance data available\n');
    end
    
    fclose(fid);
    fprintf('Minimum distance statistics saved to: %s\n', txt_filename);
    
    % 创建直方图
    if ~isempty(all_min_distances)
        fig = figure('Position', [100, 100, 800, 600]);
        histogram(all_min_distances(:, 2), 50, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k');
        xlabel('Minimum Distance (nm)');
        ylabel('Frequency');
        title('Distribution of Minimum Distances between Points');
        grid on;
        
        % 在图上添加统计信息
        text_str = sprintf('Mean = %.3f nm\nStd = %.3f nm\nN = %d', ...
            mean(all_min_distances(:, 2)), std(all_min_distances(:, 2)), ...
            size(all_min_distances, 1));
        text(0.7, 0.9, text_str, 'Units', 'normalized', 'FontSize', 10, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
        
        saveas(fig, fullfile(output_dir, 'min_distance_histogram.png'), 'png');
        savefig(fig, fullfile(output_dir, 'min_distance_histogram.fig'));
        close(fig);
        fprintf('Minimum distance histogram saved to: %s\n', fullfile(output_dir, 'min_distance_histogram.png'));
    end
end


%% ==================== 单帧图的时间顺序路径绘制 ====================
function drawTimeOrderedPath_individual(points_nm, point_colors, point_size,cylinder_radius)
    if size(points_nm, 1) < 2
        % 只有一个点就画小球
        colors_rgb = [point_colors, zeros(size(point_colors)), 1-point_colors]; % 1->红(1,0,0), 0->蓝(0,0,1)
        drawSpheres(points_nm, colors_rgb, point_size);  % 使用绘制小球的函数
        return;
    end

    % 按第四列排序（时间顺序）
    [~, time_order] = sort(points_nm(:, 4));
    ordered_points = points_nm(time_order, 1:3);
    ordered_colors = point_colors(time_order);   % 对颜色也按时间排序

    % 绘制每个点的小球
    drawSpheres(ordered_points, ordered_colors, point_size);  % 使用绘制小球的函数

% 添加立体感连线 - 只在有足够点的情况下
    if size(ordered_points, 1) >= 2
        % 构造样条曲线
        scaled_ordered_points=ordered_points;
        scaled_ordered_points(:,3)=scaled_ordered_points(:,3)*5;

        pts = scaled_ordered_points';
        curve = cscvn(pts);
        t_breaks = curve.breaks;          % 参数断点
        n_pts = size(scaled_ordered_points, 1);
        
        % 设置圆柱体连线参数
        % cylinder_radius = 0.4;  % 圆柱体半径
        samples_per_seg = 15;   % 每段采样点数
        
        % 在每个分段上绘制圆柱体
        for k = 1:(n_pts-1)
            % 在当前段参数区间上采样
            t_seg = linspace(t_breaks(k), t_breaks(k+1), samples_per_seg);
            seg_points = fnval(curve, t_seg);   % 3×M
            
            % 决定这一段线的颜色
            c1 = ordered_colors(k);
            c2 = ordered_colors(k+1);
            if c1 == 1 && c2 == 1
                % col = [1, 0.4118, 0.7059];   % 红色
                col = [1, 0, 0];   % 红色

            else
                col = [0, 0.7, 1];   % 蓝色
            end
            
            % 在采样点之间绘制小圆柱体
            for m = 1:(samples_per_seg-1)
                p1 = seg_points(:, m)';
                p2 = seg_points(:, m+1)';
                
                % 计算两点间的距离
                v = p2 - p1;
                L = norm(v);
                
                if L < 0.1
                    continue;  % 跳过太短的段
                end
                
                % 创建圆柱体
                [X, Y, Z] = cylinder(cylinder_radius, 6);  % 6边形圆柱
                Z = Z * L;  % 缩放高度
                
                % 计算旋转矩阵
                z_axis = [0, 0, 1];
                v_norm = v / L;
                
                if abs(dot(z_axis, v_norm)) < 0.999
                    rot_axis = cross(z_axis, v_norm);
                    rot_axis = rot_axis / norm(rot_axis);
                    theta = acos(dot(z_axis, v_norm));
                    
                    % 旋转矩阵
                    K = [0, -rot_axis(3), rot_axis(2);
                         rot_axis(3), 0, -rot_axis(1);
                         -rot_axis(2), rot_axis(1), 0];
                    
                    R = eye(3) + sin(theta)*K + (1-cos(theta))*(K*K);
                    
                    % 旋转点
                    points = R * [X(:)'; Y(:)'; Z(:)'];
                    X_rot = reshape(points(1,:), size(X));
                    Y_rot = reshape(points(2,:), size(Y));
                    Z_rot = reshape(points(3,:), size(Z));
                else
                    X_rot = X; Y_rot = Y; Z_rot = Z;
                end
                
                % 平移到起点
                X_rot = X_rot + p1(1);
                Y_rot = Y_rot + p1(2);
                Z_rot = Z_rot + p1(3);
                
                % 绘制圆柱体
                surf(X_rot, Y_rot, Z_rot, ...
                    'FaceColor', col, ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.9, ...
                    'FaceLighting', 'gouraud');
            end
        end
    end


    % 设置图形属性
    set(gcf,'color','black'); %窗口背景白色
    colordef black; %2D/3D图背景黑色
 set(gca, 'Box', 'on', 'TickDir', 'in', 'LineWidth', 0.8, 'FontSize', 12, 'GridLineStyle', '--');
grid on;
set(gca, 'GridColor', [0.2 0.2 0.2]);  % 白色网格 

% 标签颜色
xlabel('X (nm)', 'Color', [1 1 1]);   
ylabel('Y (nm)', 'Color', [1 1 1]);
zlabel('Z (count)', 'Color', [1 1 1]);


end

%% ==================== 总图的时间顺序路径绘制 ====================
function drawTimeOrderedPath_combined(points_nm, point_colors, point_size)
    if size(points_nm, 1) < 2
        % 只有一个点就画散点
        colors_rgb = [point_colors, zeros(size(point_colors)), 1-point_colors]; % 1->红(1,0,0), 0->蓝(0,0,1)
        scatter3(points_nm(:, 1), points_nm(:, 2), points_nm(:, 3), point_size, colors_rgb, 'filled');
        return;
    end

    % 按第四列排序（时间顺序）
    [~, time_order] = sort(points_nm(:, 4));
    ordered_points = points_nm(time_order, 1:3);
    ordered_colors = point_colors(time_order);   % 对颜色也按时间排序

    % 构造一条全局样条曲线
    pts = ordered_points';
    curve = cscvn(pts);
    t_breaks = curve.breaks;          % 参数断点，对应原始点之间的区间
    n_pts = size(ordered_points, 1);  % 原始点数

    hold on;

    % 每一对相邻的原始点 k, k+1 对应一小段参数区间 [t_k, t_{k+1}]
    samples_per_seg = 20;  % 每段上采样点数，可调

    for k = 1:(n_pts-1)
        % 在第 k 段参数区间上采样
        t_seg = linspace(t_breaks(k), t_breaks(k+1), samples_per_seg);
        seg_points = fnval(curve, t_seg);   % 3×M

        % 决定这一段线的颜色
        c1 = ordered_colors(k);
        c2 = ordered_colors(k+1);
        if c1 == 1 && c2 == 1
            col =[1, 0.4118, 0.7059];   % 红线
        else
            col = [0, 0.7, 1];   % 蓝线
        end

        plot3(seg_points(1, :), seg_points(2, :), seg_points(3, :), ...
              'Color', col, 'LineWidth', 0.5);  % 总图线宽细一些
    end

  colors_rgb2 = zeros(n_pts,3);  % [R G B]

 for i = 1:n_pts
    % 散点颜色：1->红，0->蓝
     if ordered_colors(i, 1) > 0.5
            % 红色 - 使用更鲜艳的红色
            % colors_rgb2(i,:)=[1, 0.4118, 0.7059];
           colors_rgb2(i,:)= [1, 0, 0];   % 红色
         else
            % 蓝色 - 使用更鲜艳的蓝色
            colors_rgb2(i,:) = [0, 0.7, 1];  % 电光蓝，在黑色背景中更醒目
     end
 end


    scatter3(ordered_points(:, 1), ordered_points(:, 2), ordered_points(:, 3), ...
             point_size, colors_rgb2 ,'filled','MarkerFaceAlpha',0.5);  % 总图透明度较低
        axis tight equal;

end

%% ==================== 绘制小球的函数 ====================
function drawSpheres(points, colors_rgb, point_size)
    % 根据给定点的位置，绘制小球
    hold on;
    scale_factor = 1.5;  % 设置球的大小因子（可以调整大小）

    % 设置光照效果
    lighting phong;  % 使用 Gouraud 光照模型，适合平滑表面
    material shiny;    % 设置材质为 shiny，具有高光泽
    light('Position', [1, 0, 1], 'Style', 'infinite');  % 添加一个光源

    for i = 1:size(points, 1)
        % 生成单位球体网格数据
        [x, y, z] = sphere(50);  % 10 表示球体的分辨率
        x = x * scale_factor + points(i, 1);  % 缩放并移动球体到目标位置
        y = y * scale_factor + points(i, 2);
        z = z * scale_factor + points(i, 3)*5
            % 根据二分法获取颜色
        if colors_rgb(i, 1) > 0.5
            % 红色 - 使用更鲜艳的红色
            % face_color = [1, 0.4118, 0.7059];  % 或 [1, 0.1, 0] 获得更暖的红色
            face_color = [1, 0, 0];  % 或 [1, 0.1, 0] 获得更暖的红色
        else
            % 蓝色 - 使用更鲜艳的蓝色
            face_color = [0, 0.7, 1];  % 电光蓝，在黑色背景中更醒目
        end
        
        % 绘制小球
        surf(x, y, z, 'FaceColor', face_color, 'EdgeColor', 'none');
    end
end

%% ==================== 单帧图像保存 ====================
function saveIndividualFigure3D(fig, frame_number, model_dir)
    savefig(fig, fullfile(model_dir, sprintf('frame%d_3d.fig', frame_number)));
    saveas(fig, fullfile(model_dir, sprintf('frame%d_3d.png', frame_number)), 'png');
    fprintf('  Saved 3D individual figure: frame%d_3d.fig/png\n', frame_number);
end

%% ==================== 总图像保存 ====================
function save3DFigure(fig)
    output_dir = 'combined_analysis_3d';
    savefig(fig, fullfile(output_dir, 'chromatin_conformation_3d.fig'));
    saveas(fig, fullfile(output_dir, 'chromatin_conformation_3d.png'), 'png');
    fprintf('3D chromatin conformation figure saved to %s\n', output_dir);
end

%% ==================== 保存统计数据 ====================
function saveLinearitySummary3D_timeOrder(summary_data, output_dir)
    if isempty(summary_data)
        fprintf('No summary data to save\n');
        return;
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    save(fullfile(output_dir, 'frame_summary_timeOrder.mat'), 'summary_data');
    fprintf('Summary saved to %s\n', fullfile(output_dir, 'frame_summary_timeOrder.mat'));

    % Save as a text report
    txt_filename = fullfile(output_dir, 'frame_statistics_3d.txt');
    fid = fopen(txt_filename, 'w');
    for i = 1:size(summary_data, 1)
        fprintf(fid, 'Frame %d:\n', summary_data(i, 1));
        fprintf(fid, '   Total points: %d\n', summary_data(i, 2));
        fprintf(fid, '   Red points: %d (%.2f%%)\n', summary_data(i, 3), 100*summary_data(i, 3)/summary_data(i, 2));
        fprintf(fid, '   Blue points: %d (%.2f%%)\n', summary_data(i, 4), 100*summary_data(i, 4)/summary_data(i, 2));
        fprintf(fid, '   Compaction: %.3f\n', summary_data(i, 5));
        fprintf(fid, '--------------------------------\n');
    end
    fclose(fid);
    fprintf('3D frame statistics saved to: %s\n', txt_filename);
end

%% ==================== 图形收尾 ====================
function finalizeIndividualFigure3D(fig, frame_number, compaction, num_points)
    figure(fig);
    xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (count)');
    title(sprintf('Frame %d - Compaction: %.3f\nPoints: %d', frame_number, compaction, num_points),'Color', [0 0 0]);
     
    % set(gca, 'GridAlpha', 0.3);
    % set(gcf, 'Color', 'white');
    axis tight equal;
    daspect([1 1 1]);
    % 修改Z轴刻度标签：将所有刻度值除以5
z_ticks = get(gca, 'ZTick');  % 获取当前Z轴刻度
if ~isempty(z_ticks)
    % 计算新标签（除以5）
    z_labels = arrayfun(@(x) num2str(x/5), z_ticks, 'UniformOutput', false);
    
    % 设置新的刻度标签
    set(gca, 'ZTickLabel', z_labels);
    % 禁用自动刻度更新
    % set(gca, 'ZTickMode', 'manual', 'ZTickLabelMode', 'manual');

end
end

% %% ==================== 总图收尾 ====================
function finalize3DFigure(fig)
    figure(fig);
    xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (count)');
    % set(gca, 'GridAlpha', 0.3);
    % set(gcf, 'Color', 'white');
    axis tight equal;
    daspect([1 1 1]);
     light('Position', [1, 1, 1], 'Style', 'infinite');  % 添加一个光源
    grid on;
end  