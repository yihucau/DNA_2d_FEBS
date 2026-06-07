function plotLabelsSeparately_all(loc, point_idx, tim, areaNumber)
% 计算并保存每个label的frame数据，基于前后点距离过滤离群点，并增加首尾距离过滤，不进行绘图
% 输入参数：
%   loc: 坐标数据 [x, y]
%   point_idx: 标签索引
%   tim: 时间数据
%   areaNumber: 区域编号，用于保存文件名

    % -------------------------
    % 可调参数
    % -------------------------
    threshold = 10;      % 单点前后距离阈值（离群点）
    end_to_end_threshold = 100;  % 首尾点距离阈值，超过则整条轨迹过滤

    % 找到所有点
    in_region = true(size(loc, 1), 1);
    region_points    = loc(in_region, :);
    region_point_idx = point_idx(in_region);
    region_tim       = tim(in_region);
    unique_point_idx = unique(region_point_idx);
    num_point_idx    = length(unique_point_idx);

    fprintf('总标签数量: %d\n', num_point_idx);

    % 初始化存储所有frame数据的结构体
    all_frames          = struct();
    all_frames_filtered = struct();
    filtered_count      = 0;

    for label_idx = 1:num_point_idx
        current_label = unique_point_idx(label_idx);
        label_mask    = (region_point_idx == current_label);
        points        = region_points(label_mask, :);
        times         = region_tim(label_mask);

        % 保存原始frame数据： [x, y, t, label]
        frame_data = [points, times, repmat(current_label, size(points, 1), 1)];
        all_frames.(sprintf('frame%d', current_label)) = frame_data;

        % -------------------------
        % 基于距离的离群点过滤
        % -------------------------
        if size(points, 1) >= 5
            inlier_mask = true(size(points, 1), 1);

            % 按时间排序
            [times_sorted, sort_idx] = sort(times);
            points_sorted            = points(sort_idx, :);
            N = size(points_sorted, 1);

            % ---------- 第一个点 ----------
            if N >= 2
                dist_first = norm(points_sorted(2,:) - points_sorted(1,:));
                fprintf('Label %d: k=1 dist_first=%.2f\n', current_label, dist_first);
                if dist_first > threshold
                    inlier_mask(sort_idx(1)) = false;
                end
            end

            % ---------- 中间点 ----------
            for k = 2:N-1
                d_prev = norm(points_sorted(k,:)   - points_sorted(k-1,:));
                d_next = norm(points_sorted(k,:)   - points_sorted(k+1,:));
                fprintf('Label %d: k=%d d_prev=%.2f d_next=%.2f\n', ...
                        current_label, k, d_prev, d_next);

                if d_prev > threshold && d_next > threshold
                    inlier_mask(sort_idx(k)) = false;  % 删除当前点
                end
            end

            % ---------- 最后一个点 ----------
            if N >= 2
                dist_last = norm(points_sorted(N,:) - points_sorted(N-1,:));
                fprintf('Label %d: k=%d dist_last=%.2f\n', current_label, N, dist_last);
                if dist_last > threshold
                    inlier_mask(sort_idx(N)) = false;
                end
            end

            % 保留非离群点
            points_filtered = points(inlier_mask, :);
            times_filtered  = times(inlier_mask);

            % ---------- 新增：首尾点距离过滤整个 frame ----------
            if ~isempty(points_filtered)
                % 按时间排序后，计算首尾点距离
                [times_f_sorted, sort_f_idx] = sort(times_filtered);
                points_f_sorted = points_filtered(sort_f_idx, :);

                if size(points_f_sorted, 1) >= 2
                    dist_end_to_end = norm(points_f_sorted(end,:) - points_f_sorted(1,:));
                    fprintf('Label %d: end-to-end distance = %.2f\n', ...
                            current_label, dist_end_to_end);

                    if dist_end_to_end > end_to_end_threshold
                        % 首尾距离太大，整条轨迹过滤
                        frame_data_filtered = [];
                        filter_status = sprintf('过滤（首尾距离 %.2f > %.2f）', ...
                                                dist_end_to_end, end_to_end_threshold);
                    else
                        % 正常保留
                        frame_data_filtered = [points_filtered, ...
                                               times_filtered, ...
                                               repmat(current_label, size(points_filtered,1),1)];
                        filter_status = '保留（距离过滤 + 首尾检查通过）';
                    end
                else
                    % 过滤后剩下不到两个点，意义不大，直接过滤
                    frame_data_filtered = [];
                    filter_status = '过滤（过滤后点数 <2）';
                end
            else
                frame_data_filtered = [];
                filter_status = '过滤（全部离群）';
            end

            all_frames_filtered.(sprintf('frame%d_filtered', current_label)) = frame_data_filtered;

            if ~isempty(frame_data_filtered)
                filtered_count = filtered_count + 1;
            end

        else
            all_frames_filtered.(sprintf('frame%d_filtered', current_label)) = [];
            filter_status = '过滤（点数 <5）';
        end
        % -------------------------

        % 显示进度信息
        fprintf('计算 Label %d: %d 个数据点 [%s]\n', ...
                current_label, size(points, 1), filter_status);

        if ~isempty(points)
            time_range = [min(times), max(times)];
            time_span  = time_range(2) - time_range(1);
            fprintf('  Label %d 统计: 质心(%.1f, %.1f), 时间 %.2f-%.2fs, 时长 %.2fs\n', ...
                    current_label, mean(points(:,1)), mean(points(:,2)), ...
                    time_range(1), time_range(2), time_span);
        end
    end

    % -------------------------
    % 保存原始frame数据
    % -------------------------
    save_filename = ['area', num2str(areaNumber), '_frames.mat'];
    save(save_filename, '-struct', 'all_frames');
    fprintf('\n已保存 %d 个原始frame数据到文件: %s\n', num_point_idx, save_filename);

    % -------------------------
    % 保存过滤后的frame数据
    % -------------------------
    save_filtered_filename = ['area', num2str(areaNumber), '_frames_filtered.mat'];
    save(save_filtered_filename, '-struct', 'all_frames_filtered');
    fprintf('已保存 %d 个过滤后的frame数据到文件: %s\n', filtered_count, save_filtered_filename);

end
