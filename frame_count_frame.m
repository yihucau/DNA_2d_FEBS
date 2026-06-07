% 自动查找 area*_frames_filtered.mat 文件
files = dir('area*_frames_filtered.mat');

if isempty(files)
    error('当前目录下未找到 area*_frames_filtered.mat 文件。');
elseif length(files) > 1
    fprintf('检测到多个匹配文件：\n');
    for k = 1:length(files)
        fprintf('%d: %s', k, files(k).name);
    end
    error('请确保当前目录下只有一个 area*_frames_filtered.mat 文件，或先手动选择。');
end

input_file = files(1).name;

% 自动生成输出文件名：area*_frames_filtered.mat -> area*_frames_updated.mat
output_file = strrep(input_file, '_frames_filtered.mat', '_frames_updated.mat');

% 加载数据文件到结构体，避免 eval
data = load(input_file);

% 获取所有变量名
variable_names = fieldnames(data);

% 筛选出 frame 开头的变量
frame_vars = variable_names(startsWith(variable_names, 'frame'));

if isempty(frame_vars)
    error('文件 %s 中未找到以 frame 开头的变量。', input_file);
end

% 按数字顺序排序 frame 变量
frame_nums = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), frame_vars);
[~, idx] = sort(frame_nums);
frame_vars = frame_vars(idx);

% 遍历每个 frame
for i = 1:length(frame_vars)
    frame_name = frame_vars{i};
    frame_data = data.(frame_name);
    
    % 每次处理新一帧数据时，重置计数器
    point_count_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    
    % 若原始数据不足3列，先扩展到第3列
    if size(frame_data, 2) < 3
        frame_data(:, 3) = 0;
    end
    
    % 遍历当前 frame 中的每个点
    for j = 1:size(frame_data, 1)
        x = frame_data(j, 1);
        y = frame_data(j, 2);
        
        % 创建点的唯一标识符
        point_key = sprintf('%.6f,%.6f', x, y);
        
        % 统计当前 frame 内重复出现次数
        if isKey(point_count_map, point_key)
            current_count = point_count_map(point_key) + 1;
            point_count_map(point_key) = current_count;
        else
            current_count = 1;
            point_count_map(point_key) = current_count;
        end
        
        % 填入第3列
        frame_data(j, 3) = current_count;
    end
    
    % 存回结构体
    data.(frame_name) = frame_data;
    
    fprintf('已处理 %s，共 %d 个点', frame_name, size(frame_data, 1));
end

% 保存处理后的变量
save(output_file, '-struct', 'data');

fprintf('处理完成！所有点的出现次数已填入第3列。\n');
fprintf('输入文件: %s', input_file);
fprintf('输出文件: %s', output_file);
fprintf('总共处理了 %d 个 frame', length(frame_vars));

