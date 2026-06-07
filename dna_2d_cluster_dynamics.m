% clear all 
% close all

%% 原始数据mat提取有效数据
matFiles = dir('*.mat');
% 找到名字中包含 'area' 的文件，并提取数字
filePattern = 'area(\d+)_loc_time.mat';
for i = 1:length(matFiles)
    tokens = regexp(matFiles(i).name, filePattern, 'tokens');
    if ~isempty(tokens)
        areaNumber = tokens{1}{1};  % 提取到的数字
     
        break;
    end
end
% 加载对应的 .mat 文件
fileName = ['area' areaNumber '_loc_time.mat'];
load(fileName);
eval(['loc = area' areaNumber '_loc;']);
eval(['time = tim;']);
% eval(['tid = area' areaNumber '_tid;']);

%% 计算时间差 可选 前期已经计算完成
% result_table = findPointTimeDiffs(loc, tim)
% 
%      all_time_diffs = [];
%     for i = 1:height(result_table )
%         if ~isempty(result_table .TimeDiffs{i})
%             all_time_diffs = [all_time_diffs; result_table .TimeDiffs{i}(:)];
%         end
%     end
%  save(['area', num2str(areaNumber), 'all_time_diffs.mat'], 'all_time_diffs');

%% 根据时间见间隔 在时间上剥离 内含降噪代码 在计算时间间隔的时候降噪了 后面看结构保留 详见assignPointIndexSimple 头几行
time_threshold=0.05 % 50ms 10 35 40 45 55 60 65 100
[point_idx, loc] = assignPointIndexSimple(loc, tim, time_threshold);
 
%% 分frame 空间上剥离
plotLabelsSeparately_all(loc, point_idx, tim, areaNumber)%会储两个 一个是everyframe 一个是十个点以上的

 
 