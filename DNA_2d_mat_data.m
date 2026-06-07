
clear all
close all

% 现在，loaded_data 中包含了.mat文件中的数据
% 获取当前文件夹下的所有.mat文件
matFiles = dir('*.mat');

% 确保文件夹中只有一个.mat文件
if length(matFiles) == 1
% 加载.mat文件
matFileName = matFiles(1).name;
load(matFileName);
tim=tim';
% 提取.mat文件名中的area+数字的组合
expr = 'area(\d+)';  % 使用正则表达式提取area+数字的组合
match = regexp(matFileName, expr, 'tokens', 'once');

if ~isempty(match)
areaNumber = str2double(match{1});

% 构建对应的变量名
locVarName = ['area', num2str(areaNumber), '_loc'];
ratioVarName = ['area', num2str(areaNumber), '_ratio'];
timeVarName = ['tim'];
tidVarName = ['area', num2str(areaNumber), '_tid'];;
%% 2d 5 3d 10
x=itr.gvx(:,4);
y=itr.gvy(:,4);
z=itr.dmz(:,4);
% 将向量中的每个元素除以3e-9并向下取整，得到归后的值
classified_x = floor(x / 3e-9);
classified_y = floor(y / 3e-9);
classified_z = floor(z / 3e-9);
classified_x = classified_x-min(classified_x);
classified_y = classified_y-min(classified_y);
classified_y=max(classified_y)-classified_y;
classified_z = classified_z-min(classified_z);

% 获取数据
eval([locVarName, ' = [classified_x, classified_y, classified_z];']);
eval([ratioVarName, ' = itr.dcr(:, 5);']);
eval([timeVarName, ' = tim;']);
eval([tidVarName, ' = tid;']);

% 绘制 3D 散点图
locData= [classified_x, classified_y, classified_z];
timData=tim;
figure;
scatter(locData(:, 1)*3, locData(:, 2)*3, 15, timData, 'filled');
axis equal;
set(gcf,'color','black'); %窗口背景白
colordef black; %2D/
% 添加标签
xlabel('X(nm)');
ylabel('Y(nm)');
zlabel('Z(nm)');
title(['Scatter Plot for area', num2str(areaNumber), ' with Time Color']);
savefig(['area', num2str(areaNumber), '_tim.fig'] )

% figure;
% scatter(locData(:, 1)*3, locData(:, 2)*3, 15, tid, 'filled','MarkerFaceAlpha', 0.5);
% % 设置当前颜色映射
% colormap 'jet'
% % 添加标签
% xlabel('X(nm)');
% ylabel('Y(nm)');
% axis equal;
% set(gcf,'color','black'); %窗口背景白
% colordef black; %2D/
% title(['Scatter Plot for area', num2str(areaNumber), ' TID(Tracer ID)']);
% savefig(['area', num2str(areaNumber), '_TID.fig'] )
% 
% 
% axis equal;
% xlabel('X(nm)');
% ylabel('Y(nm)');
% title(['Scatter Plot for area', num2str(areaNumber), ' TID (Tracer ID)']);

new_tid = cumsum([1, diff(tid) ~= 0])';
figure;
scatter(locData(:, 1)*3, locData(:, 2)*3, 15, new_tid, 'filled','MarkerFaceAlpha', 0.5);
% 使用离散颜色映射
num_tids = length(unique(new_tid));  % 你的TID数量
% 生成80个明显不同的颜色
rng(42);  % 设置随机种子保证可重复
custom_cmap = rand(num_tids, 3);  % 随机RGB颜色
colormap(custom_cmap);
colorbar;
xlabel('X(nm)');
ylabel('Y(nm)');
axis equal;
set(gcf,'color','black'); %窗口背景白
colordef black; %2D/
title(['Scatter Plot for area', num2str(areaNumber), ' TID(Tracer ID)']);
savefig(['area', num2str(areaNumber), '_TID_2.fig'] )
% 保存 loc, ratio 和 time 数据到一个.mat文件
save(['area', num2str(areaNumber), '_loc_time.mat'], locVarName, ratioVarName, timeVarName,tidVarName);
else
disp('File name does not match the expected pattern.');
end
else
disp('There should be only one .mat file in the folder.');
end
