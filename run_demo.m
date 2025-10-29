clc
clear 
close all
warning off
addpath(genpath('./datasets'))
addpath('./funs');
addpath('./finchpp');
addpath('./datasets');
data = {'HW.mat'}; 

k=5;
lambda = [1];

for idx = 1: length(data) 
    fprintf('-------%s\n',data{idx})
    load(data{idx})


    c = length(unique(Y));   


    for v = 1:length(lambda)
    tic;
    [y,obj,U,S0,S0_initial] = main_max2(X,c,lambda(v),k);
    total_time = toc;
    result_MEA_PKN(v,:) = [ClusteringMeasure(Y,y) total_time];
    end
    

% 绘制图形（加粗线条）
plot(obj, 'LineWidth', 2.5);
% 给图形添加标题和坐标轴标签
title('收敛函数');
xlabel('x 轴');
ylabel('y 轴');

% 显示网格线
grid on;
% 加粗坐标轴
set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold');
end
