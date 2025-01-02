clc
clear 
close all
warning off
addpath(genpath('./datasets'))
addpath('./funs');
addpath('./finchpp');
addpath('./datasets');
data = {'MSRC.mat'}; 


% lambda = [0.1 0.3 0.5 0.7 0.9 1 1.1 1.3 1.5 1.7 1.9 2];
lambda = [1];
for idx = 1: length(data) 
    fprintf('-------%s\n',data{idx})
    load(data{idx})

    c = length(unique(Y));   

%     tic;
    for v = 1:length(lambda)
    tic;
    [y,obj,U,S0,S0_initial] = main_max(X,c,lambda(v));
    total_time = toc;
    result_MEA_PKN(v,:) = [ClusteringMeasure(Y,y) total_time];
    end

end
