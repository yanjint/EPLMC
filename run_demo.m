clc
clear 
close all
warning off
addpath(genpath('./datasets'))
addpath('./funs');
addpath('./finchpp');
addpath('./datasets');
% 'ORL','MSRC','COIL20','Cite','Caltech101-7','Caltech101-20','bbcsport','3sources','MNIST.mat'
data = {'ORL.mat'}; 

% load('ORL.mat')
lambda = [0.1 0.3 0.5 0.7 0.9 1 1.1 1.3 1.5 1.7 1.9 2];
% lambda = [1];
for idx = 1: length(data) 
    fprintf('-------%s\n',data{idx})
    load(data{idx})
%     X = data;
%     Y = truelabel{1};
%     Y = truth;
%     for i=1:length(X)
%         X{i}=X{i}(1:1280,:);
%     end
%     Y = Y(1:1280,:);
%     A = generate_graph(X, k);
    c = length(unique(Y));   % 类别 classes = 7

%     tic;
    for v = 1:length(lambda)
    tic;
    [y,obj,U,S0,S0_initial] = main_max(X,c,lambda(v));
    total_time = toc;
    result_MEA_PKN(v,:) = [ClusteringMeasure(Y,y) total_time];
    end

end
