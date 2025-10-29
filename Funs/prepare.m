function [Graph, K, F, sml_ev]=prepare(X, k, dim) %k��ʾ������dim��ʾ���
Graph = generate_graph(X, k);
num = size(Graph{1}, 1);    % 210
m=length(Graph);              % 5
% ��ͼ m = 5
% �õ�graph����С����ֵ
%
for i = 1 : m
    try
        evs(i) = eigs(Graph{i}, 1,'smallestabs');
    catch
        evss = eig(Graph{i});
        evs(i) = min(evss);
    end
end
sml_ev = min(evs);


F = cell(1, m);
L = cell(1, m);
for i=1:m
    graph = Graph{i};   %   ���ڽӾ������ʽ������graph
    L{i} = diag(sum(graph)) - graph; %�õ�������˹����
    tmp = eig1(L{i}, dim+1, 0, 1);
    tmp(:, 1) = []; %��һ��û����
    tmp = tmp./repmat(sqrt(sum(tmp.^2,2)),1,dim);
    tmp(isnan(tmp))=0;
    F{i} = tmp;  %�õ���Ƕ��
%     F{i} = eig1(L{i}, dim, 0, 1);
end
K = cell(1, m);
for i = 1 :m
    K{i} = eye(num) - L{i} / 2;
end
end