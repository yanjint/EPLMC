function Graph=generate_graph(X,k)
Graph = pcan(X,k);
end

function W_cell = pcan(X_cell,k)
n = size(X_cell{1},1); %210
z=length(X_cell); %5
W_cell=cell(1,z);
for j=1:z
    X=X_cell{j}';
    D = L2_distance_1(X, X); %210*210
    [dumb, idx] = sort(D, 2); %对行进行处理
    rr = zeros(n,1);
    W = zeros(n,n);
    for i = 1:n
        di=dumb(i,2:k+2); %表示从节点的1阶开始取数(1*6维)
        rr(i) = 0.5*(k*di(k+1)-sum(di(1:k))); %
        id = idx(i,2:k+2);
        W(i,id) = (di(k+1)-di) / (k*di(k+1)-sum(di(1:k))+eps);
    end
%     W_cell{j}=(W+W')/2+eye(n);
    W = (W+W')/2; %为了实现对称，
    W_cell{j}=(W+W')/2;
end
end

function d = L2_distance_1(a,b)
if (size(a,1) == 1)
    a = [a; zeros(1,size(a,2))];
    b = [b; zeros(1,size(b,2))];
end
aa=sum(a.*a);
bb=sum(b.*b); 
ab=a'*b;
d = repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab;
d = real(d);
d = max(d,0);
end