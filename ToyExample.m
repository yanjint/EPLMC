function [result_our,obj,U,S0,S0_initial] = ToyExample( X,  c, lambda)


%获取视图数目
m = length(X);

%获取样本数
n = size(X{1},1);
for i = 1:m
    X{i} = X{i}';
end
% for j=1:length(X)
%            colnum=size(X{j},2);           % colnum = 24
%            mole = repmat(std(X{j},0,2),1,colnum);
%            mole(mole==0) = 1; 
%            X{j}=(X{j}-repmat(mean(X{j},2),1,colnum))./mole;
%  end

    for i = 1:m
        for  j = 1:n
            normItem = std(X{i}(:,j));
            if (0 == normItem)
                normItem = eps;
            end;
            X{i}(:,j) = (X{i}(:,j)-mean(X{i}(:,j)))/(normItem);
        end;
    end;



 k = 7;

options.NeighborMode = 'KNN';
options.k = k; 
options.WeightMode = 'HeatKernel';
for v = 1:length(X)
    fea_v = X{1,v};
    A(v) = {constructW(fea_v',options)}; 
end

S0_initial = A;



%% init W
W = orth(rand(n, c));

%% init La
La = cell(1,m);
for i=1:m
    AA = (A{i}' + A{i}) / 2;
    DD = diag(sum(AA));
    La{i} = DD - AA;    
end

dw = zeros(1, m);

A_rep = zeros(n);
for v = 1:m
   num = size(A{v},1);
   S10 = (A{v}+A{v}')/2;
   D10 = diag(sum(S10));
   L0 = D10 - S10;
   [F0, ~, ~] = eig1(full(L0), c+1, 0);
   F{v} = F0(:,1:c);
   A_sub{v} = F{v} * F{v}';
   A_rep = A_rep + A{v};
end

K0 = zeros(num,num);  
for v = 1:m   
    K0 = K0 + A_sub{v};
end


niter = 10000;

for it=1:niter    
    fprintf('****这是第%d轮****\n',it)
     
    for v = 1:m   
        dw(v) = 1/(2 * norm(W*W' - A{v}, 'fro'));
    end
    
    %更新A
    WW = W*W';
    idxx = cell(1,m);
    ed = cell(1,m);
    for v = 1:m
        ed{v} = L2_distance_1(X{v}, X{v});
        [~, idxx{v}] = sort(ed{v}, 2); % sort each row
    end
    
    for v = 1:m
        A{v} = zeros(num);
        for i = 1:num
            id = idxx{v}(i,2:k+2);
            di = ed{v}(i, id);
            
            numerator = di(k+1)-di+2*dw(v)*WW(i,id(:))-2*dw(v)*WW(i,id(k+1));
            bta(i) = (k*di(k+1)-sum(di(1:k))-2*k*dw(v)*WW(i,id(k+1))-2*dw(v))/2;
            denominator1 = k*di(k+1)-sum(di(1:k));
            denominator2 = 2*dw(v)*sum(WW(i,id(1:k)))-2*k*dw(v)*WW(i,id(k+1));
            A{v}(i,id) = max(numerator/(denominator1+denominator2+eps),0);
        end
         belta{v} = mean(bta);
    end
    
    La = cell(1,m);
    for i=1:m
        AA = (A{i}' + A{i}) / 2;
        DD = diag(sum(AA));
        La{i} = DD - AA;    
    end
   
    for v = 1:m   
        dw(v) = 1/(2 * norm(W*W' - A{v}, 'fro'));
    end
   
    %更新W
    A0 = zeros(n, n);
    for v = 1 : m
        A0 = A0 + dw(v)*A{v};
    end
    A_temp = 2*A0 + lambda * K0;
    [W, ~, ~] = eig1(A_temp, n, 1);
    W = W(:,2:c+1);
   
        % Calculate Obj
    res = 0;
    for v = 1 : m
        res = res + norm(W*W' - A{v}, 'fro')^2+2*trace(X{v}*La{v}*X{v}')+ belta{v} *norm(A{v},'fro')^2;
    end
    obj(it) = res - lambda*trace(W'*K0*W);
    if it > 2 && abs(obj(it) - obj(it-1)) / abs(obj(it-1)) < 10^(-4)
        break
    end

end
W = NormalizeFea(W,0);
obj = obj;
result_our = kmeans(W, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
U = W*W';
S0 = A;
end