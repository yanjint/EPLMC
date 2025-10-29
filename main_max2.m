function [result_our,obj,U,S0,S0_initial] = main_max2( X,  c, lambda,k)


for j=1:length(X)
           colnum=size(X{j},2);           
           mole = repmat(std(X{j},0,2),1,colnum);
           mole(mole==0) = 1; 
           X{j}=(X{j}-repmat(mean(X{j},2),1,colnum))./mole;
end

 

for v = 1 :length(X)
    XX{v} = X{v}';
end
options.NeighborMode = 'KNN';
options.k = k; 
options.WeightMode = 'HeatKernel';
for v = 1:length(X)
    fea_v = XX{1,v};
    U(v) = {constructW(fea_v',options)}; 
end

S0_initial = U;


m = length(X);
belta=1;

n = size(U{1},1);

%% init W
W = orth(rand(n, c));

%% init La
La = cell(1,m);
for i=1:m
    AA = (U{i}' + U{i}) / 2;
    DD = diag(sum(AA));
    La{i} = DD - AA;    
end

qw = zeros(1, m);

% U_rep = zeros(n);
for v = 1:m
   num = size(U{v},1);
   S10 = (U{v}+U{v}')/2;
   D10 = diag(sum(S10));
   L0 = D10 - S10;
   [F0, ~, ~] = eig1(full(L0), c+1, 0);
   F{v} = F0(:,1:c);
   U_sub{v} = F{v} * F{v}';

end

G = zeros(num,num);
for v = 1:m   
    G = G + U_sub{v};
end



niter = 20;

for it=1:niter    
    fprintf('****这是第%d轮****\n',it)
     %
    for v = 1:m   
        qw(v) = 1/(2 * max(norm(W*W' - U{v}, 'fro'), 1e-8));
    end
    
    %更新U
    Z = W*W'; 

idxx = cell(1,m);
ed = cell(1,m);
for v = 1:m
    ed{v} = L2_distance_1(X{v}', X{v}');
    [~, idxx{v}] = sort(ed{v}, 2); % sort each row
end
    
    
    for v = 1:m
        U{v} = zeros(num);
        for i = 1:num
            id = idxx{v}(i,2:k+2);
            di1 = ed{v}(i, id);  
            di2 = 2*qw(v)*Z(i, id);
            ad = -(di1-di2)/(2*(qw(v)+belta));
          U{v}(i,id) = EProjSimplex_new(ad);
          end
    end
    
    La = cell(1,m);
    for i=1:m
        AA = (U{i}' + U{i}) / 2;
        DD = diag(sum(AA));
        La{i} = DD - AA;    
    end
      
    %更新W
    A0 = zeros(n, n);
    for v = 1 : m
        A0 = A0 + qw(v)*U{v};
    end
    A_temp = 2*A0 + lambda * G;
    A_temp = (A_temp +A_temp')/2;
    [W, ~, ~] = eig1(A_temp, n, 1);
    W = W(:,2:c+1);


   
        
    res = 0;
    for v = 1 : m
             res = res + norm(W*W' - U{v}, 'fro')^2+2*trace(X{v}'*La{v}*X{v})+ belta *norm(U{v},'fro')^2;
    end
    obj(it) = res - lambda*trace(W*W'*G);
%     if it > 2 && abs(obj(it) - obj(it-1)) / abs(obj(it-1)) < 10^(-6)
%         break
%     end

end
W = NormalizeFea(W,0);
obj = obj;
result_our = kmeans(W, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
U = W*W';
S0 = U;
end