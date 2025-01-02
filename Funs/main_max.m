function [S,  obj] = main_max(A, X,k)
%��ȡ��ͼ��Ŀ
m = length(X);

%��ȡ������
n = size(A{1},1);

%% init
a = ones(1,m) / m;

La = cell(m,1);
for i=1:m
    AA = (A{i}' + A{i}) / 2;
    DD = diag(sum(AA));
    La{i} = DD - AA;    
end
niter = 20;
lambda = 1;

for it=1:niter    
    fprintf('****���ǵ�%d��****\n',it)
    
    %����U
    U = cell(m, 1);
    for i=1:m
        U{i} = (eye(n) + 2*La{i})^(-1)*X{i};    
    end
    
    %����S
    A_hat = zeros(n,n);
    for i =1:m
       A_hat = A_hat + (a(i)*A{i})'*(a(i)*A{i});
    end
    S = eig1(A_hat,n,1);
    
    %����a
    for i = 1 : m
        beta(i) = -trace(S'*A{i}'*A{i}*S);
    end
    a = (beta) .^-1 / sum((beta) .^-1);
    
    %����A
    for i=1:m
        [A{i},belta(i)] = update_A(U{i}', a(i), lambda, S, k);
    end
    
    La = cell(m,1);
    for i=1:m
        AA = (A{i}' + A{i}) / 2;
        DD = diag(sum(AA));
        La{i} = DD - AA;    
    end
   
    obj(it) = get_obj(X, U, S, a, A, belta);
    lambda = obj(it);
%     if it>2 && (obj(it) -obj(it-1)) / obj(it-1) < 1e-4
%           break
%     end
end
S = S;
obj = obj;
end




function obj = get_obj(X, U, S, a, A, belta)
m = length(X);

for i=1:m
    AA = (A{i}' + A{i}) / 2;
    DD = diag(sum(AA));
    La{i} = DD - AA;    
end

numerator = 0;         %����
denominator = 0;     %��ĸ
for i=1:m
%     numerator = numerator + norm(S*(a(i)*A{i}),'fro')^(2);
    numerator = numerator + a(i)*A{i};
    
    denominator = denominator + norm(X{i}-U{i},'fro')^(2) + 2*trace(U{i}'*La{i}*U{i}) + belta(i)*norm(A{i})^(2); 
    
end
numerator = norm(S*numerator,'fro')^(2);
obj = numerator / (denominator +eps);

end