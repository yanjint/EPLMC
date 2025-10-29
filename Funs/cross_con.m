function KF = cross_con(K, F, order_K)
m = length(F);
for i = 1: m
    for j = 1: m
        x = F{j};
        for k = 1: order_K
            x = K{i} * x;                %�൱�������е�U
        end
        KF{i, j} = x * diag(diag(x'*x).^-.5); %���б�ɵ�λ��
    end
end
end