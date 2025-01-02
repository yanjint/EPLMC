function [F, obj] = solve_F(K, F)
    fKf = sum(F .* (K * F))';
    ff = sum(F)';

    m_all = vec2ind(F')';

    obj(1) = sum(fKf ./ ff);

    for iter = 2:50
        for i = 1:size(F, 1)
            m = m_all(i);
            if ff(m) == 1
                % avoid generating empty cluster
                continue;
            end

            Y_A = F' * K(:, i);

            fKf_s = fKf + 2 * Y_A + K(i, i);
            fKf_s(m) = fKf(m);
            ff_k = ff + 1;
            ff_k(m) = ff(m);

            fKf_0 = fKf;
            fKf_0(m) = fKf(m) - 2 * Y_A(m) + K(i, i);
            ff_0 = ff;
            ff_0(m) = ff(m) - 1;

            delta = fKf_s ./ ff_k - fKf_0 ./ ff_0;

            [~, p] = max(delta);
            if p ~= m
                fKf([m, p]) = [fKf_0(m), fKf_s(p)];
                ff([m, p]) = [ff_0(m), ff_k(p)];

                F(i, [p, m]) = [1, 0];
                m_all(i) = p;
            end
        end
        obj(iter) = sum(fKf ./ ff);

        if iter > 2 && (obj(iter) - obj(iter - 1)) / obj(iter - 1) < 1e-9
            break;
        end
    end
end