export longhex_2d
## constructing longHex2d
function longhex_2d(; K0 = 0, K1 = 3)
    a1 = [1; 0]
    Q6 = [cos(π/3) -sin(π/3); sin(π/3) cos(π/3)]
    a2 = Q6 * a1
    a3 = Q6 * a2
    
    nX = 2K0 + 1
    X = a1 .* [reshape(collect(-K0:K0), 1, nX); reshape(collect(-K0:K0), 1, nX)]
    c = [nX nX 1 1 1 nX]
    
    L1 = 2K0
    L2 = 0
    
    for j = 1:K1
        L1 += 1
        L2 += 1

        L1vec = [reshape(collect(0:L1-1), 1, L1); reshape(collect(0:L1-1), 1, L1)]
        L2vec = [reshape(collect(0:L2-1), 1, L2); reshape(collect(0:L2-1), 1, L2)]

        X = hcat(X, repeat(X[:, c[1]]+a1, outer=(1,L2)) + repeat(a3, outer=(1,L2)) .* L2vec,
                    repeat(X[:, c[2]]+a2, outer=(1,L1)) - repeat(a1, outer=(1,L1)) .* L1vec,
                    repeat(X[:, c[3]]+a3, outer=(1,L2)) - repeat(a2, outer=(1,L2)) .* L2vec,
                    repeat(X[:, c[4]]-a1, outer=(1,L2)) - repeat(a3, outer=(1,L2)) .* L2vec,
                    repeat(X[:, c[5]]-a2, outer=(1,L1)) + repeat(a1, outer=(1,L1)) .* L1vec,
                    repeat(X[:, c[6]]-a3, outer=(1,L2)) + repeat(a2, outer=(1,L2)) .* L2vec)
        c[1] = c[6] + L2 - 1
        if j==1
            c[1] = c[1] + 1
        end
        c[2] = c[1] + L2
        c[3] = c[2] + L1
        c[4] = c[3] + L2
        c[5] = c[4] + L2
        c[6] = c[5] + L1
    end
    
    return X
end