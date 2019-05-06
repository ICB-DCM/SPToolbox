function U = createRandomRotationMatrix(dim)
while true
    M = randn(dim,dim);
    [Q,R] = qr(M);
    D = diag(R);
    D = diag(D)./abs(D);
    U = Q * D;
    if det(U) > 0
        break;
    end
end
end