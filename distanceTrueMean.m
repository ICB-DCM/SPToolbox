function [Dt,Gt] = distanceTrueMean(x,N,L)
x = reshape(x,[N,L-1]);
X = [x,-sum(x,2)];
X = reshape(X,[1,N*L]);
[D,G] = distanceDiracGaussian(X,N,L);
Dt = D;
for i=1:L-1
    dXdx(:,:,i) = [zeros(N,i-1),ones(N,1),zeros(N,L-1-i),-ones(N,1)];
end
Gt = permute(sum(bsxfun(@times, reshape(G,[N,L]),dXdx),2),[1,3,2]);
Gt = reshape(Gt,[1,N*(L-1)]);