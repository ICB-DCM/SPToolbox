function [obj_variance,grad_variance] = obj_variance(w_c,x,N,L)
w_c = w_c';
w_c = [w_c,1-sum(w_c,2)];
pvpw =  bsxfun(@times,permute(x,[1,3,2]),permute(x,[3,1,2]));
variance = sum(bsxfun(@times,permute(w_c,[3,1,2]),...
            pvpw),3);
obj_variance = sum(sum((variance - eye(N)).^2,1),2);
pVpv = 2*variance-2*eye(N);
grad_variance = permute(sum(sum(bsxfun(@times,pVpv,pvpw),1),2),[1,3,2]);
pWpw = [eye(L-1),-ones(L-1,1)];
grad_variance = permute(sum(bsxfun(@times,grad_variance,pWpw),2),[2,1]);
%
%gamma = 1e-3;
%obj_variance = obj_variance + gamma*sum(w.^2);
end