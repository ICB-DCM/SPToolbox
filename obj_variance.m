function [obj_variance,grad_variance] = obj_variance(w_c,x,N)
pvpw =  bsxfun(@times,permute(x,[1,3,2]),permute(x,[3,1,2]));
variance = sum(bsxfun(@times,permute(w_c,[3,1,2]),...
            pvpw),3);
obj_variance = sum(sum((variance - eye(N)).^2,1),2);
pVpv = 2*variance-2*eye(N);
grad_variance = permute(sum(sum(bsxfun(@times,pVpv,pvpw),1),2),[1,3,2]);
end