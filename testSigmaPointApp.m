% USAGE:
% ======
% [SP] = testSigmaPointApp(nonfun,beta,D,estruct,op_SP)
% [SP] = testSigmaPointApp(nonfun,beta,D,estruct,op_SP)
%
% INPUTS:
% =======
% nonfun ... nonlinear function with argument phi to which the sigma point
%            approximation is applied. 
%            the function handle should return the following objects
%            depending on the number of requested nargout
%            [y,sy] = nonfun(phi)
%               y ... (for nargout >= 1) 
%                   output vector with rows corresponding to time-points
%                   and columns corresponding to observables
%               sy ... (for nargout >= 2)
%                   sensitivity of output vector of dimension  
%                   [time x observables x parameters]
% beta   ... vector of common effects
% D      ... covariance matrix of the random effect b
% estruct ... struct containing all information about the employed effect
%            composition
%    .beta ... @beta(xi) defines the how the common effects are constructed
%            from hyperparameter xi
%    .delta ... @delta(xi) defines the how the random effect
%            hyperparameters are constructed from hyperparameter xi
%    .phi ... @phi(beta,b) defines the combination of random effects and
%            common effects
%    .dphidbeta ... @dphidbeta(beta,b) defines the derivative of phi w.r.t.
%            beta
%    .dphidb ... @dphidb(beta,b) defines the derivative of phi w.r.t. b
%    .sigma ... @sigma(phi) defines the parametrisation of the technical
%            noise
%    .dsigmadphi ... @dsigmadphi(phi) defines the derivative of sigma with
%            respect to phi
% op_SP  ... struct containing options on what the output struct should
%            contain
%    .nderiv ... flag indicating whether derivatives should be computed
%    .req ... vector of flags indicating whether a certain output is
%            requested
%            (1) my
%            (2) Cy
%            (3) Cxy
%            (4) mz
%            (5) Cz
%    .nsamples ... number of points to sample from for testing
%
% OUTPUTS:
% =======
% whether a certain output is requested must be specified in op_SP.req
% SP.my     ... mean of the sigma-point approximation
%            [time x observables]
% SP.Cy     ... covariance of the sigma-point approximation
%            [time x observables x observables]
% SP.Cxy    ... temporal cross-correlation of the sigma-point approximation
%            [time x time x observables x observables]
% SP.mz     ... full state vector (all observables at all times)
%            [time*observables]
% SP.Cz     ... full state covariance, this contains both covariances across
%            time and states
%            [time*observables x time*observables]
% SP.B_SP   ... locations of employed sigma points
%            [parameters x 2*parameters+1]
% SP.Y      ... evaluation of nonfun at the employed sigma points
%            [time x observables x 2*parameters+1]
% SP.my_true     ... mean of the true distribustion
%            [time x observables]
% SP.Cy_true     ... covariance of the true distribustion
%            [time x observables x observables]
% SP.Cxy_true    ... temporal cross-correlation of the true distribustion
%            [time x time x observables x observables]
% SP.mz_true     ... full state vector (all observables at all times) of the true distribustion
%            [time*observables]
% SP.Cz_true     ... full state covariance, this contains both covariances across
%            time and states of the true distribustion
%            [time*observables x time*observables]
% SP.Y_true      ... evaluation of nonfun at the employed sample points
%            [time x observables x nsamples]

function SP = testSigmaPointApp(varargin)

nonfun = varargin{1};
xi  = varargin{2};
estruct = varargin{3};
op_SP = varargin{4};
op_SP.nderiv = 0;
if(~isfield(op_SP,'nsamples'))
    op_SP.nsamples = 50000;
end

SP = getSigmaPointApp(nonfun,xi,estruct,op_SP);

% Dimensions
n_t = size(SP.Y,1);
n_y = size(SP.Y,2);

beta = estruct.beta(xi);
delta = estruct.delta(xi);

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,op_SP.type_D);

n_b = size(D,1);

SP.Y_true = NaN(n_t,n_y,op_SP.nsamples);

is = 1;
while is <= op_SP.nsamples
    try
    bsample = mvnrnd(zeros(n_b,1),D);
    SP.Y_true(:,:,is) = nonfun(estruct.phi(beta,bsample));
    is = is + 1;
    catch
    end
end

if(any([op_SP.req(1),op_SP.req(1),op_SP.req(4),op_SP.req(5)]))
    SP.my_true = nanmean(SP.Y_true,3);
    DeltaY = bsxfun(@minus,SP.Y_true,SP.my_true);
end

sigma = estruct.sigma_noise(estruct.phi(beta,SP.B_SP(:,1)));

% adapt sigma to proper size
if(op_SP.req(2))
    if(size(sigma,1) == n_t)
        if(size(sigma,2) == 1)
            C_tech = bsxfun(@times,repmat(sigma.^2,[1,n_y,n_y]),permute(eye(n_y),[3,1,2]));
        elseif(size(sigma,2) == n_y)
            C_tech = bsxfun(@times,repmat(sigma.^2,[1,1,n_y]),permute(eye(n_y),[3,1,2]));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(size(sigma,2) == n_y)
        if(size(sigma,1) == 1)
            C_tech = bsxfun(@times,repmat(sigma.^2,[n_t,1,n_y]),permute(eye(n_y),[3,1,2]));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(and(size(sigma,1)==1,size(sigma,2)==1))
        C_tech = bsxfun(@times,repmat(sigma.^2,[n_t,n_y,n_y]),permute(eye(n_y),[3,1,2]));
    else
        error('Incompatible size of sigma parametrisation!')
    end
end


% adapt sigma to proper size
if(op_SP.req(5))
    if(size(sigma,1) == n_t)
        if(size(sigma,2) == 1)
            Cz_tech = diag(reshape(repmat(sigma.^2,[1,n_y]),n_t*n_y,1));
        elseif(size(sigma,2) == n_y)
            Cz_tech = diag(reshape(repmat(sigma.^2,[1,1]),n_t*n_y,1));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(size(sigma,2) == n_y)
        if(size(sigma,1) == 1)
            Cz_tech = diag(reshape(repmat(sigma.^2,[n_t,1]),n_t*n_y,1));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(and(size(sigma,1)==1,size(sigma,2)==1))
        Cz_tech = diag(reshape(repmat(sigma.^2,[n_t,n_y]),n_t*n_y,1));
    else
        error('Incompatible size of sigma parametrisation!')
    end
end

if(op_SP.req(2))
    SP.Cy_true = 1/op_SP.nsamples*sum(bsxfun(@times,permute(DeltaY,[1,2,4,3]),permute(DeltaY,[1,4,2,3])),4) ...
            + C_tech;
end

if(op_SP.req(3))
    % Cross-covariance
    SP.Cxy_true = 1/op_SP.nsamples*sum(bsxfun(@times,permute(SP.B_SP,[3,1,4,2]),permute(DeltaY,[1,4,2,3])),4);
end

if(op_SP.req(4))
    SP.mz_true = SP.my_true(:);
end

if(op_SP.req(5))
    DeltaZ = reshape(DeltaY,[size(DeltaY,1)*size(DeltaY,2),size(DeltaY,3)]);
    SP.Cz_true = 1/op_SP.nsamples*sum(bsxfun(@times,permute(DeltaZ,[1,3,2]),permute(DeltaZ,[3,1,2])),3) + Cz_tech;  
end


% Visualization

plotY(SP,[]);

if(op_SP.req(1))
   plotmy(SP.my,SP.my_true,[]) 
end

if(op_SP.req(2))
   plotCy(SP.Cy,SP.Cy_true,[]) 
end

if(op_SP.req(3))
   plotCy(SP.Cxy,SP.Cxy_true,[]) 
end

if(op_SP.req(4))
   plotmz(SP.mz,SP.mz_true,[]) 
end

if(op_SP.req(5))
   plotCz(SP.Cz,SP.Cz_true,[]) 
end



end

