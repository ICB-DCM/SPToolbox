% USAGE:
% ======
% [SP] = getSigmaPointApp(nonfun,xi,estruct,op_SP)
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
% xi  ... hyperparameters for random and common effects
% estruct ... struct containing all information about the employed effect
%            composition
%    .beta ... @beta(xi) defines the how the common effects are constructed
%            from hyperparameter xi
%    .dbetadxi ... @dbetadxi(xi) derivative of the common effects beta with respect to the
%            the hyperparameter xi; only required if op_SP.nderiv>1
%    .delta ... @delta(xi) defines the how the random effect
%            hyperparameters are constructed from hyperparameter xi
%    .ddeltadxi ... derivative of the parametrisation delta with respect to the
%            the hyperparameters xi; only required if op_SP.nderiv>1
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
%    .type_D ... string specifying the parametrisation of the covariance
%                matrix for the random effects. either
%            'diag-matrix-logarithm' for diagonal matrix with log. parametrisation or
%            'matrix-logarithm' for full matrix with log. parametrisation
%    .approx ... string specifying the approximation method. either
%            'sp' for sigma_points
%            'sampled' for sampling based approximation
%            'dmd' for approximation by dirac mixtrue distribution
%    .n_samples ... number of dirac mixture components, only defined if op_SP.approx='dmd'
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
%
% OPTIONAL OUTPUTS (only computed when derivatives are computed):
% =======
% SP.dmydxi ... derivative of my w.r.t. xi       
%            [time x observables x parameters]
% SP.dCydxi ... derivative of Cy w.r.t. xi   
%            [time x observables x observables x parameters]
% SP.dCxydxi... derivative of Cxy w.r.t. xi
%            [time x time x observables x observables x parameters]
% SP.dmzdxi ... derivative of mz w.r.t. xi 
%            [time*observables x parameters]
% SP.dCzdxi ... derivative of Cz w.r.t. xi 
%            [time*observables x time*observables x parameters]
% SP.dB_SPdxi... derivative of B_SP w.r.t. xi 
%            [parameters x 2*parameters+1 x parameters]
% SP.dYdxi  ... derivative of Y w.r.t. xi 
%            [time x observables x 2*parameters+1 x parameters]
% REFERENCES:
% =======
% --'dmd': I. Gilitschenski and U. D. Hanebeck, "Efficient Deterministic
% Dirac Mixture Approximation of Gaussian Distributions", IEEE,2013.
% --'Julier1': S. J. Julier, J. K. Uhlmann and H. F. Durrant-Whyte, "A new
% approach for filtering nonlinear systems", IEEE American Control Conf., 
% 1995.
% --'Julier2': S. J. Julier and J. K. Uhlmann, "Unscented filtering and
% nonlinear estimation", IEEE, 2004.
% --'Menegaz': H. M. Menegaz, J. Y. Ishihara and G. A. Borges, "A new
% smallest sigma set for the Unscented Transform and its applications on
% SLAM", IEEE, 2011.
% --'Lerner': U. N. Lerner, "Hybrid Bayesian Networks for Reasoning About
% Complex Systems", Ph. D., Stanford University, 2002.
% --'Charalampidis': A. C. Charalampidis and G. P. Papavassilopoulos,
% "Development and numerical investigation of new non-linear Kalman filter
% variants", IET Control Theory&Applications, 2011.

function SP = getSigmaPointApp(varargin)

% phi = phi(beta,b), with b ~ N(0,D)
%
%    beta = beta(xi)
%    D = D(delta)
%    delta = delta(xi)
%
% Goal: Propagation of distribution in b

nonfun = varargin{1};
xi  = varargin{2};
estruct = varargin{3};
op_SP = varargin{4};
if op_SP.nderiv
    compute_derivative = 1;
else
    compute_derivative = 0;
end

beta = estruct.beta(xi);
delta = estruct.delta(xi);

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,op_SP.type_D);

n_b = size(D,1);

if(~isfield(op_SP,'approx'))
    op_SP.approx = 'sp';
end
if(strcmp(op_SP.approx,'samples'))
if(~isfield(op_SP,'samples'))
    persistent samples
    if(isempty(samples))
        samples = mvnrnd(zeros(n_b,1),eye(n_b),10000);
    end
    op_SP.samples = samples;
    end
end

%% Initialization
% Dimensions



%% Parameters of sigma point approach
L = n_b;

% Recommendations in:
%       R. van der Merwe, Sigma-point Kalman filters for probabilistic
%       inference in dynamic state-space models, Ph.D. thesis, 2004.
%
% 1) kap >= 0 to guarantee positive semi-definiteness of the covariance
%    matrix. The specific value of kap is not critical though, so a 
%    good default choice is kap= 0.
% 2) Choose 0 <= alp <= 1. 
%    alp controls the "size" of the sigma-point distribution and should
%    ideally be a small number to avoid sampling non-local effects when
%    the nonlinearities are strong.
% 3) Choose bet >= 1. bet is a non-negative weighting term which can be
%    used to incorporate knowledge of the higher order moments of the
%    distribution. For a Gaussian prior the optimal choice is bet = 2.

alp = 0.7;
bet = 2;
kap = 0;
lam = alp^2*(L+kap) - L;



if compute_derivative == 1
    n_xi = size(estruct.ddeltadxi(xi),2);
end

%% Sigma points
% Matrix root: D = S*S'
if compute_derivative == 0
    [S] = chol(D,'lower');
else
    [S,dSddelta] = chol_w_diff(D,dDddelta);
    dSdxi = permute(sum(bsxfun(@times,dSddelta,permute(estruct.ddeltadxi(xi),[3,4,1,2])),3),[1,2,4,3]);
end

switch(op_SP.approx)
    case 'dmd'
        diracD = size(D,1);
        [SPToolboxFolder,~,~] = fileparts(which('CompDMD_Location'));
        filename = sprintf('%s%i%s%i%s','DMDTrueMeaninfo\B_SP_dim',diracD,'points',op_SP.n_samples,'.csv');
        if (~exist(fullfile(SPToolboxFolder,filename),'file'))
            %dimension of dirac mixture distribution
            error('The approximation does not exist. Please run CompDMD_Location first!')
        else
            %Dirac Mixture location for normal distribution
            B_SPNorm = importdata(fullfile(SPToolboxFolder,filename));
        end
        SP.B_SP = S*B_SPNorm;
        if compute_derivative == 1
            SP.dB_SPdxi = permute(sum(bsxfun(@times,B_SPNorm,permute(dSdxi,[2,4,1,3])),1),[3,4,2,1]);
        end
        % Weights
        w_m = 1/(size(SP.B_SP,2))*ones(size(SP.B_SP,2),1);
        w_c = 1/((size(SP.B_SP,2))-1)*ones(size(SP.B_SP,2),1);
    case 'samples'
        SP.B_SP = transpose(op_SP.samples*S);
        if compute_derivative == 1
            SP.dB_SPdxi = permute(sum(bsxfun(@times,op_SP.samples,permute(dSdxi,[4,1,2,3])),2),[3,4,1,2]);
        end
        % Weights
        w_m = 1/(size(SP.B_SP,2))*ones(size(SP.B_SP,2),1);
        w_c = 1/((size(SP.B_SP,2))-1)*ones(size(SP.B_SP,2),1);
    case 'halton'    % Halton Monte Carlo sequence
        samplescdf = net(haltonset(2,'skip',100),op_SP.n_samples);
        op_SP.samples = norminv(samplescdf,0,1);
        SP.B_SP = transpose(op_SP.samples*S);
        if compute_derivative == 1
            SP.dB_SPdxi = permute(sum(bsxfun(@times,op_SP.samples,permute(dSdxi,[4,1,2,3])),2),[3,4,1,2]);
        end
        % Weights
        w_m = 1/(size(SP.B_SP,2))*ones(size(SP.B_SP,2),1);
        w_c = 1/((size(SP.B_SP,2))-1)*ones(size(SP.B_SP,2),1);
    case 'sobol'    % Sobol Monte Carlo sequence
        samplescdf = net(sobolset(2,'skip',100),op_SP.n_samples);
        op_SP.samples = norminv(samplescdf,0,1);
        SP.B_SP = transpose(op_SP.samples*S);
        if compute_derivative == 1
            SP.dB_SPdxi = permute(sum(bsxfun(@times,op_SP.samples,permute(dSdxi,[4,1,2,3])),2),[3,4,1,2]);
        end
        % Weights
        w_m = 1/(size(SP.B_SP,2))*ones(size(SP.B_SP,2),1);
        w_c = 1/((size(SP.B_SP,2))-1)*ones(size(SP.B_SP,2),1);        
    case 'Julier1'    % N = 2n+1
        SP.B_SP = [zeros(L,1),sqrt(L+kap)*S,-sqrt(L+kap)*S];
        if compute_derivative == 1
            SP.dB_SPdxi = permute([zeros(L,1,n_xi),sqrt(L+kap)*dSdxi,-sqrt(L+kap)*dSdxi],[1,3,2]);
        end
        % Weights
        w0_m = kap/(L+kap);
        wi_m = 1/(2*(L+kap));
        w_m = [w0_m;wi_m*ones(2*L,1)];
        w_c = w_m;        
    case 'Julier2'    % N = 2n+1
        w0 = -1;
        SP.B_SP = [zeros(L,1),sqrt(L/(1-w0))*S,-sqrt(L/(1-w0))*S];
        if compute_derivative == 1
            SP.dB_SPdxi = permute([zeros(L,1,n_xi),sqrt(L/(1-w0))*dSdxi,-sqrt(L/(1-w0))*dSdxi],[1,3,2]);
        end
        % Weights
        wi_m = (1-w0)/(2*L);
        w_m = [w0,wi_m*ones(1,2*L)]';
        w_c = w_m;        
    case 'Menegaz'    % N = n+1
        wl = 0.5;
        C = chol(eye(L)-repmat((1-wl)/L,L),'lower');
        wi = diag(inv(C)*wl*repmat((1-wl)/L,L)*inv(C'));
        SP.B_SP = [S*C*inv(chol(diag(wi),'lower')),-sqrt((1-wl)/L)*S*(ones(L,1)/sqrt(wl))];
        if compute_derivative == 1
            SP.dB_SPdxi =  permute([permute(sum(bsxfun(@times,inv(chol(diag(wi))),...
                permute(sum(bsxfun(@times,C,permute(dSdxi,[2,4,1,3])),1),[2,1,3,4])),1),[3,2,4,1]),...
                -sqrt((1-wl)/L)*permute(sum(bsxfun(@times,(ones(L,1)/sqrt(wl)),...
                permute(dSdxi,[2,4,1,3])),1),[3,2,4,1])],[1,3,2]);
        end
        % Weights
        w_m = [wi',wl]';
        w_c = w_m;        
    case 'Lerner'    % N = 2n^2+1
        gen1 = [perms([sqrt(3),zeros(1,L-1)])',perms([-sqrt(3),zeros(1,L-1)])'];
        gen2 = [perms([sqrt(3),sqrt(3),zeros(1,L-2)])',...
        perms([sqrt(3),-sqrt(3),zeros(1,L-2)])',...
        perms([-sqrt(3),sqrt(3),zeros(1,L-2)])',...
        perms([-sqrt(3),-sqrt(3),zeros(1,L-2)])'];
        gen = unique([gen1,gen2,zeros(L,1)]','rows','stable')';
        SP.B_SP = S*gen;
        if compute_derivative == 1
            SP.dB_SPdxi = permute(sum(bsxfun(@times,gen,permute(dSdxi,[2,4,1,3])),1),[3,4,2,1]);
        end
        % Weights
        w_m = [(4-L)/18*ones(1,2*L),1/36*ones(1,2*L^2-2*L),(L^2-7*L)/18+1]';
        w_c = w_m;
    case 'Charalampidis'    % N = phi^n
        phi = op_SP.n_samples;
        i = 1:phi;
        yi = norminv((2*i-1)/(2*phi));
        xi = sqrt(phi/sum(yi.^2))*yi;
        %Sigma point for normal distribution
        B_SPNorm = unique(nchoosek(repmat(xi,[1,L]),L),'rows','stable')';
        SP.B_SP = S*B_SPNorm;
        if compute_derivative == 1
            SP.dB_SPdxi = permute(sum(bsxfun(@times,B_SPNorm,permute(dSdxi,[2,4,1,3])),1),[3,4,2,1]);
        end
        % Weights
        w_m = (ones(1,phi^L)*1/phi^L)';
        w_c = w_m;
    otherwise
        % Sigma points
        SP.B_SP = [zeros(L,1),sqrt(L+lam)*S,-sqrt(L+lam)*S];
        if compute_derivative == 1
            SP.dB_SPdxi = permute([zeros(L,1,n_xi),sqrt(L+lam)*dSdxi,-sqrt(L+lam)*dSdxi],[1,3,2]);
        end
        % Weights
        w0_m = lam/(L+lam);
        wi_m = 1/(2*(L+lam));
        w_m = [w0_m;wi_m*ones(2*L,1)];
        
        w0_c = lam/(L+lam)+(1-alp^2+bet);
        wi_c = 1/(2*(L+lam));
        w_c = [w0_c;wi_c*ones(2*L,1)];
end

%% Propagation of sigma points
% Loop: Sigma points
% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(estruct.phi(beta,SP.B_SP(:,1)),@(phi) nonfun(phi),1e-4,1,2);

for i = 1:size(SP.B_SP,2)
    if compute_derivative == 0
        Y(:,:,i) = nonfun(estruct.phi(beta,SP.B_SP(:,i)));
    else
        [Y(:,:,i),dYdphi(:,:,:,i)] = nonfun(estruct.phi(beta,SP.B_SP(:,i)));
        dphidxi(:,:,i) = estruct.dphidbeta(beta,SP.B_SP(:,i))*estruct.dbetadxi(xi) + estruct.dphidb(beta,SP.B_SP(:,i))*SP.dB_SPdxi(:,:,i);
        dYdxi(:,:,:,i) = permute(sum(bsxfun(@times,dYdphi(:,:,:,i),permute(dphidxi(:,:,i),[4,3,1,2])),3),[1,2,4,3]);
    end
end

SP.Y = Y;
if(any(isnan(Y)))
    error('Failed to successfully integrate system at all SigmaPoints')
end

[n_t,n_y,~] = size(Y);


%% Evaluation of mean, covariance and cross-covariance
% Mean
if(any([op_SP.req(1),op_SP.req(1),op_SP.req(4),op_SP.req(5)]))
    SP.my = sum(bsxfun(@times,permute(w_m,[3,2,1])  ,Y)    ,3);
    DeltaY = bsxfun(@minus,Y,SP.my);
    if compute_derivative == 1
        SP.dmydxi = sum(bsxfun(@times,permute(w_m,[4,3,2,1]),dYdxi),4);
        DeltadYdxi = bsxfun(@minus,dYdxi,SP.dmydxi);
    end
end

% Covariance
% for measurement noise we ignore the random effects at this point
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
if compute_derivative == 1
    dsigmadphi = estruct.dsigma_noisedphi(estruct.phi(beta,SP.B_SP(:,1)));
    ndim_dsigmadphi = ndims(dsigmadphi);
    dsigmadxi = permute(sum(bsxfun(@times,dsigmadphi,permute(dphidxi(:,:,1),[3:(3+ndim_dsigmadphi-2),1,2])),ndim_dsigmadphi),[1:ndim_dsigmadphi-1,ndim_dsigmadphi+1,ndim_dsigmadphi]);
    dC_techdxi = zeros(n_t,n_y,n_y,size(dsigmadxi,ndim_dsigmadphi));
    
    if(size(dsigmadxi,1) == n_t)
        if(size(dsigmadxi,2) == 1)
            dC_techdxi = bsxfun(@times,repmat(permute(dsigmadxi,[1,2,4,3]),[1,n_y,n_y,1]),permute(eye(n_y),[3,1,2]));
        elseif(size(dsigmadxi,2) == n_y)
            dC_techdxi = bsxfun(@times,repmat(permute(dsigmadxi,[1,2,4,3]),[1,1,n_y,1]),permute(eye(n_y),[3,1,2]));
        end
    elseif(size(dsigmadxi,2) == n_y)
        if(size(dsigmadxi,1) == 1)
            dC_techdxi = bsxfun(@times,repmat(permute(dsigmadxi,[1,2,4,3]),[n_t,1,n_y,1]),permute(eye(n_y),[3,1,2]));
        end
    elseif(and(size(dsigmadxi,1)==1,size(dsigmadxi,2)==1))
        dC_techdxi = bsxfun(@times,repmat(permute(dsigmadxi,[1,2,4,3]),[n_t,n_y,n_y,1]),permute(eye(n_y),[3,1,2]));
    end
end
end


if(op_SP.req(2))
    SP.Cy = sum(bsxfun(@times,permute(w_c,[4,3,2,1]),...
            bsxfun(@times,permute(DeltaY,[1,2,4,3]),permute(DeltaY,[1,4,2,3]))),4) ...
            + C_tech;
    
    if compute_derivative == 1
        dCj =  sum(bsxfun(@times,permute(w_c,[5,4,3,2,1]),...
            bsxfun(@times,permute(DeltaY,[1,2,4,5,3]),permute(DeltadYdxi,[1,5,2,3,4]))),5);
        
        SP.dCydxi = dCj + permute(dCj,[1,3,2,4]) ...
            + dC_techdxi;
        
    end
end

if(op_SP.req(3))
    % Cross-covariance
    SP.Cxy = sum(bsxfun(@times,permute(w_c,[4,3,2,1]),...
        bsxfun(@times,permute(SP.B_SP,[3,1,4,2]),permute(DeltaY,[1,4,2,3]))),4);
    
    if compute_derivative == 1
        SP.dCxydxi = sum(bsxfun(@times,permute(w_c,[5,4,3,2,1]),...
            bsxfun(@times,permute(SP.B_SP,[3,1,4,5,2]),permute(DeltadYdxi,[1,5,2,3,4]))),5) ...
            + sum(bsxfun(@times,permute(w_c,[5,4,3,2,1]),...
                bsxfun(@times,permute(SP.dB_SPdxi,[4,1,5,2,3]),permute(DeltaY,[1,4,2,5,3]))),5);
    end
end

% Full state-time covariance
if(op_SP.req(4))
    SP.mz = SP.my(:);
    if compute_derivative == 1
        SP.dmzdxi = reshape(SP.dmydxi,[size(SP.dmydxi,1)*size(SP.dmydxi,2),size(SP.dmydxi,3)]);
    end
end

if(op_SP.req(5))
    DeltaZ = reshape(DeltaY,[size(DeltaY,1)*size(DeltaY,2),size(DeltaY,3)]);
    SP.Cz = sum(bsxfun(@times,permute(w_c,[3,2,1]),bsxfun(@times,permute(DeltaZ,[1,3,2]),permute(DeltaZ,[3,1,2]))),3) + Cz_tech;
    
    if compute_derivative == 1
        DeltadZdxi = reshape(DeltadYdxi,[size(DeltadYdxi,1)*size(DeltadYdxi,2),size(DeltadYdxi,3),size(DeltadYdxi,4)]);
        dCzdxi = sum(bsxfun(@times,permute(w_c,[4,3,2,1]),bsxfun(@times,permute(DeltaZ,[1,4,3,2]),permute(DeltadZdxi,[4,1,2,3]))),4);
        SP.dCzdxi = dCzdxi + permute(dCzdxi,[2,1,3]);
    end
end

