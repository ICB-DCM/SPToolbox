function varargout = logL_SP_w_grad(varargin)

    persistent tau
    persistent P_old
    persistent logL_old
    
    if isempty(tau)
        tau = clock;
    end
    if isempty(logL_old)
        logL_old = -inf;
    end
    
    
    %% Initialization
    theta = varargin{1};
    Data = varargin{2};
    Model = varargin{3};
    
    % Options
    options.tau_update = 0;
    options.plot = 1;
    if nargin >= 4
        options = setdefault(varargin{4},options);
    end
    
    if nargin >= 5
        extract_flag = varargin{5};
    end
    
    % Plot options
    if (etime(clock,tau) > options.tau_update) && (options.plot == 1)
        options.plot = 30;
        tau = clock;
    else
        options.plot = 0;
    end
    
    %% Evaluation of likelihood function
    % Initialization
    logL = 0;
    dlogLdtheta = zeros(length(theta),1);
    ddlogLdtheta2 = zeros(length(theta));
    
    data_type = {'SCTL','PA','SCSH'};
    
    for s = 1:length(Data)
        
        %% Assignment of global variables
        A = Model.exp{s}.A;
        B = Model.exp{s}.B;
        ind_beta = Model.exp{s}.ind_beta;
        ind_D = Model.exp{s}.ind_D;
        n_beta = length(ind_beta);
        n_D = length(ind_D);
        n_b = size(B,2);
        type_D = Model.type_D;
        
        %% Construct fixed effects and covariance matrix
        beta = theta(ind_beta);
        dbeta = eye(length(theta)); dbeta = dbeta(ind_beta,:);
        [D,invD,dD,dinvD,~,HinvD] = xi2D(theta(ind_D),type_D);
        dD_full = zeros(size(dD,1),size(dD,1),length(theta));
        dD_full(:,:,ind_D) = dD;
        
        %% Construction of time vector
        t_s = [];
        for dtype = 1:length(data_type)
            if isfield(Data{s},data_type{dtype})
                t_s = union(eval(['Data{s}.' data_type{dtype} '.time']),t_s);
            end
        end
        
        if isfield(Data{s},'SCTL')
            % Evaluation of time index set
            [~,ind_t] = ismember(Data{s}.SCTL.time,t_s);
            
            
            N = size(Data{s}.SCTL.Y,3);
            Ym = Data{s}.SCTL.Y(ind_t,:,:);
            
            kappa = [];
            
            [mY,CY,Cxy,mz,Cz,B_SP,Y,dmYdxi,dCYdxi,dCxydxi,dmzdxi,dCzdxi,dB_SPdxi,dYdxi] = getSigmaPointApp(@(phi) nonfun(phi,@(x)Model.exp{s}.model(Data{s}.SCTL.time,x,kappa)),A,B,beta,D,dbeta,dD_full);
            
            
            beta_f = beta;
            beta_f(1) = beta_f(1) + 1e-4;
            
            [mY_f,CY,Cxy,mz,Cz,B_SP,Y,dmYdxi_f,dCYdxi,dCxydxi,dmzdxi,dCzdxi,dB_SPdxi,dYdxi] = getSigmaPointApp(@(phi) nonfun(phi,@(x)Model.exp{s}.model(Data{s}.SCTL.time,x,kappa)),A,B,beta_f,D,dbeta,dD_full);
            
            mYm = mean(Ym,3);
            CYm = zeros(size(Ym,1),size(Ym,2),size(Ym,2));
            Sigma_mY = zeros(size(Ym,1),size(Ym,2));
            Sigma_CY = zeros(size(Ym,1),size(Ym,2));
            for s = 1:size(Ym,1)
                CYm(s,:) = var(Ym(s,:,:));
                Sigma_mY(s,:) = sqrt(CYm(s,:,:)/N);
                Sigma_CY(s,:) = sqrt(mean(bsxfun(@minus,Ym(s,:,:),mYm(s,:)).^4,3)...
                    -(N-3)/(N-1)*CYm(s,:).^2)/sqrt(N);
            end
               
            for j = 1:size(Y,2)
                res = (mY(ind_t,j)-mYm)./Sigma_mY(:,j) + (CY(ind_t,j)-CYm(:,j))./Sigma_CY(:,j);
                dres = squeeze(bsxfun(@times,1./Sigma_mY(:,j),dmYdxi(:,j,:))) + bsxfun(@times,1./Sigma_CY(:,j),squeeze(dCYdxi(:,j,j,:)));
                
                logL = logL - 0.5*res'*res;
                dlogLdtheta = dlogLdtheta - dres'*res;
                ddlogLdtheta2 = ddlogLdtheta2 - dres'*dres;
            end
        end
    end
    
    if extract_flag
       varargout{1} = B_SP;
       return
    end
    
    %% Output
    if nargout <= 1
        % One output
        varargout{1} =  logL;
    elseif nargout <= 2
        % Two outputs
        varargout{1} =  logL;
        varargout{2} =  dlogLdtheta;
    else
        % Two outputs
        varargout{1} =  logL;
        varargout{2} =  dlogLdtheta;
        varargout{3} =  ddlogLdtheta2;
    end
    
end

function varargout = nonfun(phi,model)
   
    if nargout == 1
        [~,~,~,y] = model(exp(phi));
        varargout{1} = y;
    elseif nargout == 2
        [~,~,~,y,~,sy] = model(exp(phi));
        varargout{1} = y;
        varargout{2} = bsxfun(@times,sy,exp(permute(phi,[3,2,1])));
    end
    
    
end
