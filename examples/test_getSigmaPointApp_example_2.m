clear all;
close all;
clc;

%% Example
% Parameterization
%% Example
% Parameterization
A = eye(3);
B = [eye(2);
     zeros(1,2)];

% xi
xi = [1;2;0.5;2;1;3];

% beta
beta_fun = @(xi) [xi(1);
                  xi(2);
                  xi(3)];
beta = beta_fun(xi);
dbetadxi = [1,0,0,0,0,0;
            0,1,0,0,0,0;
            0,0,1,0,0,0];

% D

D_fun = @(xi) [xi(4),xi(5);
               xi(5),xi(6)];
D = D_fun(xi);
dDddelta = zeros(2,2,3);
dDddelta(:,:,1) = [1,0;
                   0,0];  
dDddelta(:,:,2) = [0,1;
                   1,0];  
dDddelta(:,:,3) = [0,0;
                   0,1];
ddeltadxi = [0,0,0,1,0,0;
             0,0,0,0,1,0;
             0,0,0,0,0,1];


estruct.phi = @(beta,b) A*beta + B*b;
estruct.dphidbeta = @(beta,b) A;
estruct.dphidb = @(beta,b) B;
estruct.sigma_noise = @(phi) zeros(4,1);
estruct.dsigma_noisedphi = @(phi) zeros(4,1,3);

% phi
phi = @(beta,b) A*beta + B*b;

%% Sigma point approximation
op_SP.nderiv = 1;
op_SP.req = [1,1,1,1,1];
SP = getSigmaPointApp(@(phi) nonfun_2(phi),beta_fun(xi),D_fun(xi),estruct,op_SP,dDddelta,ddeltadxi,dbetadxi);

eps = 1e-6;
for k = 1:length(xi)
    xi_k = xi;
    xi_k(k) = xi_k(k) + eps;
    op_SP.nderiv = 0;
    op_SP.req = [1,1,1,1,1];
    SP_k{k} = getSigmaPointApp(@(phi) nonfun_2(phi),beta_fun(xi),D_fun(xi_k),estruct,op_SP);
    
%     (B_SP_k{k}(:,1) - B_SP(:,1))/eps - dB_SPdxi(:,k,1)
%     (Y_k{k}(:,:,1) - Y(:,:,1))/eps - dYdxi(:,:,k,1)
%     (my_k{k} - my)/eps - dmydxi(:,:,k)
%     for j = 1:4
%         squeeze((Cy_k{k}(j,:,:) - Cy(j,:,:))/eps - dCydxi(j,:,:,k))
%     end
%     for j = 1:4
%         squeeze(Cxy(j,:,:))
%         squeeze((Cxy_k{k}(j,:,:) - Cxy(j,:,:))/eps - dCxydxi(j,:,:,k))
%     end
end

%% Sampling approximation
n_b = size(B,2);

N = 1e3;
Bs = sqrtm(D)*randn(n_b,N);
for j = 1:N
    [Ys(:,:,j),dYs(:,:,:,j)] = nonfun_2(estruct.phi(beta,Bs(:,j)));
end

% Mean
mys = mean(Ys,3);
n_t = size(Ys,1);

% Covariance
for j = 1:size(mys,1)
    Cys(j,:,:) = cov(squeeze(Ys(j,:,:))');
end

% Cross-covariance
Cxys = zeros(size(mys,1),n_b,size(mys,2));
for k = 1:size(mys,1)
    for j = 1:N
        Cxys(k,:,:) = Cxys(k,:,:) + permute(1/N*Bs(:,j)*(Ys(k,:,j)-mys(k,:)),[3,1,2]);
    end
end

% Full state covariance
Zs = reshape(Ys,[size(Ys,1)*size(Ys,2),size(Ys,3)]);
mzs = mean(Zs,2);
Czs = cov(Zs');

%% Visualization
% State-State
I = [1,2;
     1,3;
     1,4;
     1,5;
     2,3;
     2,4;
     2,5;
     3,3;
     3,4;
     ];

figure;
for k = 1:size(mys,1)
    for l = 1:size(I,1)
        % Open subplot
        subplot(size(mys,1),size(I,1),(k-1)*size(I,1)+l);
        
        % Plot samples
        plot(squeeze(Ys(k,I(l,1),:)),squeeze(Ys(k,I(l,2),:)),'b.'); hold on;

        mysk = mys(k,I(l,:))';
        Cysk = [Cys(k,I(l,1),I(l,1)),Cys(k,I(l,1),I(l,2));
                Cys(k,I(l,1),I(l,2)),Cys(k,I(l,2),I(l,2))];
        X = getEllipse(mysk,Cysk,2);

        plot(mysk(1),mysk(2),'co','linewidth',3);
        plot(X(1,:),X(2,:),'c-','linewidth',2);
        
        % Plot approximation
        myk = SP.my(k,I(l,:))';
        Cyk = [SP.Cy(k,I(l,1),I(l,1)),SP.Cy(k,I(l,1),I(l,2));
               SP.Cy(k,I(l,1),I(l,2)),SP.Cy(k,I(l,2),I(l,2))];
        X = getEllipse(myk,Cyk,2);

        plot(myk(1),myk(2),'ro','linewidth',3);
        plot(X(1,:),X(2,:),'r-','linewidth',2);

        plot(squeeze(SP.Y(k,I(l,1),:)),squeeze(SP.Y(k,I(l,2),:)),'rx','linewidth',2);

        % Title
        title(['t = ' num2str(k) ', Dim = {' num2str(I(l,1)) ',' num2str(I(l,2)) '}']);
    end
end

% State-State for different time points
I = [1,2;
     1,3;
     1,4;
     1,5;
     2,3;
     2,4;
     2,5;
     3,3;
     3,4;
     ];
It = [1,1;
      2,2;
      3,3;
      4,4;
      ];

figure;
for k = 1:length(It)
    for l = 1:size(I,1)
        % Open subplot
        subplot(length(It),size(I,1),(k-1)*size(I,1)+l);
        
        ind1 = It(k,1)+n_t*(I(l,1)-1);
        ind2 = It(k,2)+n_t*(I(l,2)-1);

        % Plot samples
        plot(Zs(ind1,:),Zs(ind2,:),'b.'); hold on;

        mzsk = mzs([ind1,ind2]);
        Czsk = [Czs(ind1,ind1),Czs(ind1,ind2);
                Czs(ind1,ind2),Czs(ind2,ind2)];
        X = getEllipse(mzsk,Czsk,2);

        plot(mzsk(1),mzsk(2),'co','linewidth',3);
        plot(X(1,:),X(2,:),'c-','linewidth',2);
        
        % Plot approximation
        mzk = SP.mz([ind1,ind2]);
        Czk = [SP.Cz(ind1,ind1),SP.Cz(ind1,ind2);
               SP.Cz(ind1,ind2),SP.Cz(ind2,ind2)];
        X = getEllipse(mzk,Czk,2);

        plot(mzk(1),mzk(2),'ro','linewidth',3);
        plot(X(1,:),X(2,:),'r-','linewidth',2);

        % Title
        title(['t = ' num2str(k) ', Dim = {' num2str(I(l,1)) ',' num2str(I(l,2)) '}']);
    end
end

% Parameter-State
I = [1,2;
     1,3;
     1,4;
     1,5;
     2,1;
     2,2;
     2,3;
     2,4;
     2,5;
     ];

figure;
for k = 1:size(mys,1)
    for l = 1:size(I,1)
        % Open subplot
        subplot(size(mys,1),size(I,1),(k-1)*size(I,1)+l);
        
        % Plot samples
        plot(squeeze(Bs(I(l,1),:)),squeeze(Ys(k,I(l,2),:)),'b.'); hold on;

        mxysk = [0;mys(k,I(l,2))];
        Cxysk = [D(I(l,1),I(l,1)),Cxys(k,I(l,1),I(l,2));
                 Cxys(k,I(l,1),I(l,2)),Cys(k,I(l,2),I(l,2))];
        X = getEllipse(mxysk,Cxysk,2);

        plot(mxysk(1),mxysk(2),'co','linewidth',3);
        plot(X(1,:),X(2,:),'c-','linewidth',2);
        
        % Plot approximation
        mxyk = [0;SP.my(k,I(l,2))];
        Cxyk = [D(I(l,1),I(l,1)),SP.Cxy(k,I(l,1),I(l,2));
                 SP.Cxy(k,I(l,1),I(l,2)),SP.Cy(k,I(l,2),I(l,2))];
        X = getEllipse(mxyk,Cxyk,2);

        plot(mxyk(1),mxyk(2),'ro','linewidth',3);
        plot(X(1,:),X(2,:),'r-','linewidth',2);

        plot(squeeze(SP.B_SP(I(l,1),:)),squeeze(SP.Y(k,I(l,2),:)),'rx','linewidth',2);

        % Title
        title(['t = ' num2str(k) ', Dim = {' num2str(I(l,1)) ',' num2str(I(l,2)) '}']);
    end
end

% %% Testing of sensitivities
% for i = 1:length(xi)
%     xi_i = xi;
%     xi_i(i) = xi_i(i) + eps;
%     
%     [my_i,dmydxi_i,Cy_i,dCydxi_i,Cxy_i,dCxydxi_i,mz_i,dmzdxi_i,Cz_i,dCzdxi_i,B_SP_i,dB_SPdxi_i,Y_i,dYdxi_i] ...
%         = getSigmaPointApp(A,B,beta_fun(xi_i),dbetadxi,D_fun(xi_i),dDdxi,@(phi) nonfun_2(phi));
%     
%     dmzdxi_fd(:,i) = (mz_i - mz)/eps;
%     dCzdxi_fd(:,:,i) = (Cz_i - Cz)/eps;
%     
%     figure;
%     subplot(1,2,1);
%     imagesc(dCzdxi(:,:,i))
%     colorbar
% 
%     subplot(1,2,2);
%     imagesc(dCzdxi_fd(:,:,i)-dCzdxi(:,:,i))
%     colorbar
% end

