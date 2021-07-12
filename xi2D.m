function [D,invD,dD,dinvD,HD,HinvD] = xi2D(xi,type)

% Covariance matrix D
switch type
    case 'diag-matrix-logarithm'
        % Diagonal matrix logarithm parameterization
        m = length(xi);
        D = diag(exp(xi));
        invD = diag(exp(-xi));
        dD = zeros(m,m,m);
        dinvD = zeros(m,m,m);
        HD = zeros(m,m,m,m);
        HinvD = zeros(m,m,m,m);
        for i = 1:length(xi)
            dD(i,i,i) = exp(xi(i));
            dinvD(i,i,i) = -exp(-xi(i));
            HD(i,i,i,i) = exp(xi(i));
            HinvD(i,i,i,i) = exp(-xi(i));
        end

    case 'Givens-parametrization'
        % Givens Parametrization
        % Initialize and get Covariance dimension from xi vector length
        m = round((-1 + sqrt(1+8*length(xi)))/2);
        G = zeros(m,m,m*(m-1)/2);
        dG = zeros(m,m,m*(m-1)/2);
        U = diag(ones(1,m));
        dU = repmat(diag(ones(1,m)),1,1,m*(m-1)/2);
        dLambda = repmat(zeros(m,m),1,1,m);
        dD = zeros(m,m,length(xi));
        dinvD = zeros(m,m,length(xi));
        HD = zeros(m,m,length(xi),length(xi));
        HinvD = zeros(m,m,length(xi),length(xi));
        invD = zeros(m,m);
        lambda = zeros(1,m);

        % Separate desired eigenvalues lambda and rotation angles delta
        % lambda should be log-scaled while delta is \in [0,\pi]. The
        % eigenvalues are stored as lambda = log(xi_i-xi_i-1) to enforce

        % Separate desired eigenvalues lambda and rotation angles delta
        lambda = exp(xi(1:m));
        dlambda = exp(xi(1:m));

        % Angles
        delta  = xi(m+1:end);

        % Get matrix angle position indices, to know where to place cos, sin, -sin
        m1m2 = zeros(m,2);
        for k = 1:m*(m-1)/2
            for m2 = 1:m
                for m1 = 1:m2-1
                    if k == m2 - m1 + ( m1 - 1 ) * ( m - m1/2 )
                        m1m2(k,1) = m1;
                        m1m2(k,2) = m2;
                    end
                end
            end
        end

        % Build eigenvector matrix U from n(n-1)/2 single rotations matrices G
        for k = 1:m*(m-1)/2
            for i = 1:m
                for j = 1:m
                    if (i == m1m2(k,1) && j == m1m2(k,1)) || (i == m1m2(k,2) && j == m1m2(k,2))
                        G(i,j,k) = cos(delta(k));
                        dG(i,j,k) = -sin(delta(k));
                    elseif (i == m1m2(k,1) && j == m1m2(k,2))
                        G(i,j,k) = sin(delta(k));
                        dG(i,j,k) = cos(delta(k));
                    elseif (i == m1m2(k,2) && j == m1m2(k,1))
                        G(i,j,k) = -sin(delta(k));
                        dG(i,j,k) = -cos(delta(k));
                    elseif (i == j) && (i ~= m1m2(k,1)) && (m1m2(k,2))
                        G(i,j,k) = 1;
                        dG(i,j,k) = 0;
                    else
                        G(i,j,k) = 0;
                        dG(i,j,k) = 0;
                    end
                end
            end
        end

        % get the derivative of the rotation angles dU/dxi
        for k = m*(m-1)/2:-1:1
            U = squeeze(G(:,:,k)) * U;
            for l = m*(m-1)/2:-1:1
                if k == l
                    dU(:,:,k) = squeeze(dG(:,:,l)) * dU(:,:,k);
                else
                    dU(:,:,k) = squeeze(G(:,:,l)) * dU(:,:,k);
                end
            end
        end

        % get the derivative of the eigenvalue component
        for k = 1:m
            dLambda(k,k,k) = dlambda(k);
        end

        % Build D = U^t Lambda U
        D = U' * diag(lambda) * U;
        for k = 1:m*(m+1)/2
            if k <= m
                dD(1:m,1:m,k) = U' * squeeze(dLambda(:,:,k)) * U;
            else
                dD(1:m,1:m,k) = squeeze(dU(:,:,k-m))' * diag(lambda) * U;
                dD(1:m,1:m,k) = dD(1:m,1:m,k) + dD(1:m,1:m,k)';
            end
        end

    case 'Lie-generators'
        % Initialize and get Covariance dimension from xi vector length
        m = round((-1 + sqrt(1+8*length(xi)))/2);
        dLambda = zeros(m,m,m);
        ddLambda = zeros(m,m,m,m);
        dD = zeros(m,m,length(xi));
        
        % The eigenvalues and the rotations are separated: 
        % Eigenvalues are stored in Lambda
        % Rotations are stored as skew-symmetric matrix which is
        % exponentiated

        % Separate desired eigenvalues lambda and rotation angles delta
        lambda = exp(xi(1:m));
        Lambda = diag(lambda);
        for iD = 1 : m
            dlambda = zeros(m);
            dlambda(iD,iD) = 1;
            dLambda(:,:,iD) = Lambda * dlambda;
            ddLambda(:,:,iD,iD) = Lambda * dlambda;
        end

        % Get Angle entries
        rotEntries = xi(m+1:end);
        rot = zeros(m);
        entryCount = 1;

        % Set up matrix with generators of rotations
        dRot = zeros(m,m,length(xi)-m);
        drot = zeros(m,m,length(xi)-m);
        for iRot = 1 : m
            for jRot = 1 : iRot - 1
                rot(jRot, iRot) = rotEntries(entryCount);
                drot(jRot, iRot, entryCount) =  1;
                drot(iRot, jRot, entryCount) = -1;
                entryCount = entryCount + 1;
            end
        end

        % Antisymmetrize and exponentiate
        rot = rot - rot';
        Rot = expm(rot);
        % Differentiating the exponential of a linear combination of
        % generators of a semi-simple Lie algenbra means encountering
        % Baker-Campbell-Hausdorff in it's full beauty, so no analytical
        % solution. Hence, the best we can hope for is finite differences
        % with a wisely chosen step size...
        fdStep = 1e-8;
        for iD = 1 : m*(m-1)/2
            dRot(:,:,iD) = (expm(rot + fdStep*drot(:,:,iD)) - expm(rot - fdStep*drot(:,:,iD))) / (2*fdStep);
        end

        D = Rot * Lambda * Rot';
        for iD = 1 : length(xi)
            if iD <= m
                dD(:,:,iD) = Rot * dLambda(:,:,iD) * Rot';
            else
                tmp = dRot(:,:,iD-m);
                dD(:,:,iD) = tmp * Lambda * Rot' + Rot * Lambda * tmp';
            end
        end

        invD = [];
        dinvD = [];
        HD = [];
        HinvD = [];

    case 'Cholesky'
        p = length(xi);
        m = floor(sqrt(2*p));
        offDiags = xi(m+1:end);

        colIndices = zeros(p - m, 1);
        colStart = 1;
        for iCol = 2 : m
            colEnd = iCol * (iCol-1) / 2;
            colIndices(colStart : colEnd) = iCol;
            colStart = colEnd + 1;
        end

        % build Cholesky decomposition
        L = diag(exp(xi(1:m)));
        for iCol = 2 : m
            offDiagsIterSpace = max(1, 1 + (iCol-1)*(iCol-2)/2) : iCol*(iCol-1)/2;
            tmpOffDiags = offDiags(offDiagsIterSpace);
            L(1:iCol-1, iCol) = tmpOffDiags(:);
        end

        % get derivatives - diagonal elements
        dL = zeros(m, m, p);
        for iDiag = 1 : m
            dL(iDiag, iDiag, iDiag) = exp(xi(iDiag));
        end
        % get derivatives - off-diagonals
        for iOff = m + 1 : p
            iCol = colIndices(iOff - m);
            iRow = iOff - m - (iCol-1) * (iCol-2) / 2;
            dL(iRow, iCol, iOff) = 1;
        end

        % set return values
        D = L' * L;
        invD = pinv(D);

        % set return variables for derivatives
        dD = zeros(m, m, p);
        dinvD = zeros(m, m, p);
        for iPar = 1 : p
            dD(:,:,iPar) = L' * dL(:,:,iPar) + transpose(dL(:,:,iPar)) * L;
            dinvD(:,:,iPar) = - invD * dD(:,:,iPar) * invD;
        end

        HD = [];
        HinvD = [];

    case 'spherical'
        p = length(xi);
        m = floor(sqrt(2*p));
        diags = xi(1:m);
        offDiags = xi(m+1:end);

        % preallocate upper Cholesky decomposition
        L = zeros(m);

        % Assemble Cholesky parametrization column by column...
        % ... and compute derivatives of Cholesky decomposition on the fly
        dL = zeros(m, m, p);

        for iCol = 1 : m
            % norms of cholesky vectors
            iColNorm = exp(diags(iCol));
            % identify off-diagonal elements
            offDiagsIterSpace = max(1, 1 + (iCol-1)*(iCol-2)/2) : iCol*(iCol-1)/2;
            tmpOffDiags = offDiags(offDiagsIterSpace);

            % Compute the standard Cholesky decomposition from spherical
            % coordinates
            for iRow = 1 : iCol - 1
                % use logical indexing to keep things understandable
                id_sin = false(1, p - m);
                id_cos = false(1, p - m);
                id_sin(offDiagsIterSpace(1 : iRow - 1)) = true;
                id_cos(offDiagsIterSpace(iRow)) = true;
                sinterms = sin(offDiags(id_sin));
                costerm = cos(offDiags(id_cos));
                L(iRow, iCol) = prod(sinterms) * prod(costerm);

                % assemble the derivatives
                for iOff = 1 : p - m
                    d_id_sin = id_sin;
                    d_id_cos = id_cos;
                    if ~id_sin(iOff) && ~id_cos(iOff)
                        signum = 0;
                    elseif id_sin(iOff)
                        d_id_sin(iOff) = false;
                        d_id_cos(iOff) = true;
                        signum = 1;
                    elseif id_cos(iOff)
                        d_id_sin(iOff) = true;
                        d_id_cos(iOff) = false;
                        signum = -1;
                    else
                        error('something went wrong...');
                    end
                    sinterms = sin(offDiags(d_id_sin));
                    costerm = cos(offDiags(d_id_cos));
                    dL(iRow, iCol, m + iOff) = signum * prod(sinterms) * prod(costerm);
                end
            end
            % Last entry is always an exception, as only sin terms appear
            % This also works fir the first column, which has no offdiags
            L(iCol, iCol) = prod(sin(tmpOffDiags));
            % Do the same for the derivative
            for iOff = 1 : p - m
                d_id_sin = false(1, p - m);
                d_id_sin(offDiagsIterSpace) = true;
                d_id_cos = false(1, p - m);
                if any(offDiagsIterSpace == iOff)
                    signum = 1;
                    d_id_sin(iOff) = false;
                    d_id_cos(iOff) = true;
                else
                    signum = 0;
                end
                sinterms = sin(offDiags(d_id_sin));
                costerm = cos(offDiags(d_id_cos));
                dL(iCol, iCol, m + iOff) = signum * prod(sinterms) * prod(costerm);
            end

            % rescale everything with norm
            L(:, iCol) = L(:, iCol) * iColNorm;
            dL(:, iCol, :) = dL(:, iCol, :) * iColNorm;
        end

        % Assemble derivatives of diagonal elements
        for iDiag = 1 : m
            dL(:, iDiag, iDiag) = L(:, iDiag);
        end

        % set return values
        D = L' * L;
        invD = pinv(D);

        % set return variables for derivatives
        dD = zeros(m, m, p);
        dinvD = zeros(m, m, p);
        for iPar = 1 : p
            dD(:,:,iPar) = L' * dL(:,:,iPar) + transpose(dL(:,:,iPar)) * L;
            dinvD(:,:,iPar) = - invD * dD(:,:,iPar) * invD;
        end
        HD = [];
        HinvD = [];

    case 'matrix-logarithm'
        % Matrix logarithm parameterization
        m = round((-1 + sqrt(1+8*length(xi)))/2);
        logD = zeros(m,m);
        dlogD = zeros(m,m,length(xi));
        d2logD = zeros(m,m,length(xi),length(xi));
        dD = zeros(m,m,length(xi));
        dinvD = zeros(m,m,length(xi));
        HD = zeros(m,m,length(xi),length(xi));
        HinvD = zeros(m,m,length(xi),length(xi));
        invD = zeros(m,m);

        k = 1;
        for i = 1:m
            for j = 1:i
                logD(i,j) = xi(k);
                logD(j,i) = xi(k);
                dlogD(i,j,k) = 1;
                dlogD(j,i,k) = 1;
                k = k+1;
            end
        end

        for i = 1:length(xi)
            for j = 1:i
                M = expm([logD(:,:)   , dlogD(:,:,i) , dlogD(:,:,j) , d2logD(:,:,i,j) ;...
                    zeros(m)    , logD(:,:)    , zeros(m)     , dlogD(:,:,j) ;...
                    zeros(m)    , zeros(m)     , logD(:,:)    , dlogD(:,:,i) ;...
                    zeros(m)    , zeros(m)     , zeros(m)     , logD(:,:)]);

            % Optional code for the inverse matrix:
            %     invM = expm([-logD(:,:)   , -dlogD(:,:,i) , -dlogD(:,:,j) , -d2logD(:,:,i,j) ;...
            %         zeros(m)    , -logD(:,:)    , zeros(m)     , -dlogD(:,:,j) ;...
            %         zeros(m)    , zeros(m)     , -logD(:,:)    , -dlogD(:,:,i) ;...
            %         zeros(m)    , zeros(m)     , zeros(m)     , -logD(:,:)]);
            %     HD(:,:,i,j) = M(1:m,3*m+1:4*m);
            %     HD(:,:,j,i) = HD(:,:,i,j)';
            %     HinvD(:,:,i,j) = invM(1:m,3*m+1:4*m);
            %     HinvD(:,:,j,i) = HinvD(:,:,i,j)';
            end
            dD(:,:,i) = M(1:m,m+1:2*m);
            % Optional code for the derivative of the inverse matrix:
            % dinvD(:,:,i) = invM(1:m,m+1:2*m);
        end
        D = M(1:m,1:m);
        % Optional code for the inverse matrix:
        % invD = invM(1:m,1:m);

end
