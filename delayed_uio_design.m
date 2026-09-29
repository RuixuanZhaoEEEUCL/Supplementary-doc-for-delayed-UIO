function g = delayed_uio_design(A,B,Bbar,C,Q,r,opts)
%DELAYED_UIO_DESIGN Offline geometry and gains for a delayed functional UIO.
%   g = delayed_uio_design(A,B,Bbar,C,Q,r)
%   g = delayed_uio_design(A,B,Bbar,C,Q,r,opts)
%
%   Plant:  x(k+1) = A*x(k)+B*u(k)+Bbar*d(k),  y(k) = C*x(k).
%   Target: Q*x(k-r), available causally at k. No unknown-input feedthrough.
%   B=[] means no known inputs; Bbar=[] means no unknown inputs.
%   Q=eye(n) requests full state. Q=[] requests coordinates of X/I_r.
%   Q=zeros(1,n) requests the identically zero scalar target.
%
%   COLUMN BASES (orthonormal, n-by-d; zero space uses n-by-0):
%     Bbasis,Vr,Vstar,W,S,Vfix,Xbad,Xgood,Wg,I
%     Vfix spans S intersect W^perp, representing S/W (not necessarily V*).
%     Xbad and Xgood are state-space lifts of fixed spectral subspaces.
%   QUOTIENT MAPS (full row rank; rows are orthonormal):
%     Cbar, PW, PS, Pg, PI represent P_Vr,P_W,P_S,P_Wg,P_I.
%     ker(Pg)=Im(Wg). Pg is NOT the n-by-n projector onto Wg.
%   HORIZON MATRICES:
%     O,Tunknown,Tknown,Omega,J,Jleft,Hwindow
%     Ywin=[y(k-r);...;y(k)], Uwin=[u(k-r);...;u(k-1)]
%     y_xi=Hwindow*(Ywin-Tknown*Uwin)=Cbar*x(k-r).
%   GAINS AND QUOTIENT DYNAMICS:
%     Lbase          friend of W
%     Af             fixed map on S/W, represented by Vfix
%     PY,Aassign,Cassign,Lassign  Algorithm-3 assignable quotient X/S
%     L              stabilizing friend of Wg (also a friend of W)
%     Ae             (A+L*Cbar)|_{X/Wg}
%     E,F            Q=E*Pg+F*Cbar, provided the target is feasible
%     observer       Az=Ae, Bu=Pg*B, By=-Pg*L, E=E, F=F
%   DECISIONS:
%     feasible       numerical test I_r subset ker(Q), NOT design success
%     canRun         feasible AND verified stable observer realization
%     dim            named dimensions; rbar = horizon saturation index
%     spectra        fixed, bad, good, assignable, observer eigenvalues
%     checks         normalized residuals of all defining identities
%     witness        a unit direction in I_r with Q*witness nonzero if infeasible
%
%   OPTIONS (all optional, exact field names below):
%     rankTol          1e-10: cutoff rankTol*max(1,largest singular value)
%     checkTol         1e-7: algebraic/subspace verification tolerance
%     targetTol        1e-8: inclusion tolerance, normalized by norm(Q,'fro')
%     unitCircleTol    1e-8: abs(z)>=1-unitCircleTol is retained as bad
%     designObserver   true: false returns geometry/recovery without gain design
%     gainMethod       'auto': 'place', 'riccati', or 'auto'
%     assignablePoles  []: n-dim(S) discrete-time poles; [] -> .2 to .7 for place
%     riccatiTol       1e-12
%     maxRiccatiIter   30000
%     verbose          false: print dimension/feasibility summary
%     warn             true: warn on borderline poles and large window gains
%
%   'place' uses Control System Toolbox and the paper's plus-sign convention:
%     Lassign=-place(Aassign',Cassign',poles)'.
%   'riccati' uses a built-in deterministic dual Riccati iteration, without
%   any extra toolbox or statistical assumption on the unknown input. It is
%   an OPTIONAL stabilization method, not the paper's pole-assignment step.
%   'auto' chooses place if present, otherwise riccati. Specified poles are
%   never silently ignored: they require place.

    if nargin < 6
        error('delayed_uio:Arguments','Supply A,B,Bbar,C,Q,r.');
    end
    if nargin < 7 || isempty(opts), opts=struct(); end
    opts=local_options(opts);
    local_matrix(A,'A');
    n=size(A,1);
    if n<1 || size(A,2)~=n
        error('delayed_uio:Dimensions','A must be a nonempty square matrix.');
    end
    A=full(double(A));
    if isempty(B), B=zeros(n,0); end
    if isempty(Bbar), Bbar=zeros(n,0); end
    if isempty(C), C=zeros(0,n); end
    local_matrix(B,'B'); local_matrix(Bbar,'Bbar'); local_matrix(C,'C');
    B=full(double(B)); Bbar=full(double(Bbar)); C=full(double(C));
    if size(B,1)~=n || size(Bbar,1)~=n || size(C,2)~=n
        error('delayed_uio:Dimensions','B,Bbar must have n rows; C must have n columns.');
    end
    if ~isnumeric(r) || ~isscalar(r) || ~isreal(r) || ~isfinite(r) || r<0 || r~=floor(r)
        error('delayed_uio:Delay','r must be a nonnegative integer.');
    end
    r=double(r); tol=opts.rankTol; p=size(C,1); m=size(Bbar,2);
    ell=size(B,2); QisMaximal=isempty(Q);
    if ~QisMaximal
        local_matrix(Q,'Q'); Q=full(double(Q));
        if size(Q,2)~=n, error('delayed_uio:Dimensions','Q must have n columns.'); end
    end
    Bbasis=local_orth(Bbar,tol);

    % 0. Horizon kernels: V_{j+1}=ker C intersect A^{-1}(V_j+Im Bbar).
    K=local_null(C,tol); V=K; Vseq={V}; rbar=[];
    for j=0:n
        VB=local_orth([V,Bbasis],tol);
        PV=local_project(VB,tol);
        Vnext=local_orth(K*local_null(PV*A*K,tol),tol);
        if local_same(Vnext,V,opts.checkTol)
            rbar=j; break;
        end
        if size(Vnext,2)>size(V,2)
            error('delayed_uio:Rank','Controlled iteration lost monotonicity; check scaling.');
        end
        V=Vnext; Vseq{end+1}=V; %#ok<AGROW>
    end
    if isempty(rbar), error('delayed_uio:Recursion','V recursion did not stabilize.'); end
    Vstar=V; Vr=Vseq{min(r,rbar)+1}; Cbar=local_project(Vr,tol);
    pr=size(Cbar,1);

    % Horizon realization. Tunknown is built from the ORIGINAL Bbar columns.
    O=zeros((r+1)*p,n); Tu=zeros((r+1)*p,r*ell);
    Td=zeros((r+1)*p,r*m); Apow=cell(r+1,1); Apow{1}=eye(n);
    for j=0:r
        if j>0, Apow{j+1}=Apow{j}*A; end
        rows=j*p+(1:p); O(rows,:)=C*Apow{j+1};
        for i=0:j-1
            Tu(rows,i*ell+(1:ell))=C*Apow{j-i}*B;
            Td(rows,i*m+(1:m))=C*Apow{j-i}*Bbar;
        end
    end
    if any(~isfinite([O(:);Tu(:);Td(:)]))
        error('delayed_uio:HorizonOverflow','Nonfinite horizon matrices; reduce r or rescale.');
    end
    Rbasis=local_orth(Td,tol); Omega=local_project(Rbasis,tol);
    J=Omega*O*Cbar'; Jleft=local_pinv(J,tol); Hwindow=Jleft*Omega;
    if local_rank(J,tol)~=pr
        error('delayed_uio:HorizonRank', ...
            'J is not numerically injective. Horizon and recursive ranks disagree; rescale.');
    end

    % 1. Algorithm 1: conditioned subspace; use CURRENT iterate, not W^(0).
    W=zeros(n,0); Wseq={W}; stopped=false;
    for j=0:n
        Wcap=local_orth(W*local_null(Cbar*W,tol),tol);
        Wnext=local_orth([Bbasis,A*Wcap],tol);
        if local_same(Wnext,W,opts.checkTol), stopped=true; break; end
        if size(Wnext,2)<size(W,2)
            error('delayed_uio:Rank','W iteration lost monotonicity; check scaling.');
        end
        W=Wnext; Wseq{end+1}=W; %#ok<AGROW>
    end
    if ~stopped, error('delayed_uio:Recursion','W recursion did not stabilize.'); end
    PW=local_project(W,tol); % recompute from FINAL basis

    % Algorithm 1: S_{j+1}=W+(ker Cbar intersect A^{-1}S_j).
    S=eye(n); Sseq={S}; stopped=false;
    for j=0:n
        Pcurrent=local_project(S,tol);
        Tj=local_null([Pcurrent*A;Cbar],tol);
        Snext=local_orth([W,Tj],tol);
        if local_same(Snext,S,opts.checkTol), stopped=true; break; end
        if size(Snext,2)>size(S,2)
            error('delayed_uio:Rank','S iteration lost monotonicity; check scaling.');
        end
        S=Snext; Sseq{end+1}=S; %#ok<AGROW>
    end
    if ~stopped, error('delayed_uio:Recursion','S recursion did not stabilize.'); end
    PS=local_project(S,tol);

    % 2. Algorithm 2: orthonormal representatives of S/W and a W-friend.
    Vfix=local_orth(S*local_null(W'*S,tol),tol);
    if size(Vfix,2)~=size(S,2)-size(W,2)
        error('delayed_uio:QuotientDimension','S/W complement has an inconsistent dimension.');
    end
    Lbase=-PW'*(PW*A*W)*local_pinv(Cbar*W,tol);
    Abase=A+Lbase*Cbar;
    Af=Vfix'*Abase*Vfix;
    % REAL ordered Schur subspaces replace minimal-polynomial factorization.
    % Compute bad and good subspaces SEPARATELY: trailing Schur columns need
    % not be the complementary good invariant subspace in a nonnormal map.
    [Zbad,Zgood,Tbad,Tgood,lambdaFixed,nearUnit]=local_spectral(Af,opts);
    Xbad=local_orth(Vfix*Zbad,tol); Xgood=local_orth(Vfix*Zgood,tol);
    Wg=local_orth([W,Xbad],tol); Pg=local_project(Wg,tol);
    if size(Wg,2)~=size(W,2)+size(Zbad,2)
        error('delayed_uio:LiftRank','Bad-subspace lift lost rank; inspect scaling.');
    end
    Ibasis=local_orth(Wg*local_null(Cbar*Wg,tol),tol);
    PI=local_project(Ibasis,tol);
    if QisMaximal, Q=PI; end

    % 3. Algorithm 3: only the quotient X/S has arbitrarily assignable poles.
    % P_Y annihilates Cbar*S. Lift Lassign through PS', NOT Pg'.
    SY=local_orth(Cbar*S,tol); PY=local_project(SY,tol);
    Aassign=PS*Abase*PS'; Cassign=PY*Cbar*PS';
    na=size(Aassign,1); ny0=size(Cassign,1);
    Lassign=zeros(na,ny0); L=[]; Ae=[];
    design=struct('requestedMethod',opts.gainMethod,'methodUsed','not-designed', ...
        'nAssignable',na,'nEffectiveOutputs',ny0,'desiredPoles',[], ...
        'iterations',0,'riccatiResidual',NaN,'stable',false,'rho',NaN);
    if opts.designObserver
        if na==0
            if ~isempty(opts.assignablePoles)
                error('delayed_uio:PoleCount','No assignable modes: assignablePoles must be empty.');
            end
            design.methodUsed='no-assignable-modes';
        else
            [Lassign,design0]=local_gain(Aassign,Cassign,opts);
            names=fieldnames(design0);
            for k=1:numel(names), design.(names{k})=design0.(names{k}); end
        end
        L=Lbase+PS'*Lassign*PY;
        Ae=Pg*(A+L*Cbar)*Pg';
        design.rho=local_rho(Ae); design.stable=design.rho<1;
        if ~design.stable
            error('delayed_uio:Unstable', ...
                'Computed quotient is not Schur (rho=%.8g). Inspect scaling/ranks/gains.',design.rho);
        end
    end

    % 4. Functional reconstruction condition; never silently change Q.
    qnorm=norm(Q,'fro');
    if qnorm==0, targetResidual=0;
    else, targetResidual=norm(Q*Ibasis,'fro')/qnorm; end
    feasible=targetResidual<=opts.targetTol;
    Psi=[Pg;Cbar]; E=[]; F=[]; witness=zeros(n,0);
    if feasible
        EF=Q*local_pinv(Psi,tol);
        E=EF(:,1:size(Pg,1)); F=EF(:,size(Pg,1)+1:end);
    else
        [~,~,Rt]=svd(Q*Ibasis,'econ');
        witness=Ibasis*Rt(:,1); % algebraic violating direction, not yet a trajectory
    end

    % 5. Verification, with absolute residuals also retained for transparency.
    checks=struct(); absChecks=struct();
    [checks.BinW,absChecks.BinW]=local_res(PW*Bbar,norm(Bbar,'fro'));
    [checks.Wconditioned,absChecks.Wconditioned]= ...
        local_res(PW*A*(W*local_null(Cbar*W,tol)),norm(A,'fro'));
    [checks.baseFriend,absChecks.baseFriend]=local_res(PW*Abase*W,norm(Abase,'fro'));
    [checks.Sinvariant,absChecks.Sinvariant]=local_res(PS*Abase*S,norm(Abase,'fro'));
    [checks.fixedMap,absChecks.fixedMap]= ...
        local_res(PW*Abase*Vfix-PW*Vfix*Af,norm(Abase,'fro'));
    WplusV=local_orth([W,Vstar],tol);
    [checks.generalDecomposition,absChecks.generalDecomposition]= ...
        local_res(S*S'-WplusV*WplusV',1);
    [checks.badInvariant,absChecks.badInvariant]= ...
        local_res(Af*Zbad-Zbad*Tbad,norm(Af,'fro'));
    [checks.goodInvariant,absChecks.goodInvariant]= ...
        local_res(Af*Zgood-Zgood*Tgood,norm(Af,'fro'));
    [checks.unknownDecoupling,absChecks.unknownDecoupling]=local_res(Pg*Bbar,norm(Bbar,'fro'));
    [checks.outputMixing,absChecks.outputMixing]=local_res(PY*Cbar*S,norm(Cbar,'fro'));
    [checks.horizonFactorization,absChecks.horizonFactorization]= ...
        local_res(Omega*O-J*Cbar,norm(Omega*O,'fro'));
    [checks.windowState,absChecks.windowState]=local_res(Hwindow*O-Cbar,norm(Cbar,'fro'));
    [checks.windowUnknown,absChecks.windowUnknown]=local_res(Hwindow*Td,norm(Td,'fro'));
    [checks.leftInverse,absChecks.leftInverse]=local_res(Jleft*J-eye(pr),sqrt(pr));
    [checks.projectionRows,absChecks.projectionRows]= ...
        local_res(Pg*Pg'-eye(size(Pg,1)),sqrt(size(Pg,1)));
    if opts.designObserver
        AL=A+L*Cbar;
        [checks.friend,absChecks.friend]=local_res(Pg*AL*Wg,norm(AL,'fro'));
        [checks.inducedMap,absChecks.inducedMap]=local_res(Pg*AL-Ae*Pg,norm(AL,'fro'));
        [checks.fixedUnchanged,absChecks.fixedUnchanged]= ...
            local_res(Vfix'*AL*Vfix-Af,norm(Af,'fro'));
        [checks.assignableLift,absChecks.assignableLift]= ...
            local_res(PS*AL*PS'-(Aassign+Lassign*Cassign),norm(AL,'fro'));
    end
    if feasible
        [checks.recovery,absChecks.recovery]=local_res(E*Pg+F*Cbar-Q,qnorm);
    end
    vals=struct2cell(checks); vals=cell2mat(vals);
    maxResidual=max(vals);
    if maxResidual>opts.checkTol || any(~isfinite(vals))
        names=fieldnames(checks); [~,idx]=max(vals);
        error('delayed_uio:Residual', ...
            'Identity %s failed (normalized residual %.3g); rescale or revisit tolerances.', ...
            names{idx},vals(idx));
    end

    % All matrix representatives share one explicitly recorded coordinate choice.
    g=struct('A',A,'B',B,'Bbar',Bbar,'C',C,'Q',Q,'r',r,'options',opts);
    g.Bbasis=Bbasis; g.Vr=Vr; g.Vstar=Vstar; g.Vseq=Vseq; g.rbar=rbar;
    g.W=W; g.S=S; g.Vfix=Vfix; g.Wseq=Wseq; g.Sseq=Sseq;
    g.Wg=Wg; g.I=Ibasis; g.Xbad=Xbad; g.Xgood=Xgood;
    g.Cbar=Cbar; g.PW=PW; g.PS=PS; g.Pg=Pg; g.PI=PI;
    g.O=O; g.Tunknown=Td; g.Tknown=Tu; g.Rbasis=Rbasis;
    g.Omega=Omega; g.J=J; g.Jleft=Jleft; g.Hwindow=Hwindow;
    g.Lbase=Lbase; g.Abase=Abase; g.Af=Af;
    g.Zbad=Zbad; g.Zgood=Zgood; g.Tbad=Tbad; g.Tgood=Tgood;
    g.PY=PY; g.Aassign=Aassign; g.Cassign=Cassign; g.Lassign=Lassign;
    g.L=L; g.Ae=Ae; g.E=E; g.F=F; g.Psi=Psi;
    g.feasible=feasible; g.targetResidual=targetResidual;
    g.witness=witness; g.maximalTargetRequested=QisMaximal;
    g.canRun=feasible && opts.designObserver && design.stable;
    g.design=design; g.checks=checks; g.absoluteChecks=absChecks;
    g.spectra=struct('fixed',lambdaFixed,'bad',eig(Tbad),'good',eig(Tgood), ...
        'nearUnitCircle',lambdaFixed(nearUnit),'assignable',[],'observer',[]);
    if opts.designObserver
        g.spectra.assignable=eig(Aassign+Lassign*Cassign);
        g.spectra.observer=eig(Ae);
    end
    g.dim=struct('n',n,'knownInputs',ell,'unknownInputColumns',m, ...
        'unknownInputRank',size(Bbasis,2),'measuredOutputs',p,'horizonOutput',pr, ...
        'Vr',size(Vr,2),'Vstar',size(Vstar,2),'W',size(W,2),'S',size(S,2), ...
        'fixed',size(Af,1),'badFixed',size(Zbad,2),'goodFixed',size(Zgood,2), ...
        'Wg',size(Wg,2),'I',size(Ibasis,2),'maxReconstructible',size(PI,1), ...
        'observerState',size(Pg,1),'assignable',na,'targetRows',size(Q,1));
    g.diagnostics=struct('maxResidual',maxResidual,'windowNorm2',local_norm2(Hwindow), ...
        'Jcondition',local_condition(J,tol),'recoveryNorm2',local_norm2([E,F]), ...
        'boundaryClassification',any(nearUnit),'numericalOnly',true);
    g.observer=[];
    if g.canRun
        g.observer=struct('Az',Ae,'Bu',Pg*B,'By',-Pg*L,'E',E,'F',F, ...
            'Hy',Hwindow,'Hu',-Hwindow*Tu);
        g.observer.BYwindow=g.observer.By*Hwindow;
        g.observer.BUwindow=-g.observer.By*Hwindow*Tu;
        g.observer.DYwindow=F*Hwindow;
        g.observer.DUwindow=-F*Hwindow*Tu;
    end
    if opts.warn && any(nearUnit)
        warning('delayed_uio:UnitCircle', ...
            'Fixed poles near |z|=1 were conservatively included in Wg. Inspect g.spectra.');
    end
    if opts.warn && g.diagnostics.Jcondition>1e10
        warning('delayed_uio:WindowCondition','The horizon map is ill conditioned (cond %.3g).', ...
            g.diagnostics.Jcondition);
    end
    if opts.verbose
        fprintf('r=%d: dim(Vr,W,S,Wg,I)=(%d,%d,%d,%d,%d), maxInfo=%d\n', ...
            r,g.dim.Vr,g.dim.W,g.dim.S,g.dim.Wg,g.dim.I,g.dim.maxReconstructible);
        fprintf('target feasible=%d, canRun=%d, method=%s, maxResidual=%.3g\n', ...
            g.feasible,g.canRun,g.design.methodUsed,maxResidual);
    end
end

function o=local_options(o)
    if ~isstruct(o) || ~isscalar(o), error('delayed_uio:Options','opts must be a scalar structure.'); end
    defaults=struct('rankTol',1e-10,'checkTol',1e-7,'targetTol',1e-8, ...
        'unitCircleTol',1e-8,'designObserver',true,'gainMethod','auto', ...
        'assignablePoles',[],'riccatiTol',1e-12,'maxRiccatiIter',30000, ...
        'verbose',false,'warn',true);
    f=fieldnames(o); known=fieldnames(defaults);
    for j=1:numel(f)
        if ~any(strcmp(f{j},known)), error('delayed_uio:Options','Unknown option: %s',f{j}); end
    end
    for j=1:numel(known)
        if ~isfield(o,known{j}), o.(known{j})=defaults.(known{j}); end
    end
    pos={'rankTol','checkTol','targetTol','riccatiTol'};
    for j=1:numel(pos)
        t=o.(pos{j});
        if ~isnumeric(t)||~isreal(t)||~isscalar(t)||~isfinite(t)||t<=0||t>=1
            error('delayed_uio:Options','%s must lie strictly between zero and one.',pos{j});
        end
    end
    t=o.unitCircleTol;
    if ~isnumeric(t)||~isreal(t)||~isscalar(t)||~isfinite(t)||t<0||t>=0.1
        error('delayed_uio:Options','unitCircleTol must be in [0,0.1).');
    end
    bool={'designObserver','verbose','warn'};
    for j=1:numel(bool)
        t=o.(bool{j});
        if ~isscalar(t)||~(islogical(t)||(isnumeric(t)&&isreal(t)&&any(t==[0,1])))
            error('delayed_uio:Options','%s must be true or false.',bool{j});
        end
        o.(bool{j})=logical(t);
    end
    if ~(ischar(o.gainMethod)||(isstring(o.gainMethod)&&isscalar(o.gainMethod)))
        error('delayed_uio:Options','gainMethod must be auto, place or riccati.');
    end
    o.gainMethod=lower(char(o.gainMethod));
    if ~any(strcmp(o.gainMethod,{'auto','place','riccati'}))
        error('delayed_uio:Options','gainMethod must be auto, place or riccati.');
    end
    p=o.assignablePoles;
    if ~isempty(p) && (~isnumeric(p)||~isvector(p)||any(~isfinite(p(:)))||any(abs(p(:))>=1))
        error('delayed_uio:Poles','assignablePoles must be finite and strictly inside the unit circle.');
    end
    o.assignablePoles=p(:);
    t=o.maxRiccatiIter;
    if ~isnumeric(t)||~isreal(t)||~isscalar(t)||~isfinite(t)||t<1||t~=floor(t)
        error('delayed_uio:Options','maxRiccatiIter must be a positive integer.');
    end
end

function local_matrix(M,name)
    if ~isnumeric(M)||ndims(M)~=2||~isreal(M)||any(~isfinite(M(:)))
        error('delayed_uio:Matrix','%s must be a finite real numeric matrix.',name);
    end
end

function U=local_orth(M,tol)
    if isempty(M), U=zeros(size(M,1),0); return; end
    [U,S,~]=svd(M,'econ'); k=min(size(M)); d=diag(S(1:k,1:k));
    q=sum(d>tol*max(1,d(1))); U=U(:,1:q);
end

function Z=local_null(M,tol)
    n=size(M,2);
    if n==0, Z=zeros(0,0); return; end
    if size(M,1)==0, Z=eye(n); return; end
    [~,S,V]=svd(M); k=min(size(M)); d=diag(S(1:k,1:k));
    if isempty(d), q=0; else, q=sum(d>tol*max(1,d(1))); end
    Z=V(:,q+1:n);
end

function P=local_project(U,tol)
    P=local_null(U',tol)';
end

function R=local_pinv(M,tol)
    if isempty(M), R=zeros(size(M,2),size(M,1)); return; end
    [U,S,V]=svd(M,'econ'); k=min(size(M)); d=diag(S(1:k,1:k)); q=sum(d>tol*max(1,d(1)));
    R=V(:,1:q)*diag(1./d(1:q))*U(:,1:q)';
end

function q=local_rank(M,tol)
    if isempty(M), q=0; return; end
    d=svd(M); q=sum(d>tol*max(1,d(1)));
end

function ok=local_same(U,V,tol)
    ok=size(U,2)==size(V,2) && norm(U*U'-V*V','fro')<=tol;
end

function [b,g,Tb,Tg,e,near]=local_spectral(Af,opts)
    d=size(Af,1);
    if d==0
        b=zeros(0,0); g=b; Tb=b; Tg=b; e=zeros(0,1); near=false(0,1); return;
    end
    [U,T]=schur(Af,'real'); e=ordeig(T);
    near=abs(abs(e)-1)<=opts.unitCircleTol;
    bad=abs(e)>=1-opts.unitCircleTol;
    % Keep each real 2-by-2 Schur block together even at a roundoff boundary.
    j=1;
    while j<d
        if T(j+1,j)~=0
            bad(j:j+1)=any(bad(j:j+1)); near(j:j+1)=any(near(j:j+1)); j=j+2;
        else, j=j+1; end
    end
    nb=sum(bad); ng=d-nb;
    if nb==0, b=zeros(d,0); Tb=zeros(0,0);
    else, [Ub,Sb]=ordschur(U,T,bad); b=Ub(:,1:nb); Tb=Sb(1:nb,1:nb); end
    if ng==0, g=zeros(d,0); Tg=zeros(0,0);
    else, [Ug,Sg]=ordschur(U,T,~bad); g=Ug(:,1:ng); Tg=Sg(1:ng,1:ng); end
end

function [L0,info]=local_gain(A0,C0,opts)
    n=size(A0,1); p=size(C0,1); method=opts.gainMethod;
    hasPlace=exist('place','file')==2 || exist('place','file')==3;
    if strcmp(method,'auto')
        if hasPlace, method='place';
        elseif isempty(opts.assignablePoles), method='riccati';
        else, error('delayed_uio:Toolbox','Specified poles require place (Control System Toolbox).'); end
    end
    if p==0, error('delayed_uio:Assignable','Nonempty assignable quotient has no effective output.'); end
    info=struct('methodUsed',method,'desiredPoles',[],'iterations',0,'riccatiResidual',NaN);
    if strcmp(method,'place')
        if ~hasPlace, error('delayed_uio:Toolbox','place is unavailable. Use gainMethod=riccati.'); end
        poles=opts.assignablePoles;
        if isempty(poles)
            if n==1, poles=0.4; else, poles=linspace(0.2,0.7,n).'; end
        end
        if numel(poles)~=n
            error('delayed_uio:PoleCount', ...
                'Need %d poles for X/S, not n-dim(Wg) poles.',n);
        end
        unused=poles(:);
        while ~isempty(unused)
            z=unused(1); unused(1)=[];
            if abs(imag(z))>1e-12
                [dist,j]=min(abs(unused-conj(z)));
                if isempty(j)||dist>1e-10
                    error('delayed_uio:Poles','Complex poles must occur in conjugate pairs.');
                end
                unused(j)=[];
            end
        end
        L0=-place(A0',C0',poles)';
        if norm(imag(L0),'fro')>1e-9*max(1,norm(L0,'fro'))
            error('delayed_uio:ComplexGain','Unexpected complex injection for a real system.');
        end
        L0=real(L0); info.desiredPoles=poles;
    else
        if ~isempty(opts.assignablePoles)
            error('delayed_uio:Poles','riccati does not assign specified poles; choose place.');
        end
        % Dual deterministic LQR with identity weights, in predictor form.
        % L0=-A0*P*C0'/(I+C0*P*C0'); this is the plus-sign injection.
        P=eye(n); R=eye(p); W=eye(n); stopped=false;
        for it=1:opts.maxRiccatiIter
            H=R+C0*P*C0'; K=(A0*P*C0')/H;
            Pnext=A0*P*A0'+W-K*(C0*P*A0'); Pnext=(Pnext+Pnext')/2;
            if any(~isfinite(Pnext(:)))
                error('delayed_uio:Riccati','Riccati iteration overflow; check scaling/observability.');
            end
            rr=norm(Pnext-P,'fro')/max(1,norm(Pnext,'fro'));
            P=Pnext;
            if rr<=opts.riccatiTol, stopped=true; break; end
        end
        if ~stopped, error('delayed_uio:Riccati','Riccati iteration did not converge.'); end
        L0=-(A0*P*C0')/(R+C0*P*C0');
        info.iterations=it; info.riccatiResidual=rr;
    end
    if local_rho(A0+L0*C0)>=1
        error('delayed_uio:Assignable','Assignable closed-loop quotient is not Schur.');
    end
end

function rho=local_rho(A)
    if isempty(A), rho=0; else, rho=max(abs(eig(A))); end
end

function [r,a]=local_res(M,scale)
    a=norm(M,'fro'); r=a/max(1,scale);
end

function v=local_norm2(M)
    if isempty(M), v=0; else, v=norm(M,2); end
end

function k=local_condition(M,tol)
    if size(M,2)==0, k=1; return; end
    s=svd(M); q=sum(s>tol*max(1,s(1)));
    if q<size(M,2), k=Inf; else, k=s(1)/s(q); end
end
