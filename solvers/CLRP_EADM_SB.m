function [ u,v,out ] = CLRP_EADM_SB(x0,K,S,opts)
% EADM for a customized low-rank prior model for structured
% cartoon-texture image decomposition model
% with different linear operator Phi.

% min_{u,v} tau ||u||_TV + mu ||v||_* + 0.5 ||Phi(u+v)-b||^2.
% where || dot ||_TV is the anisotropic TV norm and 
% || dot ||_* is the nuclear norm .

% Input:
% x0: observed image;
% K: blurring kernel;
% S: missing mask;
% opts: model parameters setting;

% Output:
% u: cartoon component;
% v: texture component;
% out: comparative indices.

% Author: Zhiyuan Zhang
% Date: 1, December, 2020

tau  = opts.tau;   mu  = opts.mu;
beta1 = opts.beta1;  beta2 = opts.beta2;
MaxIt = opts.MaxIt;  I     = opts.I;      Tol   = opts.Tol;
[n1,n2,n3] = size(x0);

%%%%%%%%%%%%%%%%% Periodic  boundary condtion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1h = zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h = zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px  = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x
Py  = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 y
PTx = @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x
PTy = @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u  = zeros(n1,n2,n3);
v  = zeros(n1,n2,n3);
vn = zeros(n1,n2,n3);
y1 = zeros(n1,n2,n3); y2 = zeros(n1,n2,n3);
z  = zeros(n1,n2,n3);
lbd11 = zeros(n1,n2,n3); lbd12 = zeros(n1,n2,n3);
lbd2  = zeros(n1,n2,n3);

siz = size(K); center = [fix(siz(1)/2+1),fix(siz(2)/2+1)];
P   = zeros(n1,n2,n3); for i =1:n3; P(1:siz(1),1:siz(2),i) = K; end
Bm   = fft2(circshift(P,1-center));
H   = @(x) real(ifft2(Bm.*fft2(x)));       %%% Blur operator.  B x
HT  = @(x) real(ifft2(conj(Bm).*fft2(x))); %%% Transpose of blur operator.
ATf   = HT(S.*x0);
MDz = @(x) LSP(H,HT,S,x,n1,n2,n3,beta2);
MDu   = beta1*(abs(d1h).^2+abs(d2h).^2) + beta2;

Time = zeros(1,MaxIt);
SNR  = zeros(1,MaxIt);

for itr = 1:MaxIt
    
    %%% step 1: \tilde u
    tic;
    Temp= PTx(beta1*y1+lbd11)+PTy(beta1*y2+lbd12)+beta2*(z-v)+lbd2;
    un  = real(ifft2(fft2(Temp)./MDu));
    time_u = toc;
    
    %%% step 2: \tilde v
    tic;
    Temp = z-un+lbd2/beta2;
    for ii = 1:n3
        [U,D,VT] = svd(Temp(:,:,ii),'econ');
        D    = diag(D);
        ind  = find(D>mu/beta2);
        D    = diag(D(ind) - mu/beta2);
        vn(:,:,ii)   = U(:,ind) * D * VT(:,ind)';
    end
    time_v = toc;
    
    %%% step 3: \tilde y
    tic;
    dxun= Px(un);
    dyun= Py(un);
    Temp1 = dxun - lbd11/beta1;
    Temp2 = dyun - lbd12/beta1;
    nsk = sqrt(Temp1.^2 + Temp2.^2); nsk(nsk==0)=1;
    nsk = max(1-(tau/beta1)./nsk,0);
    yn1 = Temp1.*nsk;
    yn2 = Temp2.*nsk;
    time_y =toc;
    
    %%% step 4: \tilde z
    tic;
    Temp= ATf + beta2*(un+vn) - lbd2;
    [zn,~]  = pcg(MDz,Temp(:),5e-4,20);
    zn  = reshape(zn,n1,n2,n3);
    time_z = toc;
    
    gamma = 1;
    %%%% update Lagrangian multipliers
    tic;
    lbdn11= lbd11 + gamma*beta1*(yn1-dxun);
    lbdn12= lbd12 + gamma*beta1*(yn2-dyun);
    lbdn2 = lbd2  + gamma*beta2*(zn-un-vn);
    time_l = toc;
    
    %%% outputs
    Time(itr)= time_u+ time_v + max([time_y , time_z]) + time_l;
    SNR(itr) = 20*log10(norm(I(:))/norm(un(:)+vn(:)-I(:)));
    %%% stopping rule
    Stopic = compute_tol(u,v,un,vn);
    
    if opts.verbose
        fprintf('the %d th Iter. the Tol. = %.4f \n',itr,Stopic);
    end
    
    if Stopic<Tol && itr>1;
        u = un ;
        v = vn;
        fprintf('It=%d,cpu=%4.2f,\n',itr,sum(Time));
        Time = Time(1:itr);
        SNR  = SNR(1:itr);
        break;
    end
    
    u = un ; v = vn;
    y1 = yn1; y2 = yn2; z = zn;
    lbd11 = lbdn11; lbd12 = lbdn12; lbd2 = lbdn2;
end
out.It   =itr;
out.Stopic = Stopic;
out.Time = Time;
out.SNR  = SNR;
end
