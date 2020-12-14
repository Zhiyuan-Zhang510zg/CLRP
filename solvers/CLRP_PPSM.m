function [ u,v,out ] = CLRP_PPSM(x0,K,OPTK,opts,alg_opts)
% PPSM for a customized low-rank prior model for structured
% cartoon-texture image decomposition model
% with different linear operator K.

% min_{u,v} tau ||u||_TV + mu ||v||_* + 0.5 ||K(u+v)-b||^2.
% where || dot ||_TV is the anisotropic TV norm and 
% || dot ||_* is the nuclear norm .

% Input:
% x0: observed image;
% K: linear operator;
% OPTK: linear operator K tzpe: 'I', 'S', and 'B'.
% opts: model parameters setting;
% alg_opts: PPSM parameters setting;

% Output:
% u: cartoon component;
% v: texture component;
% out: comparative indices.

% Author: Zhizuan Zhang
% Date: 1, December, 2020

tau  = opts.tau;   mu  = opts.mu;
beta1 = opts.beta1;     beta2 = opts.beta2;
MaxIt = opts.MaxIt;  I     = opts.I;      Tol   = opts.Tol;
[n1,n2,n3] = size(x0);

%% the algorithm parameters
alg_gamma = alg_opts.gamma;
alg_s = alg_opts.s;
alg_r = alg_opts.r;
if alg_s*alg_r <= 2
    error('Please check the algorithm parameters!')
end

%%%%%%%%%%%%%%%%% Periodic  boundarz condtion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1h = zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h = zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px  = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x
Py  = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 z
PTx = @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x
PTy = @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u  = zeros(n1,n2,n3);
v  = zeros(n1,n2,n3);
v_tilde = zeros(n1,n2,n3);
y1 = zeros(n1,n2,n3); y2 = zeros(n1,n2,n3);
z  = zeros(n1,n2,n3);
lbd11 = zeros(n1,n2,n3); lbd12 = zeros(n1,n2,n3);
lbd2  = zeros(n1,n2,n3);

switch lower(OPTK)
    case 'i'  %%%%%%%%%%%%%%%%%%%%%%%%% K = I %%%%%%%%%%%%%%%%
        MDz   = 1 + alg_r*beta2;  ATf = x0;
    case 's'  %%%%%%%%%%%%%%%%%%%%%%%%% K = S %%%%%%%%%%%%%%%%
        MDz   = K + alg_r*beta2;  ATf = K.*x0;
    case 'b'  %%%%%%%%%%%%%%%%%%%%%%%%% K = B %%%%%%%%%%%%%%%%
        siz = size(K); center = [fix(siz(1)/2+1),fix(siz(2)/2+1)];
        P   = zeros(n1,n2,n3); for i =1:n3; P(1:siz(1),1:siz(2),i) = K; end
        Bm   = fft2(circshift(P,1-center));
        HT  = @(x) real(ifft2(conj(Bm).*fft2(x))); %%% Transpose of blur operator.
        MDz = abs(Bm).^2 + alg_r*beta2;  ATf   = HT(x0);
end
MDu   = beta1*(abs(d1h).^2+abs(d2h).^2) + beta2;

Time = zeros(1,MaxIt);
SNR  = zeros(1,MaxIt);

for itr = 1:MaxIt
    
    %%% step 1: un
    tic;
    Temp= PTx(beta1*y1 + alg_s*lbd11)+PTy(beta1*y2+alg_s*lbd12)+beta2*(z-v)+alg_s*lbd2;
    un  = real(ifft2(fft2(Temp)./MDu));
    time_u = toc;
    
    %%%% update Lagrangian multipliers
    tic;
    dun1= Px(un);
    dun2= Py(un);
    lbd11_tilde= lbd11 + beta1/alg_s*(y1-dun1);
    lbd12_tilde= lbd12 + beta1/alg_s*(y2-dun2);
    lbd2_tilde = lbd2  + beta2/alg_s*(z-un-v);
    time_l = toc;
    
    %%% step 2: \tilde v
    tic;
    Temp = v + (2*lbd2_tilde-lbd2)/alg_r/beta2;
    for ii = 1:n3
        [U,D,VT] = svd(Temp(:,:,ii),'econ');
        D    = diag(D);
        ind  = find(D>mu/beta2/alg_r);
        D    = diag(D(ind) - mu/beta2/alg_r);
        v_tilde(:,:,ii)   = U(:,ind) * D * VT(:,ind)';
    end
    time_v = toc;
    
    %%% step 3: \tilde y
    tic;
    Temp1 = y1 - (2*lbd11_tilde-lbd11)/alg_r/beta1;
    Temp2 = y2 - (2*lbd12_tilde-lbd12)/alg_r/beta1;
    nsk = sqrt(Temp1.^2 + Temp2.^2); nsk(nsk==0)=1;
    nsk = max(1-(tau/beta1/alg_r)./nsk,0);
    y1_tilde = Temp1.*nsk;
    y2_tilde = Temp2.*nsk;
    time_y = toc;
    
    %%% step 4: \tilde z
    tic;
    Temp = ATf + alg_r*beta2*z - 2*lbd2_tilde + lbd2;
    if lower(OPTK) == 'b'
        z_tilde = real(ifft2(fft2(Temp)./MDz));
    else
        z_tilde = Temp./MDz;
    end
    time_z = toc;
    
    %%% step 5: relaxation
    tic;
    vn = v - alg_gamma*(v-v_tilde);
    yn1 = y1 - alg_gamma*(y1 - y1_tilde);
    yn2 = y2 - alg_gamma*(y2 - y2_tilde);
    zn  = z - alg_gamma*(z - z_tilde);
    lbdn11 = lbd11 - alg_gamma*(lbd11 -lbd11_tilde);
    lbdn12 = lbd12 - alg_gamma*(lbd12 -lbd12_tilde);
    lbdn2  = lbd2 - alg_gamma*(lbd2 -lbd2_tilde);
    time_update = toc;
    if alg_gamma==1
        time_update = 0;
    end
    
    %%% outputs
    Time(itr)= time_u + max([time_v,time_y,time_z]) + time_l + time_update;
    SNR(itr) = 20*log10(norm(I(:))/norm(un(:)+vn(:)-I(:)));
    %%% stopping rule
    Stopic = compute_tol(u,v,un,vn);
    
    if opts.verbose
        fprintf('the %d th Iter. the Tol. = %.4f \n ',itr,Stopic);
    end
    
    if Stopic<Tol && itr >1;
        u = un;
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