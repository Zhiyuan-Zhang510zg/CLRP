function [ tol ] = compute_tol( u,v,un,vn )
% compute the tolerance for image decomposition

% input:
% u,v : the last iteration results
% un,vn: the lateset iteration results

du = u - un;
dv = v - vn;
tol1 = norm(du(:))/(1+norm(u(:)));
tol2 = norm(dv(:))/(1+norm(v(:)));
tol = max(tol1,tol2);

end

