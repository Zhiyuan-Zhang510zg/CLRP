function I = Synthesis_circles_fun(N,Width,Sigma)
% Synthesize a circles image with sigma noise
% recommend parameters
% N = 64;
% Width = 5;
% Sigma = 0.1;

Radius = [10,20,30];
n1=-N/2:N/2-1;   
n2=-N/2:N/2-1;
u=zeros(N,N);
[x,y]=meshgrid(n1,n2);
circle=x.^2+y.^2;
for i = 1 : length(Radius)
    R = Radius(i);
    r = R-Width;
    A = (circle>=r*r) .* (circle<R*R);
    u(find(A==1))=1; 
end

v = [zeros(N,1),ones(N,1)];
for k = 1: N/2-1
    v = [v,zeros(N,1),ones(N,1)];
end

I = (u+v)/2 + Sigma*rand(N);
end