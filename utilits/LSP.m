function Ax = LSP(H,HT,S,x,n1,n2,n3,beta2)
xx = reshape(x,n1,n2,n3);
Ax = HT(S.*H(xx))+beta2*xx;
Ax = Ax(:);
end