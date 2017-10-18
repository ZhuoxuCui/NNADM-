 function [Dux,Duy] = Grad(U)
  % Forward finite difference operator
  %  diff(X,N,DIM) the Nth difference function along dimension DIM. 
  Dux = [diff(U,1,2), U(:,1) - U(:,end)];
  Duy = [diff(U,1,1); U(1,:) - U(end,:)];