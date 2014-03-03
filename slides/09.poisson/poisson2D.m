function e = poisson2D(N, x0, x1, f, e)
  len = (N-1)*(N-1);
  A = spalloc(len, len, 5*len);
  h = (x1-x0)/N;
  A=spdiags( 4*ones(len,1), 0,    A) + ...
    spdiags(-1*ones(len,1),-1,    A) + ...
    spdiags(-1*ones(len,1), 1,    A) + ...
    spdiags(-1*ones(len,1),N-1,   A) + ...
    spdiags(-1*ones(len,1),-(N-1),A);
%  for j=1:N-2
%    A(j*(N-1), j*(N-1)+1) = 0.0;
%  end
%  for j=1:N-2
%    A(j*(N-1)+1, j*(N-1)) = 0.0;
%  end
  fgrid = linspace(x0, x1, N+1);
  g = h*h*feval(f, fgrid(2:end-1));
  us = A\g;
  u = reshape(us, N-1, N-1);
  clf;
  surf(fgrid(2:end-1), fgrid(2:end-1), u);
  legend('finite difference solution');
  set(gca,'fontsize',18);

  ers = feval(e, fgrid(2:end-1));
  e = norm(ers-us,'inf');
