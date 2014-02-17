function conv2
  pts = linspace(10, 100, 10);
  ers = zeros(size(pts));
  for i=1:length(pts)
    ers(i) = poisson2D(pts(i), 0, 1, @source2, @exact2);
  end
  ers
  figure; loglog(pts, ers, 'bx-', 'linewidth', 2); hold on; 
  reference = 10.^([-2 -4 -6]);
  loglog([10^1 10^2 10^3],reference,'ko-','linewidth',2);
  legend('convergence rate','second order reference');
  set(gca, 'fontsize',18);
