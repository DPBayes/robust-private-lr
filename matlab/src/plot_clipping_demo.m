SIZES = [300, 1000, 3000, 10000];
DIM = 10;
FONTSIZE = 6;

err_r2 = {};
err_rho = {};

st=[0.1:0.1:2];
[X, Y] = meshgrid(st, st);

for k=1:length(SIZES),
  [err_r2{k}, err_rho{k}] = clipping_demo(DIM, SIZES(k));
end

V = [0:20] / 20;
for k=1:length(SIZES),
  subplot(length(SIZES), 2, 2*k-1);
  c = contourf(X, Y, err_r2{k}, V);
  colorbar;
  set(gca, 'FontSize', FONTSIZE)
  xlabel('B_x / \sigma_x')
  ylabel('B_y / \sigma_y')
  title(sprintf('R^2, n=%d', SIZES(k)))
  subplot(length(SIZES), 2, 2*k);
  c = contourf(X, Y, err_rho{k});
  colorbar;
  set(gca, 'FontSize', FONTSIZE)
  xlabel('B_x / \sigma_x')
  ylabel('B_y / \sigma_y')
  title(sprintf('rho, n=%d', SIZES(k)))
end

set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0 0 15.0, 18.0])
print -depsc2 ../../../../drugsens/tex/figures/clipping_contours.eps
