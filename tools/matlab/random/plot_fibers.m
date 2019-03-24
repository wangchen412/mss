function plot_fibers(C, S, O)
  figure()
  set(gcf, 'position', [0, 0, 800, 800])
  scatter(S(:,1), S(:,2), 120, 'b', 'filled')
  xlim([-3 3])
  ylim([-3 3])

  hold on
  scatter(C(:,1), C(:,2), 120, 'r', 'filled')
  scatter(O(:,1), O(:,2), 120, 'b', 'filled')
  hold off
end
