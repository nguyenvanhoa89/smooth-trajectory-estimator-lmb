function plot_card(truth,est)
  
  num_steps = length(truth.X);
  true_card = zeros(1,num_steps);
  est_card = zeros(1,num_steps);

  for i = 1:num_steps
    true_card(i) = size(truth.X{i},2);
    est_card(i) = size(est.X{i},2);
  end
  
  figure;
  axes;
  hold on;
  grid on;
  box on;
  plot(true_card,'k-');
  plot(est_card,'Color','b','Marker','o');
  xlabel('Time Index');
  ylabel('Cardinality');

end