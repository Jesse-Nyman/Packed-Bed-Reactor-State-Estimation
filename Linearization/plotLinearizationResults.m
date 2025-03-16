%% POSSIBLE FUNCTION FOR PLOTTING; NOT CURRENTLY USED

function plotLinearizationResults(linearizeSteps,A_matrices,B_matrices,C_matrices,D_matrices)
% Just a placeholder to plot norms of matrices over time

time = linearizeSteps; % steps at which we linearized
Anorm = arrayfun(@(ii) norm(A_matrices{ii},2), 1:length(linearizeSteps));
figure('Name','A-matrix Norm Over Time','NumberTitle','off','Color','w');
plot(time,Anorm,'-o','LineWidth',2);
xlabel('Time step');
ylabel('||A||_2');
title('Evolution of A-matrix Norm Over Time');
grid on;

end
