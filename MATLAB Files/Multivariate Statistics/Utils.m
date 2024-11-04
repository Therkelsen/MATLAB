classdef Utils
    %UTILS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function [mean_values, cov_matrix, range_values] = calculate_descriptive_statistics(data, print, show_figures, data_name)
            % Calculate mean, covariance, and range separately
            mean_values = mean(data);
            cov_matrix = cov(data);
            range_values = range(data);

            % Define column names based on the number of columns in your data
            headers = arrayfun(@(x) sprintf('Dim %d', x), 1:size(data, 2), 'UniformOutput', false);
            
            % Create a table for the mean, covariance diagonal, and range with headers
            descriptive_statistics = [mean_values; diag(cov_matrix)'; range_values];
        
            % Define row names for the table
            rows = {'Mean', 'COV (Diagonal)', 'Range'};
            
            % Create a table with row names and headers
            descriptive_statistics_table = array2table(descriptive_statistics, 'RowNames', rows, 'VariableNames', headers);

            % Pass true or false through to decide whether to print the descriptive statistics.
            if print
                
                % Display mean values
                disp(array2table(mean_values, 'VariableNames', headers, 'RowNames', {'Mean'}));
                
                fprintf('\nCovariance Matrix (Diagonal):\n');
                % Display Covariance matrix diagonal
                disp(array2table(cov_matrix, 'VariableNames', headers));
                
                % Display range values
                fprintf('\nRange Values:\n');
                disp(array2table(range_values, 'VariableNames', headers, 'RowNames', {'Range'}));
            end
            % Pass true or false through to decide whether to display the figures
            if show_figures
                % Create a line plot
                figure;
                % Transpose and plot with markers
                plot(table2array(descriptive_statistics_table)', 'o-');
                
                % Add x-axis labels and title
                xlabel('Category');
                ylabel('Values');
                title(data_name + ': Mean, Covariance Diagonal, and Range for Each Category');

                % Set x-axis labels
                set(gca, 'XTick', 1:numel(headers), 'XTickLabel', headers);
                
                % Add legend for Mean, Covariance (Diagonal), and Range
                legend('Mean', 'COV (Diagonal)', 'Range');
                
                % Rotate x-axis labels for better readability (optional)
                xtickangle(45);
                
                % Display the plot
                grid on;
    
                % Create a scatter matrix with marginal histograms and boxplots
                figure;
                plotmatrix(data);
                xlabel('Variables');
                ylabel('Variables');
                title(data_name + ': Scatter Matrix with Marginal Histograms and Boxplots');
            end
        end

        function [mahalanobis_distances] = calculate_mahalanobis_distances(data, mu, S, print, show_figures, data_name)
            mahalanobis_distances = [];
            
            for j=1:size(data, 1)
                mahalanobis_distances = [mahalanobis_distances; (data(j, :) - mu) * inv(S) * (data(j, :) - mu)'];
            end
            
            % Pass true or false through to decide whether to print the descriptive statistics.
            if print
                array2table(mahalanobis_distances, 'VariableNames', {'Mahalanobis Distances'})
            end
            % Pass true or false through to decide whether to display the figures
            if show_figures
                figure
                qqplot(mahalanobis_distances, chi2rnd(ones(1, length(data(:, 1)))))
                title(data_name + ': QQplot for Mahalanobis distances vs. Χ² distances')
            end
        end

        function plot2d_CR_for_mu_ellipsis (mu_hat, SIGMA_hat, alpha, n)
            t = 0:0.01:2*pi;
            N = length(t);
            p = 2;
            [V,D] = eig(SIGMA_hat);
            a = sqrt(D(2,2))*sqrt((p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p));
            b = sqrt(D(1,1))*sqrt((p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p));
            P  = [a*cos(t); b*sin(t)];
            theta = atan2(V(2,2),V(1,2));
            T = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            P_rot = T*P + mu_hat*ones(1,N);
            plot(P_rot(1,:),P_rot(2,:),'LineWidth',3,'Color','k'),grid
        end
    end
end
