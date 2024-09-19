classdef Utils
    %UTILS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function [mean_values, cov_matrix, range_values] = calculate_descriptive_statistics(data, headers)
            % Calculate mean, covariance, and range separately
            mean_values = mean(data);
            cov_matrix = cov(data);
            range_values = range(data);
            
            % Create a table for the mean, covariance diagonal, and range with headers
            descriptive_statistics = [mean_values; diag(cov_matrix)'; range_values];
            
            % Define row names for the table
            rows = {'Mean', 'COV (Diagonal)', 'Range'};
            
            % Create a table with row names and headers
            descriptive_statistics_table = array2table(descriptive_statistics, 'RowNames', rows, 'VariableNames', headers);
            
            % Create a line plot excluding marathon, so the plot looks better.
            figure(1);
            % Transpose and plot with markers
            plot(table2array(descriptive_statistics_table)', 'o-');
            
            % Add x-axis labels and title
            xlabel('Category');
            ylabel('Values');
            title('Mean, Covariance Diagonal, and Range for Each Category');
            
            % Set x-axis labels
            set(gca, 'XTick', 1:numel(headers), 'XTickLabel', headers);
            
            % Add legend for Mean, Covariance (Diagonal), and Range
            legend('Mean', 'COV (Diagonal)', 'Range');
            
            % Rotate x-axis labels for better readability (optional)
            xtickangle(45);
            
            % Display the plot
            grid on;

            % Create a scatter matrix with marginal histograms and boxplots
            figure(2);
            plotmatrix(data);
            xlabel('Variables');
            ylabel('Variables');
            title('Scatter Matrix with Marginal Histograms and Boxplots');
        end
    end
end
