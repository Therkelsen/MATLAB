classdef Utils
    %UTILS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function descriptive_statistics = calculate_descriptive_statistics(data, headers)
            descriptive_statistics = [mean(data);
                std(data);
                range(data)];
            
            % Define row names for descriptive statistics
            rows = {'Mean', 'STD', 'Range'};
            
            % Create a table for descriptive statistics with the same column headers
            descriptive_statistics = array2table(descriptive_statistics, 'RowNames', rows, 'VariableNames', headers)
            
            % Create a line plot excluding marathon, so the plot looks better.
            figure;
            % Transpose and plot with markers
            plot(table2array(descriptive_statistics(:, 1:end))', 'o-');
            
            % Add x-axis labels and title
            xlabel('Category');
            ylabel('Values');
            title('Mean, STD, and Range for Each Category');
            
            % Add x-axis labels
            set(gca, 'XTick', 1:numel(headers), 'XTickLabel', headers);
            
            % Add legend for Mean, STD, and Range
            legend('Mean', 'STD', 'Range');
            
            % Rotate x-axis labels for better readability (optional)
            xtickangle(45);
            
            % Display the plot
            grid on;

            % Make plots of the multivariate scatter matrix plus marginal histograms and boxplots
            % Create a scatter matrix with marginal histograms and boxplots
            figure;
            plotmatrix(data);
            xlabel('Variables');
            ylabel('Variables');
            title('Scatter Matrix with Marginal Histograms and Boxplots');
            end
    end
end