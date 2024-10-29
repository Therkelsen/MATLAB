% Specify desired width and height
width = 800; % Adjust as necessary
height = 600; % Adjust as necessary

% Get all figure handles
figHandles = findall(0, 'Type', 'figure');

% Set dimensions for each figure
for i = 1:length(figHandles)
    set(figHandles(i), 'Position', [100, 100, width, height]);
end
