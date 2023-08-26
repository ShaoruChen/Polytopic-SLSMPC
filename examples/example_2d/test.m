% Create some data
x = 1:10;
y = rand(7,10);

% Create figure and subplots
figure;

% Subplot 1
subplot(1, 2, 1);
for i = 1:7
    h(i) = plot(x, y(i,:), 'LineWidth', 2); 
    hold on;
end
title('Subplot 1');

% Subplot 2
subplot(1, 2, 2);
for i = 1:7
    plot(x, y(i,:), 'LineWidth', 2); 
    hold on;
end
title('Subplot 2');

% Create legends
legend1 = legend([h(1), h(2), h(3), h(4)], 'First', 'Second', 'Third', 'Fourth', 'NumColumns', 4, 'Location', 'southoutside');
legend1.Box = 'off';

% Make a copy of the axes, which will contain the second legend
axes_copy = copyobj(gca, gcf);

% Create the second legend
legend2 = legend(axes_copy, [h(5), h(6), h(7)], 'Fifth', 'Sixth', 'Seventh', 'NumColumns', 3, 'Location', 'southoutside');
legend2.Box = 'off';

% Get the positions
pos1 = get(legend1, 'Position');
pos2 = get(legend2, 'Position');

% Adjust the position of the second legend
set(legend2, 'Position', [pos1(1) + pos1(3) - pos2(3), pos1(2), pos2(3), pos2(4)]);

% Hide the copied axes
set(axes_copy, 'Visible', 'off');
