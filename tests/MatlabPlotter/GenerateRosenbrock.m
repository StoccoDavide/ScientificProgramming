% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% MIT License                                                                     %
%                                                                                 %
% Copyright (c) 2022, Davide Stocco                                               %
%                                                                                 %
% Permission is hereby granted, free of charge, to any person obtaining a copy    %
% of this software and associated documentation files (the "Software"), to deal   %
% in the Software without restriction, including without limitation the rights    %
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       %
% copies of the Software, and to permit persons to whom the Software is           %
% furnished to do so, subject to the following conditions:                        %
%                                                                                 %
% The above copyright notice and this permission notice shall be included in all  %
% copies or substantial portions of the Software.                                 %
%                                                                                 %
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      %
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        %
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     %
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          %
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   %
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   %
% SOFTWARE.                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% clear the workspace
clc;
close all;
clear all; %#ok<CLALL>

set(0,     'DefaultFigurePosition', [5000, 5000, 560, 420]);
set(0,     'DefaultFigureWindowStyle',        'normal');
set(0,     'defaultTextInterpreter',          'latex' );
set(groot, 'defaultAxesTickLabelInterpreter', 'latex' );
set(groot, 'defaulttextinterpreter',          'latex' );
set(groot, 'defaultLegendInterpreter',        'latex' );
set(0,     'defaultAxesFontSize',             18      );

% Tests to be plotted
test = 'Rosenbrock';
path = './../Rosenbrock/';

% Solvers to be plotted
solvers = {'BroydenCombined', 'Dumped_BroydenCombined'};

% Function to be plotted
a = 1.0;
b = 10.0;
Rosenbrock_fun = @(x,y) [ a*(1-x), b*(y-x*x) ];
Rosenbrock_fun = @(x,y) norm(Rosenbrock_fun(x,y),2);

% Load the parameters

% Plot the function
figure();
hold on;
%title('Rosenbrock 2D function');
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
grid minor;

% Function plot settings
lim     = 2.5;
x_min   = -lim;
x_max   = +lim;
y_min   = -lim;
y_max   = +lim;
f_alpha = 0.5;
f_edge  = 'none';

fsurf( ...
  Rosenbrock_fun, ...
  [x_min, x_max, y_min, y_max], ...
  'FaceAlpha', f_alpha, ...
  'EdgeColor', f_edge ...
);

% Plot the iterations
x_min = NaN;
x_max = NaN;
y_min = NaN;
y_max = NaN;
for i = 1:length(solvers)

  % Catch data
  data_tmp = readtable([path, test, '_', solvers{i} ,'.csv']);
  X_tmp = table2array(data_tmp(:,1));
  X_tmp = reshape(X_tmp, [2,length(X_tmp)/2]);
  F_tmp = table2array(data_tmp(:,2));
  F_tmp = reshape(F_tmp, [2,length(F_tmp)/2]);
  D_tmp = table2array(data_tmp(:,3));
  D_tmp = reshape(D_tmp, [2,length(D_tmp)/2]);

  x_min = min(x_min, min(X_tmp(1,:)));
  x_max = max(x_max, max(X_tmp(1,:)));
  y_min = min(y_min, min(X_tmp(2,:)));
  y_max = max(y_max, max(X_tmp(2,:)));

  % plot iterations
  plot3( ...
    X_tmp(1,:), X_tmp(2,:), vecnorm(F_tmp(1:2,:),2,1), ...
    '-o', ...
    'LineWidth', 1.2, ...
    'MarkerFaceColor', 'auto' ...
  );
end
x_min = x_min * (1.0 - 0.2);
x_max = x_max * (1.0 + 0.2);
y_min = y_min * (1.0 - 0.2);
y_max = y_max * (1.0 + 0.2);
xlim([x_min, x_max]);
ylim([y_min, y_max]);

legend('$\vec{f}(x,y)$', 'BCM', 'D--BCM');

zlim([-0.5, ...
  max([Rosenbrock_fun(x_min,y_min), Rosenbrock_fun(x_min,y_max), ...
       Rosenbrock_fun(x_max,y_min), Rosenbrock_fun(x_max,y_max),])]);
hold off;

