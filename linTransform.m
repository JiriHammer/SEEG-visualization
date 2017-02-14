function y = linTransform(x, fromLimits, toLimits)
% linearly transforms values x lying in interval 'fromLimits' into another interval 'toLimits'
% eq: y = k*x + c

% (c) Jiri, Mar11, Nov16

% input
x1 = fromLimits(1);
x2 = fromLimits(end);
y1 = toLimits(1);
y2 = toLimits(end);

% linear transformation
k = (y2-y1)/(x2-x1);    % slope
c = -k*x1 + y1;         % offset
y = k*x + c;

%% old implementation
% y = (max(toLimits) - min(toLimits))/(max(fromLimits) - min(fromLimits))*(x - min(fromLimits)) + min(toLimits);
% assert(length(x) == length(y));
