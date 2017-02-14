function idx = closestval(a,b,w)
% function closestval - find closest match to a given value
%
% syntax: idx = closestval(a,b)
%         idx = closestval(a,b,w)
%
% closestval(a,b) finds the index idx of the value in vector a closest to 
% the scalar b. For several identical nearest values, idx only returns the
% first appearance in a.
% If a is a matrix, idx returns a vector with indices of the closest values
% in each column of a.
% closestval(a,b,w) with w being >0 returns the index of the closes value
% greater than or equal to b. If w<0 it works for the nearest value <=b.
% w=0 (default) picks the value with the smallest absolute difference.
% For w<0 or w>0, if a has no values <b or >b, respectively, idx is 0. 
% a, b and w have to be numeric and real. b and w should be scalars,
% otherwise the first value will be used.

% Tobias Pistohl, May 2007

if nargin<3
    w= 0;
end

if ~isreal(a) || ~isreal(b) || ~isreal(w) || ...
        ~isnumeric(a) || ~isnumeric(b) || ~isnumeric(w)
    error('This function works only on numeric real values.')
end

if ~isvector(a)
    idx= zeros(1,size(a,2));
end
   
if ~isscalar(b)
    b= b(1);
end
if ~isscalar(w)
    w= w(1);
end

w= sign(w);

switch w
    case 1
        d= a-b;
        d(find(d<0))= inf;
    case -1
        d= b-a;
        d(find(d<0))= inf;
    otherwise
        d= abs(a-b);
end

if isvector(d)
    d= d(:);
end
smallest= min(d);
for c= 1:length(smallest)
    idx(c)= find(d(:,c)==smallest(c),1);
    if isinf(smallest(c))
        idx(c)= 0;
    end
end