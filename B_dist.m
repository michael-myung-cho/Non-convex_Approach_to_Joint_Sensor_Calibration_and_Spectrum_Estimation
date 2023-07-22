function dist = B_dist(a,b)
% compute bottleneck distance between a and b
% for example
% a = [.1 .5 .7]';
% b = [.13 .56 .78]';

% a and b must be of the same length
if length(a) ~= length(b)
    fprintf('two vectors must be of the same length')
    return;
else
    sparsity = length(a);
end

a = sort(a,'ascend');
b = sort(b,'ascend');

% make a and b column vectors
if size(a,1) == 1
    a = a';
end
if size(b,1) == 1
    b = b';
end

dist = zeros(sparsity,1);
for k = 0:sparsity-1
    ashift     = circshift(a,k);
    dist(k+1)  = max(min(abs(ashift-b),1-abs(ashift-b)));
end
dist = min(dist);
    
