function dist = H_dist(a,b)
% Computes Hausdorf/bottleneck distance between column vectors a and b

% a = [.1 .5 .7]';
% b = [.13 .56 .78]';

La = length(a);
Lb = length(b);

distab = 0;
for k = 1 : La
    distab = max(distab,min(min(abs(a(k)-b),1-abs(a(k)-b))));
end

distba = 0;
for k = 1 : Lb
    distba = max(distba,min(min(abs(b(k)-a),1-abs(b(k)-a))));
end

dist = max(distab,distba);