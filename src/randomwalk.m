function [Vmed, Vmax, Vmin, Vall] = randomwalk(data, steps)

data = data(unique(convhulln(data)), :);
[q, r] = qr(rand(3));
position = q*sign(diag(diag(r)));
Vall = zeros(steps, 1);

for step = 1:steps,
    rotated = data*position';
    Vall(step) = prod(max(rotated) - min(rotated));
    
    [q, r] = qr(rand(3));
    position = q*sign(diag(diag(r)))*position;
end

Vmed = median(Vall);
Vmax = max(Vall);
Vmin = min(Vall);
