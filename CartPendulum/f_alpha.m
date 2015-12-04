function val = f_alpha(x, wmax)
% density of a squared uniform variable distributed between [-wmax, wmax]

val = zeros(size(x));

for k = 1:length(x)
    if x(k) >= 0 && x(k) <= wmax^2
        val(k) = 1/(2*wmax*sqrt(x(k)));
    end
end

end