function I = build_current_patterns(M)
%BUILD_CURRENT_PATTERNS  I^(j) = e_1 - e_{j+1}, j=1..M-1 (paper sec 5: "one
%feeding electrode, current exits in turns from the remaining ones").
%
%   I = build_current_patterns(M)   % I is [(M-1) x M], row j = I^(j)'

N = M-1;
I = zeros(N,M);
for j = 1:N
    I(j,1) = I(j,1) + 1;
    I(j,j+1) = I(j,j+1) - 1;
end
end
