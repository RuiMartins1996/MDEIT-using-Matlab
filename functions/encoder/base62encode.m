function b62 = base62encode(bytes)
    import java.math.BigInteger

    % Ensure positive number by prefixing a zero byte
    bi = BigInteger([0; bytes(:)]');

    symbols = ['0':'9', 'A':'Z', 'a':'z'];
    b62 = '';

    while bi.compareTo(BigInteger.ZERO) > 0
        divRem = bi.divideAndRemainder(BigInteger.valueOf(62));
        bi = divRem(1);                   % quotient
        rem = divRem(2).intValue();       % remainder
        b62 = [symbols(rem + 1), b62];    % MATLAB 1-based indexing
    end

    if isempty(b62)
        b62 = '0';
    end
end
