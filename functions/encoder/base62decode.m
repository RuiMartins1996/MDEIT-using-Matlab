function bytes = base62decode(b62)
    import java.math.BigInteger

    symbols = ['0':'9', 'A':'Z', 'a':'z'];
    base = BigInteger.valueOf(62);
    bi = BigInteger.ZERO;

    for i = 1:length(b62)
        idx = find(symbols == b62(i)) - 1;
        bi = bi.multiply(base).add(BigInteger.valueOf(idx));
    end

    bytes = typecast(bi.toByteArray(), 'uint8');

    % Remove leading zero added during encoding
    if bytes(1) == 0
        bytes = bytes(2:end);
    end
end