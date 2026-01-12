% function hashStr = encodeModelParameters(S)
%     jsonStr = jsonencode(S);
%     compressed = zlibencode(uint8(jsonStr));
%     hashStr = base62encode(compressed);  % Custom function below
% end


function [hashStr, shortHash] = encodeModelParameters(S)
    % Serialize the struct to JSON
    jsonStr = jsonencode(S);
    
    % Compress it
    compressed = zlibencode(uint8(jsonStr));
    
    % Full base62 encoding (can be long)
    hashStr = base62encode(compressed);

    % Generate a short hash (MD5, 16 hex chars = 64 bits of uniqueness)
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(hashStr));
    digest = typecast(md.digest, 'uint8');
    shortHash = lower(dec2hex(digest))';    % hex string
    shortHash = shortHash(:)';              % make it a row string
    shortHash = shortHash(1:16);            % take first 16 chars
end
