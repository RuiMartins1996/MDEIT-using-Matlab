function S = decodeModelParameters(hashStr)
    compressed = base62decode(hashStr);  % Custom function below
    jsonStr = char(zlibdecode(compressed));
    S = jsondecode(jsonStr);
end