% author: mferhata
% the variable FILENAME should be assigned before calling PARSE_SUBMISSION2
f = fopen (filename);
while true
    subcollection   = fgetl (f);
    if subcollection == -1
        fclose (f);
        break;
    end

    measurements    = fgetl (f);
    parsed          = [];

    shape_index     = 1;
    string_index    = 1;
    while true
        string_index    = string_index + 1;
        [m, c, e, n]    = sscanf (measurements(string_index:end), '%f');
        parsed          = [parsed; m'];
        string_index    = string_index + n; 
        
        n               = strfind (measurements(string_index:end), '(');
        if isempty (n)
            break;
        else
            string_index= string_index + n(1) - 1;
        end
    end
    switch subcollection
    case 'c+'
        cp  = parsed;
    case 'c-'
        cm  = parsed;
    case 's+'
        sp  = parsed;
    case 's-'
        sm  = parsed;
    case {'collection2_1' 'c21'}
        c21 = parsed;
    case {'collection2_2' 'c22'}
        c22 = parsed;
    case {'off' 'collection3'}
        off = parsed;
    case 'collection3'
        off = parsed;
    end
end
