function str = str_vec(vec)
%STR_VEC Output a string representation of input vector/cell VEC. Assumes
%that VEC is an 1xN or Nx1 vector/cell of numbers.
% Borrowed from:
%   http://stackoverflow.com/questions/14924181/how-to-display-print-vector-in-matlab
s = repmat('%d,', 1, length(vec));
s(end) = []; %Remove trailing comma
if isa(vec, 'cell')
    vec = cell2mat(vec);
end
str = sprintf(['(' s ')'], vec);
end

