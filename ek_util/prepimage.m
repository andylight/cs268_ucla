function Iout = prepimage(I)
% Prepare input image prior to showing via subimage. Shifts input image 
% values s.t. they are in the range [0,1].
minval = min(I(:));
maxval = max(I(:));
Iout = (I - minval) / (maxval - minval);
end
