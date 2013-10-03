function test_get_img_patch()
%TEST_GET_IMG_PATCH Summary of this function goes here
%   Detailed explanation goes here
A = [[1 0 1 0];
     [0 1 0 1];
     [1 0 1 0];
     [0 1 0 1]];
 
expect = [[-1 -1 -1 -1 -1];
          [-1 1 0 1 0];
          [-1 0 1 0 1];
          [-1 1 0 1 0];
          [-1 0 1 0 1];];
      
result = get_img_patch(A, 2, 2, 2, 2, 'fillval', -1);

if ~all(expect == result)
    disp('Doh!');
    expect
    result
else
    disp('Woohoo!');
end


end

