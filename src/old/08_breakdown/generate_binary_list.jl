
function int2bin_arr(
    x      :: Integer,
    digits :: Integer,
)


arr = zeros(Integer, digits)


for i in 1:digits
    arr[end-i+1] = x >> (i-1) & 0b1
end

return arr
end



function gen_bin_arr(
    digits :: Integer
)

    max = 1 << digits
    
    arr = zeros(Integer, max, digits)

    for i = 1 : max
        arr[i, :] = int2bin_arr(i-1, digits)
    end

    return arr
end
