# write a program that display an image with values in arbitrary range, and then scale it to 0-255
function disp(img)
    img = img .- minimum(img)
    img = img ./ maximum(img)
    img = img .* 255
    Gray.(img)
end

# write an equivalent of the matlab function meshgrid 
function meshgridForJulia(x,y)
    x = x' .* ones(length(y))
    y = ones(length(x))' .* y
    return x,y
end

