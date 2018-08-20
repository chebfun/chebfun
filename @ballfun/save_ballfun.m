function save_ballfun(f,name)
% Save the tensor of coefficient of a ballfun in a mat file
F = f.coeffs;
save(name,"F");
end