function par = read_parameter_file(sim_name)
%read_parameter_file(sim_name) reads through sim_name.par and assigns the
%values within to struct fields of the same name as the original simulation
%variables (given that the .par file contains these names as expected). It
%requires no knowledge of the order in which the parameters are stored etc.
%and thus will still work if the main code has been modified to include
%additional data in the .par and so on.

fid = fopen([sim_name '.par']);

while ~feof(fid)
    
    val = fscanf(fid, '%f', 1);
    
    if isempty(val)
        break;
    end
    
    var_name = fscanf(fid, '%s', 2);
    var_name = var_name(3:end);
    
    par.(var_name) = val;
    
end

fclose(fid);

end

