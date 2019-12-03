for i=1:4
    Model{i}.salt = Fig{i}.Damian.ModelSalt;
    OBS{i}.salt = Fig{i}.Damian.OBSSalt;
end

for i=1:4
    id = find(isnan(OBS{i}.salt));
    Model{i}.salt(id) = [];
    OBS{i}.salt(id) = [];
    
    id = find(OBS{i}.salt>35);
    Model{i}.salt(id) = [];
    OBS{i}.salt(id) = [];
end