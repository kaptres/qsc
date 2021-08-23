function FV = OBJ2FV (obj_struct)
    FV.vertices = obj_struct.vertices;
    FV.faces    = obj_struct.objects(end).data.vertices;
end
