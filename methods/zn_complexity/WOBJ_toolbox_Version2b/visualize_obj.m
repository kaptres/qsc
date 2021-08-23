% author: mferhata
function visualize_obj (obj_struct)
    vertices    = obj_struct.vertices;
    faces       = obj_struct.objects(end).data.vertices;
    set (gcf, 'color', 'w');
    p = patch ( 'vertices', vertices, 'faces', faces, ...
                'EdgeColor', 'none');
    axis equal off;
    lighting none;
    %material ([.4 .2 .5 .5]);
    %h = camlight('left');
    %for i = 1:100
    %   camorbit(2.2,1)
    %   pause(.1)
    %end
end
