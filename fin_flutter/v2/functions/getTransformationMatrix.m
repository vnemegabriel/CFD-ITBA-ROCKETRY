function [T, nodes2D] = getTransformationMatrix(nodes3D)
    % nodes3D: [4x3] coordenadas globales (X, Y, Z)

    % 1. Definir eje x' (del nodo 1 al 2)
    v12 = nodes3D(2,:) - nodes3D(1,:);
    ex = v12 / norm(v12);

    % 2. Definir eje z' 
    v13 = nodes3D(3,:) - nodes3D(1,:);
    ez = cross(ex, v13);

    ez = ez / norm(ez);

    % 3. Definir eje y' (ortogonal)
    ey = cross(ez, ex);

    % Matriz de rotación R (3x3)
    R = [ex; ey; ez];

    % Matriz de transformación T global (24x24)
    % Bloque de 6x6 por cada nodo
    T_node = zeros(6,6);
    T_node(1:3, 1:3) = R;
    T_node(4:6, 4:6) = R;

    T = blkdiag(T_node, T_node, T_node, T_node);

    % Proyectar nodos 3D a 2D locales (para usar en Bm, Bb, Bs)
    % Tomamos el nodo 1 como origen local
    origin = nodes3D(1,:);
    nodes2D = zeros(4,2);
    for i = 1:4
        rel_vec = nodes3D(i,:) - origin;
        nodes2D(i,1) = dot(rel_vec, ex);
        nodes2D(i,2) = dot(rel_vec, ey);
    end
end

% function [T, nodes2D] = getTransformationMatrix(nodes3D)
%     % Soporta N nodos dinámicamente (3 para T3, 4 para Q4)
%     numNodes = size(nodes3D, 1);
% 
%     v12 = nodes3D(2,:) - nodes3D(1,:);
%     v13 = nodes3D(3,:) - nodes3D(1,:);
%     ez = cross(v12, v13);
%     ez = ez / norm(ez);
% 
%     % FORZAR EL EJE LOCAL X A SER EL EJE GLOBAL Y (Longitudinal)
%     % Esto cura el zigzag. Garantiza que Nx y Qy extraídos signifiquen 
%     % físicamente lo mismo (longitudinal y circunferencial) para TODOS 
%     % los triángulos sin importar su orientación diagonal.
%     Y_global = [0, 1, 0];
%     ex = Y_global - dot(Y_global, ez) * ez;
%     ex = ex / norm(ex);
%     ey = cross(ez, ex);
% 
%     R = [ex; ey; ez];
% 
%     % Bloque de 6x6 por cada nodo
%     T_node = zeros(6,6);
%     T_node(1:3, 1:3) = R;
%     T_node(4:6, 4:6) = R;
% 
%     % Matriz de transformación dinámica (18x18 o 24x24)
%     T_cell = repmat({T_node}, 1, numNodes);
%     T = blkdiag(T_cell{:});
% 
%     % Proyectar nodos 3D a 2D locales
%     nodesLoc = (R * (nodes3D - nodes3D(1,:))')';
%     nodes2D = nodesLoc(:, 1:2);
% end