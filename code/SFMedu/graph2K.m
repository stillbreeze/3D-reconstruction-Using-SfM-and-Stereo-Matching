function K = graph2K(graph)

K = eye(3);
K(1,1) = graph.fx;
K(2,2) = graph.fy;
K(1, 3) = px;
K(2, 3) = py;