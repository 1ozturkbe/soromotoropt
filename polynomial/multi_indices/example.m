%% Total order (total degree) multi-indices
tomil = TotalOrderMultiIndexLattice(3, 2);
% first argument is the dimension, second argument is the polynomial degree
% (l_1 truncation)
tomil.init();

Mtotalorder = tomil.get_midx_matrix()


%% Full tensor multi-indices

tensormil = FullOrderMultiIndexLattice(3, 4);
% first argument is the dimension, second argument is the polynomial degree
% (l_\infty truncation)
tensormil.init();

Mtensor = tensormil.get_midx_matrix()