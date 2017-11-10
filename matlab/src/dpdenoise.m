function val = dpdenoise(D_x, D_y, N_int, N_ext, xy_int, xy_ext_DP, DP_xy_scale),

if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end

mod = py.importlib.import_module('dpdenoise');

xy_int_array = py.numpy.array(xy_int(:)');
xy_int_array.reshape([D_x, D_y]);

xy_ext_array = py.numpy.array(xy_ext_DP(:)');
xy_ext_array.reshape([D_x, D_y]);

data = py.dict(pyargs('D_x', int32(D_x), 'D_y', int32(D_y), ...
                      'N_int', int32(N_int), 'N_ext', int32(N_ext), ...
                      'xy_int', xy_int_array, ...
                      'xy_ext_DP', xy_ext_array, ...
                      'DP_xy_scale', DP_xy_scale));
res = py.dpdenoise.optimize_xy(data);
res = cell(res);
val = double(py.array.array('d', res{1}))';

% py.dpdenoise.test_me(int32(10), int32(1), int32(10), int32(100))
