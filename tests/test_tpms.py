import microgen
import numpy as np

def test_tpms():
    rve = microgen.Rve(dim_x=1, dim_y=1, dim_z=1, size_mesh=0.03)
    bounding_sphere_radius = 1.1*np.sqrt((0.5*rve.dx)**2 + (0.5*rve.dy)**2 + (0.5*rve.dz)**2)

    function = "cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)"
    height = 0.3

    tpms = microgen.Custom(rve, function, height)
    assert tpms.eval([0, 0, 0]) == 3 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    
    max_edges_size_at_feature_edges = 4
    tpms = microgen.Hyperboloid(max_edges_size_at_feature_edges)
    tpms.getFeatures()
    assert tpms.eval([0, 0, 0]) == -0.5**2
    assert tpms.get_bounding_sphere_squared_radius() == (1 + 0.5)**2 + 1**2
    tpms = microgen.SchwarzP(rve, height)
    assert tpms.eval([0, 0, 0]) == 3 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.SchwarzD(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 0 + 0 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.Neovius(rve, height)
    assert tpms.eval([0, 0, 0]) == (3 + 1 + 1) + (4*1*1*1) + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.SchoenIWP(rve, height)
    assert tpms.eval([0, 0, 0]) == 2*(1+1+1) - (1+1+1) + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.SchoenFRD(rve, height)
    assert tpms.eval([0, 0, 0]) == 4 - (1+1+1) + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.FischerKochS(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 0 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.PMY(rve, height)
    assert tpms.eval([0, 0, 0]) == 2 + 0 + 0 + 0 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.HoneyComb(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 1 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.Gyroid(rve, height)
    assert tpms.eval([0, 0, 0]) == 1
    # assert tpms.eval([1, 0, 0]) == 0 + 0 + 0 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius