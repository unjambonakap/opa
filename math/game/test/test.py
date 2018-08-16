import sys

from chdrft.cmds import CmdsList
from chdrft.graphics.helper import stl_to_meshdata
from chdrft.graphics.loader import stl_parser
from chdrft.main import app
from chdrft.utils.cmdify import ActionHandler
from chdrft.utils.misc import Attributize, flatten
from chdrft.utils.swig import swig
from opa.math.game.proto import common_pb2
from vispy.geometry.meshdata import MeshData
from vispy.scene.cameras import TurntableCamera
from vispy.scene.visuals import Markers
from vispy.scene.visuals import Mesh, Line
from vispy.visuals import transforms
from vispy.visuals.collections import PointCollection
import numpy as np
import numpy.random
import pprint
import vispy, vispy.app, vispy.scene, vispy.visuals
from asq.initiators import query as q

common = swig.opa_common_swig
mathgame = swig.opa_math_game_swig
import chdrft.display.vispy_utils as vispy_utils

def main():
  common = swig.opa_common_swig
  mathgame = swig.opa_math_game_swig
  pts = [[4.0, 0.5, 6.0], [4.0, 0.5, 1.0], [4.0, 0.0, 6.0], [4.5, 0.0, 6.0],
        [4.0, 0.5, 6.0], [4.0, 0.0, 6.0], [4.5, 0.5, 1.0], [4.0, 0.5, 1.0],
        [4.0, 0.5, 6.0], [4.5, 0.0, 6.0], [4.0, 0.0, 6.0], [4.0, 0.0, 1.0],
        [4.0, 0.5, 1.0], [4.0, 0.0, 1.0], [4.0, 0.0, 6.0], [4.5, 0.0, 6.0],
        [4.5, 0.5, 6.0], [4.0, 0.5, 6.0], [4.5, 0.5, 6.0], [4.5, 0.5, 1.0],
        [4.0, 0.5, 6.0], [4.5, 0.5, 1.0], [4.5, 0.5, 6.0], [4.5, 0.0, 6.0],
        [4.5, 0.0, 1.0], [4.5, 0.0, 6.0], [4.0, 0.0, 1.0], [4.5, 0.5, 1.0],
        [4.5, 0.0, 6.0], [4.5, 0.0, 1.0]]

  pts = [ [0,0,0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
  pts = ((0.500000,6.000000,0.000000),(0.500000,1.000000,0.000000),(0.000000,6.000000,0.000000),(0.500000,6.000000,0.000000),(0.500000,1.000000,0.500000),(0.500000,1.000000,0.000000))
  bbox = mathgame.compute_best_box(pts)
  print(bbox.str())
  #print(bbox.str())
  aabbox = mathgame.compute_aabb_box_norot(pts)
  print(bbox.area())
  print(bbox.vec(0))
  print(bbox.vec(1))
  print(bbox.vec(2))

  print(aabbox.area())

  tc = mathgame.FaceCollection(True)
  mesh = tc.clear().load_stl('robotics/testy.stl').triangulate().to_mesh()
  atomizer = mathgame.KAtomizer(mesh, 3)
  atomizer.debug = True

  res = atomizer.initial_assignment()
  for i in res.debug_data.init_round_debug:
    print(i.selected_candidate.priority)

  for i in range(4):
    res = atomizer.initial_assignment()

  data = list(vispy_utils.plot_mesh_dirs(mesh))
  data[0].points = None
  for cluster_data in res.cluster_data:
    print("FUUUU ", cluster_data.box.almost_empty())
    if cluster_data.box.almost_empty(): continue
    mx = tc.clear().add_box(cluster_data.box).to_mesh()
    pts = []
    pts.extend(vispy_utils.vec_pos_to_list(atomizer.get_points(cluster_data.faces)))
    data.append(vispy_utils.get_mesh_data(mx))
    data[-1].points = pts

  ctx = vispy_utils.VispyCtx()
  ctx.cur_objs = ctx.plot_meshes(*data)

  ctx.run()
main()
