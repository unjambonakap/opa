#!/usr/bin/env python

from chdrft.cmds import CmdsList
from chdrft.main import app
from chdrft.utils.cmdify import ActionHandler
from chdrft.utils.misc import Attributize
import chdrft.utils.misc as cmisc
from chdrft.utils.swig import swig
import glog
from chdrft.geo.satsim import TileGetter
import math
import numpy as np
import time
import traceback as tb
import sys
import pprint
import ephem
from astroquery.jplhorizons import Horizons
from chdrft.utils.cache import Cachable
import pandas as pd
from scipy.interpolate import interp1d

global flags, cache
flags = None
cache = None

global g_ctx


def args(parser):
  clist = CmdsList().add(test).add(earth).add(test_inter).add(create_db).add(test_pc).add(
      test_ephem)
  parser.add_argument(
      '--resources-dir', default=cmisc.path_here('../../../etu/meteor_stuff/satsim/resources'))
  parser.add_argument('--mosaic', action='store_true')
  parser.add_argument('--precision', type=float, default=2.)
  parser.add_argument('--split-tr-norm', type=float, default=1.)
  parser.add_argument('--limit-fps', type=float, default=60)
  parser.add_argument('--dump-stats', action='store_true')
  parser.add_argument('--influx-db', default='db1')
  parser.add_argument('--stl-file', type=cmisc.cwdpath)
  parser.add_argument('--noengine', action='store_true')
  ActionHandler.Prepare(parser, clist.lst, init=init)


class InfluxDbRun:

  def __init__(self, enable):
    self.enable = enable
    if not enable: return
    from influxdb import InfluxDBClient
    client = InfluxDBClient('localhost', 8086, 'root', 'root', database=flags.influx_db)
    self.client = client
    self.tags = dict(
        runid=flags.runid,
        xtime=time.time(),)

  def normalize(self, data):
    for k, v in data.items():
      nv = v
      if (isinstance(v, str) or isinstance(v, int) or isinstance(v, bool) or isinstance(v, float)):
        pass
      else:
        nv = str(v)
      yield (k, nv)

  def push(self, data={}, **kwargs):
    if not self.enable: return
    kwargs.update(data)
    final = dict(self.normalize(kwargs))
    print(final, 'FUU ', kwargs)
    cur = dict(measurement='engine_main', fields=final)
    assert self.client.write_points([cur], tags=self.tags, retention_policy=None)


def clamp(x, l, h):
  if x < l: return l
  if x > h: return h
  return x

def linearize(t, tl, vl):
  if vl[0] == vl[1]: return vl[0]
  t = max(min(t, tl[1]), tl[0])
  return vl[0] + (t - tl[0]) * (vl[1] - vl[0]) / (tl[1] - tl[0])


def sigmoid(t, vl):
  return linearize(1 / (1 + math.exp(-t)), (0, 1), vl)


def init(ctx):
  global g_ctx
  g_ctx = ctx
  swig.setup(flags.other_args[1:])
  if ctx.noengine:
    return
  ctx.engine = swig.opa_engine_swig

  class UpdaterFunc(swig.opa_engine_swig.Updater_Swig):

    def __init__(self, x):
      super().__init__()
      self.x = x

    def update(self):
      self.x()

  class InputCallerFunc(swig.opa_engine_swig.InputCaller_Swig):

    def __init__(self, x):
      super().__init__()
      self.x = x

    def update(self, arg):
      return self.x(arg)

  ctx.UpdaterFunc = UpdaterFunc
  ctx.InputCallerFunc = InputCallerFunc

  ctx.dumper = InfluxDbRun(ctx.dump_stats)


def create_db(ctx):
  ctx.dumper.client.create_database(ctx.influx_db)


def test(ctx):
  print(ctx.engine.quat_look_at_safe((1, 0, 0), (0, 0, 1)))
  #game = ctx.engine.Game()
  #game.init()
  #game.updater().add(engine.StartLoop, define_cb(lambda: print('jambon')).__disown__())
  #game.run()


class EarthEnv:

  def __init__(self, alt):
    self.alt = alt / 1000
    self.radius = 6.371

  def latlng_to_coord(self, lat, lng, nalt=None):
    if nalt is None: nalt = self.alt
    lat = lat / 360 * 2 * math.pi
    lng = lng / 360 * 2 * math.pi
    z = math.sin(lat)

    x = math.cos(lat) * math.cos(lng)
    y = math.cos(lat) * math.sin(lng)
    res = np.array([x, y, z])
    res = res * (self.radius + nalt)
    return res.tolist()


def nodes_to_jobs(lst):
  return list([dict(x=a.coord[0], y=a.coord[1], z=a.depth) for a in lst])


def limit_fps():
  try:
    game = g_ctx.game
    time.sleep(1 / flags.limit_fps)
    mouse_pos = game.input_handler().mouse().cur()
    cpos = game.camera().obj().pos_obj()
    ray = game.camera().get_click_dir(mouse_pos)
    fps = game.data_collector().fps()

    data = dict(
        mouse=mouse_pos,
        raypos=ray.pos,
        raydir=ray.dir,
        front=cpos.get_front(),
        cpos=cpos.pos(),
        fps=fps,)
    g_ctx.dumper.push(data)
    pprint.pprint(data)
  except:
    tb.print_exc()
    raise


def earth(ctx):
  tree = ctx.engine.AutoKDTree()
  tree.params.precision = ctx.precision
  tree.params.split_tr_norm = ctx.split_tr_norm
  ix = 0

  game = ctx.engine.Game()
  game.init()
  ctx.game = game

  desc = ctx.engine.init_earth_env()
  viewport = game.camera().viewport()
  desc.camera.update_viewport(viewport[0], viewport[1])
  ctx.obj = None
  ctx.iter = 0
  ctx.iter2 = 0

  tm = ctx.engine.TileResourceManager(ctx.resources_dir, ctx.mosaic)
  ctx.iter_speed = 0

  ctx.refresh_render = True
  res = []
  ephem = get_ephem('301', '2017-09-09', '2017-10-10', '1h')

  pc = ctx.engine.PointCloudSceneObj()
  pc.pc().init(1)
  p = pc.pc().get_pt(0)
  p.col = (1, 0, 0)
  p.pos = (0,0,0)
  pc.pc().refresh()
  pc.init(game)
  game.register_obj(pc, game.scene())
  pc.pc().set_r_fact(1)


  moonlat = ephem['ObsEclLat'].data
  moonlon = ephem['ObsEclLon'].data
  xt = np.linspace(0, 1, len(ephem))
  latinter = interp1d(xt, moonlat)
  loninter = interp1d(xt, moonlon)

  def set_moon_pos(x):
    x = clamp(x, 0, 1)
    nlat = latinter(x)
    nlon = loninter(x)
    coord = EarthEnv(alt=5000).latlng_to_coord(nlat, nlon)
    print(nlat, nlon, coord)
    pc.pos_obj().set_pos(coord)


  def setup_cam():

    print('ON SETUP')
    try:

      lng = 4.85000
      lng = 0
      lat0 = 45.7
      lat1 = 45.8
      lat0 = lat1 = 0
      alt0 = 38400
      alt1 = 38400
      pitch0 = math.pi / 3
      pitch1 = -math.pi / 3
      #pitch0 = pitch1 = -math.pi / 3
      max_iter = 200



      ctx.iter2 += ctx.iter_speed
      t = linearize(ctx.iter2, (-max_iter, max_iter), (0, 1))
      set_moon_pos(t)
      print(t, ctx.iter2)
      if not ctx.refresh_render: return



      alt = linearize(ctx.iter, (-max_iter, max_iter), (alt0, alt1))
      lat = linearize(ctx.iter, (-max_iter, max_iter), (lat0, lat1))
      pitch = linearize(ctx.iter, (-max_iter, max_iter), (pitch0, pitch1))

      ctx.iter += ctx.iter_speed

      env = EarthEnv(alt=alt)
      pos = env.latlng_to_coord(lat, lng)
      nrot = ctx.engine.quat_look_at_safe((-np.array(pos)).tolist(), (0, 0, 1))

      rot2 = ctx.engine.quat_from_euler((pitch, 0, 0))

      desc.camera.set_near_field(1e-1)
      desc.camera.set_far_field(1e3)
      desc.cam_pos.set_pos(pos)
      desc.cam_pos.set_rot(nrot)
      desc.cam_pos.rot_rel(rot2)
      desc.setup()

      game.camera().set_params(desc.camera.params())
      game.camera().obj().set_pos_obj(desc.cam_pos)
      game.camera().setup()
      tmp = game.camera().get_view_frustrum_points_world(desc.cam_pos)
      print(tmp)

      tmp = game.camera().hps_cam_space(desc.cam_pos)
      for hp in tmp:
        print(hp.plane.v, hp.plane.dir)

      nodes = tree.collect(desc)
      res.append(len(nodes))
      ctx.dumper.push(iter=ctx.iter, alt=alt, lat=lat, pitch=pitch, nnodes=len(nodes))
      print('GOT ', res)
      if not ctx.mosaic:
        print('COLLECTING ', lat, pos)
        with TileGetter(ctx.resources_dir) as tg:
          tg.do_jobs(nodes_to_jobs(nodes), wait=True)

      if ctx.obj is not None:
        game.remove_obj(ctx.obj)

      obj = ctx.engine.create_earth_obj(tm, nodes, game, True)
      obj.set_pos_obj(desc.obj_pos)
      #print('GOT OBJ', obj.id())
      ctx.obj = obj
    except:
      tb.print_exc()
      raise

  def handle_debug_req():
    mouse_pos = game.input_handler().mouse().cur()
    print('MOUSE POS >> ', mouse_pos)
    ray = game.camera().get_click_dir(mouse_pos)
    found, res = game.scene().isec().find_intersection(ray, True)
    if not found: return
    print(res.obj.debug_const())
    print(res.debug)
    print(res.pos, res.dist)
    print(np.array(ray.dir) * res.dist + np.array(ray.pos))
    res.obj.status().hide = 1 ^ res.obj.status().hide

  ctx.debug_cnt = 0

  def input_handler(x):
    try:
      # should do something with directorin
      x = game.input_handler()
      if not x.action().pressed(): return False
      if x.input().type() == ctx.engine.Key_LEFT: ctx.iter -= 1
      elif x.input().type() == ctx.engine.Key_RIGHT: ctx.iter += 1
      elif x.input().type() == ctx.engine.Key_DOWN: ctx.iter_speed -= 1
      elif x.input().type() == ctx.engine.Key_UP: ctx.iter_speed += 1
      elif x.input().type() == ctx.engine.Key_SPACE:
        open(f'/tmp/debug_{ctx.debug_cnt}.out', 'wb').write(tree.ctx.str())
        print(ctx.debug_cnt)
        ctx.debug_cnt += 1
      elif x.input().type() == ctx.engine.Key_W:
        return False
      elif x.input().type() == ctx.engine.Key_A:
        ctx.refresh_render = not ctx.refresh_render
      elif x.input().type() == ctx.engine.Key_D:
        handle_debug_req()
      else:
        return ctx.refresh_render
    except:
      tb.print_exc()
      raise
    return True

  f2 = ctx.UpdaterFunc(setup_cam)
  icf1 = ctx.InputCallerFunc(input_handler)
  f1 = ctx.UpdaterFunc(limit_fps)
  game.updater().add(ctx.engine.StartLoop, f1.to_updater())
  game.input_handler().cb_caller().add(ctx.engine.InputMark_Spec - 1, icf1.to_input_caller())
  game.updater().add(ctx.engine.StartLoop, f2.to_updater())

  if 0:
    setup_cam()
    setup_cam()
    return

  game.run()


def test_inter(ctx):
  m = swig.opa_math_game_swig

  p0 = (0.853582, -0.446406, -0.268551)
  p1 = (0.852168, -0.450568, -0.266081)
  p2 = (0.853156, -0.449313, -0.265034)

  return

  p1 = (1.000000, 0.000000, 0.000000)
  p2 = (0.995185, 0.098017, 0.000000)
  p3 = (0.990408, 0.097547, -0.097861)
  raypos = (1.313923, 0.000000, 0.000000)
  raydir = (-0.975810, -0.171109, 0.136079)

  p1 = (0.995185, -0.098017, 0.000000)
  p2 = (1.000000, 0.000000, 0.000000)
  p3 = (0.995200, 0.000000, -0.097861)
  raypos = (1.313923, 0.000000, 0.000000)
  raydir = (-0.990243, 0.038283, -0.133990)

  res = m.line_parallelogram_intersection(raypos, raydir, p1, p2, p3, True)
  print(res)
  return
  tr0 = ((0.774011, -0.533172), (-0.026693, -1.166229), (0.748132, -2.092785),)
  tr1 = ((1.403585, -0.173061), (-0.173060, -1.403585), (-1.403584, 0.173060),)
  for p in tr0:
    print(m.inside_polygon_v2_tr2d(tr1, p))
  for p in tr1:
    print(m.inside_polygon_v2_tr2d(tr0, p))

  print(m.is_counterclockwise_tr2d(tr0))
  print(m.is_counterclockwise_tr2d(tr1))

  print(m.tr_tr_intersection(tr0, tr1))


def test_pc(ctx):
  m = swig.opa_math_game_swig

  game = ctx.engine.Game()
  game.init()
  ctx.game = game

  f1 = ctx.UpdaterFunc(limit_fps)
  game.updater().add(ctx.engine.StartLoop, f1.to_updater())

  pc = ctx.engine.PointCloudSceneObj()
  in_fc = m.FaceCollection(True)
  in_fc.load_stl(ctx.stl_file).whiten()
  mesh = in_fc.to_mesh()
  vlist = mesh.list_vertices()
  n = len(vlist)
  #n = 10

  pc.pc().init(n)
  for i, v in enumerate(vlist):
    p = pc.pc().get_pt(i)
    p.col = (1, 0, 0)
    p.pos = mesh.get_vertex_attr(v).pos

  #for i in range(n):
  #  pc.pc().get_pt(i).col = (1,0,0)
  #  ang = 2*math.pi * i / n
  #  pc.pc().get_pt(i).pos = (math.cos(ang), math.sin(ang), 0)

  pc.pc().refresh()
  pc.init(game)

  ctx.v = 1

  def input_handler(x):
    try:
      # should do something with directorin
      x = game.input_handler()
      if not x.action().pressed(): return False
      if x.input().type() == ctx.engine.Key_LEFT: ctx.v *= 0.9
      elif x.input().type() == ctx.engine.Key_RIGHT: ctx.v /= 0.9
      else: return False

      pc.pc().set_r_fact(ctx.v)
    except:
      tb.print_exc()
      raise
    return True

  icf1 = ctx.InputCallerFunc(input_handler)
  game.input_handler().cb_caller().add(ctx.engine.InputMark_Spec - 1, icf1.to_input_caller())

  game.camera().params().near_field = 0.01

  game.camera().obj().pos_obj().set_pos((0, 0, 10))
  game.camera().obj().pos_obj().set_rot(ctx.engine.quat_look_at_safe((0, 0, -1), (1, 0, 0)))
  game.camera().setup()
  game.register_obj(pc, game.scene())

  game.run()


@Cachable.cachedf(fileless=False)
def _get_ephem(obj, start, stop, step):
  return Horizons(
      id=obj,
      location=500,
      epochs=dict(start=start, stop=stop, step=step),
      id_type='majorbody',).ephemerides().to_pandas().to_msgpack()

def get_ephem(obj, start, stop, step):
  return pd.read_msgpack(_get_ephem(obj, start, stop, step))






def test_ephem(ctx):
  x = get_ephem('301', '2017-09-09', '2017-10-10', '1h')
  print(x['datetime_str'])
  print(x['ObsEclLat'])
  print(x['ObsEclLon'])


def main():
  ctx = Attributize()
  ActionHandler.Run(ctx)


app()
