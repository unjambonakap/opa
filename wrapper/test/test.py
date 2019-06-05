#!/usr/bin/env python

from chdrft.cmds import CmdsList
from chdrft.main import app
from chdrft.utils.cmdify import ActionHandler
from chdrft.utils.misc import Attributize
import chdrft.utils.misc as cmisc
import glog
import chdrft.utils.Z as Z

global flags, cache
flags = None
cache = None


m= None

def args(parser):
  clist = CmdsList()
  ActionHandler.Prepare(parser, clist.lst, global_action=1)

def check_rs(ctx):
  n = 255
  k = 223
  #u = Z.swig_unsafe.opa_wrapper_swig.CReedSolomon(8,16,112,11,0,4,0,1)
  gamma_pw, chk_pw = 11, 112
  u = Z.swig_unsafe.opa_wrapper_swig.CReedSolomon(8,16,gamma_pw,chk_pw,0,1,0,1)
  Z.np.random.seed(1)
  import unireedsolomon as urs

  pr_gf2= m.cvar.PR_GF2
  px =m.cvar.PR_GF2.import_base(0x187)
  gf_256 = m.GF_q_u32(m.cvar.GF2, px)
  assert gf_256.is_prim_elem(gf_256.theta())
  gf_256.set_prim_elem(gf_256.theta())
  rs = m.RSCode_u32(gf_256, n, k);
  gamma = gf_256.faste(gf_256.theta(), gamma_pw)
  rs.advanced_init(gf_256.theta(), chk_pw)
  pr_gf26 = m.PolyRing_Poly_u32(gf_256)


  rs2 = urs.RSCoder(n, k, generator=2, prim=0x187)
  data = Z.os.urandom(k)
  data = b'a'*k
  tb = []
  for x in data: tb.append(gf_256.import_base(x))
  c = bytes(map(ord, rs2.encode(data[::-1])))
  cc = rs.encode(m.v_poly_u32(tb))
  ccx = bytearray(map(gf_256.export_base, list(cc)))
  #ccx[0] = 0
  #ccx[1] = 1
  ccx[-1] = 1
  ccx[-2] = 1
  ccx[-3] = 1
  cc = m.v_poly_u32(list(map(gf_256.import_base, ccx)))

  m2 = bytes(map(gf_256.export_base, rs.decode(cc)))

  res, mx = u.Decode(ccx[::-1])
  print(res)
  print(u.CorrectableErrorsInFrame())
  print(u.UncorrectableErrorsInFrame())
  print(mx)
  print(m2)



  print(pr_gf2.export_base(gf_256.getModPoly()))
  return


  print(rs.g())
  print(rs.w().str())
  print(rs.g().str())
  for x in rs.g().to_vec():
    print(pr_gf2.export_base(x))





def main():
  global m
  m = Z.swig.opa_math_common_swig
  ctx = Attributize()
  ActionHandler.Run(ctx)


app()
