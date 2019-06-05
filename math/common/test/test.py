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
  clist = CmdsList().add(test)
  ActionHandler.Prepare(parser, clist.lst, global_action=1)

def test_rs(ctx):
  Z.np.random.seed(0)

  gf_256 = m.GF_q_u32(m.cvar.GF2, 8)
  rs = m.RSCode_u32(gf_256, 255, 223);
  data = tuple(Z.np.random.randint(2, size=223*8).tolist())
  for nx in range(0, 40):
    print('Diffs ', nx)
    m1 =m.pack_vector_gfq_u32(gf_256, m.v_u32(data))
    c = rs.encode(m1)
    cx =m.unpack_vector_gfq_u32(gf_256, c)
    cx = list(cx)
    for i in range(nx):
      for j in range(1):
        cx[(32+i)*8+j] ^= 1
    c2 =m.pack_vector_gfq_u32(gf_256, cx)
    mm = rs.decode(c2)
    print(len(mm))
    mx = m.unpack_vector_gfq_u32(gf_256, mm)
    print(len(mx))
    if mx != data:
      print('REAL DIFF AT diff ', nx)
      break



def test(ctx):
  gf_256 = m.GF_q_u32(m.cvar.GF2, 8)
  #bch = m.BCHCode_P_u32.getNew(m.cvar.GF2, 255, 1, 11);
  #gf_256 = m.GF_q_u32(m.cvar.GF2, 8)
  bch = m.BCHCode_u32.getNew(m.cvar.GF2, 255, 1, 11);
  cc = m.CyclicCode_u32(m.cvar.GF2, bch.get_g(), 255)
  a=  cc.encode([1]*163)
  print(cc.decode(a))

def check_rs(ctx):
  n = 255
  k = 223
  Z.np.random.seed(1)
  import unireedsolomon as urs
  pr_gf2= m.cvar.PR_GF2
  px =m.cvar.PR_GF2.import_base(0x11b)
  rs2 = urs.RSCoder(n, k, generator=3)

  gf_256 = m.GF_q_u32(m.cvar.GF2, 8)
  print(pr_gf2.export_base(gf_256.getModPoly()))
  return


  x = gf_256._import(pr_gf2.x())



  g3 = m.cvar.PR_GF2.import_base(3)
  gf_256.set_prim_elem(g3)
  pr_gf26 = m.PolyRing_Poly_u32(gf_256)

  rs = m.RSCode_u32(gf_256, n, k);
  print(rs.g())
  print(rs.w().str())
  print(rs.g().str())
  for x in rs.g().to_vec():
    print(pr_gf2.export_base(x))

  print(rs2.g[k])


  lsb=1 # fixed, need to be consistent with pack_vector_
  for ntest in range(10):
    data = tuple(Z.np.random.randint(2, size=k*8).tolist())
    m1 =m.pack_vector_gfq_u32(gf_256, m.v_u32(data))
    data_byte = bytes(Z.Format(data).bin2byte(lsb=lsb).v)

    v = bytes([0] * (n-k)) + data_byte
    tb = []
    for x in v: tb.append(gf_256.import_base(x))
    print(len(tb))
    mprime = m.v_poly_u32(tb)
    mprime = pr_gf26._import(mprime)
    b = pr_gf26.mod(mprime, rs.g())
    u = pr_gf26.add(mprime, b)
    print(pr_gf2.export_base(u.get_safe(0)), pr_gf2.export_base(u.get_safe(32)), pr_gf2.export_base(u.get_safe(254)))

    #lsb ^= 1
    c2_byte = rs2.encode(data_byte[::-1])
    c2_byte = bytes(map(ord, c2_byte))
    print(c2_byte)
    print(data_byte[0], data_byte[-1])


    m1 =m.pack_vector_gfq_u32(gf_256, m.v_u32(data))
    c = rs.encode(m1)
    print(pr_gf2.export_base(c[0]), pr_gf2.export_base(c[-1]))
    c1 =list(m.unpack_vector_gfq_u32(gf_256, c))
    c1_byte = bytes(Z.Format(c1).bin2byte(lsb=lsb).v)
    #c1_byte = c1_byte[n-k:] + c1_byte[:n-k]

    print(c2_byte[::-1]==c1_byte)
    print(data_byte[0], data_byte[-1])
    return

    cx =m.unpack_vector_gfq_u32(gf_256, c)
    cx = list(cx)
    c2 =m.pack_vector_gfq_u32(gf_256, cx)
    mm = rs.decode(c2)




def main():
  global m
  m = Z.swig.opa_math_common_swig
  ctx = Attributize()
  ActionHandler.Run(ctx)


app()
