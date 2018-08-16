#!/usr/bin/env python

from chdrft.cmds import CmdsList
from chdrft.main import app
from chdrft.utils.cmdify import ActionHandler
from chdrft.utils.misc import Attributize
import glog
from chdrft.utils.swig import swig

global flags, cache
flags = None
cache = None


def t1():
  x = Cracker()
  checker = AsciiChecker()
  mode = ZipState()

  x.set_checker(checker)
  x.set_state(mode)
  x.load_dictionary(
      '/home/benoit/programmation/hack/data/crackstation-human-only.txt')
def t2():
  opa_init_swig(sys.argv)


def args(parser):
  clist = CmdsList().add(t1).add(t2).add(test).add(test2)
  ActionHandler.Prepare(parser, clist.lst)


def test(ctx):
  m=swig.opa_math_common_swig
  c = swig.opa_crypto_swig
  gf2=m.cvar.GF2
  pr_gf2=m.cvar.PR_GF2
  seq=[1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0]
  print(m.Poly_u32(pr_gf2, m.v_u32(seq)).to_vec())

  seq=m.v_u32(seq)
  res = c.LFSR_u32_fromSequence(gf2, seq)
  print(res.advance(m.bignum(-10)))
  for i in range(10): print(res.get_next())
  print(res.advance(m.bignum(-10)))
  for i in range(10): print(res.get_next())
  print(res.get_state())
  res.set_state(res.get_initial_state(seq))
  res.set_prev()
  res.set_prev()
  res.set_prev()
  print(res.get_poly().to_vec())
  print(res.get_state().to_vec())


def test2(ctx):
  m=swig.opa_math_common_swig
  c = swig.opa_crypto_swig
  gf2=m.cvar.GF2
  pr_gf2=m.cvar.PR_GF2
  seq=[0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0]

  print(m.Poly_u32(pr_gf2, m.v_u32(seq)).to_vec())

  seq=m.v_u32(seq)
  res = c.LFSR_u32_fromSequence(gf2, seq)

  print(res.get_poly().to_vec())
  print(res.get_state().to_vec())
def main():
  ctx = Attributize()
  ActionHandler.Run(ctx)


app()
