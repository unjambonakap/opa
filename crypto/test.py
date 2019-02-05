#!/usr/bin/env python

from chdrft.cmds import CmdsList
from chdrft.main import app
from chdrft.utils.cmdify import ActionHandler
from chdrft.utils.misc import Attributize
import glog
from chdrft.utils.swig import swig
import Crypto.Cipher.AES as AES
import Crypto.Util.Padding as Padding
import os
import random

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
  clist = CmdsList().add(t1).add(t2).add(test).add(test2).add(test_aes)
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

def test_aes(ctx):
  m=swig.opa_math_common_swig
  c = swig.opa_crypto_swig

  BS = 16
  for i in range(100):

    key = os.urandom(BS)
    a = os.urandom(random.randint(0, BS-1))
    a1_pad = Padding.pad(a, BS)
    a2_pad = c.pkcs7(a, BS)

    c1 = AES.new(key, AES.MODE_ECB).encrypt(a1_pad)
    c2 = c.Aes(key, True).encrypt_raw(a2_pad)


    print(c1)
    print(c2)
    assert c1 == c2

  for i in range(100):

    print(i)
    key = os.urandom(BS)
    a = os.urandom(random.randint(60, 100))
    a1_pad = Padding.pad(a, BS)
    a2_pad = c.pkcs7(a, BS)

    c1 = AES.new(key, AES.MODE_ECB).encrypt(a1_pad)
    c2 = c.Aes(key, True).encrypt_ecb(a2_pad)


    print(c1)
    print(c2)
    assert c1 == c2
    m2_pad = c.Aes(key, False).decrypt_ecb(c2)
    m2,ok = c.rpkcs7(m2_pad, BS)
    assert c1 == c2
    print(m2, a, m2_pad)
    assert m2 == a


  for i in range(100):
    iv = os.urandom(BS)

    print(i)
    key = os.urandom(BS)
    a = os.urandom(random.randint(60, 100))
    a1_pad = Padding.pad(a, BS)
    a2_pad = c.pkcs7(a, BS)

    c1 = AES.new(key, AES.MODE_CBC, iv=iv).encrypt(a1_pad)
    c2 = c.Aes(key, True).encrypt_cbc(a2_pad, iv)


    print(c1)
    print(c2)
    assert c1 == c2
    m2_pad = c.Aes(key, False).decrypt_cbc(c2, iv)
    m2,ok = c.rpkcs7(m2_pad, BS)
    assert c1 == c2
    print(m2, a, m2_pad)
    assert m2 == a



def main():
  ctx = Attributize()
  ActionHandler.Run(ctx)


app()
