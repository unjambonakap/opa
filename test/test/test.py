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


def args(parser):
  clist = CmdsList().add(test)
  ActionHandler.Prepare(parser, clist.lst, global_action=1)


def test(ctx):
  print('on test')

  a = Z.swig_unsafe.opa_test_swig
  x = Z.np.complex64(1 +2j)
  a.test1(x)


def main():
  ctx = Attributize()
  ActionHandler.Run(ctx)


app()
