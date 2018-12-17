#!/usr/bin/env python

from chdrft.cmds import CmdsList
from chdrft.main import app
from chdrft.utils.cmdify import ActionHandler
from chdrft.utils.misc import Attributize
import chdrft.utils.misc as cmisc
import glog
from chdrft.utils.swig import swig

global flags, cache
flags = None
cache = None


def args(parser):
  clist = CmdsList()
  ActionHandler.Prepare(parser, clist.lst, global_action=1)


def test(ctx):
  algo = swig.opa_algo_swig
  print(algo.kmp_matches('defabcdef', 'def'))


def main():
  ctx = Attributize()
  ActionHandler.Run(ctx)


app()
