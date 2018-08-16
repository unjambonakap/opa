#!/usr/bin/env python

from chdrft.cmds import CmdsList
from chdrft.main import app
from chdrft.utils.cmdify import ActionHandler
from chdrft.utils.misc import Attributize
import glog
import sys
import opa_common_swig
from chdrft.utils.swig import swig
from chdrft.display.utils import DataOp
import numpy as np

global flags, cache
flags = None
cache = None


def args(parser):
  clist = CmdsList().add(test)
  ActionHandler.Prepare(parser, clist.lst)


def test(ctx):
  print(swig)
  print(swig.opa_or_swig)
  bounds  =swig.opa_or_swig.GridSearchBounds()
  tb1=swig.opa_common_swig.vd()
  bounds.bounds.append(np.arange(37.508, 49., 0.01))
  bounds.bounds.append(np.arange(0., 1., 0.01))

  data = np.fromfile('/tmp/data1', 'float32')

  if 0:
    solver = swig.opa_or_swig.DspBinaryMatcher(data)
    res = swig.opa_or_swig.DoGridSearch(bounds, solver.get_grid_search_func())
    print(res.score)
    print(res.state)
    print(solver.grid_search_func_impl([41.5, 0.1]))
    print(solver.grid_search_func_impl([41.5, 0.7]))
    print(solver.grid_search_func_impl([41.5, 0.2]))
  else:

    r2 = DataOp.ClockRecovery(data, 45, sps_dev=0.1)
    print(r2)

def main():
  ctx = Attributize()
  ActionHandler.Run(ctx)


app()
