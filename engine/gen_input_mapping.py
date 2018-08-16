#!/usr/bin/env python

import argparse
from chdrft.utils.misc import cwdpath
import re

FLAGS = None


def build_enum_name(name, glfw_val, source):
  return f'{source}_{name},'
def build_macro(name, glfw_val, source):
  return f'OPA_INPUT_DECL({source}_{name}, {glfw_val}, {source})'


def extract_keys(content):
  return [
      (m.group(2), m.group(1), 'Key')
      for m in re.finditer('define (GLFW_KEY_(\w+))\s', content)
  ]


def extract_buttons(content):
  return [
      (m.group(2), m.group(1), 'Button')
      for m in re.finditer('define (GLFW_MOUSE_BUTTON_(\w+))\s', content)
  ]


def extract_joystick(content):
  return [
      (m.group(2), m.group(1), 'Joystick')
      for m in re.finditer('define (GLFW_JOYSTICK_(\w+))\s', content)
  ]


def go():
  content = open(FLAGS.glfw, 'r').read()
  data = []
  data += extract_keys(content)
  data += extract_buttons(content)
  data += extract_joystick(content)
  data += [('Enter', 0, 'Focus')]
  data += [('Mov', 0, 'Mouse')]
  data += [('Scroll', 0, 'Mouse')]
  if 1:
    print('''
namespace opa {
namespace engine {
enum InputEnum {
''')

    for x in data:
      print(build_enum_name(*x))

    print('''
    ENUM_END,
    };
    }
    }
    ''')
  else:
    for x in data:
      print(build_macro(*x))


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--glfw', type=cwdpath)

  global FLAGS
  FLAGS = parser.parse_args()
  go()


if __name__ == '__main__':
  main()
