#!/usr/bin/env python

import sys
import subprocess

filename = sys.argv[1]

f = subprocess.Popen(['tail','-n 1',filename], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
print f.stdout.readline()
