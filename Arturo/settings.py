#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

try:
    type_image = sys.argv[1]
except IndexError:
    print "Il manque un argument, remplacé par fabry-pérot"
    type_image = "fp0"
