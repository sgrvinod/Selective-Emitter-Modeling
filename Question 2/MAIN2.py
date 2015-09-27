#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy as np
import scipy as sci
import math as m
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import interpolate

#import python files for multilayer selective emitters and single layer emitters respectively
import selectiveemitter as selem
import singlelayer as sl

#call functions with arguments specifying emitter material
selem.generate_sel_em_emissivity("SiO2data.txt")

selem.generate_sel_em_emissivity("BNdata.txt")

sl.generate_singlelayerem("SiO2data.txt")

sl.generate_singlelayerem("BNdata.txt")

#note: program takes a 