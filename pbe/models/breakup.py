# -*- coding: utf-8 -*-
"""Breakup models

This module contains breakup models objects for population balance equation

"""


# Importing dependencies
from numpy import arange, sqrt, exp, pi


class breakupModels():
    def __init__(self) -> None:
        pass

    @staticmethod
    def beta(v1: float, v2: float) -> float:
        return 2.4 / v2 * exp(-4.5 * (2 * v1 - v2)**2 / (v2**2))
    

print(1)