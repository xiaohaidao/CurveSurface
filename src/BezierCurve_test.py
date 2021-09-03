#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .BezierCurve import Bezier

def test_BezierShowTest_test():
    arra = [[1, 0, 1],
            [1, 1, 1],
            [0, 2, 2]]


    assert Bezier(arra).UValue(0.5) == [0.75, 1.0, 1.25]
    assert Bezier(arra, has_w=True).UValue(0.5) == [0.6, 0.8]

