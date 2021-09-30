#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .BSplineCurve import BSpline

def test_BezierShowTest_test():
    ctl_points1 = [
            [0, 0, 1],
            [1, 2, 3],
            [2, 2, 1],
            [4, 0, 1]]

    knots1 = [0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
    bs = BSpline([], knots1)
    bs.degree = 2
    u1 = 5  / 2

    assert bs.BasisFuns(u1) == [1 / 8, 6 / 8, 1 / 8]
    assert bs.DersBasisFuns(u1) == [[1 / 8, 6 / 8, 1 / 8],
            [-0.5, 0.0, 0.5],
            [1.0, -2.0, 1.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]]
    assert bs.OneBasisFuns(u1) == 0.125
    assert bs.DersOneBasisFuns(u1) == [0.125, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]

