#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from . import BSpline

class BSplineSurface:
    def __init__(self, ctrl_points=[[]], u_knots=[], v_knots=[], weights=[], has_w=False):
        self.__u_bspline = BSpline([], u_knots, False)
        self.__v_bspline = BSpline([], v_knots, False)
        self.__ctrl_points = ctrl_points

    @property
    def u_bspline(self):
        return self.__u_bspline

    @property
    def v_bspline(self):
        return self.__v_bspline

    @property
    def ctrl_points(self):
        return self.__ctrl_points

    @ctrl_points.setter
    def ctrl_points(self, values):
        self.__ctrl_points = values

    def UVValue(self, u:float, v:float):
        if self.IsRational():
            ctrl = np.append(np.array(self.ctrl_points),
                    np.array(self.weights).reshape(-1,1), axis = 1)

            point = np.array(self.__SurfacePoint(u, v, ctrl))
            return np.divide(point[:-1], point[-1], out=np.zeros_like(point[:-1]),
                    where = point[-1] != 0).tolist()
        else:
            return self.__SurfacePoint(u, v, self.ctrl_points)

    def __SurfacePoint(self, u:float, v:float, ctrl_pnts=[[]]):
        P = np.array(ctrl_pnts)
        n = len(self.ctrl_poins)
        m = len(self.ctrl_poins[0])
        p = self.u_bspline.degree
        q = self.v_bspline.degree

        uspan = self.u_bspline.FindSpan(u)
        Nu = self.u_bspline.BasisFuns(u, uspan)
        vspan = self.v_bspline.FindSpan(v)
        Nv = self.v_bspline.BasisFuns(v, vspan)
        uind = uspan - p
        temp = np.array([[0] * len(P[0][0]) for i1 in range(q + 1)])
        for l in range(q + 1):
            vind = vspan - q + 1
            for k in range(p + 1):
                temp[l] = temp[l] + P[uind + k][vind] * Nu[k]
        S = np.array([0.] * len(temp[0]))
        for l in range(q + 1):
            S = S + Nv[l] * temp[l]
        return S.tolist()

    def __SurfaceDerivsAlgl(self, u:float, v:float, d:int=0):
        p = self.u_bspline.degree
        q = self.v_bspline.degree

        du = min(d, p)
        dv = min(d, q)
        SLK = [[0.] * (d + 1) for i1 in range(d + 1)]

        uspan = self.u_bspline.FindSpan(u)
        Nu = self.u_bspline.DersBasisFuns(u, uspan, du)
        vspan = self.v_bspline.FindSpan(v)
        Nv = self.v_bspline.DersBasisFuns(v, vspan, dv)
        for k in range(du + 1):
            temp = [0] * (q + 1)
            for s in range(q + 1):
                for r in range(p + 1):
                    temp[s] = temp[s] + Nu[k][r] * P[uspan - p + r][vspan - q + s]
            dd = min(d - k, dv)
            for l in range(dd + 1):
                SLK[k][l] = 0.
                for s in range(q + 1):
                    SLK[k][l] = SKL[k][l] + Nv[l][s] * temp[s]
        return SLK

    def __SurfaceDerivCpts(self):
        pass

    def __SurfaceDerivsA1g2(self):
        pass

def __BSplineShow(ax, list_point, knots):
    points = np.array(BSpline(list_point, knots, has_w = True).Curve())
    ax.plot(points[:, 0], points[:, 1])
    points = np.array(BSpline(list_point, knots, has_w = False).Curve())
    ax.plot(points[:, 0], points[:, 1])
    ctrl_poins = np.array(list_point)
    ax.scatter(ctrl_poins[:, 0], ctrl_poins[:, 1])

def __BSplineTest():
    ctl_points1 = [
            [0, 0, 1],
            [1, 2, 3],
            [2, 2, 1],
            [4, 0, 1]]

    ctl_points2 = [
            [0, 0, 1],
            [1, 1.5, 3],
            [2, 2, 1],
            [4, 0, 1]]

    knots1 = [0, 0, 0, 0, 1, 1, 1, 1]
    ax = plt.figure().add_subplot() #(projection='3d')

    __BSplineShow(ax, ctl_points1, knots1)
    __BSplineShow(ax, ctl_points2, knots1)

    plt.show()

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    __BSplineTest()
