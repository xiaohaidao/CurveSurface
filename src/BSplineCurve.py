#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

class BSpline:
    def __init__(self, ctrl_points=[], knots=[], weights=[], has_w=False):
        if has_w:
            points = np.array(ctrl_points).astype(float)
            self.__ctrl_points = points[:,:-1].tolist()
            self.__weights = points[:,-1].tolist()
        else:
            self.__ctrl_points = ctrl_points
            self.__weights = weights
        self.__knots = knots
        self.__p = 3


    @property
    def ctrl_points(self):
        return self.__ctrl_points

    @ctrl_points.setter
    def ctrl_points(self, values):
        self.__ctrl_points = values

    @property
    def knots(self):
        return self.__knots

    @knots.setter
    def knots(self, values):
        self.__knots = values

    @property
    def weights(self):
        return self.__weights

    @weights.setter
    def weights(self, values):
        self.__weights = values

    @property
    def degree(self) -> int:
        return self.__p

    @degree.setter
    def degree(self, degree:int):
        self.__p = degree

    def IsRational(self) -> bool:
        return len(self.weights) > 0

    ##
    # @param u the value range is [0,1]
    def UValue(self, u:float):
        if self.IsRational():
            ctrl = np.append(np.array(self.ctrl_points),
                    np.array(self.weights).reshape(-1,1), axis = 1)

            point = np.array(self.__CurvePoint(u, ctrl))
            return np.divide(point[:-1], point[-1], out=np.zeros_like(point[:-1]),
                    where = point[-1] != 0).tolist()
        else:
            return self.__CurvePoint(u, self.ctrl_points)

    def Curve(self, split_num = 100):
        points = []
        for i in range(split_num):
            points.append(self.UValue(i / split_num))
        return points

    def CurveDer(self, k:int, split_num = 100):
        points = []
        for i in range(split_num):
            points.append(self.DN(i / split_num, k))
        return points

    def D0(self, u:float):
        return self.DN(u, 0)

    def D1(self, u:float):
        return self.DN(u, 1)

    def D2(self, u:float):
        return self.DN(u, 2)

    def DN(self, u:float, k:int):
        if k > self.degree:
            return 0
        if self.IsRational():
            ctrl = np.append(np.array(self.ctrl_points),
                    np.array(self.weights).reshape(-1,1), axis = 1)

            point = np.array(self.__CurveDerivsAlg1(u, self.ctrl_points, k)[k])
            return np.divide(point[:-1], point[-1], out=np.zeros_like(point[:-1]),
                    where = point[-1] != 0).tolist()
        else:
            return self.__CurveDerivsAlg1(u, self.ctrl_points, k)[k]

    def InsertKnot(self, u:float):
        pass

    def __CurveKnotIns(self):
        pass

    def __SurfaceKnotIns(self):
        pass

    def __CurvePntByCornerCut(self):
        pass

    def __RefineKnotVectCurve(self):
        pass

    def __RefineKnotVectSurface(self):
        pass

    def __DecomposeCurve(self):
        pass

    def __DecomposeSurface(self):
        pass

    def __RemoveCurveKnot(self):
        pass

    def __DegreeElevateCurve(self):
        pass

    def __DegreeElevateSurface(self):
        pass

    def __DegreeReduceCurve(self):
        pass

    def __DegreeReduceSurface(self):
        pass

    def RemoveKnot(self):
        pass

    def __GenereteKnots(self):
        m = len(self.ctrl_points) + self.degree + 1
        self.knots = [0] * (m + 1)
        # split_count = self.degree + 1
        p = self.degree
        for i in range(p + 1, m - p):
            self.knots[i] = p + 1 - i
        for i in range(m - p, m + 1):
            self.knots[i] = 1
        pass

    # def ToBezier(self) -> [Bezier]:
        # pass

    # Get the range u $\in$ [U[index], U[index + 1])
    def FindSpan(self, u:float) -> int:
        U = self.knots
        # m = len(self.ctrl_points) + self.degree
        # n = m - self.degree - 1
        n = (len(U) - 1) - self.degree - 1
        if u == U[n + 1]: return n
        low = self.degree
        high = n + 1
        mid = int((low + high) / 2)
        while u < U[mid] or u >= U[mid + 1]:
            if u < U[mid]:
                high = mid
            else:
                low = mid
            mid = int((low + high) / 2)
        return mid

    def BasisFuns(self, u:float, span:int=-1) -> [float]:
        p = self.degree
        U = self.knots
        i = self.FindSpan(u) if span < 0 else span
        N = [0.] * (p + 1)
        left = N.copy()
        right = left.copy()
        N[0] = 1.
        for j in range(1, p + 1):
            left[j] = u - U[i + 1 - j]
            right[j] = U[i + j] - u
            saved = 0.
            for r in range(j):
                temp = N[r] / (right[r + 1] + left[j - r])
                N[r] = saved + right[r + 1] * temp
                saved = left[j - r] * temp
            N[j] = saved
        return N

    def DersBasisFuns(self, u:float, span:int=-1) -> [[float]]:
        p = self.degree
        U = self.knots
        n = (len(U) - 1) - self.degree - 1
        i = self.FindSpan(u) if span < 0 else span
        ndu = [[0.] * (p + 1) for i1 in range(p + 1)]
        ndu[0][0] = 1.

        left = [0.] * (p + 1)
        right = left.copy()
        for j in range(1, p + 1):
            left[j] = u - U[i + 1 - j]
            right[j] = U[i + j] - u
            saved = 0.
            for r in range(j):
                ndu[j][r] = right[r + 1] + left[j - r]
                temp = ndu[r][j - 1] / ndu[j][r]
                ndu[r][j] = saved + right[r + 1] * temp
                saved = left[j - r] * temp
            ndu[j][j] = saved
        ders = [[0.] * (p + 1) for i1 in range(n + 1)]
        for j in range(p + 1):
            ders[0][j] = ndu[j][p]
        for r in range(p + 1):
            s1, s2 = 0, 1
            a = [[0.] * (p + 1) for i1 in range(n + 1)]
            a[0][0] = 1.
            for k in range(1, n + 1):
                d = 0.
                rk, pk = (r - k), (p - k)
                if r >= k:
                    a[s2][0] = a[s1][0] / ndu[pk + 1][rk]
                    d = a[s2][0] * ndu[rk][pk]
                j1 = 1 if rk >= -1 else -rk
                j2 = (k - 1) if (r - 1 <= pk) else (p - r)
                for j in range(j1, j2 + 1):
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j]
                    d += (a[s2][j] * ndu[rk +j][pk])
                if r <= pk:
                    a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r]
                    d += (a[s2][k] * ndu[r][pk])
                ders[k][r] = d
                j, s1, s2 = s1, s2, j
        r = p
        for k in range(1, n + 1):
            for j in range(p + 1): ders[k][j] *= r
            r *= (p - k)
        return ders

    def OneBasisFuns(self, u:float, span:int=-1) -> float:
        p = self.degree
        U = self.knots
        n = (len(U) - 1) - self.degree - 1
        m = n + p + 1
        i = self.FindSpan(u) if span < 0 else span
        if (i == 0 and u == U[0]) or (i == n and u == U[m]):
            return 1.
        if u < U[i] or u >= U[i + p + 1]:
            return 0.
        N = [0.] * (p + 1)
        for j in range(p + 1):
            if u >= U[i + j] and u < U[i + j + 1]:
                N[j] = 1.
        for k in range(1, p + 1):
            saved = 0.
            if N[0] != 0.:
                saved = (u - U[i]) * N[0] / (U[i + k] - U[i])
            for j in range(p - k + 1):
                Uleft = U[i + j + 1]
                Uright = U[i + j + k + 1]
                if N[j + 1] == 0.:
                    N[j] = saved
                    saved = 0.
                else:
                    temp = N[j + 1] / (Uright - Uleft)
                    N[j] = saved + (Uright - u) * temp
                    saved = (u - Uleft) * temp
        return N[0]

    def DersOneBasisFuns(self, u:float, span:int=-1) -> [float]:
        p = self.degree
        U = self.knots
        n = (len(U) - 1) - self.degree - 1
        i = self.FindSpan(u) if span < 0 else span
        ders = [0.] * (n + 1)
        if u < U[i] or u >= U[i + p + 1]:
            return ders
        N = [[0.] * (p + 1) for i1 in range(p + 1)]
        for j in range(p + 1):
            if u >= U[i + j] and u < U[i + j + 1]:
                N[j][0] = 1.
        for k in range(1, p + 1):
            saved = 0.
            if N[0][k - 1] != 0.:
                saved = (u - U[i]) * N[0][k - 1] / (U[i + k] - U[i])
            for j in range(p - k + 1):
                Uleft = U[i + j + 1]
                Uright = U[i + j + k + 1]
                if N[j + 1][k - 1] == 0.:
                    N[j][k] = saved
                    saved = 0.
                else:
                    temp = N[j + 1][k - 1] / (Uright - Uleft) if Uright != Uleft else 0
                    N[j][k] = saved + (Uright - u) * temp
                    saved = (u - Uleft) * temp
        ders[0] = N[0][p]
        ND = [0] * (n + 1)
        for k in range(1, n + 1):
            if k > p:
                return ders
            for j in range(k + 1):
                ND[j] = N[j][p - k]
            for jj in range(1, k + 1):
                saved = 0.
                if ND[0] != 0.:
                    saved = ND[0] / (U[i + p - k + jj] - U[i])
                for j in range(k - jj + 1):
                    Uleft = U[i + j + 1]
                    Uright = U[i + j + p + jj + 1]
                    if ND[j + 1] == 0.:
                        ND[j] = (p - k + jj) * saved
                        saved = 0.
                    else:
                        temp = ND[j + 1] / (Uright - Uleft) if Uright != Uleft else 0
                        ND[j] = (p - k + jj) * (saved - temp)
                        saved = temp
            ders[k] = ND[0]
        return ders

    def __CurvePoint(self, u:float, ctrl=[]):
        p = self.degree
        span = self.FindSpan(u)
        N = np.array(self.BasisFuns(u, span))
        c = np.array([0.] * len(ctrl[0]))
        for i in range(p + 1):
            c = c + N[i] * np.array(ctrl[span - p + i])
        return c

    def __CurveDerivsAlg1(self, u:float, ctrl_pts=[], d:int=0):
        p = self.degree
        du = min(d, p)
        # CK = np.zeros((d + 1, len(ctrl_pts[0])))
        CK = np.array([[0.] * len(ctrl_pts[0]) for i1 in range(d + 1)])
        ctrl = np.array(ctrl_pts)
        span = self.FindSpan(u)
        ders = np.array(self.DersBasisFuns(u, span))
        for k in range(du + 1):
            CK[k] =0.
            for j in range(p + 1):
                CK[k] = CK[k] + ctrl[span - p + j] * ders[k][j]
        return CK.tolist()

    def __CurveDerivsAlg2(self, u:float, ctrl_pts=[], d:int=0):
        p = self.degree
        du = min(d, p)
        # CK = np.zeros((d + 1, len(ctrl_pts[0])))
        CK = np.array([[0.] * len(ctrl_pts[0]) for i1 in range(d + 1)])
        ctrl = np.array(ctrl_pts)
        span = self.FindSpan(u)
        N = self.BasisFuns(u, span)
        PK = self.__CurveDerivCpts(du, span - p, span)
        for k in range(du + 1):
            CK[k] =0.
            for j in range(p - k + 1):
                CK[k] = CK[k] + N[j][p - k] * PK[k][j]
        return CK.tolist()

    ##
    # @brief __CurveDerivCpts 计算曲线直到d阶(包括d，d<=p)的所有导曲线的控制点
    #
    # @return PK[k][i]返回k阶导曲线的第i个控制点，$k=0,1,\cdots,d$
    #    $i=r_1,\cdots,r_2-k$，r1=0且r2=n则计算所有的控制点
    def __CurveDerivCpts(self, d:int, r1:int, r2:int):
        P = self.ctrl_points
        p = self.degree
        n = (len(U) - 1) - self.degree - 1
        r = r2 - r1
        PK = [[0.] * (r + 1) for i1 in range(d + 1)]
        # PK = [[0.] * (n + 1) for i1 in range(p + 1)]
        for i in range(r + 1):
            PK[0][i] = P[r1 + i]
        for k in range(1, d + 1):
            tmp = p - k + 1
            for i in range(r - k + 1):
                PK[k][i] =  tmp * (PK[k - 1][i + 1] - PK[k - 1][i]) / (
                        U[r1 + i + p + 1] - U[r1 + i + k])
        return PK

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
