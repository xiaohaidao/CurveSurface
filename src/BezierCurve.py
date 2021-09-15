#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

class Bezier:
    def __init__(self, ctrl_points = [], weights = [], has_w = False):
        if has_w:
            points = np.array(ctrl_points).astype(float)
            self.__ctrl_points = points[:,:-1].tolist()
            self.__weights = points[:,-1].tolist()
        else:
            self.__ctrl_points = ctrl_points
            self.__weights = weights

    @property
    def ctrl_points(self):
        return self.__ctrl_points

    @ctrl_points.setter
    def ctrl_points(self, values):
        self.__ctrl_points = values

    @property
    def weights(self):
        return self.__weights

    @weights.setter
    def weights(self, values):
        self.__weights = values

    @property
    def degree(self) -> int:
        return len(ctrl_points) - 1

    def IsRational(self) -> bool:
        return len(self.weights) > 0

    ##
    # @param u the value range is [0,1]
    def UValue(self, u):
        if self.IsRational():
            ctrl = np.append(np.array(self.ctrl_points),
                    np.array(self.weights).reshape(-1,1), axis = 1)

            point = self.__Basis(ctrl, u)
            point = np.array(point)
            return np.divide(point[:-1], point[-1], out=np.zeros_like(point[:-1]),
                    where = point[-1] != 0).tolist()
        else:
            return self.__Basis(self.ctrl_points, u)

    def Curve(self, split_num = 100):
        points = []
        for i in range(split_num):
            points.append(self.UValue(i / split_num))

        return points

    def D0(self):
        pass

    def D1(self):
        pass

    def D2(self):
        pass

    def DN(self):
        pass

    ##
    # @brief __Basis get the bezier curve, the algorithm see deCasteljaul
    #
    # @param ctlpt the control point, type is [[x,y,z],...]
    # @param u the range is [0,1], type is float
    #
    # @return  curve point by u, type is [[x,y,z],...]
    def __Basis(self, ctlpt, u):
        q = np.array(ctlpt, dtype = float)
        num = len(ctlpt)
        for i in range(1, num):
            for k in range(num - i):
                q[k]  = (1 - u) * q[k] + u * q[k + 1]

        return q[0].tolist()

def __BezierShow(ax, list_point):
    points = np.array(Bezier(list_point, has_w = True).Curve())
    ax.plot(points[:, 0], points[:, 1])
    points = np.array(Bezier(list_point, has_w = False).Curve())
    ax.plot(points[:, 0], points[:, 1])

def __BezierShowTest():
    arra = [[1, 0, 1],
            [1, 1, 1],
            [0, 2, 2]]

    ctl_points1 = [
            [0, 0, 1],
            [1, 1, 1],
            [2, 0, 1]]

    ctl_points2 = [
            [0, 0, 3],
            [1, 1, 1],
            [2, 0, 3]]

    ax = plt.figure().add_subplot() #(projection='3d')

    __BezierShow(ax, arra)
    __BezierShow(ax, ctl_points1)
    __BezierShow(ax, ctl_points2)

    plt.show()

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    __BezierShowTest()
