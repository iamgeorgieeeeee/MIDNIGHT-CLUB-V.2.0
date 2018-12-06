from error_dialog import ErrorDialog

import sys
from PyQt4 import QtCore, QtGui, uic
import numpy as np
import math


global A, E, I, x1, x2, x3, y1, y2, L, cosTheta, sinTheta, AL2I, EIL3

qtCreatorFile1 = "problem.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile1)

class ProblemWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.calc_btn.clicked.connect(self.solve)
        self.clear_btn.clicked.connect(self.clear)
        #self.grasp_btn.clicked.connect(self.grasp)
        
    def clear(self):
    	self.elastic_modulus.clear()
        self.inertia.clear()
        self.area.clear()
        self.x1.clear()
        self.x2.clear()
        self.x3.clear()
        self.y1.clear()
        self.y2.clear()
        self.P.clear()
        self.w.clear()
        self.M.clear()
    
    def solve(self):
        print("entered in solve")
        E = self.elastic_modulus.text()
        I = self.inertia.text()
        A = self.area.text()
        x1 = self.x1.text()
        x2 = self.x2.text()
        x3 = self.x3.text()
        y1 = self.y1.text()
        y2 = self.y2.text()
        P = self.P.text()
        w = self.w.text()
        M = self.M.text() 
        print("finished saving inputs")
        param = [E, I, A, x1, x2, x3, y1, y2, P, w, M]
        print("before checking blanks")
        blank_inputs = self.checkBlanks(param)
        print("after checking blanks")
        if(blank_inputs == 1):
            print("blanks check passed")
            if(('-' not in E)&('-' not in I)&('-' not in A)&('-' not in x1)&('-' not in x2)&('-' not in x3)&('-' not in y1)&('-' not in y2)):
                print("no negative inputs")
                ss = self.compute(float(E), float(I), float(A), float(x1)*1000, float(x2)*1000, float(x3)*1000, float(y1)*1000, float(y2)*1000, float(P), float(w), float(M))
                print("after computation")
                n = 3
                m = 0
                for i in ss:
                    if (round(ss.item(m), n) == -0.0):
                        ss[m] = 0.0
                        m +=1
                    else:
                        m +=1

                self.x1_2.setText(str(x1))
                self.x2_2.setText(str(x2))
                self.x3_2.setText(str(x3))
                self.y1_2.setText(str(y1))
                self.y2_2.setText(str(y2))
                self.M_2.setText(str(M))
                self.w_2.setText(str(w))
                self.P_2.setText(str(P))


                self.d1.setText(str(round(ss.item(0), n)))
                self.d2.setText(str(round(ss.item(1), n)))
                self.d3.setText(str(round(ss.item(2), n)))
                self.d4.setText(str(round(ss.item(3), n)))
                self.d5.setText(str(round(ss.item(4), n)))
                self.d6.setText(str(round(ss.item(5), n)))
                self.p1.setText(str(round(ss.item(6), n)))
                self.p2.setText(str(round(ss.item(7), n)))
                self.p3.setText(str(round(ss.item(8), n)))
                self.p4.setText(str(round(ss.item(9), n)))
                self.p5.setText(str(round(ss.item(10), n)))
                self.p6.setText(str(round(ss.item(11), n)))
                self.r7.setText(str(round(ss.item(12), n)))
                self.r8.setText(str(round(ss.item(13), n)))
                self.r9.setText(str(round(ss.item(14), n)))
                self.r10.setText(str(round(ss.item(15), n)))
                self.r11.setText(str(round(ss.item(16), n)))
                self.r12.setText(str(round(ss.item(17), n)))
            else:
                print("negative valued inputs")
                self.showErrorDialog()
        else:
            print("blanks check failed")
            self.showErrorDialog()

        
    def compute(self, E, I, A, x1, x2, x3, y1, y2, P, W, M):
        print("entered compute")
    	Xb = [0, x1, x1+x2]
        Xe = [ x1, x1+x2, x1+x2+x3]
        Yb = [0, y1, y1]
        Ye = [y1, y1, y1-y2]
        w = [P, W, M]
        xd = [ Xe[0]-Xb[0], Xe[1]-Xb[1], Xe[2]-Xb[2] ]
        yd = [ Ye[0]-Yb[0], Ye[1]-Yb[1], Ye[2]-Yb[2] ]
        print("finished var assignment in compute")
        print("print xd and yd")
        print(xd)
        print(yd)
        L = [math.sqrt(((xd[0])**2)+((yd[0])**2)), math.sqrt((xd[1])**2+(yd[1])**2), math.sqrt((xd[2])**2 + (yd[2])**2)]
        print("L is calculated")
        cosTheta = [xd[0]/L[0], xd[1]/L[1], xd[2]/L[2]]
        print("cosTheta is achieved")
        sinTheta = [yd[0]/L[0], yd[1]/L[1], yd[2]/L[2]]
        print("sinTheta is achieved")
        AL2I = [(A*L[0]**2)/I, ((A*L[1]**2)/I), ((A*L[2]**2)/I)]
        EIL3 = [((E*I)/(L[0]**3)), ((E*I)/(L[1]**3)), ((E*I)/(L[2]**3))]
        print("before computeK func")
        def computeK (i):
            print("entered computeK func")
            pk = EIL3[i]*np.matrix([[((AL2I[i]*(cosTheta[i]**2))+(12*(sinTheta[i]**2))), ((AL2I[i]-12)*cosTheta[i]*sinTheta[i]), -6*sinTheta[i]*L[i]],
                  [((AL2I[i]-12)*(cosTheta[i])*(sinTheta[i])), ((AL2I[i]*(sinTheta[i]**2))+ (12*(cosTheta[i]**2))), 6*(cosTheta[i])*L[i]],
                  [(-6*(sinTheta[i])*L[i]), 6*(cosTheta[i])*L[i], 4*(L[i])**2]])

            k = np.matrix([[pk.item(0), pk.item(1), pk.item(2), -pk.item(0), -pk.item(1), pk.item(2)],
                           [pk.item(3), pk.item(4), pk.item(5), -pk.item(3), -pk.item(4), pk.item(5)],
                           [pk.item(6), pk.item(7), pk.item(8), -pk.item(6), -pk.item(7), pk.item(8)/2],
                           [-pk.item(0), -pk.item(1), -pk.item(2), pk.item(0), pk.item(1), -pk.item(2)],
                           [-pk.item(3), -pk.item(4), -pk.item(5), pk.item(3), pk.item(4), -pk.item(5)],
                           [pk.item(6), pk.item(7), pk.item(8)/2, -pk.item(6), -pk.item(7), pk.item(8)]])
            print("k will be returned after this")
            return k
            
        k1 = computeK(0)
        print(k1)
        k2 = computeK(1)
        print(k2)
        k3 = computeK(2)
        print(k3)
        print("s i calculated here")
        s = np.matrix([[k1.item(21)+k2.item(0), k1.item(22)+k2.item(1), k1.item(23)+k2.item(2), k2.item(3), k2.item(4), k2.item(5)],
                       [k1.item(27)+k2.item(6), k1.item(28)+k2.item(7), k1.item(29)+k2.item(8), k2.item(9), k2.item(10), k2.item(11)],
                       [k1.item(33)+k2.item(12), k1.item(34)+k2.item(13), k1.item(35)+k2.item(14), k2.item(15), k2.item(16), k2.item(17)],
                       [k2.item(18), k2.item(19), k2.item(20), k2.item(21)+k3.item(0), k2.item(22)+k3.item(1), k2.item(23)+k3.item(2)],
                       [k2.item(24), k2.item(25), k2.item(26), k2.item(27)+k3.item(6), k2.item(28)+k3.item(7), k2.item(29)+k3.item(8)],
                       [k2.item(30), k2.item(31), k2.item(32), k2.item(33)+k3.item(12), k2.item(34)+k3.item(13), k2.item(35)+k3.item(14)]])
        print(s)
        print("t will be calculated")
        t = np.array(s.I)
        print(t)
        print("q will be calculated here")
        q1 = [0, (w[0]*L[0]/2), (w[0]*(L[0]**2)/12), 0, (w[0]*L[0]/2), -(w[0]*(L[0]**2)/12)]
        print(q1)
        q2 = [0, 3*w[1]*L[1]/20, w[1]*(L[1]**2)/30, 0, 7*w[1]*L[1]/20, -(w[1]*(L[1]**2)/20)]
        print(q2)
        q3 = [w[2]*1000/2, 0, w[2]*L[2]*1000/8, w[2]*1000/2, 0, -(w[2]*L[2]*1000/8)]
        print(q3)
        qf = np.array([-sinTheta[0]*q1[1], (cosTheta[0]*q1[1])+q2[1], q1[5]+q2[2], q3[0], q2[4], q2[5]+q3[2]])
        q = np.array([0, 0, 0, 0, 0, 0])
        qqf = np.array(q-qf)
        print("deflection will be calculated here")
        d = np.dot(t, qqf)
        d1 = np.array([0, 0, 0, d[0], d[1], d[2]])
        d2 = d
        d3 = np.array([d[3], d[4], d[5], 0, 0, 0])
        print("f will be calculated here")
        f1 = np.dot(k1, d1)
        f2 = np.dot(k2, d2)
        f3 = np.dot(k3, d3)

        tm1 = [-sinTheta[0]*q1[1], cosTheta[0]*q1[1], q1[2], -sinTheta[0]*q1[1], cosTheta[0]*q1[1], q1[5]]
        tm2 = q2
        tm3 = q3
        print("reactions will be calculated here")
        reactions = np.array([(f1.item(3)+f2.item(0)+qf.item(0))/1000, (f1.item(4)+f2.item(1)+qf.item(1))/1000, (f1.item(5)+f2.item(2)+qf.item(2))/1000000, 
                            (f2.item(3)+f3.item(0)+qf.item(3))/1000, (f2.item(4)+f3.item(1)+qf.item(4))/1000, (f2.item(5)+f3.item(2)+qf.item(5))/1000000, 
                            (f1.item(0)+tm1[0])/1000, (f1.item(1)+tm1[1])/1000, (f1.item(2)+tm1[2])/1000000,
                            (f3.item(3)+tm3[3])/1000, (f3.item(4)+tm3[4])/1000, (f3.item(5)+tm3[5])/1000000])

        answer = np.append(d, reactions)
        print("answer is returned")
        return answer


    def showErrorDialog(self):
        self._new_window = ErrorDialog()
        self._new_window.show()

    def checkBlanks(self, param):
        t = 1
        for n in param:
            if (' ' in n) or (n == ''):
                self._new_window = ErrorDialog()
                self._new_window.show()
                t = 0
                break
        return t

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = ProblemWindow()
    window.show()
    sys.exit(app.exec_())
