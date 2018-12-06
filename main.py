import sys
from PyQt4 import QtCore, QtGui, uic
from problem import ProblemWindow


qtCreatorFile1 = "main.ui" # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile1)

class WelcomeWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(WelcomeWindow,self).__init__()
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.click_btn.clicked.connect(self.LoadProblem)
        
    def LoadProblem(self):
    	self._new_window = ProblemWindow()
    	self._new_window.show()


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = WelcomeWindow()
    window.show()
    sys.exit(app.exec_())