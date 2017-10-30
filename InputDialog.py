from PySide.QtGui import *


class InputForm(QDialog):
    def __init__(self,mainWindow = None):
        super(self.__class__, self).__init__()
        self.mainWindow = mainWindow
        self.formLayout = QFormLayout()
        self.minDpLine = QLineEdit()
        self.minDpAsymLine = QLineEdit()
        self.maxDpAsymLine = QLineEdit()
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)

        self.formLayout.addRow(self.tr("&MinDp"),self.minDpLine)
        self.formLayout.addRow(self.tr("&MinDpAsym"), self.minDpAsymLine)
        self.formLayout.addRow(self.tr("&MaxDpAsym"), self.maxDpAsymLine)

        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)

        self.formLayout.addWidget(self.buttonBox)
        self.setLayout(self.formLayout)

    def acceptInput(self):
        try:
            self.minDp = float(self.minDpLine.text())
            self.minDpAsym = float(self.minDpAsymLine.text())
            self.maxDpAsym = float(self.maxDpAsymLine.text())
            self.accept()
        except:
            self.mainWindow.show_error_dialog("Input processed_data not valid.Please input again!")

    def getData(self):
        return self.minDp,self.minDpAsym, self.maxDpAsym


