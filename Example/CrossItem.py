from pyqtgraph import QtWidgets,QtGui, QtCore


class CrossItem(QtWidgets.QGraphicsItem):

    def __init__(self,cl,cw, x, y):
        super(CrossItem, self).__init__()
        #self.setFlag(QtGui.QGraphicsItem.ItemIsMovable, False)
        #self.setFlag(QtGui.QGraphicsItem.ItemIsSelectable, False)
        #self.setFlag(QtGui.QGraphicsItem.ItemIsFocusable, True)
        #self.setAcceptsHoverEvents(True)
        self.x = x
        self.y = y
        self.cl = cl
        self.cw = cw

    def boundingRect(self):
        return QtCore.QRectF(self.x, self.y, 2*self.cl, 2*self.cl)

    def set_prop(self,cl,cw,x,y):
        self.x = x
        self.y = y
        self.cl = cl
        self.cw = cw



    def paint(self, painter, option, widget):
        painter.setPen(QtGui.QPen(QtGui.QColor("red"), 1))
        #painter.setClipRect( option.exposedRect )
        path = QtGui.QPainterPath()
        x = self.x
        y = self.y
        cw = self.cw
        cl = self.cl
        path.moveTo(x-cw/2,y+cw/2)
        path.lineTo(x-cw/2,y+cl/2)
        path.lineTo(x+cw/2,y+cl/2)
        path.lineTo(x+cw/2,y+cw/2)
        path.lineTo(x+cl/2,y+cw/2)
        path.lineTo(x+cl/2,y-cw/2)
        path.lineTo(x+cw/2,y-cw/2)
        path.lineTo(x+cw/2,y-cl/2)
        path.lineTo(x-cw/2,y-cl/2)
        path.lineTo(x-cw/2,y-cw/2)
        path.lineTo(x-cl/2,y-cw/2)
        path.lineTo(x-cl/2,y+cw/2)
        path.lineTo(x-cw/2,y+cw/2)
        painter.drawPath(path)
