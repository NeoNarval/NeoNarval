#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np



class Cursor:
	"""Ex: fig=plt.plot(np.arange(10)) // c=cursor.Cursor(fig) // plt.show(fig)"""
	def onclick(self,event):
		#if (event.inaxes == self.fig.):
		if (event.button == 1) and (event.inaxes):
# 			print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
			plt.clf()
			#plt.imshow(self.image)
			plt.plot(self.imageX,self.imageY)
			if (self.tipo=="CROSS") or ("H" in self.tipo):
# 				plt.axhline(self.ydata,color='w')
				plt.axhline(event.ydata,color='w')
			if (self.tipo=="CROSS") or ("V" in self.tipo):
# 				plt.axvline(self.xdata,color='w')
				plt.axvline(event.xdata)
			self.x=event.x
			self.y=event.y
			self.xdata=event.xdata
			self.ydata=event.ydata
			self.button=event.button
			
			plt.draw()
			#plt.close(self.fig)
			#ax=self.fig.add_subplot(111)
			#ax.plot([0,0],[0,10])
			#self.fig.figure.canvas.draw()
		elif (event.button == 3):
			plt.close()
			#plt.disconnect(self.fig)
			
	#def __init__(self,fig,tipo=None):
	def __init__(self,fig,dataX,dataY,tipo=None):
		self.x=0
		self.y=0
		self.xdata=0
		self.ydata=0
		self.button=0
		if tipo==None:
			self.tipo='CROSS'
		else:
			self.tipo=tipo
		self.fig=fig
		self.imageX=dataX#fig.get_array()
		self.imageY=dataY	
		self.cid=plt.connect('button_press_event',self.onclick)

		
	


