'''
This module is graphical API using pymatplotlib.
Specs:
-- We use matplotlib OO class level api, we do not use its high-level helper modules. Favor endured stability over simplicity. 
-- use TkAgg for interactive mode. Beaware of Tk,pyTk installation defects in various cern distributions.
-- PNG as default batch file format
-- we support http mode by sending string buf via meme type image/png. Sending a premade static plot to webserver is considered a uploading process instead of http dynamic graphical mode. Therefore covered in this module.
'''
import sys,os
import numpy,datetime
import matplotlib
from RecoLuminosity.LumiDB import CommonUtil

batchonly=False
if not os.environ.has_key('DISPLAY') or not os.environ['DISPLAY']:
    batchonly=True
    matplotlib.use('Agg',warn=False)
    from matplotlib.backends.backend_agg import FigureCanvasAgg as CanvasBackend
else:
    try:
        matplotlib.use('TkAgg',warn=False)
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as CanvasBackend
        from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
        import Tkinter as Tk
        root=Tk.Tk()
        root.wm_title("Lumi GUI in TK")
    except ImportError:
        print 'unable to import GUI backend, switch to batch only mode'
        matplotlib.use('Agg',warn=False)
        from matplotlib.backends.backend_agg import FigureCanvasAgg as CanvasBackend
        batchonly=True

from matplotlib.figure import Figure
from matplotlib.font_manager import fontManager,FontProperties
matplotlib.rcParams['lines.linewidth']=1.3
matplotlib.rcParams['grid.linewidth']=0.2
matplotlib.rcParams['xtick.labelsize']=9
matplotlib.rcParams['ytick.labelsize']=9
matplotlib.rcParams['legend.fontsize']=10
matplotlib.rcParams['axes.labelsize']=10
matplotlib.rcParams['font.weight']=550
def destroy(e) :
    sys.exit()
    
class matplotRender():
    def __init__(self,fig):
        self.__fig=fig
        self.__canvas=''
        self.colormap={}
        self.colormap['Delivered']='r'
        self.colormap['Recorded']='b'
        self.colormap['Effective']='g'
        self.colormap['Max Inst']='r'

    def plotSumX_Run(self,rawxdata,rawydata,nticks=6):
        xpoints=[]
        ypoints={}
        ytotal={}
        xidx=[]
        #print 'max rawxdata ',max(rawxdata)
        #print 'min rawxdata ',min(rawxdata)
        for x in CommonUtil.inclusiveRange(min(rawxdata),max(rawxdata),1):
            #print 'x : ',x
            xpoints.append(x)
            xidx.append(rawxdata.index(x)) #get the index of the sample points
            #print 'xidx : ',rawxdata.index(x)
        for ylabel,yvalues in rawydata.items():
            ypoints[ylabel]=[]
            for i in xidx:
                ypoints[ylabel].append(sum(yvalues[0:i])/1000.0)
            ytotal[ylabel]=sum(yvalues)/1000.0    
        ax=self.__fig.add_subplot(111)
        ax.set_xlabel(r'Run',position=(0.95,0))
        ax.set_ylabel(r'L nb$^{-1}$',position=(0,0.9))
        xticklabels=ax.get_xticklabels()
        for tx in xticklabels:
            tx.set_rotation(30)
        majorLocator=matplotlib.ticker.LinearLocator( nticks )
        majorFormatter=matplotlib.ticker.FormatStrFormatter('%d')
        minorLocator=matplotlib.ticker.LinearLocator(numticks=4*nticks)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_xbound(lower=xpoints[0],upper=xpoints[-1])
        ax.grid(True)
        keylist=ypoints.keys()
        keylist.sort()
        legendlist=[]
        for ylabel in keylist:
            cl='k'
            if self.colormap.has_key(ylabel):
                cl=self.colormap[ylabel]
            ax.plot(xpoints,ypoints[ylabel],label=ylabel,color=cl,drawstyle='steps')
            legendlist.append(ylabel+' '+'%.2f'%(ytotal[ylabel])+' '+'nb$^{-1}$')
        #font=FontProperties(size='medium',weight='demibold')
        #annotations
        trans=matplotlib.transforms.BlendedGenericTransform(ax.transData,ax.transAxes)
        ax.text(xpoints[0],1.025,str(xpoints[0]),transform=trans,horizontalalignment='center',size='x-small',color='green',bbox=dict(facecolor='white'))
        ax.text(xpoints[-1],1.025,str(xpoints[-1]),transform=trans,horizontalalignment='center',size='x-small',color='green',bbox=dict(facecolor='white'))
        ax.legend(tuple(legendlist),loc='upper left')
        self.__fig.subplots_adjust(bottom=0.1,left=0.1)
        
    def plotSumX_Fill(self,rawxdata,rawydata,rawfillDict,nticks=6):
        #print 'plotSumX_Fill rawxdata ',rawxdata
        ytotal={}
        ypoints={}
        xpoints=rawfillDict.keys()
        xpoints.sort()
        beginfo=''
        endinfo=''
        for ylabel,yvalue in rawydata.items():
            ypoints[ylabel]=[]
            ytotal[ylabel]=sum(rawydata[ylabel])/1000.0
            for idx,fill in enumerate(xpoints):
                runlist=rawfillDict[fill]
                if idx==0:
                    beginfo=str(fill)+':'+str(runlist[0])
                if idx==len(xpoints)-1:
                    endinfo=str(fill)+':'+str(runlist[-1])
                xidx=rawxdata.index(max(runlist))
                ypoints[ylabel].append(sum(yvalue[0:xidx])/1000.0)
        ax=self.__fig.add_subplot(111)
        ax.set_xlabel(r'LHC Fill Number',position=(0.84,0))
        ax.set_ylabel(r'L nb$^{-1}$',position=(0,0.9))
        ax.set_xbound(lower=xpoints[0],upper=xpoints[-1])
        xticklabels=ax.get_xticklabels()
        majorLocator=matplotlib.ticker.LinearLocator( nticks )
        majorFormatter=matplotlib.ticker.FormatStrFormatter('%d')
        #minorLocator=matplotlib.ticker.MultipleLocator(sampleinterval)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        #ax.xaxis.set_minor_locator(minorLocator)
        ax.grid(True)
        keylist=ypoints.keys()
        keylist.sort()
        legendlist=[]
        for ylabel in keylist:
            cl='k'
            if self.colormap.has_key(ylabel):
                cl=self.colormap[ylabel]
            ax.plot(xpoints,ypoints[ylabel],label=ylabel,color=cl,drawstyle='steps')
            legendlist.append(ylabel+' '+'%.2f'%(ytotal[ylabel])+' '+'nb$^{-1}$')
        #font=FontProperties(size='medium',weight='demibold')
        #annotations
        trans=matplotlib.transforms.BlendedGenericTransform(ax.transData,ax.transAxes)
        ax.text(xpoints[0],1.025,beginfo,transform=trans,horizontalalignment='left',size='x-small',color='green',bbox=dict(facecolor='white'))
        ax.text(xpoints[-1],1.025,endinfo,transform=trans,horizontalalignment='left',size='x-small',color='green',bbox=dict(facecolor='white'))
        ax.legend(tuple(legendlist),loc='upper left')
        self.__fig.subplots_adjust(bottom=0.1,left=0.1)
        
    def plotSumX_Time(self,rawxdata,rawydata,minTime,maxTime,nticks=6):
        '''
        input:
           rawxdata runDict{runnumber:[delivered,recorded,recorded_hltpath]}
           rawydata {label:[rundata]}
        '''
        xpoints=[]
        ypoints={}
        ytotal={}
        xidx=[]
        runs=rawxdata.keys()
        runs.sort()
        for idx,run in enumerate(runs):
            xpoints.append(matplotlib.dates.date2num(rawxdata[run][0]))
            xidx.append(idx)
        for ylabel,yvalue in rawydata.items():
            ypoints[ylabel]=[]
            for i in xidx:
                ypoints[ylabel].append(sum(yvalue[0:i])/1000.0)
            ytotal[ylabel]=sum(yvalue)/1000.0
        ax=self.__fig.add_subplot(111)
        dateFmt=matplotlib.dates.DateFormatter('%d/%m')
        majorLoc=matplotlib.ticker.LinearLocator(numticks=nticks)
        ax.xaxis.set_major_locator(majorLoc)
        minorLoc=matplotlib.ticker.LinearLocator(numticks=nticks*4)
        ax.xaxis.set_major_formatter(dateFmt)
        ax.set_xlabel(r'Date',position=(0.84,0))
        ax.set_ylabel(r'L nb$^{-1}$',position=(0,0.9))
        ax.xaxis.set_minor_locator(minorLoc)
        ax.set_xbound(lower=xpoints[0],upper=xpoints[-1])
        xticklabels=ax.get_xticklabels()
        for tx in xticklabels:
            tx.set_horizontalalignment('left')
        ax.grid(True)
        keylist=ypoints.keys()
        keylist.sort()
        legendlist=[]
        for ylabel in keylist:
            cl='k'
            if self.colormap.has_key(ylabel):
                cl=self.colormap[ylabel]
            ax.plot(xpoints,ypoints[ylabel],label=ylabel,color=cl,drawstyle='steps')
            legendlist.append(ylabel+' '+'%.2f'%(ytotal[ylabel])+' '+'nb$^{-1}$')
        #annotations
        trans=matplotlib.transforms.BlendedGenericTransform(ax.transData,ax.transAxes)
        #print 'run boundary ',runs[0],runs[-1]
        #print 'xpoints boundary ',xpoints[0],xpoints[-1]
        ax.text(xpoints[0],1.025,str(runs[0]),transform=trans,horizontalalignment='left',size='x-small',color='green',bbox=dict(facecolor='white'))        
        ax.text(xpoints[-1],1.025,str(runs[-1]),transform=trans,horizontalalignment='left',size='x-small',color='green',bbox=dict(facecolor='white'))
        
        ax.legend(tuple(legendlist),loc='upper left')
        self.__fig.autofmt_xdate(bottom=0.18,rotation=0)
        self.__fig.subplots_adjust(bottom=0.1,left=0.1)
        
    def plotPerdayX_Time(self,rawxdata,rawydata,minTime,maxTime,nticks=6):
        xpoints=[]
        ypoints={}
        ymax={}
        xidx=[]
        runs=rawxdata.keys()
        runs.sort()
        minDay=minTime.toordinal()
        maxDay=maxTime.toordinal()
        daydict={}
        for day in CommonUtil.inclusiveRange(minDay,maxDay,1):
            daydict[day]=[]#run index list
        for run in runs:
            runstartday=rawxdata[run][0].toordinal()
            if CommonUtil.findInList(daydict.keys() ,runstartday)!=-1:
                daydict[runstartday].append(runs.index(run))
        xpoints=daydict.keys()
        xpoints.sort()
        for ylabel,yvalue in rawydata.items():
            ypoints[ylabel]=[]
            ymax[ylabel]=[]
            for day,runindices in daydict.items():
                sumlumi=0.0
                maxlumi=0.0
                if len(runindices)!=0:
                    sumlumi=sum(yvalue[min(runindices):max(runindices)+1])/1000.0
                ypoints[ylabel].append(sumlumi)
            ymax[ylabel]=max(ypoints[ylabel])
        ax=self.__fig.add_subplot(111)
        dateFmt=matplotlib.dates.DateFormatter('%d/%m')
        majorLoc=matplotlib.ticker.LinearLocator(numticks=nticks)
        minorLoc=matplotlib.ticker.LinearLocator(numticks=nticks*4)
        ax.xaxis.set_major_formatter(dateFmt)
        ax.set_xlabel(r'Date',position=(0.84,0))
        ax.set_ylabel(r'L nb$^{-1}$',position=(0,0.9))
        ax.xaxis.set_major_locator(majorLoc)
        ax.xaxis.set_minor_locator(minorLoc)
        xticklabels=ax.get_xticklabels()
        for tx in xticklabels:
            tx.set_horizontalalignment('right')
        ax.grid(True)
        keylist=ypoints.keys()
        keylist.sort()
        legendlist=[]
        for ylabel in keylist:
            cl='k'
            if self.colormap.has_key(ylabel):
                cl=self.colormap[ylabel]
            ax.plot(xpoints,ypoints[ylabel],label=ylabel,color=cl,drawstyle='steps')
            legendlist.append(ylabel+' Max '+'%.2f'%(ymax[ylabel])+' '+'nb$^{-1}$')
        ax.legend(tuple(legendlist),loc='upper left')
        #ax.set_xlim(left=minDay,right=maxDay)
        ax.set_xbound(lower=xpoints[0],upper=xpoints[-1])
        self.__fig.autofmt_xdate(bottom=0.18,rotation=0)
        self.__fig.subplots_adjust(bottom=0.18,left=0.3)

    def plotPeakPerday_Time(self,daydict,minDay,maxDay,nticks=6):
        '''
        Input: daydict={}#{day:[run,lsnum,instlumi]}
        '''
        xpoints=[]
        ypoints=[]
        legendlist=[]
        days=daydict.keys()
        days.sort()
        beginfo=str(daydict[days[0]][0])+':'+str(daydict[days[0]][1])
        endinfo=str(daydict[days[-1]][0])+':'+str(daydict[days[-1]][1])
        maxinfo=''
        ymax=0.0
        for day in CommonUtil.inclusiveRange(minDay,maxDay,1):
            xpoints.append(day)
            if not daydict.has_key(day):
                ypoints.append(0.0)
            else:
                daymaxdata=daydict[day]
                ypoints.append(daymaxdata[2])
                if daydict[day][2]>ymax:
                    ymax=daydict[day][2]
                    runmax=daydict[day][0]
                    lsmax=daydict[day][1]
                    maxinfo=str(runmax)+':'+str(lsmax)
        ax=self.__fig.add_subplot(111)
        dateFmt=matplotlib.dates.DateFormatter('%d/%m')
        majorLoc=matplotlib.ticker.LinearLocator(numticks=nticks)
        minorLoc=matplotlib.ticker.LinearLocator(numticks=nticks*4)
        ax.xaxis.set_major_formatter(dateFmt)
        ax.set_xlabel(r'Date',position=(0.84,0))
        ax.set_ylabel(r'L $\mu$b$^{-1}$s$^{-1}$',position=(0,0.9))
        ax.xaxis.set_major_locator(majorLoc)
        ax.xaxis.set_minor_locator(minorLoc)
        xticklabels=ax.get_xticklabels()
        for tx in xticklabels:
            tx.set_horizontalalignment('right')
        ax.grid(True)
        cl=self.colormap['Max Inst']
        #print 'xpoints ',xpoints
        #print 'ypoints ',ypoints
        #print 'maxinfo ',maxinfo
        #print 'beginfo ',beginfo
        #print 'endinfo ',endinfo
        ax.plot(xpoints,ypoints,label='Max Inst',color=cl,drawstyle='steps')
        legendlist.append('Max Inst %.2f'%(ymax)+' '+'$\mu$b$^{-1}$s$^{-1}$')
        ax.legend(tuple(legendlist),loc='upper left')
        ax.set_xbound(lower=minDay,upper=maxDay)
        #annotations
        trans=matplotlib.transforms.BlendedGenericTransform(ax.transData,ax.transAxes)
        ax.text(xpoints[0],1.025,beginfo,transform=trans,horizontalalignment='left',size='x-small',color='green',bbox=dict(facecolor='white'))
        ax.text(xpoints[-1],1.025,endinfo,transform=trans,horizontalalignment='left',size='x-small',color='green',bbox=dict(facecolor='white'))
        self.__fig.autofmt_xdate(bottom=0.18,rotation=0)
        self.__fig.subplots_adjust(bottom=0.1,left=0.1)

    def plotInst_RunLS(self,rawxdata,rawydata,nticks=6):
        '''
        Input: rawxdata [run,fill,norbit,starttime,stoptime,totalls,ncmsls]
               rawydata {label:[instlumi]}
        '''
        runnum=rawxdata[0]
        fill=rawxdata[1]
        norbit=rawxdata[2]
        starttime=rawxdata[3]
        stoptime=rawxdata[4]
        totalls=rawxdata[-2]
        ncmsls=rawxdata[-1]
        peakinst=max(rawydata['Delivered'])
        lslength=float(norbit)*3564*25.0e-9
        totaldelivered=sum(rawydata['Delivered'])*lslength
        totalrecorded=sum(rawydata['Recorded'])*lslength
        xpoints=range(1,totalls+1)        
        #print len(xpoints)
        ypoints={}
        ymax={}
        for ylabel,yvalue in rawydata.items():
            ypoints[ylabel]=yvalue
            ymax[ylabel]=max(yvalue)
        left=0.15
        width=0.7
        bottom=0.1
        height=0.65
        bottom_h=bottom+height
        rect_scatter=[left,bottom,width,height]
        rect_table=[left,bottom_h,width,0.25]
        
        nullfmt=matplotlib.ticker.NullFormatter()
        nullloc=matplotlib.ticker.NullLocator()
        axtab=self.__fig.add_axes(rect_table,frameon=False)
        axtab.set_axis_off()
        axtab.xaxis.set_major_formatter(nullfmt)
        axtab.yaxis.set_major_formatter(nullfmt)
        axtab.xaxis.set_major_locator(nullloc)
        axtab.yaxis.set_major_locator(nullloc)

        ax=self.__fig.add_axes(rect_scatter)
        
        majorLoc=matplotlib.ticker.LinearLocator(numticks=nticks)
        minorLoc=matplotlib.ticker.LinearLocator(numticks=nticks*4)
        ax.set_xlabel(r'LS',position=(0.96,0))
        ax.set_ylabel(r'L $\mu$b$^{-1}$s$^{-1}$',position=(0,0.9))
        ax.xaxis.set_major_locator(majorLoc)
        ax.xaxis.set_minor_locator(minorLoc)
        ax.set_xbound(lower=xpoints[0],upper=xpoints[-1])
        xticklabels=ax.get_xticklabels()
        for tx in xticklabels:
            tx.set_horizontalalignment('right')
        ax.grid(True)
        keylist=ypoints.keys()
        keylist.sort()
        legendlist=[]

        for ylabel in keylist:
            cl='k'
            if self.colormap.has_key(ylabel):
                cl=self.colormap[ylabel]
            ax.plot(xpoints,ypoints[ylabel],'.',label=ylabel,color=cl)
            legendlist.append(ylabel)      
        #ax.axhline(0,color='green',linewidth=0.2)
        ax.axvline(xpoints[ncmsls-1],color='green',linewidth=0.2)
  
        colLabels=('run','fill','max inst(/$\mu$b/s)','delivered(/$\mu$b)','recorded(/$\mu$b)')
        cellText=[[str(runnum),str(fill),'%.3f'%(peakinst),'%.3f'%totaldelivered,'%.3f'%(totalrecorded)]]
       
        sumtable=axtab.table(cellText=cellText,colLabels=colLabels,colWidths=[0.12,0.1,0.27,0.27,0.27],cellLoc='center',loc='center')
        trans=matplotlib.transforms.BlendedGenericTransform(ax.transData,ax.transAxes)        
        axtab.add_table(sumtable)
        
        ax.text(xpoints[0],1.02,starttime[0:17],transform=trans,horizontalalignment='left',size='x-small',color='green',bbox=dict(facecolor='white'))   
        ax.text(xpoints[ncmsls-1],1.02,stoptime[0:17],transform=trans,horizontalalignment='left',size='x-small',color='green',bbox=dict(facecolor='white'))        
        ax.legend(tuple(legendlist),loc='upper right',numpoints=1)

    def drawHTTPstring(self):
        self.__canvas=CanvasBackend(self.__fig)    
        cherrypy.response.headers['Content-Type']='image/png'
        buf=StringIO()
        self.__canvas.print_png(buf)
        return buf.getvalue()
    
    def drawPNG(self,filename):
        self.__canvas=CanvasBackend(self.__fig)    
        self.__canvas.print_figure(filename)
    
    def drawInteractive(self):
        if batchonly:
            print 'interactive mode is not available for your setup, exit'
            sys.exit()    
        self.__canvas=CanvasBackend(self.__fig,master=root)
        self.__canvas.show()
        self.__canvas.get_tk_widget().pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
        toolbar=NavigationToolbar2TkAgg(self.__canvas,root)
        toolbar.update()
        self.__canvas._tkcanvas.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
        button = Tk.Button(master=root,text='Quit',command=sys.exit)
        button.pack(side=Tk.BOTTOM)
        Tk.mainloop()
if __name__=='__main__':
    fig=Figure(figsize=(8,8),dpi=100)
    a=fig.add_subplot(111)
    t=numpy.arange(0.0,3.0,0.01)
    s=numpy.sin(2*numpy.pi*t)
    a.plot(t,s)
    m=matplotRender(fig)
    m.drawPNG('testmatplotrender.png')
    m.drawInteractive()
    #print drawHTTPstring()   
