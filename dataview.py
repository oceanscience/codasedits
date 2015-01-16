#!/usr/bin/env python
"""
Base classes for shipboard ADCP data viewer:

-  *Main*
-  *GuiApp*
HISTORY:  -revised Jan 15 2015, by Diana Cardoso (DC) (Bedford Institute of Oceanography, NS,Canada) for user defined start time
"""
import wx
import sys

import numpy as np
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import  FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import  NavigationToolbar2WxAgg
from matplotlib.ticker import ScalarFormatter
formatter = ScalarFormatter(useOffset = False)

from pycurrents.system import logutils
log = logutils.getLogger(__file__)
import logging

import pycurrents.adcpgui.layout_aui as layout_aui
from pycurrents.adcpgui.cplotter import heading, jitter
from pycurrents.adcpgui.cplotter import CData, CPlotter

from pycurrents.adcpgui.cplotter import CDataDiff
from pycurrents.adcp.adcp_specs import codas_disableparams

from pycurrents.adcp import dataplotter
import pycurrents.adcpgui.tools as tools

pad = 0.05
saturate = ['0','100']
refbins = [2, 10]

display_opt = {}

def reset_display_opt():
    global display_opt

    display_opt.clear()
    display_opt['use_bins'] = False     # plot against bins or depth
    display_opt['use_cbarjet'] = False
    display_opt['axes'] = ['u', 'v', 'pg', 'amp', 'w','e','fvel','pvel']
    display_opt['show_spd'] = False
    display_opt['nbins'] = 128
    display_opt['velrange'] = [-.6,.6]
    display_opt['refbins'] = refbins
    display_opt['depth'] = None
    display_opt['mask'] = None


class WxTextCtrlHandler(logging.Handler):
    """Redirect stdout to a textctrl """
    def __init__(self, ctrl):
        logging.Handler.__init__(self)
        self.ctrl = ctrl

    def emit(self, record):
        s = self.format(record) + '\n'
        wx.CallAfter(self.ctrl.AppendText, s)

class Main(wx.Frame, layout_aui.Layout):

    # This base class inherits from layout_aui.Layout; subclasses
    # that have their own Layout must inherit from that before
    # inheriting from dataview.Main.  I put the Layout inheritance
    # here in the base class for clarity, since __init__ calls
    # its methods.

    minsize = (490,180)
    color = '#F6F2E1'

    def __init__(self, app, **kw):

        self._init_attributes()
        self.__dict__.update(kw)
        wx.Frame.__init__(self, None, -1, self.title, size=self.size)
        self.SetMinSize(self.minsize)
        self.SetBackgroundColour(self.color)
        self.app = app
        if hasattr(self.app.CD.data, 'resid_stats'):
            self.plot_choices.insert(-1,'resid_stats_fwd')

        self.ax_panel = []
        self.init_layout(display_opt)
        self._init_defaults()       #set plot defaults
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def _init_attributes(self):
        """
        Initialize attributes that may differ among subclasses.
        """
        ## It is not clear what needs to be here rather than being
        # a class attribute, but at least masking evidently does need
        # to be here. When it was a class attribute, the version in
        # edit.Main was not being used in the layout.
        self.plot_choices = ['u', 'v', 'fvel', 'pvel', 'amp', 'pg',
                        'w', 'e', 'pflag', heading]
        self.masking = ['no flags','codas', 1]
        self.show_saturate = False


    def OnClose(self, event):
        self.app.clear_asc_files()
        self._mgr.UnInit()
        self.app.OnQuit()
        #self.Destroy()

    def _init_defaults(self):
        self.set_plot_defaults()

    def set_plot_defaults(self):
        self.tb_start.SetValue('%4.2f' % self.startdd)
        self.tb_step.SetValue('%4.2f' % self.ddstep)
        self.velrng_min.SetValue(str(display_opt['velrange'][0]))
        self.velrng_max.SetValue(str(display_opt['velrange'][1]))
        self.refbinstart.SetValue(str(display_opt['refbins'][0]))
        self.refbinend.SetValue(str(display_opt['refbins'][1]))
        self.rb.SetSelection(self.masking[-1])
        self.saturate_vel.Show(self.show_saturate)
        for i, v in enumerate(display_opt['axes']):
            for p in self.plot_choices:
                self.ax_panel[i].Append(str(p))
            self.ax_panel[i].SetStringSelection(v)

    def OnDrawNow(self, event):
        self.draw_figures()

    def OnDrawPrev(self, event):
        
        self.draw_data(sgn = -1)

    def OnDrawNext(self, event):
        self.draw_data()

    def OnSaturate(self, event):
        pass

    def OnShowSpd(self, event):
        # "show speed" button toggled
        display_opt['show_spd'] = event.IsChecked()
        self.app.show_allspd(display_opt['axes'], display_opt['show_spd'])
        self.app.plotfig.canvas.draw()

    def OnUseBins(self, event):
        display_opt['use_bins'] = event.IsChecked()
        self.draw_figures()

    def OnUseJet(self, event):
        display_opt['use_cbarjet'] = event.IsChecked()
        self.draw_figures()

    def OnMasking(self, event):
        index = event.GetSelection()
        self.masking[-1] = index
        self.draw_figures()

    def OnPageChanging(self, event):
        pass

    def update_panel(self, event, num):
        # this is for automatic re-draw if panel variable changes
        ss = str(event.GetString())
        display_opt['axes'][num] = ss
        self.app.draw_ax(num, ss, speed=display_opt['show_spd'])
        self.app.plotfig.canvas.draw()

    def fix_refbins(self):
        if self.app.CD.data is not None:
            ssint = max(display_opt['refbins'][0], 1)
            self.refbinstart.SetValue(str(ssint))

            ssint = min(display_opt['refbins'][1], self.app.CD.data.nbins)
            self.refbinend.SetValue(str(ssint))

    def draw_data(self, sgn = 1):
        self.app.clear_asc_files()
        self.startdd += sgn*self.ddstep - sgn*pad
        self.tb_start.SetValue('%7.3f' %self.startdd)
        self.draw_figures()

    def draw_figures(self, newdata=False):
        self.app.set_lims(display_opt['velrange'])
        self.app.update_timerange(self.startdd, self.ddstep)
        self.printddrange() #only edit has override
        if newdata is True:
            self.app.CD.get_data(newdata)
        try:
            self.app.draw_cplots( self.masking[-1],
                                  display_opt['show_spd'],
                                  display_opt['use_bins'],
                                  display_opt['use_cbarjet'])
        except wx.PyDeadObjectError: # in case that pcolorframe has been closed
            pass
        try:
            self.fix_refbins()
            self.app.draw_topo()
        except wx.PyDeadObjectError: # in case that topoframe has been closed
            pass

    def printddrange(self):  
        dd1 = self.startdd
        dd2 = dd1 + self.ddstep
        log.info('extracting data: %7.3f - %7.3f' % (dd1,dd2))


class Plot(wx.Panel):
    def __init__(self, parent, id = -1, dpi = None,
                 add_toolbar=True, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)
        self.figure = Figure(dpi=dpi, figsize=(2,2))
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)

        if add_toolbar:
            self.toolbar = NavigationToolbar2WxAgg(self.canvas)
            self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)

        if add_toolbar:
            sizer.Add(self.toolbar, 0 , wx.LEFT | wx.EXPAND)

        self.SetSizer(sizer)


def nearest_notch(num, incr=.10, fcn='floor'):
    """
    Returns nearest notch (incr) below or above (floor or ceil) data
    """
    func = getattr(np, fcn)
    newfloor = func(num*(1./incr))
    return newfloor*incr


class GuiApp(object):

    use_diff = False
    use_thresholds = False

    def __init__(self, numaxes, db_params, ddaystep, ddaystart,plot_title, 
                 set_bad=None, delta_t=.02083):  #.5/24  i.e. 30min  EDIT Jan15, 2014: DC, added ddaystart 
        wxapp = wx.GetApp()
        if wxapp is None:
            wxapp = wx.PySimpleApp()
            wxapp.SetExitOnFrameDelete(True)
            self.made_app = True
        else:
            self.made_app = False
        self.wxapp = wxapp

        self.numaxes = numaxes
        self.db_params = db_params
        self.ddstep = ddaystep
        self.startdd = ddaystart  #EDIT Jan15, 2014: DC, added this line
        self.delta_t = delta_t
        self.plot_title = plot_title
        self.set_bad = set_bad
        reset_display_opt()
        self.init()

    def init(self):
#EDIT Jan15, 2014: DC ****START
        self.CDlist = []
        if self.startdd==999:
            for p in self.db_params: #Bunch; one per sonar (might be 'compare')
                CD = self.get_CD(p)
                CD.ddstep = self.ddstep
                self.CDlist.append((p.sonar,CD))

            self.CD = self.CDlist[0][1]
            if (len(self.CDlist) > 1 and
               (self.CDlist[1][1].startdd_all < self.CD.startdd_all)):
                self.CD = self.CDlist[1][1]
            self.startdd = nearest_notch(self.CD.startdd_all) #round down, up#self.startdd = nearest_notch(self.CD.startdd_all) #round down, up
        else:  
            for p in self.db_params: #Bunch; one per sonar (might be 'compare')
                CD = self.get_CD(p)
                CD.ddstep = self.ddstep
                CD.startdd = self.startdd
                self.CDlist.append((p.sonar,CD))

            self.CD = self.CDlist[0][1]
            if (len(self.CDlist) > 1 and
               (self.CDlist[1][1].startdd_all < self.CD.startdd_all)):
                self.CD = self.CDlist[1][1]    
#EDIT Jan15, 2014: DC ****END               
        self.enddd = nearest_notch(self.CD.enddd_all, incr=0.1, fcn='ceil')
        self.yearbase = self.CD.yearbase

        # get figure first
        self.add_pcolorframe()

        if self.use_diff:
            CDdiff = CDataDiff(self.CDlist)
            CDdiff.ddstep = self.ddstep
            self.CDlist.append(('diff', CDdiff))
            self.color_stsbar = self.pcolorframe.CreateStatusBar()
            self.color_stsbar.SetStatusText(CDdiff.title)

        self.CP = self.get_CP()

        self.mainframe = self.Main_frame(self,
                                    numaxes = self.numaxes,
                                    dday = self.CD.startdd_all,
                                    startdd = self.startdd,
                                    enddd = self.enddd,
                                    ddstep = self.ddstep,
                                    yearbase = self.yearbase,
                                    plot_title = self.plot_title)

        for l in log.handlers:
            l.setLevel(logging.ERROR)
        # set lower level
        handler = WxTextCtrlHandler(self.mainframe.log_display)
        log.addHandler(handler)
        self.handler = handler
        log.setLevel(logging.DEBUG)

        self.add_topoframe()
        self.update_timerange(self.startdd, self.ddstep)

        self.draw_cplots(self.mainframe.masking[-1])
        self.draw_topo()
        self.topoframe.Show()
        self.pcolorframe.Show()
        self.mainframe.Show()

        return True  #### try removing this

    def get_CD(self, param):
        return CData(**param)

    def get_CP(self):
        # find number of axes before making panel
        numaxes = int(self.numaxes)
        display_opt['axes'] = display_opt['axes'][:numaxes]
        self.axdict = tools.make_axes(fig=self.plotfig, numax=numaxes)
        return CPlotter(self.plotfig, self.axdict)

    def clear_asc_files(self):
        # this is a placeholder, used by edit only
        pass

    #-------- topography --------
    def add_topoframe(self):
        self.topoframe = wx.Frame(None, -1, 'topo', size=(400,400))
        self.topofig = Plot(self.topoframe).figure
        self.topofig.text(.5, .95, self.plot_title, ha='center')


    #----------- color plots -------------------
    def add_pcolorframe(self):
        self.pcolorframe = wx.Frame(None, -1, 'panels', size=(800,500))
        self.plotfig = Plot(self.pcolorframe).figure
        self.plotfig.text(.5, .95, self.plot_title, ha='center')

    def update_timerange(self, startdd, ddstep):
        """ Propagate new timerange to CD object """
        for k, v in self.CDlist:
            v.startdd = startdd
            v.ddstep = ddstep

    #TODO separate this method in 2 stage, so panels can be redraw
    # without geting a new set of data, ex. change in colorbar
    def draw_cplots(self, masking, speed=False,
                          use_bins=False, use_cbarjet=False):
        """ Redraws the colorplots or clear the axes """
        thresholds = {}
        jittercutoff = None
        if self.use_thresholds:
            T = self.mainframe.thresholds
            for k in T:
                if not T[k]['enabled']:
                    T[k]['value'] = codas_disableparams[k]
            thresholds = dict([(k,v['value']) for k, v in T.iteritems()])
            jittercutoff = thresholds['jitter_cutoff']

        for k, v in self.CDlist:  # one entry per instrument
            if v.data is not None:
                v.set_grid(display_opt['depth'], use_bins)
                v.fix_flags(masking, thresholds)

        for axnum in range(len(self.axdict['pcolor'])):
            self.draw_ax(axnum, display_opt['axes'][axnum],
                         speed, jittercutoff, use_cbarjet)
        self.plotfig.canvas.draw()

    def get_ax_attr(self, name):
        return name, '', self.CD

    def draw_ax(self, num, name, speed=False, jittercutoff=None, use_cbarjet=False):
        name, sonar, CDdata = self.get_ax_attr(name)
        if CDdata.data is not None:
            self.CP.draw_ax(CDdata, num, name, speed, jittercutoff,
                            use_cbarjet, sonar, set_bad=self.set_bad)
        else:
            self.CP.draw_ax(CDdata, num, name='clear', sonar=sonar)

    def set_lims(self, velrange):
        for name in ['u','v','fvel','pvel']:
            self.CP.clims[name] = velrange

    def show_allspd(self, names, speed=False):
        for axnum in range(len(self.axdict['pcolor'])):
            axname = names[axnum]
            name, sonar, CDdata = self.get_ax_attr(axname)
            if CDdata.data is None:
                continue
            self.CP.show_spd(self.CD, axnum, name, speed=speed)

    def draw_topo(self):
        """
         Redraws the topo figure
        """
        # TODO: this should (try to) respect refbins
        # should start at 1?
        # should rewrite to values used upon 'show now'
        if self.CD.data is not None:
            self.topofig.clear()
            self.topax = self.topofig.add_subplot(111)
            topz= self.CD.data.dep[display_opt['refbins'][0]-1]
            endz= self.CD.data.dep[display_opt['refbins'][-1]-1] #include this
            #add dz to include last bin
            cellm=np.diff(self.CD.data.dep)[0]
            deltaz= (endz-topz) + cellm

            dataplotter.vecplot(self.CD.data,
                                 refbins=[0,],     # we are specifying it now
                                 startz = topz,
                                 deltaz=deltaz,    # was 30,
                                 deltat = self.delta_t,  #30 min
                                 offset=cellm/2.,
                                 topo_kw=dict(reset_source=True),
                                 ax=self.topax)
        else:
            log.debug( 'no data found')
            self.topofig.clear()
            self.topax = self.topofig.add_subplot(111)
            self.topax.text(.5, .5, 'no data found', color='r', size=12,
                          transform=self.topax.transAxes, ha='center')
        self.topofig.canvas.draw()

    def run(self):
        if self.made_app:
            self._running_mainloop = True
            self.wxapp.MainLoop()
        else:
            # Presumably, if the app was already made, we are
            # running in ipython and don't need or want to block with
            # MainLoop().
            self._running_mainloop = False

    def OnQuit(self):
        log.removeHandler(self.handler)
        for w in wx.GetTopLevelWindows():
            w.Destroy()

