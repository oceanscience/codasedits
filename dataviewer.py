#!/usr/bin/env python
"""
Script for ADCP CODAS data viewer.

**Usage**

    *basic mode: view CODAS database*
        call from from a processing directory: dataviewer.py
        call from anywhere                   : dataviewer.py path/to/proc/dir

    *editing mode: `gautoedit`*
        call from processing "edit" dir      : dataviewer.py -e

    *singleping viewer (UHDAS only, correct configuration files in place)*
        call from processing "edit" dir      : dataviewer.py -p
        call from processing "edit" dir      : dataviewer.py -p [options]

    *compare 2 databases*
        call from anywhere:    dataviewer -c path1 path2

**Arguments**
    *path*
        path to processing directory or sonar name

**Options**
    -s              steps, duration (days) to view in panels
    -n              number of panels
    -v              view codas adcp data
    -e              edit codas adcp data
    -p              view codas and singleping data
    -c              view (compare) data from 2 sonars
    -t [--title]    title for panel and topo plots (overrides default)
    --dbname        path to database (up to but not including 'dir.blk')

    # for ping only, these are needed if 'dbinfo' is insufficient (or not exist)
    --cruisename    (ping only) base name for configuration file
    --beamangle     (ping only) beamangle
    --configtype    (ping only) read 'python' or 'matlab' config files
    --sonar         (ping only) such as 'os38bb', 'wh300', to read raw data

**Examples**::

   dataviewer.py path/to/cruiseID/sonar -n 3 -t sonar
   dataviewer.py -c path/to/sonar1 path/to/sonar2
HISTORY:  -revised Jan 15 2015, by Diana Cardoso (DC) (Bedford Institute of Oceanography, NS,Canada) for user defined start time
"""
import sys, os

from pycurrents.adcp.uhdasfile import guess_dbname
from pycurrents.system.misc import Bunch, Cachefile
from optparse import OptionParser

def get_dbparam(path, options):
    """
    For given *path* returns Bunch with database parameters:
        * sonar, crisename, beamangle, configtype, dbpathname

    *options* is the options object from the OptionParser
    """
    db_param = Bunch()
    db_param.yearbase = None
    db_param.beamangle = None
    db_param.badbeam = None

    dbinfo = get_dbinfo(path)
    if dbinfo is not None:
        for name in ('beamangle', 'configtype', 'cruisename', 
                     'sonar', 'badbeam'):
            if name in dbinfo.cachedict:
                db_param[name] = dbinfo.cachedict[name]
    else:
        db_param.cruisename, db_param.sonar = get_sonar(path)

    db_param.update(dict([(k, v) for (k, v) in
                    options.__dict__.items() if v is not None]))
    db_param.dbpathname = path

    if db_param.beamangle is not None:
        db_param.beamangle = int(db_param.beamangle)

    return db_param

def get_sonar(path):
    """ For a given path retuns cruisename and full sonar path """
    sonar_path = os.path.split(path)
    if sonar_path[1] == 'edit':
        sonar = os.path.split(sonar_path[0])
    else:
        sonar = sonar_path
    cruisename = os.path.split(sonar[0])[1]
    return cruisename, sonar_path[1]

def get_dbpath(path):
    """ For a given path retuns database path """
    try:
        dbpath = guess_dbname(path)
    except:
        print 'Could not find a dbname starting in %s' % (path)
        sys.exit()
    return dbpath

def get_dbinfo(dbpath):
    """ For a given database path if exist reads dbinfo.txt and returns """
    dbinfo = os.path.join(os.path.split(os.path.split(dbpath)[0])[0],
                      'dbinfo.txt')
    if os.path.exists(dbinfo):
        dbinfo = Cachefile(dbinfo)
        dbinfo.read()
    else:
        dbinfo = None
    return dbinfo

def main(arglist = None):
    parser = OptionParser(__doc__)
#EDIT Jan15, 2014: DC ****START
    parser.add_option("-d", "--start", dest="ddaystart",
                                help="start (days) to view in panels",
                                default=None)
#EDIT Jan15, 2014: DC ****End                                
    parser.add_option("-s", "--step", dest="ddaystep",
                                help="duration (days) to view in panels",
                                default=None)
    parser.add_option("-n", "--numpanels", dest="numaxes",
                                help="number of panels",
                                default=None)
    parser.add_option("-v", "--view", action="store_true", dest="view",
                                help="view codas adcp data",
                                default=False)
    parser.add_option("-e", "--edit", action="store_true", dest="edit",
                                help="edit codas adcp data",
                                default=False)
    parser.add_option("-p", "--ping", action="store_true", dest="ping",
                                help="view codas and singleping data",
                                default=False)
    parser.add_option("-c", "--compare", action="store_true", dest="compare",
                                help="view (compare) data from 2 sonars",
                                default=False)
    parser.add_option("--whitebg", action="store_true", dest="whitebg",
                                help="show pcolor masked (bad) values as white",
                                default=False)
    parser.add_option("--dbname", dest="dbpathname",
                                help='\n'.join(["path to database," +
                                "(up to but ", "not including 'dir.blk')"]))
    parser.add_option("-t", "--title", dest="plot_title",
                                help="title for panel and topo plots\n" + 
                                     "   (overrides default which is dbname)")
    parser.add_option("--cruisename", dest="cruisename",
                                help="cruise ID or title for plots " +
                                "(or use 'cruiseid'")
    parser.add_option("-m", "--minutes", dest="minutes",
                                help="number of minutes to average in topo plot" +
                                "(default is 30)")
    parser.add_option("--beamangle", dest="beamangle")
    parser.add_option("--configtype", dest="configtype",
                              help="read 'python' or 'matlab' config files" )
    parser.add_option("--sonar", dest="sonar",
                              help="such as 'os38bb', 'wh300'" )

    (options, args) = parser.parse_args(args=arglist)

    if (options.edit or options.ping):
        if os.path.split(os.getcwd())[1] != 'edit':
            print 'Must be in edit dir'
            sys.exit()
        elif len(args)== 1:
            print 'You can not point to dir, must be in it !'
            sys.exit()
    elif options.compare and len(args) < 2:
        print 'Must specify 2 sonar paths'
        sys.exit()

    if len(args)== 1:
        options.path=args[0]
    else:
        options.path = os.getcwd()

    mode = 'view'
    if options.edit:
        mode = 'edit'
    elif options.ping:
        mode = 'ping'
    elif options.compare:
        mode = 'compare'

    numaxes = options.numaxes
    if numaxes is None:
        numaxes = 4
        if mode == 'compare':
            numaxes = 6
#EDIT Jan15, 2014: DC ****START
    if options.ddaystart is None:
        ddaystart = 999
    else:
        ddaystart = float(options.ddaystart)
#EDIT Jan15, 2014: DC ****END       
    if options.ddaystep is None:
        ddaystep = 0.8
    else:
        ddaystep = float(options.ddaystep)

    if options.minutes is None:
        options.minutes = 30

    if options.whitebg:
        set_bad='w'
    else:
        set_bad=[.85,.85,.85]

    db_params = []
    if mode == 'compare':
        for s in (args[0], args[1]):
            db_param = Bunch()
            db_param.yearbase = None
            db_param.beamangle = None
            db_param.sonar = s
            path = os.path.join(options.path, s)
            db_param.dbpathname = get_dbpath(path)
            db_params.append(db_param)
        ## NOTE: This title correctly reflects what is done in 
        ##       compare.py but is completely hardwired here.
        plot_title = 'diff: %s - %s' % (db_params[1].sonar, 
                                        db_params[0].sonar)
    else:
        if options.dbpathname is not None:
            dbpath = options.dbpathname
        else:
            dbpath = get_dbpath(options.path)
        db_params.append(get_dbparam(dbpath, options))
        plot_title = db_params[0].sonar
    if options.plot_title is not None:
        plot_title = options.plot_title

    modemod = "pycurrents.adcpgui." + mode
    mod = __import__(modemod, fromlist=['GuiApp'])
    gg = mod.GuiApp(numaxes, db_params, ddaystep, ddaystart, plot_title, 
                    set_bad=set_bad, delta_t=int(options.minutes)/(24*60.)) #EDIT Jan15, 2014: DC, added ddaystart 
    return gg

if __name__ == '__main__':
    gg = main()
    gg.run()
