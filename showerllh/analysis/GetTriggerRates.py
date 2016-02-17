#!/usr/bin/env python
import sys, os, tables, datetime
import numpy as np
try:
    import MySQLdb as mysql
except:
    try:
        import pymysql as mysql
    except:
      raise RuntimeError("No access to mysql modules")

filtercomp = tables.Filters(complevel = 9)

def datetime_to_mjd(dt):
    MJD0 = datetime.datetime(year=1858, month=11, day=17)
    diff = dt - MJD0
    return diff.days + (diff.seconds + diff.microseconds / 1e6) / 86400.0


class triggerratedesc(tables.IsDescription):
    RunNumber = tables.Int32Col()
    MJDRun = tables.Float64Col()
    MJDStart = tables.Float64Col()
    MJDEnd = tables.Float64Col()
    SMT3 = tables.Float64Col(shape = (2,))
    SMT8 = tables.Float64Col(shape = (2,))
    ST = tables.Float64Col(shape = (2,))
    VT = tables.Float64Col(shape = (2,))
    SLOP = tables.Float64Col(shape = (2,))
    FRT = tables.Float64Col(shape = (2,))

h5_f = tables.openFile("triggerrates.h5", "w")
h5_t = h5_f.createTable("/", "TriggerRates", triggerratedesc, "Trigger rates for different triggers for a given run", filters = filtercomp)
h5_r = h5_t.row

c = mysql.connect(host = "cygnus.icecube.wisc.edu", user = "icecube", passwd = "skua", db = "live")
db = c.cursor()

def GetGoodRunList(db, StartRun, EndRun, SnapshotID):
    db.execute("""SELECT r.runNumber, 
                         r.tStart, 
                         r.tStop, 
                         ssr1.good_i3
                      FROM livedata_run AS r 
                      JOIN livedata_snapshotrun 
                      AS ssr1 ON ssr1.id = 
                      (
                            SELECT id
                            FROM livedata_snapshotrun AS ssr2
                            WHERE ssr2.run_id = r.id
                              AND ssr2.snapshot_id <= %s
                            ORDER BY snapshot_id DESC
                           LIMIT 1
                        )
                        JOIN livedata_snapshot AS ss
                        ON ss.id = ssr1.snapshot_id
                        WHERE r.runNumber >= %s
                        AND r.runNumber <= %s
                        AND ssr1.good_i3 = %s;""" % (SnapshotID, StartRun, EndRun, 1) )
        
    data = db.fetchall()
    return data

def GetTriggerRates(db, runnumber):
    db.execute("""SELECT trig.trigger_type_id, 
                         trig.rate, 
                         trig.rate_err
                  FROM livedata_trigger AS trig 
                  JOIN livedata_run AS run 
                  ON trig.r_id = run.id 
                  WHERE run.runNumber = %d;""" % runnumber)
    data = db.fetchall()
    return data

grl = GetGoodRunList(db, int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))

for r in grl:
    trigdata = np.array(GetTriggerRates(db, r[0]))
    h5_r["RunNumber"] = r[0]
    h5_r["MJDRun"] = (datetime_to_mjd(r[2]) + datetime_to_mjd(r[1])) / 2.
    h5_r["MJDStart"] = datetime_to_mjd(r[1])
    h5_r["MJDEnd"] = datetime_to_mjd(r[2])
    for t in trigdata:
        if r[0] == 118175: print t
        if t[0] == "32":
          h5_r["SMT8"] = np.array([float(t[1]), float(t[2])])
        if t[0] == "38":
           h5_r["SMT3"] = np.array([float(t[1]), float(t[2])])
        if t[0] == "36":
           h5_r["ST"] = np.array([float(t[1]), float(t[2])])
        if t[0] == "90":
           h5_r["SLOP"] = np.array([float(t[1]), float(t[2])])
        if t[0] == "91":
           h5_r["VT"] = np.array([float(t[1]), float(t[2])])
        if t[0] == "140" or t[0] == "39":
           h5_r["FRT"] = np.array([float(t[1]), float(t[2])])
    h5_r.append()
h5_t.flush()
h5_f.close()


