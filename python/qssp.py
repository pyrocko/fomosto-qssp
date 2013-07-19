# -*- coding: utf-8 -*-

import numpy as num
import logging, os, shutil, sys, glob

from multiprocessing import Pool

from tempfile import mkdtemp
from subprocess import Popen, PIPE
from os.path import join as pjoin

from guts import *
from guts_array import *
from pyrocko import trace, util, cake
from pyrocko.gf import store, builder, meta

from pyrocko.moment_tensor import MomentTensor, symmat6

from cStringIO import StringIO

logger = logging.getLogger('fomosto.qssp')

# how to call the programs
program_bins = {
    'qssp': 'qssp',
}

qssp_components = {
    1: 'ae an az gr sd te tn ue un uz ve vn vz'.split(), 
    2: 'ar at ap gr sd tt tp ur ut up vr vt vp'.split(),
}

def str_float_vals(vals):
    return ' '.join( [ '%e' % val for val in vals ] )

class QSSPModel(Object):
    data = Array.T(shape=(None,6), dtype=num.float)

    def string_for_config(self):
        srows = []
        for i, row in enumerate(self.data):
            srows.append( '%i %s' % (i+1, str_float_vals(row)) )
        
        return '\n'.join(srows)

    def get_nlines(self):
        return self.data.shape[0]

default_model_data = num.loadtxt(StringIO('''
  -86.00000     0.27410     0.00000  0.600000E-08  100.0      0.0
  -85.00000     0.27410     0.00000  0.700000E-08  100.0      0.0
  -84.00000     0.27534     0.00000  0.800000E-08  100.0      0.0
  -83.00000     0.27680     0.00000  0.100000E-07  100.0      0.0
  -82.00000     0.27825     0.00000  0.110000E-07  100.0      0.0
  -81.00000     0.27969     0.00000  0.130000E-07  100.0      0.0
  -80.00000     0.28112     0.00000  0.160000E-07  100.0      0.0
  -79.00000     0.28255     0.00000  0.180000E-07  100.0      0.0
  -78.00000     0.28396     0.00000  0.220000E-07  100.0      0.0
  -77.00000     0.28538     0.00000  0.250000E-07  100.0      0.0
  -76.00000     0.28678     0.00000  0.300000E-07  100.0      0.0
  -75.00000     0.28818     0.00000  0.350000E-07  100.0      0.0
  -74.00000     0.28957     0.00000  0.410000E-07  100.0      0.0
  -73.00000     0.29096     0.00000  0.470000E-07  100.0      0.0
  -72.00000     0.29233     0.00000  0.550000E-07  100.0      0.0
  -71.00000     0.29370     0.00000  0.640000E-07  100.0      0.0
  -70.00000     0.29561     0.00000  0.740000E-07  100.0      0.0
  -69.00000     0.29751     0.00000  0.860000E-07  100.0      0.0
  -68.00000     0.29940     0.00000  0.990000E-07  100.0      0.0
  -67.00000     0.30127     0.00000  0.114000E-06  100.0      0.0
  -66.00000     0.30313     0.00000  0.130000E-06  100.0      0.0
  -65.00000     0.30498     0.00000  0.149000E-06  100.0      0.0
  -64.00000     0.30682     0.00000  0.171000E-06  100.0      0.0
  -63.00000     0.30865     0.00000  0.195000E-06  100.0      0.0
  -62.00000     0.31047     0.00000  0.223000E-06  100.0      0.0
  -61.00000     0.31227     0.00000  0.254000E-06  100.0      0.0
  -60.00000     0.31407     0.00000  0.288000E-06  100.0      0.0
  -59.00000     0.31586     0.00000  0.327000E-06  100.0      0.0
  -58.00000     0.31763     0.00000  0.371000E-06  100.0      0.0
  -57.00000     0.31940     0.00000  0.420000E-06  100.0      0.0
  -56.00000     0.32116     0.00000  0.475000E-06  100.0      0.0
  -55.00000     0.32290     0.00000  0.537000E-06  100.0      0.0
  -54.00000     0.32464     0.00000  0.605000E-06  100.0      0.0
  -53.00000     0.32637     0.00000  0.682000E-06  100.0      0.0
  -52.00000     0.32809     0.00000  0.767000E-06  100.0      0.0
  -51.00000     0.32980     0.00000  0.862000E-06  100.0      0.0
  -50.00000     0.32980     0.00000  0.978000E-06  100.0      0.0
  -49.00000     0.32980     0.00000  0.110900E-05  100.0      0.0
  -48.00000     0.32980     0.00000  0.125800E-05  100.0      0.0
  -47.00000     0.32980     0.00000  0.142800E-05  100.0      0.0
  -46.00000     0.32809     0.00000  0.163800E-05  100.0      0.0
  -45.00000     0.32637     0.00000  0.188100E-05  100.0      0.0
  -44.00000     0.32464     0.00000  0.216400E-05  100.0      0.0
  -43.00000     0.32290     0.00000  0.249400E-05  100.0      0.0
  -42.00000     0.32116     0.00000  0.287800E-05  100.0      0.0
  -41.00000     0.31940     0.00000  0.332600E-05  100.0      0.0
  -40.00000     0.31763     0.00000  0.385100E-05  100.0      0.0
  -39.00000     0.31586     0.00000  0.446600E-05  100.0      0.0
  -38.00000     0.31407     0.00000  0.518700E-05  100.0      0.0
  -37.00000     0.31227     0.00000  0.603500E-05  100.0      0.0
  -36.00000     0.31047     0.00000  0.703400E-05  100.0      0.0
  -35.00000     0.30865     0.00000  0.821400E-05  100.0      0.0
  -34.00000     0.30682     0.00000  0.960900E-05  100.0      0.0
  -33.00000     0.30498     0.00000  0.112620E-04  100.0      0.0
  -32.00000     0.30313     0.00000  0.132250E-04  100.0      0.0
  -31.00000     0.30247     0.00000  0.154290E-04  100.0      0.0
  -30.00000     0.30180     0.00000  0.180120E-04  100.0      0.0
  -29.00000     0.30114     0.00000  0.210420E-04  100.0      0.0
  -28.00000     0.30047     0.00000  0.245990E-04  100.0      0.0
  -27.00000     0.29980     0.00000  0.287770E-04  100.0      0.0
  -26.00000     0.29913     0.00000  0.336880E-04  100.0      0.0
  -25.00000     0.29846     0.00000  0.394660E-04  100.0      0.0
  -24.00000     0.29778     0.00000  0.462670E-04  100.0      0.0
  -23.00000     0.29711     0.00000  0.542800E-04  100.0      0.0
  -22.00000     0.29643     0.00000  0.637270E-04  100.0      0.0
  -21.00000     0.29575     0.00000  0.748740E-04  100.0      0.0
  -20.00000     0.29507     0.00000  0.880350E-04  100.0      0.0
  -19.00000     0.29507     0.00000  0.103071E-03  100.0      0.0
  -18.00000     0.29507     0.00000  0.120676E-03  100.0      0.0
  -17.00000     0.29507     0.00000  0.141288E-03  100.0      0.0
  -16.00000     0.29507     0.00000  0.165420E-03  100.0      0.0
  -15.00000     0.29507     0.00000  0.193674E-03  100.0      0.0
  -14.00000     0.29507     0.00000  0.226753E-03  100.0      0.0
  -13.00000     0.29507     0.00000  0.265483E-03  100.0      0.0
  -12.00000     0.29507     0.00000  0.310828E-03  100.0      0.0
  -11.00000     0.29507     0.00000  0.363918E-03  100.0      0.0
  -10.00000     0.29946     0.00000  0.412707E-03  100.0      0.0
   -9.00000     0.30379     0.00000  0.466348E-03  100.0      0.0
   -8.00000     0.30806     0.00000  0.525168E-03  100.0      0.0
   -7.00000     0.31227     0.00000  0.589501E-03  100.0      0.0
   -6.00000     0.31643     0.00000  0.659697E-03  100.0      0.0
   -5.00000     0.32053     0.00000  0.736116E-03  100.0      0.0
   -4.00000     0.32458     0.00000  0.819129E-03  100.0      0.0
   -3.00000     0.32858     0.00000  0.909122E-03  100.0      0.0
   -2.00000     0.33253     0.00000  0.100649E-02  100.0      0.0
   -1.00000     0.33643     0.00000  0.111164E-02  100.0      0.0
    0.00000     0.34029     0.00000  0.122500E-02  100.0      0.0
    0.00000     5.80000     3.36000  0.272000E+01   1340.0    600.0
   20.00000     5.80000     3.36000  0.272000E+01   1340.0    600.0
   20.00000     6.50000     3.75000  0.292000E+01   1340.0    600.0
   35.00000     6.50000     3.75000  0.292000E+01   1340.0    600.0
   35.00000     8.04000     4.47000  0.331980E+01   1340.0    600.0
   77.50000     8.04500     4.48500  0.334550E+01   1340.0    600.0
  120.00000     8.05000     4.50000  0.337130E+01   1340.0    600.0
  165.00000     8.17500     4.50900  0.339850E+01    250.0    100.0
  210.00000     8.30000     4.51800  0.342580E+01    250.0    100.0
  210.00000     8.30000     4.52200  0.342580E+01    360.0    150.0
  260.00000     8.48250     4.60900  0.345610E+01    360.0    150.0
  310.00000     8.66500     4.69600  0.348640E+01    360.0    150.0
  360.00000     8.84750     4.78300  0.351670E+01    360.0    150.0
  410.00000     9.03000     4.87000  0.354700E+01    360.0    150.0
  410.00000     9.36000     5.07000  0.375570E+01    360.0    150.0
  460.00000     9.52800     5.17600  0.381750E+01    360.0    150.0
  510.00000     9.69600     5.28200  0.387930E+01    360.0    150.0
  560.00000     9.86400     5.38800  0.394100E+01    360.0    150.0
  610.00000    10.03200     5.49400  0.400280E+01    360.0    150.0
  660.00000    10.20000     5.60000  0.406460E+01    360.0    150.0
  660.00000    10.79000     5.95000  0.437140E+01    780.0    300.0
  710.00000    10.92290     6.07970  0.440100E+01    780.0    300.0
  760.00000    11.05580     6.20950  0.443050E+01    780.0    300.0
  809.50000    11.14400     6.24740  0.445960E+01    780.0    300.0
  859.00000    11.23000     6.28410  0.448850E+01    780.0    300.0
  908.50000    11.31400     6.31990  0.451730E+01    780.0    300.0
  958.00000    11.39600     6.35460  0.454590E+01    780.0    300.0
 1007.50000    11.47610     6.38830  0.457440E+01    780.0    300.0
 1057.00000    11.55430     6.42110  0.460280E+01    780.0    300.0
 1106.50000    11.63080     6.45300  0.463100E+01    780.0    300.0
 1156.00000    11.70560     6.48410  0.465910E+01    780.0    300.0
 1205.50000    11.77870     6.51430  0.468700E+01    780.0    300.0
 1255.00000    11.85040     6.54380  0.471480E+01    780.0    300.0
 1304.50000    11.92050     6.57250  0.474240E+01    780.0    300.0
 1354.00000    11.98930     6.60060  0.476990E+01    780.0    300.0
 1403.50000    12.05680     6.62800  0.479730E+01    780.0    300.0
 1453.00000    12.12310     6.65470  0.482450E+01    780.0    300.0
 1502.50000    12.18810     6.68090  0.485150E+01    780.0    300.0
 1552.00000    12.25210     6.70660  0.487850E+01    780.0    300.0
 1601.50000    12.31510     6.73170  0.490520E+01    780.0    300.0
 1651.00000    12.37720     6.75640  0.493190E+01    780.0    300.0
 1700.50000    12.43830     6.78070  0.495840E+01    780.0    300.0
 1750.00000    12.49870     6.80460  0.498470E+01    780.0    300.0
 1799.50000    12.55840     6.82820  0.501090E+01    780.0    300.0
 1849.00000    12.61740     6.85140  0.503700E+01    780.0    300.0
 1898.50000    12.67590     6.87450  0.506290E+01    780.0    300.0
 1948.00000    12.73390     6.89720  0.508870E+01    780.0    300.0
 1997.50000    12.79150     6.91990  0.511430E+01    780.0    300.0
 2047.00000    12.84870     6.94230  0.513980E+01    780.0    300.0
 2096.50000    12.90570     6.96470  0.516520E+01    780.0    300.0
 2146.00000    12.96250     6.98700  0.519040E+01    780.0    300.0
 2195.50000    13.01920     7.00930  0.521540E+01    780.0    300.0
 2245.00000    13.07580     7.03160  0.524030E+01    780.0    300.0
 2294.50000    13.13250     7.05400  0.526510E+01    780.0    300.0
 2344.00000    13.18920     7.07650  0.528980E+01    780.0    300.0
 2393.50000    13.24620     7.09910  0.531420E+01    780.0    300.0
 2443.00000    13.30340     7.12180  0.533860E+01    780.0    300.0
 2492.50000    13.36100     7.14490  0.536280E+01    780.0    300.0
 2542.00000    13.41900     7.16810  0.538690E+01    780.0    300.0
 2591.50000    13.47740     7.19170  0.541080E+01    780.0    300.0
 2641.00000    13.53640     7.21560  0.543450E+01    780.0    300.0
 2690.50000    13.59610     7.23980  0.545820E+01    780.0    300.0
 2740.00000    13.65640     7.26450  0.548170E+01    780.0    300.0
 2740.00000    13.65640     7.26450  0.548170E+01    780.0    300.0
 2789.67000    13.66790     7.27680  0.550510E+01    780.0    300.0
 2839.33000    13.67930     7.28920  0.552840E+01    780.0    300.0
 2889.00000    13.69080     7.30150  0.555150E+01    780.0    300.0
 2889.00000     8.00880     0.00000  0.991450E+01  57822.0      0.0
 2939.33000     8.09630     0.00000  0.999420E+01  57822.0      0.0
 2989.66000     8.18210     0.00000  0.100722E+02  57822.0      0.0
 3039.99000     8.26620     0.00000  0.101485E+02  57822.0      0.0
 3090.32000     8.34860     0.00000  0.102233E+02  57822.0      0.0
 3140.66000     8.42930     0.00000  0.102964E+02  57822.0      0.0
 3190.99000     8.50830     0.00000  0.103679E+02  57822.0      0.0
 3241.32000     8.58560     0.00000  0.104378E+02  57822.0      0.0
 3291.65000     8.66110     0.00000  0.105062E+02  57822.0      0.0
 3341.98000     8.73500     0.00000  0.105731E+02  57822.0      0.0
 3392.31000     8.80720     0.00000  0.106385E+02  57822.0      0.0
 3442.64000     8.87760     0.00000  0.107023E+02  57822.0      0.0
 3492.97000     8.94640     0.00000  0.107647E+02  57822.0      0.0
 3543.30000     9.01340     0.00000  0.108257E+02  57822.0      0.0
 3593.64000     9.07870     0.00000  0.108852E+02  57822.0      0.0
 3643.97000     9.14240     0.00000  0.109434E+02  57822.0      0.0
 3694.30000     9.20430     0.00000  0.110001E+02  57822.0      0.0
 3744.63000     9.26450     0.00000  0.110555E+02  57822.0      0.0
 3794.96000     9.32300     0.00000  0.111095E+02  57822.0      0.0
 3845.29000     9.37980     0.00000  0.111623E+02  57822.0      0.0
 3895.62000     9.43490     0.00000  0.112137E+02  57822.0      0.0
 3945.95000     9.48830     0.00000  0.112639E+02  57822.0      0.0
 3996.28000     9.54000     0.00000  0.113127E+02  57822.0      0.0
 4046.62000     9.59000     0.00000  0.113604E+02  57822.0      0.0
 4096.95000     9.63830     0.00000  0.114069E+02  57822.0      0.0
 4147.28000     9.68480     0.00000  0.114521E+02  57822.0      0.0
 4197.61000     9.72970     0.00000  0.114962E+02  57822.0      0.0
 4247.94000     9.77280     0.00000  0.115391E+02  57822.0      0.0
 4298.27000     9.81430     0.00000  0.115809E+02  57822.0      0.0
 4348.60000     9.85400     0.00000  0.116216E+02  57822.0      0.0
 4398.93000     9.89200     0.00000  0.116612E+02  57822.0      0.0
 4449.26000     9.92840     0.00000  0.116998E+02  57822.0      0.0
 4499.60000     9.96300     0.00000  0.117373E+02  57822.0      0.0
 4549.93000     9.99590     0.00000  0.117737E+02  57822.0      0.0
 4600.26000    10.02710     0.00000  0.118092E+02  57822.0      0.0
 4650.59000    10.05660     0.00000  0.118437E+02  57822.0      0.0
 4700.92000    10.08440     0.00000  0.118772E+02  57822.0      0.0
 4751.25000    10.11050     0.00000  0.119098E+02  57822.0      0.0
 4801.58000    10.13490     0.00000  0.119414E+02  57822.0      0.0
 4851.91000    10.15760     0.00000  0.119722E+02  57822.0      0.0
 4902.24000    10.17850     0.00000  0.120021E+02  57822.0      0.0
 4952.58000    10.19780     0.00000  0.120311E+02  57822.0      0.0
 5002.91000    10.21540     0.00000  0.120593E+02  57822.0      0.0
 5053.24000    10.23120     0.00000  0.120867E+02  57822.0      0.0
 5103.57000    10.24540     0.00000  0.121133E+02  57822.0      0.0
 5153.90000    10.25780     0.00000  0.121391E+02  57822.0      0.0
 5153.90000    11.09140     3.43850  0.127037E+02    780.0    300.0
 5204.61000    11.10360     3.44880  0.127289E+02    780.0    300.0
 5255.32000    11.11530     3.45870  0.127530E+02    780.0    300.0
 5306.04000    11.12650     3.46810  0.127760E+02    780.0    300.0
 5356.75000    11.13710     3.47700  0.127980E+02    780.0    300.0
 5407.46000    11.14720     3.48560  0.128188E+02    780.0    300.0
 5458.17000    11.15680     3.49370  0.128387E+02    780.0    300.0
 5508.89000    11.16590     3.50130  0.128574E+02    780.0    300.0
 5559.60000    11.17450     3.50850  0.128751E+02    780.0    300.0
 5610.31000    11.18250     3.51530  0.128917E+02    780.0    300.0
 5661.02000    11.19010     3.52170  0.129072E+02    780.0    300.0
 5711.74000    11.19710     3.52760  0.129217E+02    780.0    300.0
 5762.45000    11.20360     3.53300  0.129351E+02    780.0    300.0
 5813.16000    11.20950     3.53810  0.129474E+02    780.0    300.0
 5863.87000    11.21500     3.54270  0.129586E+02    780.0    300.0
 5914.59000    11.21990     3.54680  0.129688E+02    780.0    300.0
 5965.30000    11.22430     3.55050  0.129779E+02    780.0    300.0
 6016.01000    11.22820     3.55380  0.129859E+02    780.0    300.0
 6066.72000    11.23160     3.55670  0.129929E+02    780.0    300.0
 6117.44000    11.23450     3.55910  0.129988E+02    780.0    300.0
 6168.15000    11.23680     3.56100  0.130036E+02    780.0    300.0
 6218.86000    11.23860     3.56260  0.130074E+02    780.0    300.0
 6269.57000    11.23990     3.56370  0.130100E+02    780.0    300.0
 6320.29000    11.24070     3.56430  0.130117E+02    780.0    300.0
 6371.00000    11.24090     3.56450  0.130122E+02    780.0    300.0
'''.lstrip()))


class QSSPSource(Object):
    lat = Float.T(default=0.0) 
    lon = Float.T(default=0.0)
    depth = Float.T(default=10.0)
    torigin = Float.T(default=0.0)
    trise = Float.T(default=1.0)

    def string_for_config(self):
        return '%(lat)15e %(lon)15e %(depth)15e %(torigin)15e %(trise)15e' % \
                self.__dict__

class QSSPSourceMT(QSSPSource):
    munit = Float.T(default=1.0)
    mrr = Float.T(default=1.0)
    mtt = Float.T(default=1.0)
    mpp = Float.T(default=1.0)
    mrt = Float.T(default=0.0)
    mrp= Float.T(default=0.0)
    mtp = Float.T(default=0.0)

    def string_for_config(self):
        return '%(munit)15e %(mrr)15e %(mtt)15e %(mpp)15e ' \
               '%(mrt)15e %(mrp)15e %(mtp)15e ' % self.__dict__ + \
               QSSPSource.string_for_config(self)

class QSSPSourceDC(QSSPSource):
    moment = Float.T(default=1.0e9)
    strike = Float.T(default=0.0)
    dip = Float.T(default=90.0)
    rake = Float.T(default=0.0)

    def string_for_config(self):
        return '%(moment)15e %(strike)15e %(dip)15e %(rake)15e ' % \
                self.__dict__ + QSSPSource.string_for_config(self)

class QSSPReceiver(Object):
    lat = Float.T(default=10.0)
    lon = Float.T(default=0.0)
    name = String.T(default='')
    tstart = Float.T(default=0.0)

    def string_for_config(self):
        return "%(lat)15e %(lon)15e '%(name)s' %(tstart)e" % self.__dict__

class QSSPGreen(Object):
    depth = Float.T(default=10.0)
    filename = String.T(default='GF_10km')
    calculate = Bool.T(default=True)

    def string_for_config(self):
        return "%(depth)15e '%(filename)s' %(calculate)i" % self.__dict__

class QSSPConfig(Object):

    time_window = Float.T(default=900.0)
    frequency_max = Float.T(optional=True)
    slowness_max = Float.T(default=0.4)
    antialiasing_factor = Float.T(default=0.1)

    lowpass_order = Int.T(default=0, optional=True)
    lowpass_corner = Float.T(default=1.0, optional=True)
    output_slowness_min = Float.T(default=0.0, optional=True)
    output_slowness_max = Float.T(optional=True)

    spheroidal_modes = Bool.T(default=True)
    toroidal_modes = Bool.T(default=True)

    cutoff_harmonic_degree_sd = Int.T(default=0)
    crit_frequency_sge = Float.T(default=0.0)
    crit_harmonic_degree_sge = Int.T(default=0)

    include_physical_dispersion = Bool.T(default=False)

    source_patch_radius = Float.T(default=0.0) 

    model = QSSPModel.T(default=QSSPModel(data=default_model_data))

    def items(self):
        return dict( self.T.inamevals(self) )

class QSSPConfigFull(QSSPConfig):
    receiver_depth = Float.T(default=0.0)
    sampling_interval = Float.T(default=5.0)

    output_filename = String.T(default='receivers')
    output_format = Int.T(default=1)
    output_time_window = Float.T(optional=True)

    gf_directory = String.T(default='QSSP_Green')
    greens_functions = List.T(QSSPGreen.T())

    sources = List.T(QSSPSource.T())
    receivers = List.T(QSSPReceiver.T())

    @staticmethod
    def example():
        conf = QSSPConfigFull()
        conf.sources.append( QSSPSourceMT() )
        lats = [ 20. ]
        conf.receivers.extend( QSSPReceiver(lat=lat) for lat in lats )
        conf.greens_functions.append( QSSPGreen() )
        return conf

    @property
    def components(self):
        return qssp_components[self.output_format]
    
    def get_output_filenames(self, rundir):
        return [ pjoin(rundir, self.output_filename+'.'+c) for c in self.components ]

    def ensure_gf_directory(self):
        util.ensuredir(self.gf_directory)

    def string_for_config(self):

        def aggregate(l):
            return len(l), '\n'.join( x.string_for_config() for x in l )

        assert len(self.greens_functions) > 0
        assert len(self.sources) > 0
        assert len(self.receivers) > 0
        
        d = self.__dict__.copy()

        if self.output_time_window is None:
            d['output_time_window'] = self.time_window

        if self.output_slowness_max is None:
            d['output_slowness_max'] = self.slowness_max

        if self.frequency_max is None:
            d['frequency_max'] = 0.5/self.sampling_interval

        d['gf_directory'] = os.path.abspath(self.gf_directory) + '/'

        d['n_receiver_lines'], d['receiver_lines'] = aggregate(self.receivers)
        d['n_source_lines'], d['source_lines'] = aggregate(self.sources)
        d['n_gf_lines'], d['gf_lines'] = aggregate(self.greens_functions)
        d['n_model_lines'] = self.model.get_nlines()
        d['model_lines'] = self.model.string_for_config()


        if len(self.sources) == 0 or isinstance(self.sources[0], QSSPSourceMT):
            d['point_source_type'] = 1
        else:
            d['point_source_type'] = 2

        template = '''# autogenerated QSSP input by qssp.py
#
# This is the input file of FORTRAN77 program "qssp2010" for calculating
# synthetic seismograms of a self-gravitating, spherically symmetric,
# isotropic and viscoelastic earth.
#
# by
# Rongjiang  Wang <wang@gfz-potsdam.de>
# Helmholtz-Centre Potsdam
# GFZ German Reseach Centre for Geosciences
# Telegrafenberg, D-14473 Potsdam, Germany
#
# Last modified: Potsdam, July, 2010
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# If not specified, SI Unit System is used overall!
#
# Coordinate systems:
# spherical (r,t,p) with r = radial,
#                        t = co-latitude,
#                        p = east longitude.
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#	UNIFORM RECEIVER DEPTH
#	======================
# 1. uniform receiver depth [km]
#-------------------------------------------------------------------------------------------
    %(receiver_depth)e
#-------------------------------------------------------------------------------------------
#
#	TIME (FREQUENCY) SAMPLING
#	=========================
# 1. time window [sec], sampling interval [sec]
# 2. max. frequency [Hz] of Green's functions
# 3. max. slowness [s/km] of Green's functions
# 4. anti-aliasing factor (> 0 & < 1), if it is <= 0 or >= 1/e (~ 0.4), then
#    default value of 1/e is used (e.g., 0.1 = alias phases will be suppressed
#    to 10%% of their original amplitude)
#
#    Note: The computation effort increases linearly the time window and
#          quadratically with the cut-off frequency.
#-------------------------------------------------------------------------------------------
    %(time_window)e   %(sampling_interval)e
    %(frequency_max)e 
    %(slowness_max)e
    %(antialiasing_factor)e
#-------------------------------------------------------------------------------------------
#
#	SELF-GRAVITATING EFFECT
#	=======================
# 1. the critical frequency [Hz] and the critical harmonic degree, below which
#    the self-gravitating effect should be included
#-------------------------------------------------------------------------------------------
    %(crit_frequency_sge)e     %(crit_harmonic_degree_sge)i
#-------------------------------------------------------------------------------------------
#
#	WAVE TYPES
#	==========
# 1. selection (1/0 = yes/no) of speroidal modes (P-SV waves), selection of toroidal modes
#    (SH waves), and cutoff harmonic degree for static deformation
#-------------------------------------------------------------------------------------------
    %(spheroidal_modes)i     %(toroidal_modes)i    %(cutoff_harmonic_degree_sd)i
#-------------------------------------------------------------------------------------------
#	GREEN'S FUNCTION FILES
#	======================
# 1. number of discrete source depths, estimated radius of each source patch [km] and
#    directory for Green's functions
# 2. list of the source depths [km], the respective file names of the Green's
#    functions (spectra) and the switch number (0/1) (0 = do not calculate
#    this Green's function because it exists already, 1 = calculate or update
#    this Green's function. Note: update is required if any of the above
#    parameters is changed)
#-------------------------------------------------------------------------------------------
   %(n_gf_lines)i   %(source_patch_radius)e   '%(gf_directory)s'
   %(gf_lines)s
#-------------------------------------------------------------------------------------------
#
#	MULTI-EVENT SOURCE PARAMETERS
#	=============================
# 1. number of discrete point sources and selection of the source data format
#    (1 or 2)
# 2. list of the multi-event sources
#    Format 1:
#    M-Unit   Mrr  Mtt  Mpp  Mrt  Mrp  Mtp  Lat   Lon   Depth  T_origin T_rise
#    [Nm]                                   [deg] [deg] [km]   [sec]    [sec]
#    Format 2:
#    Moment   Strike    Dip       Rake      Lat   Lon   Depth  T_origin T_rise
#    [Nm]     [deg]     [deg]     [deg]     [deg] [deg] [km]   [sec]    [sec]
#-------------------------------------------------------------------------------------------
  %(n_source_lines)i     %(point_source_type)i 
%(source_lines)s
#-------------------------------------------------------------------------------------------
#
#	RECEIVER PARAMETERS
#	===================
# 1. output file name and and selection of output format:
#       1 = cartesian: vertical(z)/north(n)/east(e);
#       2 = spherical: radial(r)/theta(t)/phi(p)
#      (Note: if output format 2 is selected, the epicenter (T_origin = 0)
# 2. output time window [sec] (<= Green's function time window)
# 3. selection of order of Butterworth low-pass filter (if <= 0, then no filtering), corner
#    frequency (smaller than the cut-off frequency defined above)
# 4. lower and upper slowness cut-off [s/km] (slowness band-pass filter)
# 5. number of receiver
# 6. list of the station parameters
#    Format:
#    Lat     Lon    Name     Time_reduction
#    [deg]   [deg]           [sec]
#    (Note: Time_reduction = start time of the time window)
#-------------------------------------------------------------------------------------------
  '%(output_filename)s'  %(output_format)i
  %(output_time_window)e 
  %(lowpass_order)i       %(lowpass_corner)f
  %(output_slowness_min)e    %(output_slowness_max)e
  %(n_receiver_lines)i
%(receiver_lines)s
#-------------------------------------------------------------------------------------------
#
#	                LAYERED EARTH MODEL (IASP91)
#                   ============================
# 1. number of data lines of the layered model and selection for including
#    the physical dispersion according Kamamori & Anderson (1977)
#-------------------------------------------------------------------------------------------
    %(n_model_lines)i    %(include_physical_dispersion)i
#-------------------------------------------------------------------------------------------
#
#	MULTILAYERED MODEL PARAMETERS (source site)
#	===========================================
# no depth[km] vp[km/s] vs[km/s] ro[g/cm^3]      qp     qs
#-------------------------------------------------------------------------------------------
%(model_lines)s
#---------------------------------end of all inputs-----------------------------------------
'''

        return template % d

        
class QSSPError(Exception):
    pass
        
class QSSPRunner:
    
    def __init__(self, tmp=None, keep_tmp=False):
            
        self.tempdir = mkdtemp(prefix='qssprun-', dir=tmp)
        self.keep_tmp = keep_tmp
        self.program = program_bins['qssp']
        self.config = None
    
    def run(self, config):
        self.config = config
        config.ensure_gf_directory()
        
        input_fn = pjoin(self.tempdir, 'input')
                
        f = open(input_fn, 'w')
        input_str = config.string_for_config()
        
        logger.debug('===== begin qssp input =====\n'
                     '%s===== end qssp input =====' % input_str)
        
        f.write( input_str )
        f.close()
        program = self.program
        
        old_wd = os.getcwd()

        os.chdir(self.tempdir)
        
        try:
            proc = Popen(program, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        except OSError:
            os.chdir(old_wd)
            raise QSSPError('could not start qssp: "%s"' % program)
        
        (output_str, error_str) = proc.communicate('input\n')
       
        logger.debug('===== begin qssp output =====\n'
                     '%s===== end qssp output =====' % output_str)
        if error_str:
            logger.error('===== begin qssp error =====\n'
                         '%s===== end qssp error =====' % error_str)

        errmess = []
        if proc.returncode != 0:
            errmess.append('qssp had a non-zero exit state: %i' % proc.returncode)
        if error_str:
            errmess.append('qssp emitted something via stderr')
        if output_str.lower().find('error') != -1:
            errmess.append("the string 'error' appeared in qssp output")

        if errmess:
            os.chdir(old_wd)
            raise QSSPError('''
===== begin qssp input =====
%s===== end qssp input =====
===== begin qssp output =====
%s===== end qssp output =====
===== begin qssp error =====
%s===== end qssp error =====
%s
qssp has been invoked as "%s"'''.lstrip() %
            (input_str, output_str, error_str, '\n'.join(errmess), program))
        
        self.qssp_output = output_str
        self.qssp_error = error_str
        
        os.chdir(old_wd)
        
    def get_traces(self):
       
        fns = self.config.get_output_filenames(self.tempdir)
        traces = []
        for comp, fn in zip(self.config.components, fns):
            data = num.loadtxt(fn, skiprows=1, dtype=num.float)
            nsamples, ntraces = data.shape
            ntraces -= 1
            tmin = data[0,0]
            deltat = (data[-1,0] - data[0,0])/(nsamples-1)
            for itrace in xrange(ntraces):
                tr = trace.Trace( '', '%04i' % itrace, '', comp,
                        tmin=tmin, deltat=deltat, ydata=data[:,itrace+1] )
                
                traces.append(tr)
        
        return traces

    def __del__(self):
        if not self.keep_tmp:
            shutil.rmtree(self.tempdir)

class QSSPGFBuilder(builder.GFBuilder):
    def __init__(self, store_dir, step, block_size=None, tmp=None ):
        self.gfmapping = [
            (MomentTensor( m=symmat6(1,0,0,1,0,0) ), {'un': (0, -1), 'ue': (3, -1), 'uz': (5, -1) }),
            (MomentTensor( m=symmat6(0,0,0,0,1,1) ), {'un': (1, -1), 'ue': (4, -1), 'uz': (6, -1) }),
            (MomentTensor( m=symmat6(0,0,1,0,0,0) ), {'un': (2, -1),                'uz': (7, -1) }),
            (MomentTensor( m=symmat6(0,1,0,0,0,0) ), {'un': (8, -1),                'uz': (9, -1) }),
        ]

        self.step = step
        self.store = store.Store(store_dir, 'w')

        if self.step == 0:
            block_size = (1,1,self.store.meta.ndistances)
        else:
            if block_size is None:
                block_size = (1,1,51)

        if len(self.store.meta.ns) == 2:
            block_size = block_size[1:]

        builder.GFBuilder.__init__(self, self.store.meta, block_size=block_size)
        self.qssp_config = self.store.get_extra('qssp')

        self.tmp = tmp
        if self.tmp is not None:
            util.ensuredir(self.tmp)

        conf = QSSPConfigFull(**self.qssp_config.items())
        util.ensuredir(conf.gf_directory)
        
    def work_block(self, index):
        if len(self.store.meta.ns) == 2:
            (sz, firstx), (sz, lastx), (ns, nx) = \
                    self.get_block_extents(index)

            rz = self.store.meta.receiver_depth
        else:
            (rz, sz, firstx), (rz, sz, lastx), (nr, ns, nx) = \
                    self.get_block_extents(index)

        gf_filename = 'GF_%gkm_%gkm' % (sz/km, rz/km) 

        conf = QSSPConfigFull(**self.qssp_config.items())
        gf_path = os.path.join(conf.gf_directory, '?_' + gf_filename)

        if self.step == 0 and len(glob.glob(gf_path)) == 7:
            logger.info('Skipping step %i, block %i / %i (GF already exists)' %
                (self.step, index+1, self.nblocks))
            return 

        logger.info('Starting step %i, block %i / %i' % 
                (self.step, index+1, self.nblocks))

        runner = QSSPRunner(tmp=self.tmp)
        
        conf.receiver_depth = rz/km
        conf.sampling_interval = 1.0/self.gf_set.sample_rate

        dx = self.gf_set.distance_delta

        if self.step == 0:
            distances = [ firstx ]
        else:
            distances = num.linspace(firstx, 
                    firstx + (nx-1)*dx, nx)

        conf.receivers = [ 
                QSSPReceiver( lat=90-d*cake.m2d, lon=180., tstart=-200.) for 
                    d in distances ]


        gfs = [ QSSPGreen( filename=gf_filename,
            depth= sz/km, calculate=(self.step==0) ) ]

        conf.greens_functions = gfs


        if self.step == 0:
            conf.sources = [ QSSPSourceMT(lat=89.99999, lon=0.0) ]
            runner.run(conf)

        else:
            for mt, gfmap in self.gfmapping[:[3,4][self.gf_set.ncomponents==10]]:
                m = mt.m_up_south_east()

                conf.sources = [ QSSPSourceMT(lat=90-0.001*dx*cake.m2d, lon=0.0, 
                    mrr = m[0,0], mtt = m[1,1], mpp = m[2,2],
                    mrt = m[0,1], mrp = m[0,2], mtp = m[1,2]) ]

                runner.run(conf)
        

                rawtraces = runner.get_traces()

                self.store.lock()

                for itr, tr in enumerate(rawtraces):
                    if tr.channel in gfmap:
                        rec = conf.receivers[itr % nx]
                        x = distances[itr % nx] 
                        ig, factor = gfmap[tr.channel]
                        gf_tr = store.GFTrace(tr.get_ydata() * factor,
                                int(round((tr.tmin+rec.tstart) / tr.deltat)), tr.deltat)

                        if len(self.store.meta.ns) == 2:
                            self.store.put((sz,x,ig), gf_tr)
                        else:
                            self.store.put((rz,sz,x,ig), gf_tr)

                self.store.unlock()
   
        logger.info('Done with step %i, block %i / %i' % 
                (self.step, index+1, self.nblocks))


km = 1000.

class Input(Object):
    modelling_code_id = String.T(default='qssp')
    meta = meta.GFSet.T(default = meta.GFSetTypeA.D(
            id = 'my-qssp-gf-store',
            ncomponents = 10,
            sample_rate = 0.2,
            receiver_depth = 0*km,
            source_depth_min = 10*km,
            source_depth_max = 20*km,
            source_depth_delta = 10*km,
            distance_min = 100*km,
            distance_max = 1000*km,
            distance_delta = 10*km))

    qssp = QSSPConfig.T(default = QSSPConfig.D())

def init(store_dir):
    input_fn = os.path.join(store_dir, 'input')

    if os.path.exists(input_fn):
        raise store.StoreError('File "%s" already exists' % input_fn)

    inp = Input()
    util.ensuredirs(input_fn)
    dump(inp, filename=input_fn)

def __work_block(args):
    store_dir, step, iblock = args
    builder = QSSPGFBuilder(store_dir, step)
    builder.work_block(iblock)

def build(store_dir, force=False, nworkers=None):
    input_fn = os.path.join(store_dir, 'input')
    
    inp = load(filename=input_fn)
    store.Store.create(store_dir, meta=inp.meta, extra={'qssp': inp.qssp}, force=force)

    for step in (0,1):
        builder = QSSPGFBuilder(store_dir, step)
        iblocks = builder.all_block_indices()
         
        if nworkers is None or nworkers > 1:
            p = Pool(nworkers)
            p.map(__work_block, [ (store_dir, step, iblock) for iblock in iblocks ])

        else:
            map(__work_block, [ (store_dir, step, iblock) for iblock in iblocks ])

