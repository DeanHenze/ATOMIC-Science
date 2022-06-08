# IORG is a python module for calculating the Organization Index 
# (Tompkins and Semie, 2017) as used in RCEMIP (Wing et al., 2020).
# see https://myweb.fsu.edu/awing/iorgREADME.html


def ecdf(x_inputs):
    """
    Calculates the empirical cumulative distribution function of x_inputs
    Arguments
        x_inputs: 1D array, that the ECDF is calculated from
    Returns
        y_values: 1D array, the ECDF
        x_values: 1D array, the values ECDF is evaluated at using x_inputs
    """

    import numpy as np

    x_sorted = np.sort(x_inputs)
    x_values = np.linspace(start=min(x_sorted),
                           stop=max(x_sorted),
                           num=len(x_sorted))

    y_values = np.zeros(np.shape(x_values))
    for i in range(len(x_values)):
        temp  = x_inputs[x_inputs <= x_values[i]]
        value = np.float64(len(temp))/np.float64(len(x_inputs))
        y_values[i] = value

    return(x_values,y_values)

def iorg_crm(w,wthreshold=0):
    """
    Calculates the Organization Index for a cartesian domain
    Arguments
        w: 2-D array of demensions (y,x) of OLR or 500hPa vertical
            acceleration RCEMIP (Wing et al., 2018) uses OLR
        wthreshold: integer, 0 for w being an array of OLR, 1 for
            vertical acceleration
    Returns
        iorg:
    """

    from sklearn.neighbors import NearestNeighbors
    from scipy.ndimage import measurements
    from scipy.ndimage import label
    import numpy as np

    # create a masked array of convective entities
    if wthreshold != 1:
        if wthreshold != 0:
            print('incorrect threshold flag, using default (OLR<173)')
        wmask = (w<173)*1
    if wthreshold == 1:
        wmask = (w>0.5)*1

    # duplicates domain in all directions to account for convection on borders
    wmaskbig = np.concatenate([wmask,wmask,wmask],axis=0)
    wmaskbig = np.concatenate([wmaskbig,wmaskbig,wmaskbig],axis=1)

    # four point convective connection
    sss = [[0,1,0],
           [1,1,1],
           [0,1,0]]

    # finds connected convective entities, the number of clusters,
      # and the centroid of each cluster
    Conn = label(wmaskbig,structure=sss)[0]
    nentities = label(wmaskbig)[1]
    centroids = measurements.center_of_mass(wmaskbig,Conn,range(1,nentities+1))

    nnd,IDX,num = [],[],0
    if nentities > 1:
        # finds the nearest neighbor of each convective cluster
        classifier = NearestNeighbors(n_neighbors=1)
        for i in range(0,nentities):
            if centroids[i][0] >= np.shape(w)[0] \
               and centroids[i][0] <= (np.shape(w)[0]*2)-1 \
               and centroids[i][1] >= np.shape(w)[1] \
               and centroids[i][1] <= (np.shape(w)[1]*2)-1:
                num += 1
                classifier.fit(np.array(centroids[0:i]+centroids[i+1:]))
                m,n = classifier.kneighbors(np.reshape(centroids[i],(1,2)))
                IDX.append(n)
                nnd.append(m)

        if len(nnd) > 1:
            IDX = np.squeeze(np.array(IDX))
            nnd = np.squeeze(np.array(nnd))
        if len(nnd) <= 1:
            IDX = np.array(IDX)
            nnd = np.array(nnd)

        if len(nnd) > 0:
            # calculates ECDF of the nearest neighbors idstances
            x_values_ecdf,y_values_ecdf = ecdf(nnd)

            # calculates the poisson distribution of w
            lam = num/(np.shape(w)[0]*np.shape(w)[1])
            dd = np.array(x_values_ecdf)
            poisson = 1-np.exp(-1*lam*np.pi*dd**2)
            cdf_theory = poisson

            # calculates the area under the plot of the poisson vs ECDF
            iorg = np.trapz(y_values_ecdf,poisson)

        else:
            cdf_nnd = 0
            d = 0
            cdf_theory = 0
            iorg = 0
            y_values_ecdf,x_values_ecdf,poisson,iorg = 0,0,0,0

    else:
        cdf_nnd = 0
        d = 0
        cdf_theory = 0
        iorg = 0
        y_values_ecdf,x_values_ecdf,poisson,iorg = 0,0,0,0

    return(iorg)

def iorg_gcm(w,wthreshold=0):
    """
    Calculates the Organization Index for a speherical domain
    Arguments
        w: 2-D array of demensions (y,x) of OLR or 500hPa vertical
            acceleration RCEMIP (Wing et al., 2018) uses OLR
        wthreshold: integer, 0 for w being an array of OLR, 1 for
            vertical acceleration
    Returns
        iorg:
    """

    from sklearn.neighbors import NearestNeighbors
    from scipy.ndimage import measurements
    from scipy.ndimage import label
    import numpy as np

    # create a masked array of convective entities
    if wthreshold != 1:
        if wthreshold != 0:
            print('incorrect threshold flag, using default (OLR<173)')
        wmask = (w<173)*1
    if wthreshold == 1:
        wmask = (w>0.5)*1

    # duplicates domain in all directions to account for convection on borders
    wmaskbig = np.concatenate([wmask,wmask,wmask],axis=1)

    # four point convective connection
    sss = [[0,1,0],
           [1,1,1],
           [0,1,0]]

    # finds connected convective entities, the number of clusters,
      # and the centroid of each cluster
    Conn = label(wmaskbig,structure=sss)[0]
    nentities = label(wmaskbig)[1]
    centroids = measurements.center_of_mass(wmaskbig,Conn,range(1,nentities+1))

    nnd,IDX,num = [],[],0
    if nentities > 1:
        # finds the nearest neighbor of each convective cluster
        classifier = NearestNeighbors(n_neighbors=1)
        for i in range(0,nentities):
            if centroids[i][1] >= np.shape(w)[1] \
               and centroids[i][1] <= (np.shape(w)[1]*2)-1:
                num += 1
                classifier.fit(np.array(centroids[0:i]+centroids[i+1:]))
                m,n = classifier.kneighbors(np.reshape(centroids[i],(1,2)))
                IDX.append(n)
                nnd.append(m)

        if len(nnd) > 1:
            IDX = np.squeeze(np.array(IDX))
            nnd = np.squeeze(np.array(nnd))
        if len(nnd) <= 1:
            IDX = np.array(IDX)
            nnd = np.array(nnd)

        if len(nnd) > 0:
            # calculates ECDF of the nearest neighbors idstances
            x_values_ecdf,y_values_ecdf = ecdf(nnd)

            # calculates the poisson distribution of w
            lam = num/(np.shape(w)[0]*np.shape(w)[1])
            dd = np.array(x_values_ecdf)
            poisson = 1-np.exp(-1*lam*np.pi*dd**2)
            cdf_theory = poisson

            # calculates the area under the plot of the poisson vs ECDF
            iorg = np.trapz(y_values_ecdf,poisson)

        else:
            cdf_nnd = 0
            d = 0
            cdf_theory = 0
            iorg = 0
            y_values_ecdf,x_values_ecdf,poisson,iorg = 0,0,0,0
    else:
        cdf_nnd = 0
        d = 0
        cdf_theory = 0
        iorg = 0
        y_values_ecdf,x_values_ecdf,poisson,iorg = 0,0,0,0

    return(iorg)

def calc_iorg(w_f,geometry,threshold=0,tbeg=-999,tend=-999):
    """
    Loops through a time range to calculate Iorg
    Arguments
        w_f: 3-D array (time, y, x) of the convection data
        geometry: string to determine which IORG calculation to use,
            acceptable input is 'cartesian' or 'spherical'
        threshold: optional integer (default=0) determining which
            type of convection data is being used, 0 for olr or 1 for
            vertical velocity
        tbeg: optional integer (default=-999) determining which time
            integer to start at, -999 uses 25 days from the end
        tend: optional integer (default=-999) determining which time
            integer to end at, -999 uses the end of the data
    Returns
        iorg_f: 1-D array of Iorg data from time index tbeg to tend
    """

    import numpy as np

    if geometry != 'cartesian' and geometry != 'spherical':
        print('inapproriate geometry entry, choices are ``cartesian`` or ``spherical``')
        return(np.nan)
    else:
        if tend == -999:
            tend = np.shape(w_f)[0]
        if tbeg == -999:
            tbeg = 1800

        # loop through time range
        cdf_nnd_f,d_f,cdf_theory_f,iorg_f = [],[],[],[]
        for v in range(tbeg,tend):
            if v%24 == 0:
                print('   processing day %d of %i\r'%(v/24,int(tend/24)), end='')
            # calculate Iorg
            if geometry == 'cartesian':
                result = iorg_crm(w_f[v,:,:],threshold)
            if geometry == 'spherical':
                result = iorg_gcm(w_f[v,:,:],threshold)
            iorg_f.append(result)
        iorg_f = np.array(iorg_f)

        return(iorg_f)

###

