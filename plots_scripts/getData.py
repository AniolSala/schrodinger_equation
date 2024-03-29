import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import warnings
import os


class DataParser:
    """
    Class thought to manipulate data with a csv type-like:

    #time-value  #space-value    #sol-1      #sol-2     #etc
      t0            x0            y0,0          z0,0
      t0            x1            y0,1          z1,1
      ...           ...           ...          ...
      t0            xN-1          y0,N-1        z0,N-1
      t1            x0            y1,0          z1,0
      t1            x1            y1,1          z1,1
      t1            x2            y1,2          z1,2
      ...           ...           ...          ...
      t1            xN-1          y1,N-1        z1,N-1

    where y and z are some solutions for a given time (there can be
    an arbitrary number of solutions).

    Each subgroup of the solutions corresponding to the same time value
    is defined as a `frame`. Thus, each frame has N values.

    """

    def __init__(self, filename, dataDir="data", verbose=False):
        self._errorMsg = ""
        self.verbose = verbose
        self._dataDirName = dataDir

        self.filePath = filename
        self.data = self.getData()
        self.params = self.getParams()

    @property
    def filePath(self):
        return self._filename

    @filePath.setter
    def filePath(self, filename):
        if os.path.isfile(filename):
            self._filename = os.path.abspath(filename)
            self.printmsg("Found file {}".format(
                self.getShortDir(self._filename)))
        else:
            self.printmsg(
                "File not found in current directory. Searching for a `data` directory...", end=" ")
            dataDir = self.getDataDir()
            filepath = os.path.join(dataDir, filename)
            if not os.path.isfile(filepath):
                self._filename = None
                raise ValueError(f"File `{filename}` not found")
            else:
                self.printmsg("Done.")
                self.printmsg(f"Working with file {self.getShortDir(filepath)}")
                self._filename = filepath

    def getData(self):
        """
        Return the data written on the requested file in a numpy
        array.in a numpy array.

        """
        return np.loadtxt(self._filename).T

    def getParams(self):
        """
        Returns the grid parameters of the solution:
            dt = Time interval
            dx = Space interval
            nt = Number of time steps written in the solution
            nx = Number of space steps written in the solution

        NOTE: This parameters are determined from the solution
        written in the file. If not all the steps used to calcu-
        late the solution are written on the file containing the
        data, then dt and dx may be different of the time and space
        intervals used to calculate the solution.

        """

        self.printmsg("Getting parameters...", end=' ')
        ncols = len(self.data[0, :])

        # nx
        nx = 0
        for t in self.data[0, :]:
            if t != self.data[0, 0]:
                break
            nx += 1

        if ncols % nx != 0:
            errormsg = """Invalid number of columns (each frame has {} columns
            but the file has {} columns)""".format(nx, ncols)
            raise AttributeError(errormsg)

        # nt
        nt = ncols // nx

        # dt and dx
        dt, dx = self.data[0, nx] -\
            self.data[0, 0], self.data[1, 1] - self.data[1, 0]
        self.printmsg("Done")
        return [nt, nx, dt, dx]

    def getPlotsList(self):
        """
        Returns a list, each element of this list being the solution
        for a given time (`frame`). Thus, for example, if we define the
        variable

            plotsList = DataParser("filename").getPlotsList()

        then plotsList[0] will be a numpy array containing the solution
        for time t = 0.

        """

        if not self.params:
            return []
        nt, nx, _, _ = self.params
        return [self.data[:, i * nx:(i + 1) * nx] for i in range(nt)]

    def getAreas(self, col=None, freq=1):
        """
        Returns the area of the solution for each frame using the
        simps method of scipy.integrate.

        By default it uses the values of the last column of the data,
        unluess `col` is specified.

        To compute the area all values are taken unless `freq` is
        specified.

        """

        lastCol = data.shape[0] - 1 if col is None else col

        areas = [simps(pl[lastCol, ::freq], pl[lastCol, ::freq])
                 for pl in self.getPlotsList()]
        return np.array(areas)

    def getDataDir(self):
        """
        Useful when executing script from `build` directory or `plots`
        directory.

        When executing the script in a different directory that does not
        contain the requested data, search for the file in a directory
        called self._dataDirName (by dafaule, self._dataDirName = "data")

        When no data directory is found, raise warning and return the
        current directory.

        """

        currentDir = os.getcwd()
        # Search for `data` directory in previous directories
        prevDirs = [
            os.path.join(currentDir, ".." * 0, self._dataDirName),
            os.path.join(currentDir, ".." * 1, self._dataDirName),
            os.path.join(currentDir, ".." * 2, self._dataDirName),
        ]
        # Search for `data` directory in subdirectories
        subDirs = [os.path.join(d, self._dataDirName)
                   for d in os.listdir(currentDir) if os.path.isdir(d)]
        # Check if we found some `data` directory
        for dty in prevDirs + subDirs:
            if os.path.isdir(dty):
                return os.path.abspath(dty)

        warnings.warn("Couldn't find directory `data`", Warning)
        return currentDir

    def getShortDir(self, dty):
        """Return abbreviated absolute path. Just for readability"""

        absPath = dty if dty == os.path.abspath(dty) else os.path.abspath(dty)
        pDirList = [d for d in absPath.split(os.sep) if d]  # Empty dirs no
        if len(pDirList) >= 5:
            abbName = os.sep.join(pDirList[:3]) + os.sep + "..." + os.sep +\
                os.sep.join(pDirList[-2:])
        else:
            abbName = absPath
        return abbName

    def printmsg(self, msg, end='\n'):
        if self.verbose:
            print(msg, end=end)


if __name__ == "__main__":
    pdata = DataParser("schrEq1D_normmode_1_box.txt", True)
    # params = pdata.getParams()
