###########################################
# written by Mikiko Ito on Feb. 13th, 2023
# Tube Variation Correction (TVC)
# based on "TubeVariationCorrection_DAClinearity.py"
# add a DAC value flexibility
# on iteration method
############################################
import os  #import path, listdir, mkdir
from time import sleep
import struct
import numpy as np
import csv
import matplotlib as m
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from datetime import date

class TVC():
    def __init__(self):
        self.DEBUG = False

        self.LOGfile_indxCurr = 'LOG_indxCurr.csv'
        self.LOGfile_intst = 'LOG_intst.csv'
        # detector parameters
        self.CONST_Npixel_x = 2560        # Tube array가 이동하는 방향
        self.CONST_Npixel_y = 2048        # Tube array 방향
        self.CONST_ActiveArea_x_max = 2350      # collimator로 짤리는 영역 y_max

        self.CONST_SizePixel = 0.124      #[mm]
        # tube array parameters
        self.CONST_Ntube = 7
        self.CONST_PitchTube = 30.0       #[mm]
        self.CONST_SizeStep = 30.0        #[mm]
        self.CONST_SID = 400.0            #[mm]
        self.PositionTube = None    # need to be set by setPosTube()
        self.initVariables()

        # Variables which should be defined in advance
        self.waitingTime = 100
        self.targetIntensity = 3700                     # Need to be defined
        self.limitVariation = 0.05                      # Target +/-10%
        self.list_DAC_LSB = [9.13]*self.CONST_Ntube     # intensity increase per 1 DAC
        self.CONST_Nfiles = 25                          # 11 Dummy + (shot + dummy)X7

    def initVariables(self):
        """initialize status variables"""
        self.status_CALfinished = False
        self.status_running = False
        self.cnt_iter = 0

    def _createDir(self, dirPath):
        if not os.path.exists(dirPath):
            os.mkdir(dirPath)
            print("{dir} was created.".format(dir=dirPath))

    def setPathCALdirectory(self, path_directory):

        self.Directory = path_directory
        self.DirectoryCAL = str(path_directory) + '/cal/'
        self.DirectoryArchive = str(path_directory) + '/archive/'
        self.DirectoryLog = str(path_directory) + '/log/'
        print(path_directory)
        print("self.DirectoryCAL: ", self.DirectoryCAL)
        print("self.DirectoryArchive: ", self.DirectoryArchive)
        print("self.DirectoryLog: ", self.DirectoryLog)
        self._createDirCalArchiveLog()

    def _createDirCalArchiveLog(self):
        """
        # generating Cal directory
        # self.DirectoryCAL = 'D:/Data/tubeVarCAL/test1/'
        # self.DirectoryArchive = 'D:/Data/tubeVarCAL/archive/'
        # self.DirectoryLog = 'D:/Data/tubeVarCAL/log/'
        """
        self._createDir(self.DirectoryCAL)
        self._createDir(self.DirectoryArchive)
        self._createDir(self.DirectoryLog)


    def _calculateTubeCenter(self): #--> output: self.r_tubes
        pitch_idx = self.CONST_PitchTube / self.CONST_SizePixel
        indx_s = int(int(self.CONST_Npixel_y / 2) - int(self.CONST_Ntube / 2) * pitch_idx - 0.5 * pitch_idx * (self.CONST_Ntube % 2 - 1))
        self.r_tubes = [int(indx_s + i*pitch_idx) for i in range(self.CONST_Ntube)]
        # Need to check the order of the tube center according to iTube
        #self.r_tubes.reverse()
        print("self.r_tubes = ", self.r_tubes)

    def setDEBUG_ON(self):
        self.DEBUG = True

    def setDEBUG_OFF(self):
        self.DEBUG = False

    def setTarget(self, val): #val: intensity
        self.Target = int(val)

    def setPosLine(self, val): # val = position of Tube array [mm]   0-150 mm
        self._calculateTubeCenter()  # calculate Tube centers
        self.c_tibes = int(int(val)/self.CONST_SizePixel)
        print(self.r_tubes)
        print(self.c_tibes)

    def setTubeVoltage(self, val): # val = tube voltage (ex, 60: 60kV)
        self.tVol = int(val)

    def setTubeCurrent(self, val): # val = tube current
        self.tCurr = int(val)

    def isCALfinished(self):
        '''
        return a status of CAL process (self.status_CALfinished)
        if the variance of all tube intensities are within limitation from the target intensity,
        :return: True --> CAL will be finished
        else: return False
        '''
        return self.status_CALfinished

    def isProcessRunning(self):
        '''
        return a status of run process (self.status_running)
        if run function is staill running (in process)
        :return: True
        if run function is not running (finished or not started)
        :return False
        '''
        return self.status_running


    def setCurrentIndex(self, list_indxCurr):
        self._addDateIterINFO(list_indxCurr)
        self._writeCSV(self.DirectoryLog + self.LOGfile_indxCurr, list_indxCurr)
        self.list_indxCurr = list_indxCurr[2:]


    def _checkALLFilesSaved(self, directory):
        cnt = 0
        while cnt<self.waitingTime:
            sleep(1)
            cnt += 1
            self.fileList = os.listdir(directory)
            if len(self.fileList) >= self.CONST_Nfiles: return True

        if cnt >= 10: raise Exception("E01: we cannot find {N} files in {dir}".format(N=self.Ntube, dir=self.DirectoryCAL))
        return False

    def _deleteDummyFiles(self):
        if not len(self.fileList) >= self.CONST_Nfiles:
            raise Exception("Warning!!! len(self.fileList) != 25, please check # of files in the CAL directory")
        # select first 25 files in fileList (18 Dummy files + 7 Data files)
        self.fileList = self.fileList[:self.CONST_Nfiles]
        list_datafile = []

        for i, f in enumerate(self.fileList):
            if i<11:
                self._moveFileArchive(f)
            else:
                if i%2 == 0: self._moveFileArchive(f)
                else:
                    list_datafile.append(f)

        self.fileList = list_datafile
        if len(self.fileList) == self.CONST_Ntube: return True
        return False


    def _readData(self, file_path):
        with open(file_path, 'rb') as input_file:  # 16-bit unsigned
            fileContent = input_file.read()
            TotalN_pixels = (len(fileContent)) // 2
            data_1D = struct.unpack("H" * TotalN_pixels, fileContent)  # "H" = 16-bit unsigned
            data_1D = np.asarray(data_1D)
            data_2D = np.reshape(data_1D, (self.Npixel_y, self.Npixel_x))
        return data_2D

    def _showImage(self, data_2D, xx=[], yy=[], style='-r', tit='DATA name', xl='x_index', yl='y_index', saveOption=False):
        plt.subplots(figsize=(14, 10))
        plt.imshow(data_2D)
        if len(xx)>0 : plt.plot(xx, yy, style)
        plt.colorbar()
        plt.clim(2700, 3250)
        plt.xlabel(xl)
        plt.ylabel(yl)
        plt.title(tit)
        if saveOption: plt.savefig('./image2D_output.png')
        plt.show()

    def _showImage_rect(self, data_2D, x_min, x_max, y_min, y_max):
        plt.subplots(figsize=(14, 10))
        plt.imshow(data_2D)
        plt.gca().add_patch(Rectangle((y_min, x_min), int(y_max-y_min), int(x_max-x_min), linewidth=2, edgecolor='r', facecolor='none'))
        plt.colorbar()
        plt.xlabel("x_index")
        plt.ylabel("y_index")
        plt.title("_getIntensity()")
        plt.show()


    def _getIntensity(self, iTube, data2D):
        i_min = int(max(0, (self.r_tubes[iTube] - 50)))
        i_max = int(min((self.r_tubes[iTube] + 50), self.CONST_Npixel_y))
        j_min = int(max(0, (self.c_tibes - 50)))
        j_max = int(min((self.c_tibes + 50), self.CONST_ActiveArea_x_max))

        if j_max<j_min:
            j_min = j_max

        if self.DEBUG: self._showImage_rect(data2D, i_min, i_max, j_min, j_max)

        dataROI = data2D[i_min:i_max, j_min:j_max]
        dataROI.flatten().astype("int")
        dataROI = dataROI[~np.isnan(dataROI)]
        dataROI = dataROI[dataROI>0]
        print("x_min, x_max, y_min, y_max: ", i_min, i_max, j_min, j_max, len(dataROI))
        if len(dataROI)>0: iI = int(np.mean(dataROI))
        else: iI = 0
        if self.DEBUG: print(iTube, "x_min, x_max: ", i_min, i_max, "  -- y_min, y_max: ", j_min, j_max, " --- mean of intensity in ROI: ", iI)
        return iI, i_min, i_max, j_min, j_max

    def _addDateIterINFO(self, list):
        today = date.today()
        d = today.strftime("%Y-%m-%d")
        if len(list) == self.CONST_Ntube:
            list.insert(0, self.cnt_iter)
            list.insert(0, d)
        else:
            raise Exception("E00: The number of Current index is not the same as the number of tubes.")
        return list

    def _writeCSV(self, fileName, list_result):
        with open(fileName, 'a', newline='') as fd:
            writer = csv.writer(fd)
            writer.writerow(list_result)

    def _moveFilesArchive(self):
        for filename in self.fileList:
            self._moveFileArchive(filename)

    def _moveFileArchive(self, filename):
        src = os.path.join(self.DirectoryCAL, filename)
        dst = os.path.join(self.DirectoryArchive, filename)
        os.rename(src, dst)

    def _getListIntensity(self, directory):
        today = date.today()
        d = today.strftime("%Y-%m-%d")
        list_intst = []
        # need to set position of Line(tube array) using setPosLine()
        # Case_01 : check if all files were saved after line-mode exposure
        if self._checkALLFilesSaved(directory):  # len(fileList) == 25: Dummy Files + data Files
            if self._deleteDummyFiles():  # delete Dummy file --> len(fileList) == 7 : TVC starts
                for iTube, f in enumerate(self.fileList):
                    print(iTube, self.DirectoryCAL + f)
                    img = self._readData(directory + f)
                    intensity, x_min, x_max, y_min, y_max = self._getIntensity(iTube, img)
                    if self.DEBUG: self._showImage_rect(img, x_min, x_max, y_min, y_max)
                    list_intst.append(intensity)
                self._addDateIterINFO(list_intst)
                self._writeCSV(self.DirectoryLog + self.LOGfile_intst, list_intst)
                if (self.ArchiveON): self._moveFilesArchive()
            else:
                raise Exception("False from self._deleteDummyFiles(): please check # of files which should be 7 in CAL directory")
        else:
            raise Exception("False from self._checkALLFilesSaved(): please check # of files which should be 25 in CAL directory")

        return list_intst[2:]

    def _calculateNewTarget(self, list_intst):
        list_intensity = [i for i in list_intst]
        list_intensity.sort()
        new_target = sum(list_intensity[1:-1])/len(list_intensity[1:-1])
        print("---- set NEW target to self.Target: ", self.Target, '  -->  ', new_target)
        self.setTarget(new_target)



    def _calculateNewIndxCurr(self):
        self.cnt_iter += 1
        print(self.cnt_iter, "--- Calculate New DAC index for Current ---")

        newIndxCurr = []
        for i, intst in enumerate(self.list_intst):
            newIndx = int(self.list_indxCurr[i] - (intst - self.Target)/self.list_DAC_LSB[i])
            print(i, intst, "new index: ", self.list_indxCurr[i], (intst - self.Target), int((intst - self.Target)/self.list_DAC_LSB[i]), newIndx)
            newIndxCurr.append(newIndx)
        self._addDateIterINFO(newIndxCurr)
        self._writeCSV(self.DirectoryLog + self.LOGfile_indxCurr, newIndxCurr)
        print("NEW DAC index: ", newIndxCurr)
        return newIndxCurr[2:]

    def _calculateVariance(self):
        self._showHistVariance()
        for intst in self.list_intst:
            if abs(self.Target - intst)>self.Target*self.limitVariation: return False
        return True

    def _showHistVariance(self):
        x, y = np.arange(self.CONST_Ntube), self.list_intst
        id = ['Tube_{num}'.format(num = n) for n in range(self.CONST_Ntube)]
        xmin, xmax = -0.8, self.CONST_Ntube - 0.2
        ymin, ymax = self.Target*(1 - self.limitVariation), self.Target*(1 + self.limitVariation)


        plt.subplots(figsize=(10, 8))
        plt.fill([xmin, xmin, xmax, xmax], [ymin, ymax, ymax, ymin], color='lightgray', alpha=0.5)
        plt.bar(x, y)
        plt.hlines(ymax, xmin=xmin, xmax=xmax, colors='r', linestyles='dashdot')
        plt.hlines(self.Target, xmin=xmin, xmax=xmax, colors='r', linestyles='solid')
        plt.hlines(ymin, xmin=xmin, xmax=xmax, colors='r', linestyles='dashdot')
        plt.xlim(xmin, xmax)
        plt.xticks(x, id)
        plt.ylabel('X-ray image intensity')
        plt.title("Variation of X-ray tube output at iter#_{num}".format(num=self.cnt_iter))
        plt.show()

    def checkUniformity(self, directory, MODE_rename):
        self._calculateTubeCenter()
        fileList = os.listdir(directory)
        list_intst = []

        if MODE_rename:
            for i in range(len(fileList)):
                iLine, iTube = i // self.CONST_Ntube, 6-(i % self.CONST_Ntube)
                PosLine = int(iLine * self.CONST_SizeStep)  # [mm]
                self.setPosLine(PosLine)
                fname = str(i) + ".raw"
                data2D = self._readData(directory + fname)

                iIntst, x_min, x_max, y_min, y_max = self._getIntensity(iTube, data2D)
                #if self.DEBUG: self._showImage_rect(data2D, x_min, x_max, y_min, y_max)
                if iIntst>0:
                    list_intst.append(iIntst)
                print(i, iLine, iTube, " PosLine: ", PosLine, " -- iIntst: ", iIntst)

        else:
            for i, f in enumerate(fileList):
                iLine, iTube = i//self.CONST_Ntube, i%self.CONST_Ntube
                PosLine = int(iLine*self.CONST_SizeStep) # [mm]
                self.setPosLine(PosLine)
                data2D = self._readData(directory + f)

                iIntst, x_min, x_max, y_min, y_max = self._getIntensity(iTube, data2D)
                if self.DEBUG: self._showImage_rect(data2D, x_min, x_max, y_min, y_max)
                list_intst.append(iIntst)
                print(i, iLine, iTube, " PosLine: ", PosLine, " -- iIntst: ", iIntst)


        m = int(len(fileList)//self.CONST_Ntube)
        if len(fileList) == self.CONST_Ntube*m :
            list_intst = np.asarray(list_intst)
            list_intst = np.reshape(list_intst, (-1, self.CONST_Ntube))
            list_intst = np.fliplr(list_intst)
            self._showImage(list_intst, tit="Uniformity Check", xl='tube Number', yl= 'Step Number')

    def getDACLinearity(self):
        directory = "D:/Data/Calibration_tube/DAC/"
        folderNames = ["DAC_N100", "DAC_N80", "DAC_N60", "DAC_N40", "DAC_N20", "DAC_0",
                       "DAC_P20", "DAC_P40", "DAC_P60", "DAC_P80", "DAC_P100"]
        DAC_min = -100
        DAC_max = 100
        DAC_interval = 20
        i = 0
        for dac in range(DAC_min, DAC_max+1, DAC_interval):
            DirectoryDAC = directory + folderNames[i] + '/'
            print(dac, DirectoryDAC)
            self.list_intst = self._getListIntensity(DirectoryDAC)
            i+=1

    def _overWriteCSV(self, fileName, list_result):
        with open(fileName, 'w') as fd:
            writer = csv.writer(fd)
            writer.writerow(list_result)

    def _readCSV(self, file_path):
        with open(file_path, 'r') as fd:
            list_data = list(csv.reader(fd))
            list_data = [int(i) for i in list_data[0]]
        return list_data

    def saveDACindex(self, list_IndxCurr):
        file_path = self.Directory + '/list_indxDAC.csv'
        self._overWriteCSV(file_path, list_IndxCurr)

    def readDACindex(self):
        file_path = self.Directory + '/list_indxDAC.csv'
        return self._readCSV(file_path)


    def run(self, n_iter):
        self.ArchiveON = False
        self.status_running = True
        #self.list_indxCurr = self._setCurrentIndex(list_indxCurr)
        self.list_intst = self._getListIntensity(self.DirectoryCAL)
        if int(n_iter) == 0: self._calculateNewTarget(self.list_intst)
        self.status_CALfinished = self._calculateVariance()
        print("status_CALfinished: ", self.status_CALfinished)
        print("--- list of intensity: ", self.list_intst, '\n--- list of new DAC index for TubeCurr.: ', self.list_indxCurr)
        list_newIndxCurr = self._calculateNewIndxCurr()
        self.status_running = False
        print("Get NEW DAC index: ", list_newIndxCurr)
        return list_newIndxCurr

if __name__ == "__main__":

    #input variables
    LinePosition, tubeVolt, tubeCurr = 150.0, 60, 0.5
    cnt_iter = 0
    list_indxCurr = [0, 0, 0, 0, 0, 0, 0]

    tvc = TVC()

    tvc.setPathCALdirectory('D:/Data/Calibration_tube')  # generate cal, archive, log directories
    tvc.setPosLine(LinePosition)
    tvc.setTubeVoltage(tubeVolt)
    tvc.setTubeCurrent(tubeCurr)
    tvc.saveDACindex(list_indxCurr)  # 처음에 CAL 시작할 때 시작 DAC index list 파일 생성
    while not tvc.isCALfinished() or cnt_iter<3:
        print("\n  --- new iteration: {n_iter} --- ".format(n_iter=cnt_iter))
        list_indxCurr = tvc.readDACindex()
        print("list_indxDAC: ", list_indxCurr)
        tvc.setCurrentIndex(list_indxCurr)
        list_newDACindx = tvc.run(cnt_iter)  # get new DAC index
        tvc.saveDACindex(list_newDACindx)
        cnt_iter+=1
