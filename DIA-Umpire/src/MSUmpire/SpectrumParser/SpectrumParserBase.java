/*
 * This file is part of DIA-Umpire.
 *
 * DIA-Umpire is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * DIA-Umpire is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with DIA-Umpire. If not, see <https://www.gnu.org/licenses/>.
 *
 */
package MSUmpire.SpectrumParser;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.eclipse.collections.impl.map.mutable.primitive.IntFloatHashMap;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;

/**
 * Parent class of spectrum parser
 * @author Chih-Chiang Tsou
 */
public abstract class SpectrumParserBase {
    //public ScanCollection scanCollection = null;
    public String filename=null;
    public InstrumentParameter parameter=null;
    public int TotalScan=0;
    public int NoCPUs = 4;
    public SpectralDataType.DataType datatype= SpectralDataType.DataType.DDA;
    public DIA_Setting dIA_Setting=null;

    public TreeMap<Integer, Integer> MsLevelList=null;
    protected TreeMap<Float, Integer> ElutionTimeToScanNoMap=null;
//    protected HashMap<Integer, Float> ScanToElutionTime=null;
    protected IntFloatHashMap ScanToElutionTime=null;
    public int NoMS1Scans = 0;
    
    public SpectrumParserBase(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, DIA_Setting dIA_Setting, int NoCPUs){
        this.filename = filename;
        this.ElutionTimeToScanNoMap = new TreeMap<>();
        //scanCollection.ElutionTimeToScanNoMap = new TreeMap<>();
//        this.ScanToElutionTime = new HashMap<>();
        this.ScanToElutionTime = new IntFloatHashMap();
        this.MsLevelList = new TreeMap<>();
        this.dIA_Setting = dIA_Setting;
        this.parameter = parameter;
        this.datatype = datatype;
        this.NoCPUs = NoCPUs;
         if (datatype != SpectralDataType.DataType.DDA) {
            if (datatype != SpectralDataType.DataType.DIA_V_Window) {
                dIA_Setting.DIAWindows = new TreeMap<>();
            }
        }        
    }
    
    public static SpectrumParserBase GetInstance(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, DIA_Setting dIA_Setting, int NoCPUs) throws Exception{
        if(filename.toLowerCase().endsWith("mzxml")){
            return new mzXMLParser(filename, parameter, datatype, dIA_Setting, NoCPUs);
        }
        else if (filename.toLowerCase().endsWith("mzml") || filename.toLowerCase().endsWith("raw"))
            return new mzXMLParser(mzXMLParser.to_mzXML(filename, NoCPUs),
                    parameter, datatype, dIA_Setting, NoCPUs);
        LogManager.getRootLogger().error("The spectral lfile: "+filename +" is not supported.");
        return null;
    }
    
    public float GetMS1CycleTime() {
        return (ElutionTimeToScanNoMap.lastKey() - ElutionTimeToScanNoMap.firstKey()) / NoMS1Scans;
    }

//    public HashMap<Integer, Float> GetScanElutionTimeMap() {
    public IntFloatHashMap GetScanElutionTimeMap() {
        return ScanToElutionTime;
    }
    
    protected void FSElutionIndexWrite() throws IOException {
        try {
            LogManager.getRootLogger().debug("Writing RTidx to file:" + FilenameUtils.removeExtension(filename) + ".RTidxFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".RTidxFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(ElutionTimeToScanNoMap);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }

        try {
            LogManager.getRootLogger().debug("Writing Scanidx to file:" + FilenameUtils.removeExtension(filename) + ".ScanidxFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".ScanidxFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(MsLevelList);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }

        try {
            LogManager.getRootLogger().debug("Writing ScanRT to file:" + FilenameUtils.removeExtension(filename) + ".ScanRTFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".ScanRTFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(ScanToElutionTime);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }

        if (datatype !=  SpectralDataType.DataType.DDA) {
            try {
                LogManager.getRootLogger().debug("Writing DIAWindows to file:" + FilenameUtils.removeExtension(filename) + ".DIAWindowsFS..");
                FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".DIAWindowsFS", false);
                FSTObjectOutput oos = new FSTObjectOutput(fout);
                oos.writeObject(dIA_Setting.DIAWindows);
                oos.close();
                fout.close();
            } catch (Exception ex) {
                LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
        }
         if (datatype ==  SpectralDataType.DataType.WiSIM) {
            try {
                LogManager.getRootLogger().debug("Writing MS1 windows to file:" + FilenameUtils.removeExtension(filename) + ".MS1WindowsFS..");
                FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".MS1WindowsFS", false);
                FSTObjectOutput oos = new FSTObjectOutput(fout);
                oos.writeObject(dIA_Setting.MS1Windows);
                oos.close();
                fout.close();
            } catch (Exception ex) {
                LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
        }
    }
    
    protected boolean FSElutionIndexRead() throws NumberFormatException, FileNotFoundException, IOException {
        if (!new File(FilenameUtils.removeExtension(filename) + ".RTidxFS").exists()) {
            return false;
        }
        if (!new File(FilenameUtils.removeExtension(filename) + ".ScanRTFS").exists()) {
            return false;
        }
        if (!new File(FilenameUtils.removeExtension(filename) + ".ScanidxFS").exists()) {
            return false;
        }
        try {
            LogManager.getRootLogger().debug("Reading RTidx:" + FilenameUtils.removeExtension(filename) + ".RTidxFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".RTidxFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            ElutionTimeToScanNoMap = (TreeMap<Float, Integer>) in.readObject();
            //scanCollection.ElutionTimeToScanNoMap = ElutionTimeToScanNoMap;
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            //i.printStackTrace();
            LogManager.getRootLogger().debug("RTidxFS serialization file failed");
            //LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }

        try {
            LogManager.getRootLogger().debug("Reading ScanRT:" + FilenameUtils.removeExtension(filename) + ".ScanRTFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".ScanRTFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
//            ScanToElutionTime = (HashMap<Integer, Float>) in.readObject();
            ScanToElutionTime = (IntFloatHashMap) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            //i.printStackTrace();
            LogManager.getRootLogger().debug("ScanRTFS serialization file failed");
            //LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            throw new RuntimeException(ex);
        }

        try {
            LogManager.getRootLogger().debug("Reading Scanidx:" + FilenameUtils.removeExtension(filename) + ".ScanidxFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".ScanidxFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            MsLevelList = (TreeMap<Integer, Integer>) in.readObject();
            for (int value : MsLevelList.values()) {
                if (value == 1) {
                    NoMS1Scans++;
                }
            }
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            //i.printStackTrace();
            LogManager.getRootLogger().debug("ScanidxFS serialization file failed");
            //LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }

        if (datatype != SpectralDataType.DataType.DDA) {
            if (!new File(FilenameUtils.removeExtension(filename) + ".DIAWindowsFS").exists()) {
                return false;
            }
            try {
                LogManager.getRootLogger().debug("Reading DIAWindows:" + FilenameUtils.removeExtension(filename) + ".DIAWindowsFS...");
                FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".DIAWindowsFS");
                FSTObjectInput in = new FSTObjectInput(fileIn);
                dIA_Setting.DIAWindows = (TreeMap<XYData, ArrayList<Integer>>) in.readObject();
                in.close();
                fileIn.close();
            } catch (Exception ex) {
                //i.printStackTrace();
                LogManager.getRootLogger().debug("DIAWindowsFS serialization file failed");
                //LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                return false;
            }
        }
        
         if (datatype == SpectralDataType.DataType.WiSIM) {
            if (!new File(FilenameUtils.removeExtension(filename) + ".MS1WindowsFS").exists()) {
                return false;
            }
            try {
                LogManager.getRootLogger().debug("Reading MS1 Windows:" + FilenameUtils.removeExtension(filename) + ".MS1WindowsFS...");
                FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".MS1WindowsFS");
                FSTObjectInput in = new FSTObjectInput(fileIn);
                dIA_Setting.MS1Windows = (TreeMap<XYData, ArrayList<Integer>>) in.readObject();
                in.close();
                fileIn.close();
            } catch (Exception ex) {
                //i.printStackTrace();
                LogManager.getRootLogger().debug("MS1WindowsFS serialization file failed");
                //LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                return false;
            }
        }
        return true;
    }

    public int GetScanNoByRT(float RT) {
        int ScanNo = 0;
        if (RT <= ElutionTimeToScanNoMap.firstKey()) {
            ScanNo = ElutionTimeToScanNoMap.firstEntry().getValue();
        } else if (RT >= ElutionTimeToScanNoMap.lastKey()) {
            ScanNo = ElutionTimeToScanNoMap.lastEntry().getValue();
        } else {
            ScanNo = ElutionTimeToScanNoMap.lowerEntry(RT).getValue();
        }
        return ScanNo;
    }

    protected ScanCollection InitializeScanCollection() {
        ScanCollection scanCollection = new ScanCollection(parameter.Resolution);
        scanCollection.Filename = filename;
        scanCollection.ElutionTimeToScanNoMap = ElutionTimeToScanNoMap;
        return scanCollection;
    }
    
    protected int GetEndScan(float endTime) {
        int EndScanNo;
        if (endTime <= ElutionTimeToScanNoMap.firstKey()) {
            EndScanNo = ElutionTimeToScanNoMap.firstEntry().getValue();
        } else {
            if (endTime >= ElutionTimeToScanNoMap.lastKey()) {
                EndScanNo = ElutionTimeToScanNoMap.lastEntry().getValue();
            } else {
                EndScanNo = ElutionTimeToScanNoMap.higherEntry(endTime).getValue();
            }
        }
        return EndScanNo;
    }

    protected int GetStartScan(float startTime) {
        int StartScanNo;
        if (startTime <= ElutionTimeToScanNoMap.firstKey()) {
            StartScanNo = ElutionTimeToScanNoMap.firstEntry().getValue();
        } else {
            if (startTime >= ElutionTimeToScanNoMap.lastKey()) {
                StartScanNo = ElutionTimeToScanNoMap.lastEntry().getValue();
            } else {
                StartScanNo = ElutionTimeToScanNoMap.lowerEntry(startTime).getValue();
            }
        }
        return StartScanNo;
    }
    //Get all the DIA MS2 scans according to a isolation window range
     public ScanCollection GetScanCollectionDIAMS2(XYData DIAWindow, boolean IncludePeak,float startRT, float endRT){
        if (dIA_Setting == null) {
            LogManager.getRootLogger().error("This is not DIA data" + filename);
            return null;
        }
        return GetScanDIAMS2(DIAWindow, IncludePeak, startRT, endRT);
    }
    public abstract ScanCollection GetScanDIAMS2(XYData DIAWindow, boolean IncludePeak, float startRT, float endRT);
    
    public ScanCollection GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak){
        return GetAllScanCollectionByMSLabel(MS1Included, MS2Included, MS1Peak, MS2Peak, 0f, 999999f);
    }

    public abstract ScanCollection GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak, float startTime, float endTime);
        
    
    //Get all the DIA MS1 scans according to MS1 m/z range, this was only for WiSIM data
    public ScanCollection GetScanCollectionMS1Window(XYData MS1Window, boolean IncludePeak) {
        if (dIA_Setting == null) {
            LogManager.getRootLogger().error("This is not DIA data" + filename);
            return null;
        }
        return GetScanCollectionMS1Window(MS1Window, IncludePeak, 0f, 999999f);
    }
    public abstract ScanCollection GetScanCollectionMS1Window(XYData MS1Window, boolean IncludePeak, float startTime, float endTime) ;
    
}
