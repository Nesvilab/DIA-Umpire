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
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import java.io.*;
import java.util.ArrayList;
import java.util.concurrent.*;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.logging.log4j.LogManager;
import org.xml.sax.SAXException;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

/*
 * Using multi-threading to parse mzXML files Load all scans at once in the
 * constructor To change this template, choose Tools | Templates and open the
 * template in the editor.
 * modified from UCSD mzML spectrum parser
 */
/**
 *
 * @author Chih-Chiang Tsou
 */
public final class mzMLParser extends SpectrumParserBase{

    public ScanCollection scanCollection = null;
    
    public mzMLParser(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, DIA_Setting dIA_Setting, int NoCPUs) throws FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException {
        super(filename, parameter, datatype, dIA_Setting, NoCPUs);        
        scanCollection = InitializeScanCollection(); 
        if(!FSElutionIndexRead()){
            ParseAllScans();
        }
    }
    
    //Check if the the file has been parsed
    private void CheckStatus() {
        if (scanCollection.ScanHashMap.isEmpty()) {
            ParseAllScans();
        }
    }
    
    private void ParseAllScans() {
        LogManager.getRootLogger().info("Using MzMLUnmarshaller....");
        MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(new File(filename));
        MzMLObjectIterator<Spectrum> itr = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);        
        ArrayList<mzMLSpecConverter> ScanList=new ArrayList<>();        
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        LogManager.getRootLogger().info("Starting to convert from jmzspec");
        while (itr.hasNext()) {
            Spectrum jmzSpec = itr.next();
            mzMLSpecConverter converter=new mzMLSpecConverter(jmzSpec,parameter);
            executorPool.execute(converter);    
            ScanList.add(converter);
        }
        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            LogManager.getRootLogger().info("interrupted..");
        }
        
        LogManager.getRootLogger().info("Building elution time and scan index");
        for (mzMLSpecConverter converter : ScanList) {
            ScanData spec = converter.spec;
            scanCollection.AddScan(spec);
            ElutionTimeToScanNoMap.put(spec.RetentionTime, spec.ScanNum);
            ScanToElutionTime.put(spec.ScanNum, spec.RetentionTime);
            MsLevelList.put(spec.ScanNum, spec.MsLevel);

            if(spec.MsLevel==1){
                 NoMS1Scans++;     
            }
            
            if (datatype != SpectralDataType.DataType.DDA && spec.MsLevel == 2){
                if (datatype == SpectralDataType.DataType.DIA_V_Window) {
                    for (XYData window : dIA_Setting.DIAWindows.keySet()) {
                        if (window.getX() <= spec.isolationWindowTargetMz && window.getY() >= spec.isolationWindowTargetMz) {
                            dIA_Setting.DIAWindows.get(window).add(spec.ScanNum);
                            break;
                        }
                    }
                }                
                else {
                    if (spec.isolationWindowLoffset > 0f && spec.isolationWindowRoffset > 0f) {
                        if (!dIA_Setting.DIAWindows.containsKey(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset))) {
                            ArrayList<Integer> scanList2 = new ArrayList<>();
                            dIA_Setting.DIAWindows.put(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset), scanList2);
                        }
                        dIA_Setting.DIAWindows.get(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset)).add(spec.ScanNum);
                    } else {
                        if (datatype == SpectralDataType.DataType.DIA_F_Window) {
                            spec.isolationWindowLoffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.2f;
                            spec.isolationWindowRoffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.8f;
                            if (!dIA_Setting.DIAWindows.containsKey(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset))) {
                                ArrayList<Integer> scanList2 = new ArrayList<>();
                                dIA_Setting.DIAWindows.put(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset), scanList2);
                            }
                            dIA_Setting.DIAWindows.get(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset)).add(spec.ScanNum);
                        }
                    }
                }
            }
        } 
        try {        
            FSElutionIndexWrite();
        } catch (IOException ex) {            
        }
    }
        
    @Override
    public ScanCollection GetScanDIAMS2(XYData DIAWindow, boolean IncludePeak, float startTime, float endTime) {
        if (dIA_Setting == null) {
            LogManager.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
        CheckStatus();
        ScanCollection swathScanCollection = new ScanCollection(parameter.Resolution);

        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);
        EndScanNo = GetEndScan(endTime);
        for (int scannum : dIA_Setting.DIAWindows.get(DIAWindow)) {
            if (scannum >= StartScanNo && scannum <= EndScanNo) {
                ScanData scan = scanCollection.ScanHashMap.get(scannum);
                swathScanCollection.AddScan(scan);
                swathScanCollection.ElutionTimeToScanNoMap.put(scan.RetentionTime, scan.ScanNum);
            }
        }
        return swathScanCollection;
    }

    @Override
    public ScanCollection GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak, float startTime, float endTime) {
       CheckStatus();
        ScanCollection newscanCollection = InitializeScanCollection();
        LogManager.getRootLogger().debug("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB (" + NoCPUs + " threads)");

        ArrayList<Integer> IncludedScans = new ArrayList<>();
        
        for(int ScanNum : MsLevelList.keySet()){
            if(MsLevelList.get(ScanNum)==1 && MS1Included){
                IncludedScans.add(ScanNum);
            }
            if(MsLevelList.get(ScanNum)==2 && MS2Included){
                IncludedScans.add(ScanNum);
            }
        }
        
        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
        
        for(int scannum : IncludedScans){
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                ScanData scan=scanCollection.ScanHashMap.get(scannum);                
                newscanCollection.AddScan(scan);
                newscanCollection.ElutionTimeToScanNoMap.put(scan.RetentionTime, scan.ScanNum);
            }
        }        
        System.gc();
        LogManager.getRootLogger().debug("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");
        return newscanCollection;
    }
    
    @Override
    public ScanCollection GetScanCollectionMS1Window(XYData MS1Window, boolean IncludePeak, float startTime, float endTime)  {
        if (dIA_Setting == null) {
            LogManager.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
        CheckStatus();
        ScanCollection MS1WindowScanCollection = new ScanCollection(parameter.Resolution);
       
        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
        ArrayList<Integer> IncludedScans=new ArrayList<>();
        for(int scannum : dIA_Setting.MS1Windows.get(MS1Window)){
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                IncludedScans.add(scannum);
            }
        }
                
        for(int scannum : IncludedScans){
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                ScanData scan=scanCollection.ScanHashMap.get(scannum);                
                MS1WindowScanCollection.AddScan(scan);
                MS1WindowScanCollection.ElutionTimeToScanNoMap.put(scan.RetentionTime, scan.ScanNum);
            }
        }        
        System.gc();
        LogManager.getRootLogger().debug("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");
        return MS1WindowScanCollection;
    }
    
}
