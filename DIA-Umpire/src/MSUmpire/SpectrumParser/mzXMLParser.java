/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.SpectrumParser;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.*;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.eclipse.collections.impl.list.mutable.primitive.IntArrayList;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;
import umich.ms.fileio.filetypes.mzml.MZMLFile;

/* * 
 * mzXML parser
 */
/**
 *
 * @author Chih-Chiang Tsou
 */
public final class mzXMLParser  extends SpectrumParserBase{

    public TreeMap<Integer, Long> ScanIndex=null;
    public mzXMLParser(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, DIA_Setting dIA_Setting, int NoCPUs) throws Exception {
        super(filename,parameter,datatype,dIA_Setting,NoCPUs);
        ReadElutionAndScanIndex();
    }

    //Parser elution time index and scan index and save them as binary files
    private void ReadElutionAndScanIndex() throws Exception {
        if (!FSScanPosRead()) {
            InferOffsets();
//            ParseIndex();
            WriteIndexSerialization();
        }
        if (!FSElutionIndexRead()) {
            ParseElutionIndex();
            FSElutionIndexWrite();
        }
    }

     //Wirte seralization file for scan index
    private void WriteIndexSerialization() {
        FSScanPosWrite();
    }
    
     //Wirte seralization file for scan index
    private void FSScanPosWrite() {
        try {
            Logger.getRootLogger().debug("Writing ScanPos to file:" + FilenameUtils.removeExtension(filename) + ".ScanPosFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".ScanPosFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(ScanIndex);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    
    //Read seralization file for scan index
    private boolean FSScanPosRead() {
        if (!new File(FilenameUtils.removeExtension(filename) + ".ScanPosFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().debug("Reading ScanPos:" + FilenameUtils.removeExtension(filename) + ".ScanPosFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".ScanPosFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            ScanIndex = (TreeMap<Integer, Long>) in.readObject();
            TotalScan = ScanIndex.size();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            Logger.getRootLogger().debug("ScanIndex serialization file failed");
            //Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    /**
     * Infer scan tag offsets.
     * @throws IOException
     */
    private void InferOffsets() throws IOException {
        this.ScanIndex = new TreeMap<>();
        try (final InputStream is = Files.newInputStream(Paths.get(this.filename));) {
            final LongArrayList offsets = Infer_mzXML_offsets.infer_mzXML_scanno(is);
            this.TotalScan = offsets.size();
            for (int i=0;i<offsets.size()-1;++i){
                this.ScanIndex.put(i+1,offsets.get(i));
            }
            this.ScanIndex.put(Integer.MAX_VALUE, offsets.getLast());
        }
    }

    //Parse scan index at the bottom of mzXML file
    private void ParseIndex() throws FileNotFoundException, IOException {
        TotalScan = 0;
        ScanIndex = new TreeMap<>();
        try (RandomAccessFile fileHandler = new RandomAccessFile(filename, "r")) {
            StringBuilder sb = new StringBuilder();

            String CurrentLine = "";
            long currentLastPt = fileHandler.length() - 1;
            boolean indexexist=false;
            int linecount=0;
            while (!(CurrentLine.trim().startsWith("<index name=") | CurrentLine.trim().startsWith("</msRun>"))) {
                //Read backward
                for (long filePointer = currentLastPt; filePointer != -1; filePointer--) {
                    fileHandler.seek(filePointer);
                    int readByte = fileHandler.readByte();
                    if (readByte == 0xA) {
                        if (filePointer == currentLastPt) {
                            continue;
                        } else {
                            currentLastPt = filePointer;
                            break;
                        }
                    } else if (readByte == 0xD) {
                        if (filePointer == currentLastPt - 1) {
                            continue;
                        } else {
                            currentLastPt = filePointer;
                            break;
                        }
                    }
                    sb.append((char) readByte);
                }
                linecount++;
                CurrentLine = sb.reverse().toString();
                sb = new StringBuilder();

                if (CurrentLine.trim().startsWith("</index>")) {
                    indexexist = true;
                }

                if (!indexexist && linecount > 10) {
                    fileHandler.close();
                    Logger.getRootLogger().debug("File : " + filename + " doesn't have index. the processing will stop.");
                    System.exit(1);
                }

                if (CurrentLine.trim().startsWith("<offset id")) {
                    int scanNo = Integer.parseInt(CurrentLine.substring(CurrentLine.indexOf("<offset id=\"") + 12).split("\"")[0]);
                    long index = (long) Long.parseLong(CurrentLine.substring(CurrentLine.indexOf(">") + 1, CurrentLine.indexOf("</offset>")));
                    if (index < 0) {
                        index = index + 2147483647l + 2147483648l;
                    }
                    if (ScanIndex.containsKey(scanNo + 1) && ScanIndex.get(scanNo + 1) == index) {
                        Logger.getRootLogger().debug("File : " + filename + " index is not correct, ScanNo:" + scanNo + " and " + scanNo + 1 + " have same index");
                        Logger.getRootLogger().debug("Please use indexmzXML from  TPP package to fix incorrect index of the mzXML file.");
                        Logger.getRootLogger().debug("command: indexmzXML filename.mzXML");
                        System.exit(1);
                    }
                    ScanIndex.put(scanNo, index);
                } else if (CurrentLine.trim().startsWith("<indexOffset>")) {
                    long IndexEnd = (long) Long.parseLong(CurrentLine.substring(CurrentLine.indexOf("<indexOffset>") + 13, CurrentLine.indexOf("</indexOffset>")));
                    if (IndexEnd < 0) {
                        IndexEnd = IndexEnd + 2147483647l + 2147483648l;
                    }
                    ScanIndex.put(Integer.MAX_VALUE, IndexEnd);
                }
            }
            TotalScan = ScanIndex.size();
            sb = null;
            fileHandler.close();
        }
    }

     //Parse elution time-scan number mapping
    //For DIA data, isolation window ranges are parsed in this method
    private void ParseElutionIndex() throws Exception {
        
        if(ScanIndex==null | ScanIndex.isEmpty()){
            return;
        }

        try (RandomAccessFile fileHandler = new RandomAccessFile(filename, "r")) {
            Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();
            Long currentIdx = iter.next().getValue();
            while (iter.hasNext()) {
                long startposition = currentIdx;
                long nexposition = iter.next().getValue();
                currentIdx = nexposition;
                fileHandler.seek(startposition);
                
                byte[] bufr = new byte[(int) (nexposition - startposition)];
                fileHandler.readFully(bufr, 0, (int) (nexposition - startposition));

                String temp = new String(bufr);

                float rt = 0f;
                int scanno = 0;
                int mslevel = 0;
                //float precursorF=0f;
                if (!temp.contains("<scan")) {
                    fileHandler.close();
                    return;
                }

                if (temp.contains("<scan num=") && (temp.contains("retentionTime=\"PT"))) {
                    String substr = temp.substring(temp.indexOf("<scan num=") + 11);
                    scanno = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    
                    rt = Float.parseFloat(temp.substring(temp.indexOf("retentionTime=\"PT") + 17).split("S\"")[0]);
                    rt=rt/60f;
                    mslevel = Integer.parseInt(temp.substring(temp.indexOf("msLevel=") + 9, temp.indexOf("msLevel=") + 10));
                    if (temp.contains("scanType=\"calibration\"")) {
                        mslevel = -1;
                    }
                    if (mslevel == 1) {
                        NoMS1Scans++;                        
                        if (temp.contains("scanType=\"SIM\"") && datatype == SpectralDataType.DataType.WiSIM) {
                            int startidx = temp.indexOf("lowMz=\"") + 7;
                            int stopidx = startidx + 1;
                            for (int i = startidx + 1; i < temp.length(); i++) {
                                if (temp.charAt(i) == '\"') {
                                    stopidx = i;
                                    break;
                                }
                            }
                            float lowmz = Float.parseFloat(temp.substring(startidx, stopidx));
                            startidx = temp.indexOf("highMz=\"") + 8;
                            stopidx = startidx + 1;
                            for (int i = startidx + 1; i < temp.length(); i++) {
                                if (temp.charAt(i) == '\"') {
                                    stopidx = i;
                                    break;
                                }
                            }
                            float highmz = Float.parseFloat(temp.substring(startidx, stopidx));
                            for (XYData MS1win : dIA_Setting.MS1Windows.keySet()) {
                                if (MS1win.getX() <= lowmz && MS1win.getY() >= highmz) {
                                    dIA_Setting.MS1Windows.get(MS1win).add(scanno);
                                }
                            }
                        }                        
                    }
                    //If it is DIA data, parse isolation window ranges 
                    if (datatype != SpectralDataType.DataType.DDA) {
                        if (mslevel == 2) {
                            if (datatype == SpectralDataType.DataType.MSX) {
                                substr = temp;
                                while (substr.contains("</precursorMz>")) {
                                    int stopidx = substr.indexOf("</precursorMz>");
                                    int startidx = 0;
                                    for (int i = stopidx; i > 0; i--) {
                                        if (substr.charAt(i) == '>') {
                                            startidx = i + 1;
                                            break;
                                        }
                                    }
                                    float precursormz = Float.parseFloat(substr.substring(startidx, stopidx));

                                    startidx = substr.indexOf("windowWideness=\"") + 16;
                                    stopidx = startidx + 1;
                                    for (int i = startidx + 1; i < substr.length(); i++) {
                                        if (substr.charAt(i) == '\"') {
                                            stopidx = i;
                                            break;
                                        }
                                    }
                                    float windowwideness = Float.parseFloat(substr.substring(startidx, stopidx));                                    
                                    //Assuming the precursor m/z is at the center of isolation window, it's for Thermo MSX data
                                    float Loffset = windowwideness / 2f;
                                    float Roffset = windowwideness / 2f;

                                    if (!dIA_Setting.DIAWindows.containsKey(new XYData(precursormz - Loffset, precursormz + Roffset))) {
                                        ArrayList<Integer> scanList = new ArrayList<>();
                                        dIA_Setting.DIAWindows.put(new XYData(precursormz - Loffset, precursormz + Roffset), scanList);
                                    }
                                    dIA_Setting.DIAWindows.get(new XYData(precursormz - Loffset, precursormz + Roffset)).add(scanno);
                                    substr = substr.substring(substr.indexOf("</precursorMz>") + 14);
                                }
                            } else if (datatype == SpectralDataType.DataType.DIA_F_Window || datatype == SpectralDataType.DataType.pSMART || datatype == SpectralDataType.DataType.WiSIM) {
                                int stopidx = temp.indexOf("</precursorMz>");
                                if (stopidx == -1) {
                                    Logger.getRootLogger().error("Parsing </precursorMz> failed. scan number :" + scanno);                                    
                                    System.exit(3);
                                }
                                int startidx = 0;
                                for (int i = stopidx; i > 0; i--) {
                                    if (temp.charAt(i) == '>') {
                                        startidx = i + 1;
                                        break;
                                    }
                                }
                                float precursormz = Float.parseFloat(temp.substring(startidx, stopidx));
                                //By default, assuming it's 5600 data, 
                                //and assume the precursor m/z is at 0.25 * window size Da to the lower bound of isolation window
                                float Loffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.2f;
                                float Roffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.8f;
                                
                                //If the scan contains "windowWideness", then it is a Thermo data, overwrite the isolation window ranges
                                if (temp.contains("windowWideness=\"")) {
                                    startidx = temp.indexOf("windowWideness=\"") + 16;
                                    stopidx = startidx + 1;
                                    for (int i = startidx + 1; i < temp.length(); i++) {
                                        if (temp.charAt(i) == '\"') {
                                            stopidx = i;
                                            break;
                                        }
                                    }
                                    float windowwideness = Float.parseFloat(temp.substring(startidx, stopidx));
                                     //Again assume the precursor m/z is at the center of isolation window, because it is a Thermo data
                                    Loffset = windowwideness / 2f;
                                    Roffset = windowwideness / 2f;
                                }

                                if (!dIA_Setting.DIAWindows.containsKey(new XYData(precursormz - Loffset, precursormz + Roffset))) {
                                    ArrayList<Integer> scanList = new ArrayList<>();
                                    dIA_Setting.DIAWindows.put(new XYData(precursormz - Loffset, precursormz + Roffset), scanList);
                                }
                                dIA_Setting.DIAWindows.get(new XYData(precursormz - Loffset, precursormz + Roffset)).add(scanno);
                            } else if (datatype == SpectralDataType.DataType.DIA_V_Window) {
                                //if the DIA data is variable window size setting, then use the pre-defined setting
                                int stopidx = temp.indexOf("</precursorMz>");
                                int startidx = 0;
                                for (int i = stopidx; i > 0; i--) {
                                    if (temp.charAt(i) == '>') {
                                        startidx = i + 1;
                                        break;
                                    }
                                }
                                float precursormz = Float.parseFloat(temp.substring(startidx, stopidx));
                                for (XYData window : dIA_Setting.DIAWindows.keySet()) {
                                    if (window.getX() <= precursormz && window.getY() >= precursormz) {
                                        dIA_Setting.DIAWindows.get(window).add(scanno);
                                        break;
                                    }
                                }
                            } else if (datatype == SpectralDataType.DataType.MSe) {
                                float mzlowF = 0f;
                                float mzhighF = 10000f;
                                if (!dIA_Setting.DIAWindows.containsKey(new XYData(mzlowF, mzhighF))) {
                                    ArrayList<Integer> scanList = new ArrayList<>();
                                    dIA_Setting.DIAWindows.put(new XYData(mzlowF, mzhighF), scanList);
                                }
                                dIA_Setting.DIAWindows.get(new XYData(mzlowF, mzhighF)).add(scanno);
                            }
                        }
                    }
                } else {
                    Logger.getRootLogger().error("index of mzXML error");
                    System.exit(1);
                }
                ElutionTimeToScanNoMap.put(rt, scanno);
                ScanToElutionTime.put(scanno, rt);
                MsLevelList.put(scanno, mslevel);
            }
            ScanToElutionTime.compact();
            fileHandler.close();
        }
    }
     
     //Get all the DIA MS2 scans according to a isolation window range
    @Override
    public ScanCollection GetScanDIAMS2(XYData DIAWindow, boolean IncludePeak, float startTime, float endTime) {
        if (dIA_Setting == null) {
            Logger.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
        ScanCollection swathScanCollection = new ScanCollection(parameter.Resolution);
        List<MzXMLthreadUnit> ScanList = new ArrayList<>();

        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
//        ArrayList<Integer> IncludedScans=new ArrayList<>();
        final BitSet IncludedScans=new BitSet();
        for(int scannum :dIA_Setting.DIAWindows.get(DIAWindow)){
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                IncludedScans.set(scannum, true);
            }
        }
        ScanList=ParseScans(IncludedScans);        
        for (MzXMLthreadUnit result : ScanList) {
            swathScanCollection.AddScan(result.scan);
            swathScanCollection.ElutionTimeToScanNoMap.put(result.scan.RetentionTime, result.scan.ScanNum);
        }        
        ScanList.clear();
        ScanList = null;
        return swathScanCollection;
    }
        
     //Get all the DIA MS1 scans according to MS1 m/z range, this was only for WiSIM data
    @Override
    public ScanCollection GetScanCollectionMS1Window(XYData MS1Window, boolean IncludePeak, float startTime, float endTime)  {
        if (dIA_Setting == null) {
            Logger.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
        ScanCollection MS1WindowScanCollection = new ScanCollection(parameter.Resolution);
       
        List<MzXMLthreadUnit> ScanList = null;

        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
//        ArrayList<Integer> IncludedScans=new ArrayList<>();
        final BitSet IncludedScans=new BitSet();
        for(int scannum : dIA_Setting.MS1Windows.get(MS1Window)){
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                IncludedScans.set(scannum, true);
            }
        }
        
        ScanList=ParseScans(IncludedScans);
        
        for (MzXMLthreadUnit result : ScanList) {
            MS1WindowScanCollection.AddScan(result.scan);
            MS1WindowScanCollection.ElutionTimeToScanNoMap.put(result.scan.RetentionTime, result.scan.ScanNum);
        }
        ScanList.clear();
        ScanList = null;
        
        return MS1WindowScanCollection;
    }
    static private int step = -1;
    //Parse scans given a list of scan numbers
//    private List<MzXMLthreadUnit>  ParseScans(ArrayList<Integer> IncludedScans){
    private List<MzXMLthreadUnit> ParseScans(final BitSet IncludedScans){
         List<MzXMLthreadUnit> ScanList=new ArrayList<>();
         ArrayList<ForkJoinTask<?>> futures = new ArrayList<>();
        final ForkJoinPool fjp= new ForkJoinPool(NoCPUs);
        Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();        
        Entry<Integer, Long> ent = iter.next();
        long currentIdx = ent.getValue();
        int nextScanNo = ent.getKey();
        final RandomAccessFile fileHandler;
        try{fileHandler = new RandomAccessFile(filename, "r");}
        catch(FileNotFoundException e){throw new RuntimeException(e);}
        byte[] buffer = new byte[1<<10];
        if(step==-1)
            step=fjp.getParallelism()*32;
        while (iter.hasNext()) {
            ent = iter.next();
            long startposition = currentIdx;
            long nexposition = ent.getValue();
            int currentScanNo = nextScanNo;
            nextScanNo = ent.getKey();
            currentIdx = nexposition;

            if (IncludedScans.get(currentScanNo)) {
                try {
                    final int bufsize =  (int) (nexposition - startposition);
                    if(buffer.length<bufsize)
                        buffer = new byte[Math.max(bufsize,buffer.length<<1)];
//                    byte[] buffer = new byte[bufsize];
//                    RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
                    fileHandler.seek(startposition);
                    fileHandler.read(buffer, 0, bufsize);
//                    fileHandler.close();
//                    String xmltext = new String(buffer);
                    String xmltext = new String(buffer,0,bufsize,StandardCharsets.ISO_8859_1);
                    if (ent.getKey() == Integer.MAX_VALUE) {
                        xmltext = xmltext.replaceAll("</msRun>", "");
                    }
                    boolean ReadPeak = true;
                    final MzXMLthreadUnit unit = new MzXMLthreadUnit(xmltext, parameter, datatype, ReadPeak);
                    futures.add(fjp.submit(unit));
                    ScanList.add(unit);

                    if ((ScanList.size() % step) == 0) {
                        futures.get(futures.size()-step).get();
                        if (iter.hasNext() && fjp.getActiveThreadCount() < fjp.getParallelism()) {
                            step *= 2;
//                            System.out.println("MzXMLthreadUnit: fjp.getActiveThreadCount()\t" + fjp.getActiveThreadCount()+"\t"+step);
                        }
                    }
                } catch (Exception ex) {
                    Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                }
            }
        }
        try {fileHandler.close();}
        catch (IOException ex) {throw new RuntimeException(ex);}
        fjp.shutdown();
        try {fjp.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);}
        catch (InterruptedException ex) {throw new RuntimeException(ex);}
//        for (MzXMLthreadUnit unit : ScanList) {
//            executorPool.execute(unit);
//        }
//        executorPool.shutdown();
//
//        try {
//            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//        } catch (InterruptedException e) {
//            Logger.getRootLogger().info("interrupted..");
//        }
        return ScanList;
    }
    
    @Override
    public ScanCollection GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak, float startTime, float endTime) {
        ScanCollection scanCollection = InitializeScanCollection();
        Logger.getRootLogger().debug("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB (" + NoCPUs + " threads)");

//        ArrayList<Integer> IncludedScans = new ArrayList<>();
        final IntArrayList IncludedScans = new IntArrayList();
        
        for(int ScanNum : MsLevelList.keySet()){
            if(MsLevelList.get(ScanNum)==1 && MS1Included){
                IncludedScans.add(ScanNum);
            }
            if(MsLevelList.get(ScanNum)==2 && MS2Included){
                IncludedScans.add(ScanNum);
            }
        }
         
        List<MzXMLthreadUnit> ScanList = null;

        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
        
//        ArrayList<Integer> temp=new ArrayList<>();
        final BitSet temp = new BitSet();
//        for(int scannum : IncludedScans){
        for(int i=0; i<IncludedScans.size(); ++i){
            final int scannum = IncludedScans.get(i);
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                temp.set(scannum, true);
            }
        }
        
        ScanList=ParseScans(temp);
        
        for (MzXMLthreadUnit result : ScanList) {
            scanCollection.AddScan(result.scan);
        }
        ScanList.clear();
        ScanList = null;
        
        System.gc();
        Logger.getRootLogger().debug("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");
        return scanCollection;
    }
 
}
