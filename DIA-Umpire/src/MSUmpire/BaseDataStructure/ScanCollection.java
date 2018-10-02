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
package MSUmpire.BaseDataStructure;

import java.util.*;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ScanCollection {

    public TreeMap<Integer, ScanData> ScanHashMap;
    public String Filename;
    private int NumScan;
    private int NumScanLevel1;
    private int NumScanLevel2;
    private int StartScan = 1000000;
    private int EndScan = 0;
    private int Resolution;
    private float MinPrecursorInt = Float.MAX_VALUE;
    public TreeMap<Float, Integer> ElutionTimeToScanNoMap;

    public ScanCollection(int Resolution) {
        Comparator<Float> cmp = new Comparator<Float>() {

            @Override
            public int compare(Float o1, Float o2) {
                return (int) (o1 - o2);
            }
        };
        ScanHashMap = new TreeMap<>();
        this.Resolution = Resolution;
        NumScan = 0;
        NumScanLevel1 = 0;
        NumScanLevel2 = 0;
        ElutionTimeToScanNoMap = new TreeMap<>();
    }

    private ArrayList<Integer> ms1ScanIndex = new ArrayList<>();    
    private ArrayList<Integer> ms2ScanIndex = new ArrayList<>();
    
    
    public ArrayList<Integer> GetScanNoArray(int mslevel) {
        if (mslevel == 1) {
            return ms1ScanIndex;
        }
        if (mslevel == 2) {
            return ms2ScanIndex;
        }
        return null;
    }

    public ArrayList<Integer> GetMS2DescendingArray() {
        return ms2ScanIndex;
    }

    public void AddScan(ScanData scan) {
        if (!ScanHashMap.containsKey(scan.ScanNum)) {
            ScanHashMap.put(scan.ScanNum, scan);

            if (scan.MsLevel == 1) {
                NumScanLevel1++;
                ms1ScanIndex.add(scan.ScanNum);
            }
            if (scan.MsLevel == 2) {
                NumScanLevel2++;
                ms2ScanIndex.add((scan.ScanNum));
            }
            NumScan++;
            if (scan.ScanNum >= EndScan) {
                EndScan = scan.ScanNum;
            }
            if (scan.ScanNum <= StartScan) {
                StartScan = scan.ScanNum;
            }
        }
    }

    public ScanData GetParentMSScan(int ScanNo) {        
        Integer preScanNo = null;
        ScanData PreScan = null;
        while ((preScanNo = ScanHashMap.lowerKey(ScanNo)) != null) {
            PreScan = ScanHashMap.get(preScanNo);
            if (PreScan.MsLevel == 1) {
                break;
            }
            ScanNo = preScanNo;
        }
        return PreScan;
    }

    public ScanData GetScan(int ScanNO) {
        return ScanHashMap.get(ScanNO);
    }

    public boolean ScanAdded(int ScanNo) {
        return ScanHashMap.containsKey(ScanNo);
    }

    public void CentoridingAllScans(int Resolution, float MiniIntF) {
        for (ScanData scan : ScanHashMap.values()) {
            scan.Centroiding(Resolution, MiniIntF);
        }
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
    
    public ScanCollection GetSubCollectionByElutionTimeAndMZ(float startTime, float endTime, float startmz, float endmz, int msLevel, boolean IsAddCalibrationScan) {
        ScanCollection scanCollection = new ScanCollection(Resolution);
        scanCollection.ElutionTimeToScanNoMap = ElutionTimeToScanNoMap;
        scanCollection.Filename = Filename;

        if (endTime == -1) {
            endTime = 9999999f;
        }
        
        //Find the start scan num and end scan num
        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetScanNoByRT(startTime);
        EndScanNo = GetScanNoByRT(endTime);

        NavigableMap<Integer, ScanData> SubScaNavigableMap = ScanHashMap.subMap(StartScanNo, true, EndScanNo, true);
        for (ScanData scan : SubScaNavigableMap.values()) {
            if (endmz == -1) {
                if (((msLevel == 0 || scan.MsLevel == msLevel) && (IsAddCalibrationScan == true || scan.Scantype != "calibration")) && scan.PointCount() > 0 && scan.TotIonCurrent() > 0) {
                    scanCollection.AddScan(scan);
                }
            } else //filter mz
            {
                if (((msLevel == 0 || scan.MsLevel == msLevel) && (IsAddCalibrationScan == true || scan.Scantype != "calibration")) && scan.PointCount() > 0 && scan.TotIonCurrent() > 0) {
                    scanCollection.AddScan(scan.GetNewSubScanBymzRange(startmz, endmz));
                }
            }
        }
        return scanCollection;
    }
    private XYPointCollection _tic = null;

    public XYPointCollection GetTIC() {
        if (_tic == null) {
            _tic = new XYPointCollection();
            for (ScanData scan : ScanHashMap.values()) {
                _tic.AddPoint(scan.RetentionTime, scan.TotIonCurrent());
            }
        }
        return _tic;
    }
    private XYPointCollection _basepeak = null;

    public XYPointCollection GetBasePeak() {
        if (_basepeak == null) {
            _tic = new XYPointCollection();
            for (ScanData scan : ScanHashMap.values()) {
                float TIC = 0;
                for (int i = 0; i < scan.PointCount(); i++) {
                    float intensity = scan.Data.get(i).getY();
                    if (intensity > TIC) {
                        TIC = intensity;
                    }
                }
                _tic.AddPoint(scan.RetentionTime, TIC);
            }
        }
        return _tic;
    }

    public XYPointCollection GetXIC(float startMZ, float endMZ) {
        XYPointCollection xic = new XYPointCollection();
        for (ScanData scan : ScanHashMap.values()) {
            float intensity = 0f;
            XYPointCollection submz = scan.GetSubSetByXRange(startMZ, endMZ);
            for (int i = 0; i < submz.PointCount(); i++) {
                intensity += submz.Data.get(i).getY();
            }
            xic.AddPoint(scan.RetentionTime, intensity);
        }
        return xic;
    }

    public float GetElutionTimeByScanNo(int scanNo) {
        return GetScan(scanNo).RetentionTime;
    }

    //Remove peaks whose the intensity low than the threshold
    public void RemoveBackground(int mslevel, float background) {
        for (ScanData scan : ScanHashMap.values()) {
            if(scan.MsLevel==mslevel){
                scan.background=background;
                scan.RemoveSignalBelowBG();
            }
        }
    }
}
