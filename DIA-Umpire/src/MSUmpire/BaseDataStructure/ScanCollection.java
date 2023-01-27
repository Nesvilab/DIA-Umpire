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
