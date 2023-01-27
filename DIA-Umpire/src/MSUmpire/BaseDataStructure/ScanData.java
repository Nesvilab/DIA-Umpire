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

import MSUmpire.SpectralProcessingModule.BackgroundDetector;
import MSUmpire.SpectralProcessingModule.ScanPeakGroup;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ScanData extends XYPointCollection{
    public int ScanNum;
    public int MsLevel;
    public float RetentionTime;
    public float StartMz;
    public float EndMz;
    public float BasePeakMz;
    public float BasePeakIntensity;
    private float _totIonCurrent = 0f;
    public float PrecursorMz;
    public int PrecursorCharge;
    public String ActivationMethod;
    public float PrecursorIntensity;
    public String Scantype;
    public int precision;
    public String compressionType;
    public boolean centroided;
    public int precursorScanNum;
    public int PeaksCountString;
    public float background = 0f;
    public String MGFTitle;
    public ScanData TopPeakScan;
    public float windowWideness;
    public String scanType;
    public float isolationWindowTargetMz;
    public float isolationWindowLoffset;
    public float isolationWindowRoffset;

    public void Centroiding(int Resolution, float MinMZ) {
        CentroidingbyLocalMaximum(Resolution, MinMZ);
        centroided = true;
    }

    public float PrecursorMass(){
        return (float) (PrecursorCharge * (PrecursorMz - ElementaryIon.proton.getTheoreticMass()));
    }
    
    public XYData GetHighestPeakInMzWindow(float targetmz, float PPM) {
        float lowmz = InstrumentParameter.GetMzByPPM(targetmz, 1, PPM);
        int startidx = GetLowerIndexOfX(lowmz);
        XYData closetPeak = null;
        for (int idx = startidx; idx < Data.size(); idx++) {
            XYData peak = Data.get(idx);
            if (InstrumentParameter.CalcPPM(targetmz, peak.getX()) <= PPM) {
                if (closetPeak == null || peak.getY() > closetPeak.getY()) {
                    closetPeak = peak;
                }
            } else if (peak.getX() > targetmz) {
                break;
            }
        }
        return closetPeak;
    }
    
    public void GenerateTopPeakScanData(int toppeaks) {
        SortedXYCollectionClass Intsorted = new SortedXYCollectionClass();
        for (int i = 0; i < Data.size; i++) {
            Intsorted.add(new XYData(Data.get(i).getY(), Data.get(i).getX()));
        }
        TopPeakScan = new ScanData();
        for (int i = Intsorted.size() - 1; TopPeakScan.PointCount() < toppeaks && i >= 0; i--) {
            XYData peak = (XYData) Intsorted.get(i);
            TopPeakScan.AddPoint(new XYData(peak.getY(), peak.getX()));
        }
    }

    public void Normalization() {
        if (MaxY != 0) {
            for (int i = 0; i < PointCount(); i++) {
                XYData pt = Data.get(i);
                pt.setY(pt.getY() / MaxY);
            }
        }
    }
   

    public void RemoveSignalBelowBG() {
        SortedXYCollectionClass newData = new SortedXYCollectionClass();
        for (int i = 0; i < Data.size(); i++) {
            if (Data.get(i).getY() > background) {
                newData.add(Data.get(i));
            }
        }
        Data.clear();
        Data = newData;
        Data.Finalize();
    }

    public float TotIonCurrent() {
        if (_totIonCurrent == 0f) {
            for (int i = 0; i < PointCount(); i++) {
                _totIonCurrent += Data.get(i).getY();
            }
        }
        return _totIonCurrent;
    }

    public void SetTotIonCurrent(float totioncurrent) {
        _totIonCurrent = totioncurrent;
    }

    public ScanData CloneScanData() {
        ScanData newscanData = new ScanData();
        for (XYData pt : Data) {
            newscanData.AddPoint(pt);
        }
        newscanData.ScanNum = ScanNum;
        newscanData.MsLevel = MsLevel;
        newscanData.RetentionTime = RetentionTime;
        newscanData.StartMz = StartMz;
        newscanData.EndMz = EndMz;
        newscanData.BasePeakMz = BasePeakMz;
        newscanData.BasePeakIntensity = BasePeakIntensity;
        newscanData.SetTotIonCurrent(_totIonCurrent);
        newscanData.PrecursorMz = PrecursorMz;
        newscanData.PrecursorCharge = PrecursorCharge;
        newscanData.ActivationMethod = ActivationMethod;
        newscanData.PrecursorIntensity = PrecursorIntensity;
        newscanData.Scantype = Scantype;
        newscanData.precision = precision;
        newscanData.compressionType = compressionType;
        newscanData.centroided = centroided;
        return newscanData;
    }

    public ScanData GetNewSubScanBymzRange(float startmz, float endmz) {
        ScanData newScanData = CloneScanData();
        newScanData.Data = GetSubSetByXRange(startmz, endmz).Data;
        return newScanData;
    }

    public float GetTopNIntensity(int N){
     
         ArrayList<Float> IntList = new ArrayList<>();
        for (int i = 0; i < Data.size(); i++) {
            IntList.add(Data.get(i).getY());
        }
        Collections.sort(IntList);
        if(IntList.size()>10){
            return IntList.get(IntList.size()-N);
        }
        return -1;        
    }
    
    public void Preprocessing(InstrumentParameter parameter) {
         if (!centroided) {
            Centroiding(parameter.Resolution, parameter.MinMZ);
        }
        if (parameter.EstimateBG) {
            BackgroundDetector detector=new BackgroundDetector(this);
            //detector.DetermineConstantBackground();
            detector.AdjacentPeakHistogram();
        } else {
            if (MsLevel == 1) {
                background = parameter.MinMSIntensity;
            }
            if (MsLevel == 2) {
                background = parameter.MinMSMSIntensity;
            }
        }

        if (parameter.Denoise) {
            RemoveSignalBelowBG();
        }
        if (parameter.Deisotoping && MsLevel == 1) {
             ScanPeakGroup scanpeak= new ScanPeakGroup(this, parameter);
             scanpeak.Deisotoping();
        }
    }
}
