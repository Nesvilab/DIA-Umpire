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
package MSUmpire.SpectralProcessingModule;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SortedXYCollectionClass;
import MSUmpire.BaseDataStructure.XYData;
import java.util.ArrayList;

/**
 * Isotope peak groups for a single scan
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ScanPeakGroup {

    public ScanData Scan;
    public ArrayList<IsotopePeakGroup> peakGroupList;
    private InstrumentParameter parameter;
    int startcharge = -1;
    int endcharge = -1;
    int maxNoPeaks = -1;
    int minNoPeaks = -1;
    float PPMThreshold=0f;

    public ScanPeakGroup(ScanData scan, InstrumentParameter parameter) {
        this.Scan = scan;
        this.parameter = parameter;
        peakGroupList = new ArrayList<>();
        if (Scan.MsLevel == 1) {
            startcharge = parameter.StartCharge;
            endcharge = parameter.EndCharge;
            maxNoPeaks = parameter.MaxNoPeakCluster;
            minNoPeaks = parameter.MinNoPeakCluster;
            PPMThreshold=parameter.MS1PPM;
        } else if (Scan.MsLevel == 2) {
            startcharge = parameter.MS2StartCharge;
            endcharge = parameter.MS2EndCharge;
            maxNoPeaks = parameter.MaxMS2NoPeakCluster;
            minNoPeaks = parameter.MinMS2NoPeakCluster;
            PPMThreshold=parameter.MS2PPM;
        }       
    }
    
    public void Deisotoping(){
         PeakGroupDetection();
         RemovedGroupedIsotopicPeaks();
    }

    private void RemovedGroupedIsotopicPeaks() {
        for (IsotopePeakGroup peakGroup : peakGroupList) {
            for (int i = 1; i < peakGroup.PeakGroupList.size(); i++) {
                if (Scan.Data.contains(peakGroup.PeakGroupList.get(i))) {
                    Scan.Data.remove(peakGroup.PeakGroupList.get(i));
                }
            }
        }
        peakGroupList.clear();
        peakGroupList = null;
        Scan.Data.Finalize();
    }

    private void RemoveNonGroupedPeak() {
        Scan.Data = null;
        Scan.Data = new SortedXYCollectionClass();
        for (IsotopePeakGroup peakGroup : peakGroupList) {
            for (XYData pt : peakGroup.PeakGroupList) {
                if (!Scan.Data.contains(pt)) {
                    Scan.AddPoint(pt);
                }
            }
        }
        peakGroupList.clear();
        peakGroupList = null;
        Scan.Data.Finalize();
    }

    public void PeakGroupDetection() {
        for (int i = startcharge; i <= endcharge; i++) {
            IsotopeGroupDetectionForCharge(i);
        }
    }

    private void IsotopeGroupDetectionForCharge(int charge) {
        boolean[] grouped = new boolean[Scan.PointCount()];
        for (int i = 0; i < Scan.PointCount(); i++) {
            grouped[i] = false;
        }
        for (int i = 0; i < Scan.PointCount(); i++) {
            if (!grouped[i]) {
                DetectPeakGroupFromPeakIndex(i, charge, grouped);
            }
        }
        grouped = null;
    }

    private void DetectPeakGroupFromPeakIndex(int startidx, int charge, boolean[] grouped) {
        grouped[startidx] = true;
        float startmz = Scan.Data.get(startidx).getX();

        boolean[] PeakFound = new boolean[maxNoPeaks];
        int[] PeakFoundIdx = new int[maxNoPeaks];
        for (int i = 0; i < maxNoPeaks; i++) {
            PeakFound[i] = false;
        }
        PeakFound[0] = true;
        PeakFoundIdx[0] = startidx;

        for (int peak = 1; peak < maxNoPeaks; peak++) {
            float mzdistance = peak * (1f / charge);
            float minippm = 10000f;
            int closetidx = -1;

            for (int i = startidx + 1; i < Scan.PointCount(); i++) {
                float ppm = InstrumentParameter.CalcPPM(Scan.Data.get(i).getX(), startmz + mzdistance);
                if (ppm < PPMThreshold + (PPMThreshold * peak * 0.5)) {
                    if (ppm < minippm && Scan.Data.get(PeakFoundIdx[peak-1]).getY()>Scan.Data.get(i).getY()) {
                        closetidx = i;
                        minippm = ppm;
                    }
                } else {
                    if (Scan.Data.get(i).getX() > startmz + mzdistance) {
                        break;
                    }
                }
            }
            if (minippm < 10000f) {                
                PeakFound[peak] = true;
                PeakFoundIdx[peak] = closetidx;
                startidx = closetidx;
            } else {
                break;
            }
        }

        IsotopePeakGroup peakGroup = new IsotopePeakGroup(charge);
        peakGroup.AddPeak(Scan.Data.get(PeakFoundIdx[0]));

        //float PrecursorInt = Scan.Data.get(PeakFoundIdx[0]).Y;
        //float ListInt = PrecursorInt;
        //float PrecursorMZ = Scan.Data.get(PeakFoundIdx[0]).X;
        grouped[PeakFoundIdx[0]] = true;
        //float MW = PrecursorMZ * charge - (charge * 1.00727f);

        for (int i = 1; i < maxNoPeaks; i++) {
            if (!PeakFound[i]) {
                break;
            }
            peakGroup.AddPeak(Scan.Data.get(PeakFoundIdx[i]));
            grouped[PeakFoundIdx[i]] = true;
        }
        if (peakGroup.PeakGroupList.size() >= minNoPeaks) {
            peakGroupList.add(peakGroup);
        }
        PeakFound = null;
        PeakFoundIdx = null;
    }
}
