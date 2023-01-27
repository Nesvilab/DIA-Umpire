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
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.MathPackage.PearsonCorr;
import MSUmpire.PeakDataStructure.PeakCurve;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Calculate peak profile peak correlation given two peak curves
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCurveCorrCalc {

    public static float CalPeakCorr_Overlap(PeakCurve peakA, PeakCurve peakB, int Astart, int Aend, int Bstart, int Bend, int NoPeakPerMin) throws IOException {
        return CalPeakCorr_Overlap(peakA, peakB, Astart, Aend, Bstart, Bend, false, NoPeakPerMin);
    }

    public static float CalPeakCorr(PeakCurve peakA, PeakCurve peakB, int NoPointPerMin) throws IOException {        
        PearsonCorr corr = new PearsonCorr();
        float startRT = Math.max(peakA.StartRT(), peakB.StartRT());
        float endRT = Math.min(peakA.EndRT(), peakB.EndRT());
        XYPointCollection PeakACollection = peakA.GetSmoothPeakCollection(startRT, endRT);
        XYPointCollection PeakBCollection = peakB.GetSmoothPeakCollection(startRT, endRT);
        float corre = 0f;
        
        //double corre2 = 0f;
        if (PeakACollection.Data.size() > 0 && PeakBCollection.Data.size() > 0) {
            corre = corr.CalcCorr(PeakACollection, PeakBCollection, NoPointPerMin);   
        }
        return corre;
    }
    public static float CalPeakCorr_Overlap(PeakCurve peakA, PeakCurve peakB, int Astart, int Aend, int Bstart, int Bend, boolean output, int NoPeakPerMin) throws IOException {
        PearsonCorr corr = new PearsonCorr();
        float startRT = Math.max(peakA.GetPeakRegionList().get(Astart).getX(), peakB.GetPeakRegionList().get(Bstart).getX());
        float endRT = Math.min(peakA.GetPeakRegionList().get(Aend).getZ(), peakB.GetPeakRegionList().get(Bend).getZ());
        XYPointCollection PeakACollection = peakA.GetSmoothPeakCollection(startRT, endRT);
        XYPointCollection PeakBCollection = peakB.GetSmoothPeakCollection(startRT, endRT);
        float corre = 0f;
        if (PeakACollection.Data.size() > 0 && PeakBCollection.Data.size() > 0) {
            corre = corr.CalcCorr(PeakACollection, PeakBCollection, NoPeakPerMin);
            if (output) {
                FileWriter writer = new FileWriter("PeakA.csv");
                for (int i = 0; i < PeakACollection.PointCount(); i++) {
                    writer.write(PeakACollection.Data.get(i).getX() + "," + PeakACollection.Data.get(i).getY() + "\n");
                }
                writer.close();
                writer = new FileWriter("PeakB.csv");
                for (int i = 0; i < PeakBCollection.PointCount(); i++) {
                    writer.write(PeakBCollection.Data.get(i).getX() + "," + PeakBCollection.Data.get(i).getY() + "\n");
                }
                writer.close();
            }
        }
        return corre;
    }
}
