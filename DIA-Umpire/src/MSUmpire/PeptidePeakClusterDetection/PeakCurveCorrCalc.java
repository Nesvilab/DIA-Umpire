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
