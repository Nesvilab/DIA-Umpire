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

import MSUmpire.BaseDataStructure.ScanData;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Background detection algorithm
 * @author Chih-Chiang Tsou
 */
public class BackgroundDetector {

    public ScanData Scan;
    private float Ratio = 2f;

    public BackgroundDetector(ScanData Scan) {
        this.Scan = Scan;
    }

    public void DetermineConstantBackground() {
        if (Scan.Data.isEmpty()) {
            return;
        }
        ArrayList<Float> IntList = new ArrayList<>();
        IntList.add(1f);
        for (int i = 0; i < Scan.Data.size(); i++) {
            IntList.add(Scan.Data.get(i).getY());
        }
        Collections.sort(IntList);
        int idx = 1;
        while (idx < IntList.size()) {
            if (IntList.get(idx).equals(IntList.get(idx - 1))) {
            } else {
                Scan.background = IntList.get(idx - 1);
            }
            idx++;
        }

        IntList.clear();
        IntList = null;
    }

    public void AdjacentPeakHistogram() {
        if (Scan.PointCount() < 10) {
            return;
        }
        ArrayList<Float> IntList = new ArrayList<>();
        for (int i = 0; i < Scan.Data.size(); i++) {
            IntList.add(Scan.Data.get(i).getY());
        }        
        Collections.sort(IntList);
        float upper = IntList.get((int) (IntList.size() * 0.7f));
        float lower = IntList.get(0);

        if(upper<=lower+0.001){
            return;
        }
        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        int count4 = 0;
        int noise = 0;

        float bk = 0f;
        float interval = (upper - lower) / 20f;
        
        for (bk = lower; bk < upper; bk += interval) {
            count1 = 0;
            count2 = 0;
            count3 = 0;
            count4 = 0;
            noise = 0;
            int preidx = -1;
            for (int i = 1; i < Scan.Data.size(); i++) {
                if (Scan.Data.get(i).getY() > bk) {
                    if (preidx != -1) {
                        float dist = Scan.Data.get(i).getX() - Scan.Data.get(preidx).getX();
                        //writer.write(dist + "\t");
                        if (dist > 0.95 && dist < 1.05 && Scan.Data.get(preidx).getY() > Scan.Data.get(i).getY()) {
                            count1++;
                        } else if (dist > 0.45 && dist < 0.55 && Scan.Data.get(preidx).getY() > Scan.Data.get(i).getY()) {
                            count2++;
                        } else if (dist > 0.3 && dist < 0.36 && Scan.Data.get(preidx).getY() > Scan.Data.get(i).getY()) {
                            count3++;
                        } else if (dist > 0.24 && dist < 0.26 && Scan.Data.get(preidx).getY() > Scan.Data.get(i).getY()) {
                            count4++;
                        } else if (dist < 0.23f) {
                            noise++;
                        }                        
                    }
                    preidx=i;
                }                
            }
            if (noise < (count1 + count2 + count3 + count4) * Ratio) {
                break;
            }
        }
        if (bk > 0f) {
            Scan.background = bk;
            Scan.RemoveSignalBelowBG();
        }
    }
}
