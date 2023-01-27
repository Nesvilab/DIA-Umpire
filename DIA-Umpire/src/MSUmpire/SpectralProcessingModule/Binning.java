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

import MSUmpire.BaseDataStructure.XYPointCollection;
import java.util.ArrayList;

/**
 * Bin peaks in a spectrum for the purpose of calculating spectral similarity 
 * @author Chih-Chiang Tsou
 */
public class Binning {
    
    public XYPointCollection Binning(XYPointCollection scan, float threshold, ArrayList<Integer> removeList) {

        float Binsize = 0.1f;
        XYPointCollection returnCollection = new XYPointCollection();
        int ArraySize = (int) Math.ceil((scan.Data.get(scan.PointCount() - 1).getX() - scan.Data.get(0).getX()) / Binsize) + 2;

        float[] mzarray = new float[ArraySize];
        float[] valueindex = new float[ArraySize];

        for (int i = 0; i < ArraySize; i++) {
            valueindex[i] = scan.Data.get(0).getX() + Binsize * i;
        }
        int arrayidx = 1;
        for (int i = 0; i < scan.PointCount(); i++) {
            if (scan.Data.get(i).getY() > threshold) {
                
                while (scan.Data.get(i).getX() > valueindex[arrayidx]) {
                    arrayidx++;
                }
                float intensity = scan.Data.get(i).getY();
                float mz = scan.Data.get(i).getX();
                float intenlow = intensity * (Binsize - (mz - valueindex[arrayidx - 1])) / Binsize;
                float intenup = intensity * (Binsize - (valueindex[arrayidx] - mz)) / Binsize;

                if (intenlow > mzarray[arrayidx - 1]) {
                    mzarray[arrayidx - 1] = intenlow;
                }
                if (intenup > mzarray[arrayidx]) {
                    mzarray[arrayidx] = intenup;
                }
            }
        }
        for (int i = 0; i < mzarray.length; i++) {
            if (removeList == null || removeList.isEmpty() || !removeList.contains(i)) {
                if (mzarray[i] > threshold) {
                    returnCollection.AddPoint(i, mzarray[i]);
                }
            }
        }
        return returnCollection;
    }
}
