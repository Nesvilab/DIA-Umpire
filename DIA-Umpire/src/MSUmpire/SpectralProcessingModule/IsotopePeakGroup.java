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

import MSUmpire.BaseDataStructure.XYData;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.util.ArrayList;

/**
 * Isotope peak group class
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class IsotopePeakGroup {

    public int Charge;
    public ArrayList<XYData> PeakGroupList = new ArrayList<>();

    public IsotopePeakGroup(int charge) {
        this.Charge = charge;
    }

    public float NeutralMass(){
        return Charge * (float)((PrecursorMz() - ElementaryIon.proton.getTheoreticMass()));
    }
    
    public float PrecursorMz() {
        return PeakGroupList.get(0).getX();
    }

    public XYData GetPeakXYPointByPeakidx(int pkidx) {
        return PeakGroupList.get(pkidx);
    }

    public void AddPeak(XYData peakPoint) {
        PeakGroupList.add(peakPoint);
    }
    public float[] GetIsoPattern() {
        float[] Pattern = new float[3];
        for (int i = 0; i < 3; i++) {
            if (i == 0) {
                Pattern[i] = 1;
            }
            Pattern[i] = PeakGroupList.get(i).getY() / PeakGroupList.get(0).getY();
        }
        return Pattern;
    }
        
}
