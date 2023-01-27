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
package MSUmpire.PSMDataStructure;

import com.compomics.util.experiment.identification.matches.ModificationMatch;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Fragments from spectral library for a peptide ion
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepFragmentLib implements Serializable{
    private static final long serialVersionUID = 8372548274L;

    public String ModSequence;
    public String Sequence;
    public String ModificationString;
    public int Charge;
    public float MaxProbability = 0f;
    public float PrecursorMz;
    public HashMap<String, FragmentPeakGroup> FragmentGroups = new HashMap<>();
    public ArrayList<ModificationMatch> Modifications = new ArrayList<>();
    public ArrayList<Float> RetentionTime=new ArrayList<>();
    public float MS1Score;

    public String GetKey() {
        return ModSequence + "_" + Charge;
    }
    
    public HashMap<String, FragmentPeakGroup> CloneFragmentGroup(){
        HashMap<String, FragmentPeakGroup> NewFragmentGroups=new HashMap<>();
        for(FragmentPeakGroup frag : FragmentGroups.values()){
            FragmentPeakGroup newfrag=frag.clone();
            NewFragmentGroups.put(newfrag.GetFragKey(), newfrag);
        }
        return NewFragmentGroups;
    }

    public void AddFragments(ArrayList<FragmentPeak> FragmentPeaks) {
        float topintensity = 0f;
        for (FragmentPeak frag : FragmentPeaks) {
            if (frag.intensity > topintensity) {
                topintensity = frag.intensity;
            }
        }
        for (FragmentPeak fragment : FragmentPeaks) {
            if (!FragmentGroups.containsKey(fragment.GetFragKey())) {
                FragmentPeakGroup frag = new FragmentPeakGroup();
                frag.IonType = fragment.IonType;
                frag.FragMZ = fragment.FragMZ;  
                frag.Charge = fragment.Charge;
                FragmentGroups.put(fragment.GetFragKey(), frag);
            }
            FragmentGroups.get(fragment.GetFragKey()).CorrGroup.add(fragment.corr);
            FragmentGroups.get(fragment.GetFragKey()).IntensityGroup.add(fragment.intensity / topintensity);
            FragmentGroups.get(fragment.GetFragKey()).PPMGroup.add(fragment.ppm);
            FragmentGroups.get(fragment.GetFragKey()).ApexDeltaGroup.add(fragment.ApexDelta);
            FragmentGroups.get(fragment.GetFragKey()).RTOverlapPGroup.add(fragment.RTOverlapP);
        }
    }
}
