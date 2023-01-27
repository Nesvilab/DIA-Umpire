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

import java.io.Serializable;
import java.util.ArrayList;

/**
 * Precursor-fragment group class
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FragmentPeakGroup extends PeptideFragment implements Serializable{
    private static final long serialVersionUID = 736478436L;

    public ArrayList<Float> IntensityGroup = new ArrayList<>();
    public ArrayList<Float> CorrGroup = new ArrayList<>();
    public ArrayList<Float> PPMGroup = new ArrayList<>();
    public ArrayList<Float> ApexDeltaGroup = new ArrayList<>();
    public ArrayList<Float> RTOverlapPGroup = new ArrayList<>();
    
    private float AvgInt=-1f;
    
    public void ClearGroups(){
        GetAvgInt();
        IntensityGroup=null;
        CorrGroup=null;
        PPMGroup=null;
        ApexDeltaGroup=null;
        RTOverlapPGroup=null;                
    }
    
    public float GetAvgInt(){
        if(AvgInt==-1){
            AvgInt=0;
            for(float intensity : IntensityGroup){
                AvgInt+=intensity;
            }
            AvgInt/=IntensityGroup.size();
            AvgInt=(float) Math.sqrt(AvgInt);
        }        
        return AvgInt;
    }
    
    public String GetCorrString() {
        String output = "";
        for (float corr : CorrGroup) {
            output += corr + ";";
        }
        return output;
    }

    public String GetPPMString() {
        String output = "";
        for (float ppm : PPMGroup) {
            output += ppm + ";";
        }
        return output;
    }

    public String GetIntString() {
        String output = "";
        for (float intensity : IntensityGroup) {
            output += intensity + ";";
        }
        return output;
    }
    
    public String GetApexDeltaString() {
        String output = "";
        for (float delta : ApexDeltaGroup) {
            output += delta + ";";
        }
        return output;
    }
    
    public String GetRTOverlapString() {
        String output = "";
        for (float overlap : RTOverlapPGroup) {
            output += overlap + ";";
        }
        return output;
    }
   
    
    @Override
    public FragmentPeakGroup clone(){
        FragmentPeakGroup fragmentPeakGroup=new FragmentPeakGroup();
        fragmentPeakGroup.AvgInt=AvgInt;
        fragmentPeakGroup.Charge=Charge;
        fragmentPeakGroup.CorrGroup=CorrGroup;
        fragmentPeakGroup.IntensityGroup=IntensityGroup;
        fragmentPeakGroup.ApexDeltaGroup=ApexDeltaGroup;
        fragmentPeakGroup.RTOverlapPGroup=RTOverlapPGroup;
        fragmentPeakGroup.PPMGroup=PPMGroup;
        fragmentPeakGroup.ObservedMZ=ObservedMZ;
        fragmentPeakGroup.IonType=IonType;
        return fragmentPeakGroup;
    }
}
