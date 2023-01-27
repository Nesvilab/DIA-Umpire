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
import com.compomics.util.general.IsotopicDistribution;
import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.MolecularFormula;
import java.io.Serializable;
import java.util.ArrayList;

/**
 * PSM data structure 
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PSM implements Serializable{
    private static final long serialVersionUID = 1098763624628L;

    public float Probability;
    public String Sequence;
    public int Charge;
    public String ModSeq = "";
    public String TPPModSeq = "";
    public int Rank;
    public int ScanNo;
    public ArrayList<String> ParentProtIDs;
    public String RawDataName;
    public String SpecNumber;
    public float RetentionTime=-1;
    public float NeighborMaxRetentionTime;
    public ArrayList<ModificationMatch> Modifications;
    public String PreAA;
    public String NextAA;
    public int MissedCleavage;
    public float MassError;
    public float hyperscore;
    public float nextscore;
    public float bscore;
    public float yscore;
    public float cscore;
    public float zscore;
    public float ascore;
    public float xscore;
    public float expect;
    public PepIonID pepIonID;

    public String GetRawNameString() {
        return SpecNumber.substring(0, SpecNumber.indexOf("."));
    }

    public PSM() {
        Modifications = new ArrayList<>();
        ParentProtIDs = new ArrayList<>();
    }

    public String GetModificationString() {
        String ModificationString = "";
        for (ModificationMatch mod : Modifications) {
            ModificationString += mod.getTheoreticPtm() + "(" + mod.getModificationSite() + ");";
        }
        return ModificationString;
    }
    
    public String GetPepKey() {
        return ModSeq + "_" + Charge;
    }

    public void AddParentProtein(String protein) {
        if (!ParentProtIDs.contains(protein)) {
            ParentProtIDs.add(protein);
        }
    }

    public boolean IsDecoy(String decoytag) {
        boolean decoy = true;
        for (String pro : ParentProtIDs) {
            if (!(pro.startsWith(decoytag)|pro.endsWith(decoytag))) {
                decoy = false;
            }
        }
        return decoy;
    }

    public float NeutralPrecursorMz() {
        return (NeutralPepMass + Charge * 1.00727f) / Charge;
    }

    public float GetObsrIsotopicMz(int pkdix) {
        //pkdix starts from 0
        return ObserPrecursorMz() + (float) (pkdix) / Charge;
    }

    public float NeutralPepMass;
    public float ObserPrecursorMass;

    public float ObserPrecursorMz() {
        return (ObserPrecursorMass + Charge * 1.00727f) / Charge;
    }

    private MolecularFormula GetMolecularFormula() {
        MolecularFormula formula = new MolecularFormula(GetAASequenceImpl());
        //formula.addMolecularFormula(GetModMolecularFormula());
        return formula;
    }
    
    AASequenceImpl AAimple;
    private AASequenceImpl GetAASequenceImpl() {
        if (AAimple == null) {
            AAimple = new AASequenceImpl(Sequence);
        }
        return AAimple;
    }
    
    

    float[] TheoIso;
    public float[] IsotopicDistrubtionRatio(int NoOfIsoPeaks) {
        if (TheoIso == null) {
            IsotopicDistribution calc = GetIsotopicDistribution();
            Double[] isopeak = calc.getPercMax();
            float firstPattern = (float) (double) isopeak[0];

            TheoIso = new float[NoOfIsoPeaks];
            for (int i = 0; i < NoOfIsoPeaks; i++) {
                TheoIso[i] = (float) (double) (isopeak[i] / firstPattern);
            }
        }
        return TheoIso;
    }

    IsotopicDistribution calc;
    private IsotopicDistribution GetIsotopicDistribution() {
        if (calc == null) {
            calc = new IsotopicDistribution(GetMolecularFormula());
            calc.calculate();
        }
        //AAimple.addModification(new ModificationImplementation(Modifications , Modifications, new HashMap(), 0));
        return calc;
    }
}
