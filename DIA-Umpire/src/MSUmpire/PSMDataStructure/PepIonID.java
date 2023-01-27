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

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PeakCluster;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.IonFactory;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.general.IsotopicDistribution;
import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.MolecularFormula;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Peptide ion ID data structure
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepIonID implements Serializable {
    private static final long serialVersionUID = 33203837556783L;

    public int Index;
    public String ModSequence;
    public float Weight;
    public boolean Is_NonDegenerate;
    public int Charge;
    public float GroupWeight;
    public transient float FilteringWeight;
    public String Sequence;
    public float MaxProbability = -1f;
    public boolean GlycoMS1Valid=false;
    public boolean GlycoMS2Valid=false;
    public String MS1ClusIndex = "";
    public String MS2ClusIndex = "";
    public int NoNsite;
    public int NoNsiteFragObs=-1;
    public int IsDecoy=-1;//added 0828, needs to be deleted for older serialization
        
    private transient HashSet<String> PSMSpecMap;
    private ArrayList<PSM> PSMList;
    public ArrayList<String> ParentProtID_PepXML = new ArrayList<>();
    public transient ArrayList<String> ParentProtString_ProtXML = new ArrayList<>();
    public ArrayList<ProtID> ParentProtID_ProtXML = new ArrayList<>();
    public ArrayList<PeakCluster> MS1PeakClusters = new ArrayList<>();
    public ArrayList<PeakCluster> MS2UnfragPeakClusters = new ArrayList<>();
    public float PeakArea[];
    public float PeakHeight[];
    public float PeakClusterScore;
    private float RT = -1f;
    public float PeakRT = -1f;
    public ArrayList<Float> PredictRT = new ArrayList<>();
    
    private float RTSD = -1f;
    private ArrayList<Ion> Fragments;
    public ArrayList<ModificationMatch> Modifications = new ArrayList<>();
    public ArrayList<FragmentPeak> FragmentPeaks = new ArrayList<>();    
    private Peptide peptide = null;
    private float mz = -1f;
    public float ObservedMz;
    public String TPPModSeq;
    public float MS1AlignmentProbability = -1f;
    public float MS2AlignmentProbability = -1f;
    public float UScoreProbability_MS1= -1f;
    public float UScoreProbability_MS2= -1f;

    public Peptide GetPepFactory() {
        if (peptide == null) {
            ReFreshPepFactory();
        }
        return peptide;
    }
       
    public HashSet<String> GetPSMSpecMap(){
        if(PSMSpecMap==null){
            PSMSpecMap=new HashSet<>();
            for(PSM psm : PSMList){
                PSMSpecMap.add(psm.SpecNumber);
            }
        }
        return PSMSpecMap;
    }
    
    transient float intensity=-1f;
    public float GetTotalPeakHeight(){
        if(intensity<=0f){
            intensity=0f;
            for(float apx : PeakHeight){
                intensity+=apx;
            }
        }
        return intensity;
    }
    
    transient String skylinemod="";
            
    
    public float GetMS1PPM(){
        float ppm = InstrumentParameter.CalcSignedPPM(ObservedMass(), CalcNeutralPepMass());
        return ppm;
    }
    
    transient int fragcount=-1;
    public int GetFragCount(){
        if (fragcount == 0) {
            HashSet<String> FragMap=new HashSet<>();
            for (FragmentPeak peak : FragmentPeaks) {
                FragMap.add(peak.IonType);
            }
            fragcount=FragMap.size();
        }
        return fragcount;
    }
    
    public float GetMonoIsMS1() {
        if (PeakHeight != null) {
            return PeakHeight[0];
        }
        return 0f;
    }
    
    public float GetPeakMz(int isotopicpeak){
        return ObservedMz+ (float)((isotopicpeak * ElementaryIon.proton.getTheoreticMass() / Charge));
    }
    
    public float GetRTRange(){
        float range=0f;
        for(PeakCluster cluster : MS1PeakClusters){
            float temp=cluster.endRT-cluster.startRT;
            if(temp>range){
                range=temp;                
            }
        }
        for(PeakCluster cluster : MS2UnfragPeakClusters){
            float temp=cluster.endRT-cluster.startRT;
            if(temp>range){
                range=temp;                
            }
        }
        return  range;
    }
    
    public int getNMissedCleavages() {
        int cpt = 0;
        for (int i = 0; i < Sequence.length() - 1; i++) {
            if ((Sequence.charAt(i) == 'K' || Sequence.charAt(i) == 'R') && Sequence.charAt(i + 1) != 'P') {
                cpt++;
            }
        }
        return cpt;
    }
    
    public void ClearPepFragFactory(){
        peptide=null;
        Fragments=null;
    }
 
    public float TargetedProbability(){
        return Math.max(UScoreProbability_MS1, UScoreProbability_MS2);
    }

    public void ReFreshPepFactory() {
        peptide = new Peptide(Sequence, Modifications);
        mz = -1f;
    }

    public String GetMS1ClusIndex() {
        if ("".equals(MS1ClusIndex)) {
            for (PeakCluster cluster : MS1PeakClusters) {
                MS1ClusIndex += cluster.Index + ";";
            }
        }
        return MS1ClusIndex;
    }

    public String GetMS2ClusIndex() {
        if ("".equals(MS2ClusIndex)) {
            for (PeakCluster cluster : MS2UnfragPeakClusters) {
                MS2ClusIndex += cluster.Index + ";";
            }
        }
        return MS2ClusIndex;
    }

    public PepIonID ClonePepIonID() {
        PepIonID newpeptide = new PepIonID();
        newpeptide.Sequence = Sequence;
        newpeptide.Charge = Charge;
        newpeptide.Is_NonDegenerate = Is_NonDegenerate;
        newpeptide.ModSequence = ModSequence;
        newpeptide.Modifications = (ArrayList<ModificationMatch>)Modifications.clone();
        if (PeakArea != null) {
            newpeptide.PeakArea = new float[PeakArea.length];
            newpeptide.PeakHeight = new float[PeakHeight.length];
        }
        return newpeptide;
    }

    
    public ArrayList<Ion> GetFragments() {
        if (Fragments == null) {
            SetFragments();
        }
        return Fragments;
    }
    
    public float GetPepAbundanceByTopCorrFragAcrossSample(ArrayList<String> PepFrag) {
        float totalabundance = 0f;
        int count=0;
        if (PepFrag!=null && FragmentPeaks != null) {
            for (FragmentPeak frag : FragmentPeaks) {
                if (PepFrag.contains(frag.GetFragKey())) {
                    totalabundance += frag.intensity;
                    count++;
                }
            }
        }
        return totalabundance;
    }
    
    public float GetPepAbundanceByTopFragments(int topN) {        
        ArrayList<FragmentPeak> includelist = new ArrayList<>();
        for (int i = 0; i < topN; i++) {
            FragmentPeak bestfrag=null;
            for (FragmentPeak fragment : FragmentPeaks) {                
                if (fragment.corr>0.7f && !includelist.contains(fragment) && (bestfrag == null || fragment.intensity > bestfrag.intensity)) {
                    bestfrag = fragment;
                }
            }
            if (bestfrag != null) {
                includelist.add(bestfrag);
            }
        }
        float totalabundance = 0f;
        for (FragmentPeak frag : includelist ) {
            totalabundance+=frag.intensity;
        }
        return totalabundance;
    }
    
    public String GetKey() {
        return ModSequence + "_" + Charge;
    }
    
    public boolean IsDecoy(String decoytag) {
        if (IsDecoy == -1) {
            IsDecoy = 1;
            for (String pro : ParentProtID_PepXML) {
                if (!(pro.startsWith(decoytag) | pro.endsWith(decoytag))) {
                    IsDecoy = 0;
                    break;
                }
            }
        }
        return IsDecoy==1;
    }

    public float LookUpFragmentMZ(PeptideFragmentIon.IonType fragmentIonType, int num) {        
        for (Ion frag : GetFragments()) {
            if (frag.getType() == fragmentIonType && "".equals(frag.getNeutralLossesAsString()) && ((PeptideFragmentIon) frag).getNumber() == num) {
                return (float) (frag.getTheoreticMass() + ElementaryIon.proton.getTheoreticMass());
            }
        }
        return 0f;
    }
    
    private void SetFragments() {
        peptide = GetPepFactory();
        HashMap<Integer, HashMap<Integer, ArrayList<Ion>>> allfragment = IonFactory.getInstance().getFragmentIons(peptide);
        Fragments = new ArrayList<>();
        Fragments.addAll(allfragment.get(Ion.IonType.PEPTIDE_FRAGMENT_ION.index).get(PeptideFragmentIon.B_ION));
        Fragments.addAll(allfragment.get(Ion.IonType.PEPTIDE_FRAGMENT_ION.index).get(PeptideFragmentIon.Y_ION));
    }

    public void SetRTSD(float value){
        if(RTSD<value){
            RTSD=value;
        }
    }
    
    public float CalcRT() {
        if(PSMList.isEmpty()){
            return -1;
        }
        
        RT = 0f;
        for (PSM psm : PSMList) {
            RT += psm.NeighborMaxRetentionTime;
        }
        RT /= PSMList.size();

        RTSD = 0f;
        for (PSM psm : PSMList) {
            RTSD += (psm.NeighborMaxRetentionTime - RT) * (psm.NeighborMaxRetentionTime - RT);
        }

        RTSD /= PSMList.size();
        RTSD = (float) Math.sqrt(RTSD);
        return RT;
    }

    public void SetRT(float value) {
        RT = value;
    }

    public float GetRTSD() {
        if (RTSD == -1) {
            CalcRT();
        }
        return RTSD;
    }

    public float GetRT() {
        if (PeakRT != -1) {
            return PeakRT;
        }
        return GetIDRT();
    }
    
    
    public float GetAvgPredictRT(){
        float avg=0f;
        for(float rt : PredictRT){
            avg+=rt;
        }
        return avg/=PredictRT.size();
    }

    
    public String PredictRTString(){
        String rtout="";
        for(float rt : PredictRT){
            rtout+=rt+";";
        }
        return rtout;
    }

    public float GetIDRT() {
        if (RT == -1) {
            CalcRT();
        }
        return RT;
    }

    public void CreateQuantInstance(int NoPeakClusters) {
        PeakArea = new float[NoPeakClusters];
        PeakHeight = new float[NoPeakClusters];
    }

    public PepIonID() {
        PSMList = new ArrayList<>();
    }

    public void AddPSM(PSM psm) {
        if (!GetPSMSpecMap().contains(psm.SpecNumber)) {
            if (psm.Probability > MaxProbability) {
                MaxProbability = psm.Probability;
                ObservedMz = psm.ObserPrecursorMz();
            }
            PSMList.add(psm);
            GetPSMSpecMap().add(psm.SpecNumber);
            psm.pepIonID = this;
            for (String protein : psm.ParentProtIDs) {
                if (!ParentProtID_PepXML.contains(protein)) {
                    ParentProtID_PepXML.add(protein);
                }
            }

            if (Modifications.isEmpty()) {
                AddModifications(psm);
            }
            if (Modifications.size() != psm.Modifications.size()) {
                LogManager.getRootLogger().error("Modification size doesn't match.........:" + psm.SpecNumber);
            }
        }
    }

    public void AddModifications(PSM psm) {
        for (ModificationMatch mod : psm.Modifications) {
            Modifications.add(mod);
        }
    }

    public int GetSpectralCount(){
        return PSMList.size();
    }
    
    public ArrayList<PSM> GetPSMList() {
        return PSMList;
    }

    public void SetInfobyPSM(PSM psm) {
        this.Charge = psm.Charge;
        this.ModSequence = psm.ModSeq;
        this.TPPModSeq = psm.TPPModSeq;
        this.Sequence = psm.Sequence;
        this.MaxProbability = psm.Probability;
        this.ObservedMz = psm.ObserPrecursorMz();
    }

    public void SetMz(float mz) {
        this.mz = mz;
    }

    public float ObservedMass() {
        return Charge * (float)((ObservedMz - ElementaryIon.proton.getTheoreticMass()));
    }

    public float NeutralPrecursorMz() {
        if (mz == -1f) {
            mz = (float)((CalcNeutralPepMass() + Charge * ElementaryIon.proton.getTheoreticMass())) / Charge;
        }
        return mz;
    }

    public float CalcNeutralPepMass() {
        return GetPepFactory().getMass().floatValue();
    }

    public String GetModificationString() {
        String ModificationString = "";
        for (ModificationMatch mod : Modifications) {
            ModificationString += mod.getTheoreticPtm() + "(" + mod.getModificationSite() + ");";
        }
        return ModificationString;
    }

    public MolecularFormula GetMolecularFormula() {
        MolecularFormula formula = new MolecularFormula(GetAASequenceImpl());        
        return formula;
    }
    AASequenceImpl AAimple = null;

    public AASequenceImpl GetAASequenceImpl() {
        if (AAimple == null) {
            AAimple = new AASequenceImpl(Sequence);
        }
        return AAimple;
    }
    IsotopicDistribution calc = null;

    public float[] IsotopicDistrubtion(int NoOfIsoPeaks) {
        IsotopicDistribution calc = GetIsotopicDistribution();
        Double[] isopeak = calc.getPercMax();
        float firstPattern = (float) (double) isopeak[0];

        float[] TheoIso = new float[NoOfIsoPeaks];
        for (int i = 0; i < NoOfIsoPeaks; i++) {
            TheoIso[i] = (float) (double) (isopeak[i] / firstPattern);
        }
        return TheoIso;
    }

    private IsotopicDistribution GetIsotopicDistribution() {
        if (calc == null) {
            calc = new IsotopicDistribution(GetMolecularFormula());
            //calc.calculate();
        }
        //AAimple.addModification(new ModificationImplementation(Modifications , Modifications, new HashMap(), 0));
        return calc;
    }

    public String ParentProteins() {
        String proteins="";
        for(ProtID prot : ParentProtID_ProtXML){
            proteins+=prot.getAccNo()+";";
        }
        return proteins;
    }

    public float GetMS1() {
        if (PeakHeight != null) {
            return PeakHeight[0];
        }
        return 0f;
    }

    public void RemoveRedundantFrag() {
        HashMap<String, FragmentPeak> UniqueFrag = new HashMap<>();
        for (FragmentPeak frag : FragmentPeaks) {
            if (!UniqueFrag.containsKey(frag.GetFragKey())) {
                UniqueFrag.put(frag.GetFragKey(), frag);
            } else if (UniqueFrag.get(frag.GetFragKey()).intensity * UniqueFrag.get(frag.GetFragKey()).corr < frag.intensity * frag.corr) {
                UniqueFrag.put(frag.GetFragKey(), frag);
            }
        }
        FragmentPeaks.clear();
        for (FragmentPeak frag : UniqueFrag.values()) {
            FragmentPeaks.add(frag);
        }
    }

}
